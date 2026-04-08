import math
import time
import torch
import tensorly as tl
from utils_torch import find_biggest_flattenings_ps, sv_dimension

tl.set_backend("pytorch")
torch.set_default_dtype(torch.float64)


def multiSPM(
    T,
    rank=None,
    maxiter=5000,
    ntries=5,
    gradtol=1e-15,
    ranksel=1e-4,
    ftol=1e-2,
    flats=None,
    symmetries=None,
):
    """
    multiSPM - Decompose any tensor using the Multi Subspace Power Method

    This function computes the CP decomposition of a symmetric, partially
    symmetric, or asymmetric tensor using the Multi Subspace Power Method
    (multiSPM). It supports automatic or user-defined rank selection,
    symmetry handling, and flattening strategies.

    ** Usage **
        lambdas, factors, symvec, stat = multiSPM(T, ...)

    ** INPUT **
         T : torch.Tensor
             Tensor of dimension d_1 x d_2 x ... x d_n
      rank : int, optional
             Tensor rank. If not provided, it is estimated from singular values.
   maxiter : int, optional
             Maximum number of power iterations (default: 5000).
    ntries : int, optional
             Maximum number of random restarts (default: 5).
   gradtol : float, optional
             Convergence tolerance for the power method (default: 1e-15).
   ranksel : float, callable, or str, optional
             Rule for rank selection based on singular values. Can be:
               - scalar threshold,
               - function handle,
               - "return_sv" to return singular values only.
             (default: 1e-4)
      ftol : float, optional
             Tolerance for restarting the power method (default: 1e-2).
     flats : optional
             Flattenings to use. Can be specified as index groups or boolean
             masks. If None, flattenings maximizing admissible rank are used.
 symmetries : optional
             Tensor symmetry structure. Can be provided as:
               - list of index groups,
               - vector labeling symmetric modes,
               - list of multiplicities per symmetry group.
             (default: asymmetric tensor)

    ** OUTPUT **
    lambdas : torch.Tensor
              Scaling factors of the decomposition.
    factors : list of torch.Tensor
              Factor matrices for each unique mode.
     symvec : torch.Tensor
              Vector encoding the symmetry structure.
       stat : dict
              Statistics of the decomposition:
                - extracttime : time for preprocessing
                - powertime   : time in power iterations
                - deflatetime : time in deflation
                - avgiter     : average iterations per component
                - nrr         : number of restarts
                - totaltime   : total runtime

    ** NOTE **
    - This implementation requires the 'utils.py' file to be available in the
      Python path (for flattening selection and rank estimation helpers).
    - The method follows the multiSPM framework using tensor flattenings and
      recursive power iterations.

    ** Reference **
    K. Wang, J. M. Pereira, J. Kileel, A. Seigal,
    "Multi-subspace power method for decomposing all tensors".
    """
    timer = time.perf_counter()

    if not torch.is_tensor(T):
        T = torch.tensor(T, dtype=torch.float64)
    else:
        T = T.clone()

    device = T.device
    dtype = T.dtype

    dims = list(T.shape)
    order = len(dims)

    symvec, nudims, udims, nsyms, symorder = _parse_symmetries(symmetries, order, dims, device)
    uflats, max_rank, _ = _parse_flats(flats, order, symvec, udims, nsyms, device)

    cs_nsyms = torch.cat([torch.tensor([0], dtype=torch.long, device=device), torch.cumsum(nsyms, dim=0)])
    nflats = uflats.shape[0]
    sym_breaking = nflats == 1

    unord_flats = torch.zeros((nflats, order), dtype=torch.bool, device=device)
    for i in range(nudims):
        for j in range(nflats):
            cnt = int(uflats[j, i].item())
            if cnt > 0:
                unord_flats[j, cs_nsyms[i]:cs_nsyms[i] + cnt] = True

    if rank is not None and rank > max_rank:
        raise ValueError("Requested rank exceeds maximum admissible rank from the chosen flattenings")

    m1 = int(torch.sum(unord_flats[0]).item())

    if sym_breaking:
        comp_flat = torch.zeros(order, dtype=torch.bool, device=device)
        zero_groups = torch.where(uflats[0] == 0)[0]
        for i in zero_groups.tolist():
            comp_flat[cs_nsyms[i]:cs_nsyms[i + 1]] = True

        dim_order = (
            torch.where(unord_flats[0])[0].tolist()
            + torch.where((~comp_flat) & (~unord_flats[0]))[0].tolist()
            + torch.where(comp_flat)[0].tolist()
        )
    else:
        ufb = (uflats > 0)

        reord_sym = (
            torch.where(ufb[0] & ufb[1])[0].tolist()
            + torch.where(ufb[0] & (~ufb[1]))[0].tolist()
            + torch.where((~ufb[0]) & ufb[1])[0].tolist()
            + torch.where((~ufb[0]) & (~ufb[1]))[0].tolist()
        )

        symorder = [reord_sym[i] for i in symorder]
        uflats = uflats[:, reord_sym]
        ufb = ufb[:, reord_sym]
        udims = udims[reord_sym]
        nsyms = nsyms[reord_sym]

        m2 = int(torch.sum(nsyms[ufb[0]]).item())
        mc = int(torch.sum(nsyms[ufb[0] & ufb[1]]).item())
        mu = int(torch.sum(nsyms[ufb[0] | ufb[1]]).item())

        unord_flats = torch.repeat_interleave(uflats, nsyms, dim=1) > 0

        dim_order = (
            torch.where(unord_flats[0] & unord_flats[1])[0].tolist()
            + torch.where(unord_flats[0] & (~unord_flats[1]))[0].tolist()
            + torch.where((~unord_flats[0]) & unord_flats[1])[0].tolist()
            + torch.where((~unord_flats[0]) & (~unord_flats[1]))[0].tolist()
        )

    dim_order_new = dim_order.copy()
    for idx, pos in enumerate(symorder):
        dim_order_new[pos] = dim_order[idx]
    dim_order = dim_order_new

    T = T.permute(dim_order)
    dims = [dims[i] for i in dim_order]

    matT_1 = reshapeF(T, math.prod(dims[:m1]), -1)
    U1, S, Vh1 = torch.linalg.svd(matT_1, full_matrices=False)
    V1 = Vh1.transpose(0, 1)

    if rank is None:
        rank = rank_selector(S, ranksel)
        if rank is None:
            return S

    U1 = U1[:, :rank]
    U1_copy = U1.clone()
    S = S[:rank]
    C = torch.diag(1.0 / (S + 1e-16))
    C_copy = C.clone()
    V1 = V1[:, :rank]

    if not sym_breaking:
        perm2 = list(range(0, mc)) + list(range(m1, mu)) + list(range(mc, m1)) + list(range(mu, order))
        T2 = T.permute(perm2)
        left_dims2 = [dims[i] for i in list(range(0, mc)) + list(range(m1, mu))]
        matT_2 = reshapeF(T2, math.prod(left_dims2), -1)
        U2, S2, Vh2 = torch.linalg.svd(matT_2, full_matrices=False)
        V2 = Vh2.transpose(0, 1)
        U2 = U2[:, :rank]
        V2 = V2[:, :rank]

    factors = [torch.zeros((int(dim.item()), rank), dtype=dtype, device=device) for dim in udims]
    lambdas = torch.zeros(rank, dtype=dtype, device=device)

    lap = time.perf_counter() - timer
    stat = {
        "extracttime": lap,
        "powertime": 0.0,
        "deflatetime": 0.0,
        "avgiter": 0.0,
        "nrr": 0,
    }

    Aks = [torch.zeros(int(dim.item()), dtype=dtype, device=device) for dim in udims]

    def power_method(U, nsyms_local, r_local):
        Ak_inds = [i for i in range(len(nsyms_local)) if int(nsyms_local[i]) > 0]
        best_f = None
        best_Aks = None

        for tries in range(1, ntries + 1):
            for jj in Ak_inds:
                x = torch.randn_like(Aks[jj])
                Aks[jj] = x / (torch.linalg.norm(x) + 1e-16)

            f = torch.tensor(0.0, dtype=dtype, device=device)

            for iter_idx in range(1, maxiter + 1):
                active_Aks = [Aks[jj] for jj in Ak_inds]
                active_nsyms = [int(nsyms_local[jj]) for jj in Ak_inds]
                new_Aks, f, max_shift = power_method_iteration(U, active_Aks, active_nsyms, r_local)
                for loc, jj in enumerate(Ak_inds):
                    Aks[jj] = new_Aks[loc]
                if float(max_shift) < gradtol:
                    break

            stat["avgiter"] += iter_idx

            if 1.0 - float(f) < ftol:
                break
            elif tries == 1 or float(f) > float(best_f):
                stat["nrr"] += 1
                best_f = f.clone()
                best_Aks = [a.clone() for a in Aks]
            else:
                stat["nrr"] += 1
                for j in range(len(Aks)):
                    Aks[j] = best_Aks[j].clone()

        return f if best_f is None else max(f, best_f)

    for k in range(rank - 1, -1, -1):
        t0 = time.perf_counter()
        f = power_method(U1, uflats[0], k + 1)
        stat["powertime"] += time.perf_counter() - t0

        if sym_breaking:
            if torch.any(~(uflats[0] > 0)):
                comp_mult = ((nsyms - uflats[0]) * (uflats[0] > 0).to(nsyms.dtype)).tolist()
                Akpow = tensor_product(Aks, comp_mult)
                U_half = reshapeF(
                    Akpow.reshape(1, -1) @ reshapeF(V1, len(Akpow), -1),
                    -1,
                    rank,
                )
                M = U_half @ U_half.T
                evals, evecs = torch.linalg.eigh(M)
                v = evecs[:, -1]
                f_c = power_method(v, nsyms * (uflats[0] == 0).to(nsyms.dtype), 1)
        elif mc > 0 and m2 > mc:
            mask = (ufb[0] & ufb[1]).to(nsyms.dtype)
            Akpow = tensor_product(Aks, (nsyms * mask).tolist())
            U_half = reshapeF(
                Akpow.reshape(1, -1) @ reshapeF(U2, len(Akpow), -1),
                -1,
                rank,
            )
            M = U_half @ U_half.T
            evals, evecs = torch.linalg.eigh(M)
            v = evecs[:, -1]

            mask = ((~ufb[0]) & ufb[1]).to(nsyms.dtype)
            f_c = power_method(v, nsyms * mask, 1)

            if order > mu:
                Akpow = tensor_product(Aks, (nsyms * mask).tolist())
                U_half = reshapeF(
                    Akpow.reshape(1, -1) @ reshapeF(V1, len(Akpow), -1),
                    -1,
                    rank,
                )
                M = U_half @ U_half.T
                evals, evecs = torch.linalg.eigh(M)
                v = evecs[:, -1]

                mask = ((~ufb[0]) & (~ufb[1])).to(nsyms.dtype)
                f_c = power_method(v, nsyms * mask, 1)
        else:
            mask = (ufb[0] & (~ufb[1])).to(nsyms.dtype)
            Akpow = tensor_product(Aks, (nsyms * mask).tolist())
            U_half = reshapeF(
                Akpow.reshape(1, -1) @ reshapeF(V2, len(Akpow), -1),
                -1,
                rank,
            )
            M = U_half @ U_half.T
            evals, evecs = torch.linalg.eigh(M)
            v = evecs[:, -1]

            mask = ((~ufb[0]) & (~ufb[1])).to(nsyms.dtype)
            f_c = power_method(v, nsyms * mask, 1)

            if m2 > mc:
                Akpow = tensor_product(Aks, (nsyms * mask).tolist())
                U_half = torch.einsum("ijr,j->ir", reshapeF(V1, -1, len(Akpow), rank), Akpow)
                M = U_half @ U_half.T
                evals, evecs = torch.linalg.eigh(M)
                v = evecs[:, -1]

                mask = (ufb[1] & (~ufb[0])).to(nsyms.dtype)
                f_c = power_method(v, nsyms * mask, 1)

        alpha = (tensor_product(Aks, uflats[0].tolist()).reshape(1, -1) @ U1_copy).transpose(0, 1).reshape(-1)
        beta = (tensor_product(Aks, (nsyms - uflats[0]).tolist()).reshape(1, -1) @ V1).transpose(0, 1).reshape(-1)

        Ctbeta = (beta.reshape(1, -1) @ C_copy).transpose(0, 1).reshape(-1)
        lambdas[k] = (
            torch.linalg.norm(alpha) * torch.linalg.norm(beta)
            / (torch.dot(alpha, Ctbeta) + 1e-16)
        )

        for i in range(nudims):
            factors[i][:, k] = Aks[i]

        if k > 0:
            x = get_hh_reflector((beta.reshape(1, -1) @ C).transpose(0, 1).reshape(-1))
            if torch.linalg.norm(x) > 1e-16:
                C = RHR(C, x)
                U1 = RHR(U1, x)

        timenow = time.perf_counter() - timer
        stat["deflatetime"] += timenow - lap
        lap = timenow

    stat["avgiter"] = stat["avgiter"] / rank
    stat["totaltime"] = time.perf_counter() - timer

    if not sym_breaking:
        reordered = [None] * nudims
        for i, j in enumerate(reord_sym):
            reordered[j] = factors[i]
        factors = reordered

    return lambdas, factors, symvec, stat


def reshapeF(x, *shape):
    """
    reshapeF - Reshape tensor using MATLAB-style (column-major) ordering.

    Permutes dimensions before and after reshape so tensor flattening matches
    MATLAB conventions used in the algorithm.
    """
    if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
        shape = tuple(shape[0])
    return x.permute(*reversed(range(x.ndim))).contiguous().view(*reversed(shape)).permute(*reversed(range(len(shape))))


def rank_selector(S, rank_sel):
    """
    rank_selector - Determine rank from singular values.

    Selects rank via a user rule, thresholding, or returns None if only
    singular values are requested.
    """
    if callable(rank_sel):
        return rank_sel(S)
    if isinstance(rank_sel, (int, float)):
        typical = torch.sum(S**2) / torch.sum(torch.abs(S))
        return int(torch.sum(torch.abs(S) > rank_sel * typical).item())
    if rank_sel == "return_sv":
        return None
    raise ValueError("Rank selector option not implemented yet")


def tensor_product(Aks, nsyms=None):
    """
    tensor_product - Compute repeated tensor/Kronecker product of vectors.

    Forms a vectorized product of Aks, repeating each vector according to nsyms.
    Used for contractions in the power method.
    """
    if nsyms is None:
        nsyms = [1] * len(Aks)

    device = Aks[0].device if len(Aks) > 0 else torch.device("cpu")
    dtype = Aks[0].dtype if len(Aks) > 0 else torch.float64

    Akpow = torch.tensor([1.0], dtype=dtype, device=device)
    for i in range(len(Aks)):
        for _ in range(int(nsyms[i])):
            Akpow = reshapeF(Akpow.reshape(-1, 1) @ Aks[i].reshape(1, -1), -1)
    return Akpow


def power_method_iteration(U, Aks, nsyms, r):
    """
    power_method_iteration - Perform one recursive power method update.

    Updates factor vectors via tensor contractions and normalization,
    returning updated vectors, objective value, and maximum change.
    """
    m = len(Aks)

    if m == 1:
        Ak_old = Aks[0]

        if int(nsyms[0]) == 1:
            UAk = reshapeF(U, -1, r)
        else:
            n = Ak_old.shape[0]
            Akpow = tensor_product(Aks, [int(nsyms[0]) - 1])
            UAk = reshapeF(
                Akpow.reshape(1, -1) @ reshapeF(U, -1, n * r),
                n,
                r,
            )

        Ak_new = UAk @ (Ak_old.reshape(1, -1) @ UAk).reshape(-1)
        f = torch.dot(Ak_new.reshape(-1), Ak_old.reshape(-1))

        if int(nsyms[0]) > 1:
            fval = float(f)
            ns = int(nsyms[0])
            if fval <= 2.0 / 3.0:
                c = math.sqrt(1.0 - 1.0 / ns) * (1.0 - fval / 2.0)
            else:
                c = math.sqrt(1.0 - 1.0 / ns) * math.sqrt(max(2.0 * fval * max(1.0 - fval, 0.0), 0.0))
            Ak_new = Ak_new + c * Ak_old

        Ak_new = Ak_new / (torch.linalg.norm(Ak_new) + 1e-16)
        max_shift = torch.max(torch.abs(Ak_new - Ak_old))
        Aks[0] = Ak_new
        return Aks, f, max_shift

    m2 = math.ceil(m / 2)

    Akpow = tensor_product(Aks[:m2], nsyms[:m2])
    U_half = reshapeF(
        Akpow.reshape(1, -1) @ reshapeF(U, len(Akpow), -1),
        -1,
        r,
    )
    Aks_right, _, max_shift = power_method_iteration(U_half, Aks[m2:], nsyms[m2:], r)
    Aks[m2:] = Aks_right

    Akpow = tensor_product(Aks[m2:], nsyms[m2:])
    U_half = torch.einsum("ijr,j->ir", reshapeF(U, -1, len(Akpow), r), Akpow)
    Aks_left, f, max_shift_ = power_method_iteration(U_half, Aks[:m2], nsyms[:m2], r)
    Aks[:m2] = Aks_left

    max_shift = torch.maximum(max_shift, max_shift_)
    return Aks, f, max_shift


def get_hh_reflector(y):
    """
    get_hh_reflector - Construct Householder reflector vector.

    Returns normalized vector defining a reflection used in deflation.
    """
    y = y.clone()
    norm_y = torch.linalg.norm(y)
    if norm_y < 1e-16:
        return y
    s = torch.sign(y[-1])
    if s == 0:
        s = torch.tensor(1.0, dtype=y.dtype, device=y.device)
    y[-1] = y[-1] + norm_y * s
    denom = torch.sqrt(torch.abs(y[-1]) * norm_y)
    if denom < 1e-16:
        return y
    y = y / denom
    return y


def RHR(A, x):
    """
    RHR - Apply right Householder reflection.

    Updates matrix A by reflecting with x from the right (deflation step).
    """
    return A[:, :-1] - (A @ x).reshape(-1, 1) @ x[:-1].reshape(1, -1)


def LHR(A, x):
    """
    LHR - Apply left Householder reflection.

    Updates matrix A by reflecting with x from the left.
    """
    return A[:-1, :] - x[:-1].reshape(-1, 1) @ (x.reshape(1, -1) @ A)


def process_uflats(nsyms, uflats, ranks):
    """
    process_uflats - Select valid flattenings.

    Chooses a symmetry-breaking flattening or a valid pair and returns
    the associated maximum admissible rank.
    """
    ranks = torch.as_tensor(ranks, dtype=torch.int64, device=uflats.device)

    if ranks.numel() > 1:
        sort_idx = torch.argsort(ranks, descending=True)
        ranks = ranks[sort_idx]
        uflats = uflats[sort_idx]

    for j in range(uflats.shape[0]):
        if torch.any((uflats[j] > 0) & ((nsyms - uflats[j]) > 0)):
            return uflats[j:j + 1], int(ranks[j].item()), [j]

        for i in range(j):
            if not torch.all(uflats[i] + uflats[j] == nsyms):
                if int(ranks[i].item()) > int(ranks[j].item()):
                    uf_ind = [i, j]
                    max_rank = int(ranks[j].item())
                else:
                    uf_ind = [j, i]
                    max_rank = int(ranks[i].item())
                return uflats[uf_ind], max_rank, uf_ind

    raise ValueError("Provide either a symmetry breaking flattening, or two flattenings that are not complementary")


def _parse_symmetries(symmetries, order, dims, device):
    """
    _parse_symmetries - Normalize symmetry specification.

    Converts user input into symmetry groups, dimensions, multiplicities,
    and ordering used internally.
    """   
    if symmetries is None or (isinstance(symmetries, (list, tuple)) and len(symmetries) == 0):
        symvec = torch.arange(1, order + 1, dtype=torch.long, device=device)
        nudims = order
    elif isinstance(symmetries, (list, tuple)) and len(symmetries) > 0 and isinstance(symmetries[0], (list, tuple, torch.Tensor)):
        symvec = torch.zeros(order, dtype=torch.long, device=device)
        lsym = len(symmetries)
        for i, grp in enumerate(symmetries, start=1):
            grp_t = torch.as_tensor(grp, dtype=torch.long, device=device)
            if grp_t.numel() == 0:
                continue
            if int(torch.sum(symvec[grp_t] != 0).item()) != 0:
                raise ValueError("Symmetry groups overlap")
            symvec[grp_t] = i
        nudims = lsym + int(torch.sum(symvec == 0).item())
        zero_idx = torch.where(symvec == 0)[0]
        for k, idx in enumerate(zero_idx, start=lsym + 1):
            symvec[idx] = k
    else:
        sym_t = torch.as_tensor(symmetries, dtype=torch.long, device=device)
        if sym_t.numel() == order:
            vals, inverse = torch.unique(sym_t, sorted=True, return_inverse=True)
            symvec = inverse + 1
            nudims = len(vals)
        else:
            nudims = int(sym_t.numel())
            symvec = torch.repeat_interleave(torch.arange(1, nudims + 1, device=device), sym_t)

    udims = torch.zeros(nudims, dtype=torch.long, device=device)
    nsyms = torch.zeros(nudims, dtype=torch.long, device=device)
    symorder = []

    dims_t = torch.as_tensor(dims, dtype=torch.long, device=device)

    for i in range(1, nudims + 1):
        syms_i = torch.where(symvec == i)[0]
        symorder.extend(syms_i.tolist())
        udims[i - 1] = dims_t[syms_i[0]]
        if not torch.all(dims_t[syms_i] == udims[i - 1]):
            raise ValueError("Modes in the same symmetry group must have identical dimensions")
        nsyms[i - 1] = syms_i.numel()

    return symvec, nudims, udims, nsyms, symorder


def _parse_flats(flats, order, symvec, udims, nsyms, device):
    """
    _parse_flats - Process flattenings.

    Builds or converts flattenings, computes their rank bounds,
    and selects a valid configuration.
    """
    if flats is None:
        uflats, flat_ranks = find_biggest_flattenings_ps(udims.tolist(), nsyms.tolist(), 3, device=device)
    else:
        if isinstance(flats, (list, tuple)) and len(flats) > 0 and isinstance(flats[0], (list, tuple, torch.Tensor)) and not (
            isinstance(flats[0], (list, tuple)) and len(flats[0]) == order and all(isinstance(v, bool) for v in flats[0])
        ):
            cell_flats = flats
            flats_bool = torch.zeros((len(cell_flats), order), dtype=torch.bool, device=device)
            for i, grp in enumerate(cell_flats):
                grp_t = torch.as_tensor(grp, dtype=torch.long, device=device)
                if grp_t.numel() > 0:
                    flats_bool[i, grp_t] = True
            flats = flats_bool
        else:
            flats = torch.as_tensor(flats, dtype=torch.bool, device=device)

        flat_ranks = []
        uflats_rows = []
        for i in range(flats.shape[0] - 1, -1, -1):
            row = flats[i].to(torch.long)
            urow = torch.zeros(len(udims), dtype=torch.long, device=device)
            for j in range(order):
                urow[symvec[j] - 1] += row[j]
            uflats_rows.append(urow)

            left_rank = sv_dimension(udims.tolist(), urow.tolist())
            left_rank = left_rank - int(torch.sum(udims[urow > 0]).item())
            right_rank = sv_dimension(udims.tolist(), (nsyms - urow).tolist())
            flat_ranks.append(min(left_rank, right_rank))

        uflats = torch.stack(uflats_rows, dim=0)
        flat_ranks = torch.tensor(flat_ranks, dtype=torch.int64, device=device)

    return process_uflats(nsyms, uflats, flat_ranks)


if __name__ == "__main__":
    torch.manual_seed(0)

    def cp_to_tensor(lambdas, factors, symvec):
        order = len(symvec)
        full_factors = [factors[int(symvec[i]) - 1] for i in range(order)]
        shape = [full_factors[i].shape[0] for i in range(order)]
        T = torch.zeros(*shape, dtype=full_factors[0].dtype, device=full_factors[0].device)

        r = len(lambdas)
        for k in range(r):
            outer = lambdas[k] * full_factors[0][:, k]
            for mode in range(1, order):
                outer = torch.einsum("...,j->...j", outer, full_factors[mode][:, k])
            T = T + outer

        return T

    def relerr(T, T_hat):
        return torch.linalg.norm(reshapeF(T, -1) - reshapeF(T_hat, -1)) / torch.linalg.norm(reshapeF(T, -1))

    print("\n## multiSPM tests ##")

    print("\n## Example 1 - Asymmetric noiseless ##")
    d1, d2, d3, d4 = 18, 16, 14, 12
    r = 80


    A = torch.randn(d1, r, dtype=torch.float64)
    B = torch.randn(d2, r, dtype=torch.float64)
    C = torch.randn(d3, r, dtype=torch.float64)
    D = torch.randn(d4, r, dtype=torch.float64)
    lambdas_true = torch.ones(r, dtype=torch.float64)

    symvec_true = torch.tensor([1, 2, 3, 4], dtype=torch.long)
    T = cp_to_tensor(lambdas_true, [A, B, C, D], symvec_true)

    start = time.perf_counter()
    lambdas_hat, factors_hat, symvec_hat, stats = multiSPM(
        T,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-6,
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat, symvec_hat)
    print("Relative reconstruction error:", relerr(T, T_hat))

    print("\n## Example 2 - Asymmetric noiseless (forced flattenings) ##")
    start = time.perf_counter()
    lambdas_hat, factors_hat, symvec_hat, stats = multiSPM(
        T,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-6,
        flats=([0, 1], [0, 2]),
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat, symvec_hat)
    print("Relative reconstruction error:", relerr(T, T_hat))

    print("\n## Example 3 - Partially symmetric noiseless (first two modes symmetric) ##")
    m, n, p = 18, 15, 13
    r = 100

    A = torch.randn(m, r, dtype=torch.float64)
    B = torch.randn(n, r, dtype=torch.float64)
    C = torch.randn(p, r, dtype=torch.float64)
    lambdas_true = torch.ones(r, dtype=torch.float64)

    symvec_true = torch.tensor([1, 1, 2, 3], dtype=torch.long)
    T = cp_to_tensor(lambdas_true, [A, B, C], symvec_true)

    start = time.perf_counter()
    lambdas_hat, factors_hat, symvec_hat, stats = multiSPM(
        T,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-6,
        symmetries=[1, 1, 2, 3],
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat, symvec_hat)
    print("Relative reconstruction error:", relerr(T, T_hat))

    print("\n## Example 4 - Partially symmetric noiseless (two symmetry groups) ##")
    m, n = 16, 14
    r = 100

    A = torch.randn(m, r, dtype=torch.float64)
    B = torch.randn(n, r, dtype=torch.float64)
    lambdas_true = torch.ones(r, dtype=torch.float64)

    symvec_true = torch.tensor([1, 1, 2, 2], dtype=torch.long)
    T = cp_to_tensor(lambdas_true, [A, B], symvec_true)

    start = time.perf_counter()
    lambdas_hat, factors_hat, symvec_hat, stats = multiSPM(
        T,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-6,
        symmetries=[2, 2],
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat, symvec_hat)
    print("Relative reconstruction error:", relerr(T, T_hat))

    print("\n## Example 5 - Fully symmetric noiseless ##")
    d = 18
    n = 4
    r = 150

    A = torch.randn(d, r, dtype=torch.float64)
    lambdas_true = torch.ones(r, dtype=torch.float64)

    symvec_true = torch.tensor([1] * n, dtype=torch.long)
    T = cp_to_tensor(lambdas_true, [A], symvec_true)

    start = time.perf_counter()
    lambdas_hat, factors_hat, symvec_hat, stats = multiSPM(
        T,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-6,
        symmetries=[n],
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat, symvec_hat)
    print("Relative reconstruction error:", relerr(T, T_hat))

    print("\n## Example 6 - Noisy partially symmetric ##")
    m, n, p = 18, 15, 13
    r = 100

    A = torch.randn(m, r, dtype=torch.float64)
    B = torch.randn(n, r, dtype=torch.float64)
    C = torch.randn(p, r, dtype=torch.float64)
    lambdas_true = torch.ones(r, dtype=torch.float64)

    symvec_true = torch.tensor([1, 1, 2, 3], dtype=torch.long)
    T = cp_to_tensor(lambdas_true, [A, B, C], symvec_true)
    T_noisy = T + 0.05 * torch.randn_like(T)

    start = time.perf_counter()
    lambdas_hat, factors_hat, symvec_hat, stats = multiSPM(
        T_noisy,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-2,
        symmetries=[1, 1, 2, 3],
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat, symvec_hat)
    print("Relative reconstruction error vs clean:", relerr(T, T_hat))
    print("Relative reconstruction error vs noisy:", relerr(T_noisy, T_hat))

    print("\n## Example 7 - Singular values only ##")
    sv = multiSPM(
        T,
        rank=None,
        ranksel="return_sv",
        symmetries=[1, 1, 2, 3],
    )
    print("Returned singular values shape:", tuple(sv.shape))
    print("Leading singular values:", sv[: min(10, len(sv))])