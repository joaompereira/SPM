import math
import time
import torch
import tensorly as tl
from torch.linalg import norm
from utils_torch import find_best_flatpair

tl.set_backend("pytorch")


def MSPM_asym(
    T,
    rank=None,
    maxiter=5000,
    ntries=5,
    gradtol=1e-15,
    ranksel=1e-4,
    ftol=1e-2,
    flats=None,
):
    """
    MSPM_asym - Decompose an asymmetric tensor using the Multi Subspace Power Method

    This function computes the CP decomposition of an asymmetric tensor
    using the Multi Subspace Power Method (MSPM). It is a specialized
    version of multiSPM tailored for asymmetric tensors.

    ** Usage **
        lambdas, factors, stat = MSPM_asym(T, ...)

    ** INPUT **
         T : torch.Tensor
             Asymmetric tensor of dimension d_1 x d_2 x ... x d_n
      rank : int, optional
             Tensor rank. If not provided, it is estimated from singular values.
   maxiter : int, optional
             Maximum number of iterations for the power method
             (default: 5000).
    ntries : int, optional
             Maximum number of attempts for the power method
             (default: 5).
   gradtol : float, optional
             Gradient tolerance for convergence (default: 1e-15).
   ranksel : float, callable, or str, optional
             Tolerance for rank selection using singular values.
             You may provide a function that takes the singular values
             as input and decides the rank, or set it to "return_sv".
             In the latter case, the function terminates early and returns
             the singular values of the first flattening.
             (default: 1e-4).
      ftol : float, optional
             Function value tolerance for restarting the power method
             (default: 1e-2).
     flats : optional
             Flattenings to be considered for the decomposition. This can be
             provided as a pair of index lists or as a boolean matrix where
             each row indicates a flattening. For example, for a tensor
             T(i,j,k,l), the flattening T(ij,kl) can be indicated by
             [True, True, False, False] or ([0,1], [2,3]).
             (default: None, which automatically selects the flattenings
             that allow the highest rank).

    ** OUTPUT **
    lambdas : torch.Tensor
              Scaling factors for the decomposition.
    factors : list of torch.Tensor
              Factor matrices of the decomposition.
       stat : dict
              Dictionary containing various statistics of the decomposition process:
                - extracttime : time spent on extracting tensor properties
                - powertime   : time spent on the power method
                - deflatetime : time spent on deflation
                - avgiter     : average number of iterations per rank
                - nrr         : number of restarts during the power method
                - totaltime   : total runtime of the decomposition

    ** NOTE **
    - This implementation requires the 'utils.py' file to be available in the
      Python path.
    - The method follows the MSPM framework using tensor flattenings,
      recursive power iterations, and deflation.

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

    flats = _parse_flats(flats, order, device)

    if flats is None:
        flats, max_rank = find_best_flatpair(dims)
        flats = flats.to(device=device)
    else:
        flat_ranks = []
        for i in [1, 0]:
            mask = flats[i]
            left_dims = [dims[j] for j in range(order) if bool(mask[j].item())]
            right_dims = [dims[j] for j in range(order) if not bool(mask[j].item())]
            left_rank = math.prod(left_dims) - sum(left_dims)
            right_rank = math.prod(right_dims)
            flat_ranks.append(min(left_rank, right_rank))
        max_rank = max(flat_ranks)
        assert max_rank > 1

    assert rank is None or rank <= max_rank

    mask1 = flats[0]
    mask2 = flats[1]

    dim_order = (
        [i for i in range(order) if bool((mask1[i] & mask2[i]).item())]
        + [i for i in range(order) if bool((mask1[i] & (~mask2[i])).item())]
        + [i for i in range(order) if bool(((~mask1[i]) & mask2[i]).item())]
        + [i for i in range(order) if bool(((~mask1[i]) & (~mask2[i])).item())]
    )

    T = T.permute(dim_order)
    dims = [dims[i] for i in dim_order]

    mask1 = mask1[dim_order]
    mask2 = mask2[dim_order]

    m1 = int(torch.sum(mask1).item())
    m2 = int(torch.sum(mask2).item())
    mc = int(torch.sum(mask1 & mask2).item())
    mu = int(torch.sum(mask1 | mask2).item())

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
    C = torch.diag(1.0 / S)
    C_copy = C.clone()
    V1 = V1[:, :rank]

    perm2 = list(range(0, mc)) + list(range(m1, mu)) + list(range(mc, m1)) + list(range(mu, order))
    T2 = T.permute(perm2)

    dims_t = dims
    left_shape_2 = [dims_t[i] for i in list(range(0, mc)) + list(range(m1, mu))]
    matT_2 = reshapeF(T2, math.prod(left_shape_2), -1)
    U2, S2, Vh2 = torch.linalg.svd(matT_2, full_matrices=False)
    V2 = Vh2.transpose(0, 1)
    U2 = U2[:, :rank]
    V2 = V2[:, :rank]

    factors = [torch.zeros((dim, rank), dtype=dtype, device=device) for dim in dims]
    lambdas = torch.zeros(rank, dtype=dtype, device=device)

    lap = time.perf_counter() - timer
    stat = {
        "extracttime": lap,
        "powertime": 0.0,
        "deflatetime": 0.0,
        "avgiter": 0.0,
        "nrr": 0,
    }

    Aks = [torch.zeros(dim, dtype=dtype, device=device) for dim in dims]

    def power_method(U, Aks_local, r_local):
        mk = len(Aks_local)
        if mk == 0:
            return Aks_local, 1.0
        best_f = float('-inf')
        best_Aks = None

        for tries in range(1, ntries + 1):
            for j in range(mk):
                x = torch.randn_like(Aks_local[j])
                Aks_local[j] = x / torch.linalg.norm(x)

            for iter_idx in range(1, maxiter + 1):
                Aks_local, f, max_shift = power_method_iteration(U, Aks_local, r_local)
                if max_shift < gradtol:
                    break

            stat["avgiter"] += iter_idx

            if 1 - f < ftol:
                return Aks_local, f
            if tries == 1 or f > best_f:
                stat["nrr"] += 1
                best_f = f
                best_Aks = [a.clone() for a in Aks_local]
            else:
                stat["nrr"] += 1
                Aks_local = [a.clone() for a in best_Aks]

        return Aks_local, best_f

    for k in range(rank - 1, -1, -1):
        t0 = time.perf_counter()
        Aks[:m1], f = power_method(U1, Aks[:m1], k + 1)
        stat["powertime"] += time.perf_counter() - t0

        if mc > 0 and m2 > mc:
            if mc > 0:
                Akpow = tensor_product(Aks[:mc])
                U_mat = reshapeF(U2, len(Akpow), -1)
                U_half = (Akpow.reshape(1, -1) @ U_mat).squeeze(0)
                U_half = reshapeF(U_half, -1, rank)
                M = U_half @ U_half.transpose(0, 1)
                evals, evecs = torch.linalg.eigh(M)
                v = evecs[:, -1]
                if len(Aks[m1:mu]) > 0:
                    Aks[m1:mu], _ = power_method(v, Aks[m1:mu], 1)

            if order > mu:
                Akpow = tensor_product(Aks[m1:mu])
                U_mat = reshapeF(V1, len(Akpow), -1)
                U_half = (Akpow.reshape(1, -1) @ U_mat).squeeze(0)
                U_half = reshapeF(U_half, -1, rank)
                M = U_half @ U_half.transpose(0, 1)
                evals, evecs = torch.linalg.eigh(M)
                v = evecs[:, -1]
                if len(Aks[mu:order]) > 0:
                    Aks[mu:order], _ = power_method(v, Aks[mu:order], 1)
        else:
            Akpow = tensor_product(Aks[mc:m1])
            U_mat = reshapeF(V2, len(Akpow), -1)
            U_half = (Akpow.reshape(1, -1) @ U_mat).squeeze(0)
            U_half = reshapeF(U_half, -1, rank)
            M = U_half @ U_half.transpose(0, 1)
            evals, evecs = torch.linalg.eigh(M)
            v = evecs[:, -1]
            if len(Aks[mu:order]) > 0:
                Aks[mu:order], _ = power_method(v, Aks[mu:order], 1)

            if m2 > mc:
                if len(Aks[mu:order]) > 0:
                    Akpow = tensor_product(Aks[mu:order])
                    U_reshaped = reshapeF(V1, -1, len(Akpow), rank)
                    U_half = torch.einsum("ijr,j->ir", U_reshaped, Akpow)
                    M = U_half @ U_half.transpose(0, 1)
                    evals, evecs = torch.linalg.eigh(M)
                    v = evecs[:, -1]
                if len(Aks[m1:mu]) > 0:
                    Aks[m1:mu], _ = power_method(v, Aks[m1:mu], 1)

        alpha = (tensor_product(Aks[:m1]).reshape(1, -1) @ U1_copy).transpose(0, 1).reshape(-1)
        beta = (tensor_product(Aks[m1:order]).reshape(1, -1) @ V1).transpose(0, 1).reshape(-1)

        Ctbeta = (beta.reshape(1, -1) @ C_copy).transpose(0, 1).reshape(-1)
        lambdas[k] = 1.0 / torch.dot(alpha, Ctbeta)

        for i in range(order):
            factors[i][:, k] = Aks[i]

        if k > 0:
            x = get_hh_reflector((beta.reshape(1, -1) @ C).transpose(0, 1).reshape(-1))
            if torch.linalg.norm(x) > 0:
                C = RHR(C, x)
                U1 = RHR(U1, x)

        timenow = time.perf_counter() - timer
        stat["deflatetime"] += timenow - lap
        lap = timenow

    stat["avgiter"] = stat["avgiter"] / rank
    stat["totaltime"] = time.perf_counter() - timer

    reordered_factors = [None] * order
    for i, orig_pos in enumerate(dim_order):
        reordered_factors[orig_pos] = factors[i]
    factors = reordered_factors

    return lambdas, factors, stat


def reshapeF(x, *shape):
    """
    reshapeF - Reshape a tensor using MATLAB-style ordering.

    This helper reshapes a tensor so that flattening and matricization
    follow the same column-major convention used in the MATLAB version.
    """
    if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
        shape = tuple(shape[0])
    return x.permute(*reversed(range(x.ndim))).contiguous().view(*reversed(shape)).permute(*reversed(range(len(shape))))


def rank_selector(S, rank_sel):
    """
    tensor_product - Form the tensor product of factor vectors.

    This helper computes the Kronecker product of the vectors in Aks
    and returns the corresponding vectorized rank-1 term.
    """
    if callable(rank_sel):
        return rank_sel(S)
    if isinstance(rank_sel, (int, float)):
        typical = torch.sum(S**2) / torch.sum(torch.abs(S))
        return int(torch.sum(torch.abs(S) > rank_sel * typical).item())
    if rank_sel == "return_sv":
        return None
    raise ValueError("Rank selector option not implemented yet")


def tensor_product(Aks):
    """
    tensor_product - Form the tensor product of factor vectors.

    This helper computes the Kronecker product of the vectors in Aks
    and returns the corresponding vectorized rank-1 term.
    """
    if len(Aks) == 0:
        return torch.tensor(1.0, dtype=torch.float64)

    Akpow = Aks[0].reshape(-1, 1)
    for k in range(1, len(Aks)):
        Akpow = torch.kron(Aks[k].reshape(-1, 1), Akpow)
    Akpow = Akpow.flatten()
    return Akpow


def pagemtimes_vector(U, v):
    """
    pagemtimes_vector - Multiply stacked matrices by a vector.

    This helper applies a pagewise matrix-vector product across the
    third dimension of U.
    """
    return torch.einsum("abr,b->ar", U, v)


def power_method_iteration(U, Aks, r):
    """
    power_method_iteration - Perform one iteration of the power method.

    This function updates the current factor vectors by recursively
    contracting the working subspace and normalizing the result.
    """
    m = len(Aks)
    if m == 0:
        return Aks, 1.0, 0.0

    if m == 1:
        U = U.reshape(-1, r)

        Ak_old = Aks[0]

        tmp = U.T @ Ak_old
        Ak_new = U @ tmp

        f = torch.dot(Ak_new, Ak_old)

        Ak_new = Ak_new / torch.linalg.norm(Ak_new)

        max_shift = torch.max(torch.abs(Ak_new - Ak_old))

        Aks[0] = Ak_new

        return Aks, f, max_shift

    m2 = math.ceil(m / 2)

    Akpow = tensor_product(Aks[:m2])
    U_mat = reshapeF(U, len(Akpow), -1)
    U_half = (Akpow.reshape(1, -1) @ U_mat).squeeze(0)
    U_half = reshapeF(U_half, -1, r)
    right_Aks, _, max_shift = power_method_iteration(U_half, Aks[m2:], r)
    Aks[m2:] = right_Aks

    Akpow = tensor_product(Aks[m2:])
    U_reshaped = reshapeF(U, -1, len(Akpow), r)
    U_half = torch.einsum("ijr,j->ir", U_reshaped, Akpow)
    left_Aks, f, max_shift_ = power_method_iteration(U_half, Aks[:m2], r)
    Aks[:m2] = left_Aks

    max_shift = torch.maximum(max_shift, max_shift_)
    return Aks, f, max_shift


def get_hh_reflector(y):
    """
    get_hh_reflector - Get vector for Householder reflection.

    The last column of the corresponding Householder reflection is a
    multiple of the input vector y.
    """
    y = y.clone()
    norm_y = torch.linalg.norm(y)
    if norm_y == 0:
        return y
    s = torch.sign(y[-1])
    if s == 0:
        s = torch.tensor(1.0, dtype=y.dtype, device=y.device)
    y[-1] = y[-1] + norm_y * s
    y = y / torch.sqrt(torch.abs(y[-1]) * norm_y)
    return y


def RHR(A, x):
    """
    RHR - Apply Householder reflection from the right.

    This helper applies the Householder reflection defined by x to A
    from the right.
    """
    return A[:, :-1] - (A @ x).reshape(-1, 1) @ x[:-1].reshape(1, -1)


def LHR(A, x):
    """
    LHR - Apply Householder reflection from the left.

    This helper applies the Householder reflection defined by x to A
    from the left.
    """
    return A[:-1, :] - x[:-1].reshape(-1, 1) @ (x.reshape(1, -1) @ A)


def _parse_flats(flats, order, device):
    """
    _parse_flats - Process user-provided flattenings.

    This helper converts the flattening input into the boolean format
    used internally by MSPM_asym.
    """
    if flats is None:
        return None

    if isinstance(flats, (list, tuple)) and len(flats) == 2 and not torch.is_tensor(flats):
        out = torch.zeros((2, order), dtype=torch.bool, device=device)

        idx0 = flats[0]
        idx1 = flats[1]

        if len(idx0) > 0 and isinstance(idx0[0], bool):
            out[0] = torch.tensor(idx0, dtype=torch.bool, device=device)
        else:
            out[0, torch.tensor(idx0, dtype=torch.long, device=device)] = True

        if len(idx1) > 0 and isinstance(idx1[0], bool):
            out[1] = torch.tensor(idx1, dtype=torch.bool, device=device)
        else:
            out[1, torch.tensor(idx1, dtype=torch.long, device=device)] = True

        return out
    

if __name__ == "__main__":

    def cp_to_tensor(lambdas, factors):
        r = len(lambdas)
        shape = [U.shape[0] for U in factors]
        T = torch.zeros(*shape, dtype=factors[0].dtype, device=factors[0].device)

        for k in range(r):
            outer = lambdas[k] * factors[0][:, k]
            for mode in range(1, len(factors)):
                outer = torch.einsum("... , j -> ...j", outer, factors[mode][:, k])
            T = T + outer

        return T

    print("\n## Asymmetric MSPM tests ##")

    d1, d2, d3 = 40, 35, 30
    r = 20

    A = torch.randn(d1, r, dtype=torch.float64)
    B = torch.randn(d2, r, dtype=torch.float64)
    C = torch.randn(d3, r, dtype=torch.float64)

    lambdas_true = torch.ones(r, dtype=torch.float64)

    T = cp_to_tensor(lambdas_true, [A, B, C])


    print("\n## Example 1 - Noiseless (default flattenings) ##")

    start = time.perf_counter()
    lambdas_hat, factors_hat, stats = MSPM_asym(
        T,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-6,
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat)

    T_flat = reshapeF(T, -1)
    T_hat_flat = reshapeF(T_hat, -1)

    err = norm(T_flat - T_hat_flat) / norm(T_flat)
    print("Relative reconstruction error:", err)
    
    print("\n## Example 1 - Noiseless (forced flattenings) ##")

    start = time.perf_counter()
    lambdas_hat2, factors_hat2, stats = MSPM_asym(
        T,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-6,
        flats=([0,1],[0,2]),
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat2 = cp_to_tensor(lambdas_hat2, factors_hat2)

    err2 = norm(T.reshape(-1) - T_hat2.reshape(-1)) / norm(T.reshape(-1))
    print("Relative reconstruction error (forced flats):", err2)


    print("\n## Example 2 - Noisy ##")

    T_noisy = T + 0.1 * torch.randn(d1, d2, d3, dtype=torch.float64)

    start = time.perf_counter()
    lambdas_hat, factors_hat, stats = MSPM_asym(
        T_noisy,
        rank=r,
        maxiter=5000,
        ntries=5,
        gradtol=1e-12,
        ftol=1e-2,
    )
    print("Time taken:", time.perf_counter() - start)

    T_hat = cp_to_tensor(lambdas_hat, factors_hat)

    print(
        "Relative reconstruction error vs clean:",
        norm(T.reshape(-1) - T_hat.reshape(-1)) / norm(T.reshape(-1)),
    )

    print(
        "Relative reconstruction error vs noisy:",
        norm(T_noisy.reshape(-1) - T_hat.reshape(-1)) / norm(T_noisy.reshape(-1)),
    )