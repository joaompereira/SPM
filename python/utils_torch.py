import math
import torch
import tensorly as tl
tl.set_backend('pytorch')
from tensorly.tenalg import khatri_rao

torch.set_default_dtype(torch.float64)

def sv_dimension(dims, mults):
    """
    sv_dimension - Compute symmetric vector space dimension.

    Returns the dimension associated with mode sizes and multiplicities.
    """
    if isinstance(dims, torch.Tensor):
        dims = dims.tolist()
    if isinstance(mults, torch.Tensor):
        mults = mults.tolist()

    out = 1
    for d, m in zip(dims, mults):
        m = int(m)
        if m == 0:
            continue
        num = 1
        den = 1
        for k in range(1, m + 1):
            num *= (d + k - 1)
            den *= k
        out *= num // den
    return int(out)


def find_biggest_flattenings_ps(udims, nsyms, m=3, device=None):
    """
    find_biggest_flattenings_ps - Find top flattenings for partially symmetric tensors.

    Enumerates candidate flattenings and returns those with largest rank bounds.
    """
    if device is None:
        device = torch.device("cpu")

    udims = list(udims)
    nsyms = list(nsyms)
    nudims = len(udims)

    all_rows = []
    all_ranks = []

    def rec(idx, current):
        if idx == nudims:
            if all(v == 0 for v in current) or all(current[i] == nsyms[i] for i in range(nudims)):
                return

            left_rank = sv_dimension(udims, current) - sum(udims[i] for i in range(nudims) if current[i] > 0)
            right_mults = [nsyms[i] - current[i] for i in range(nudims)]
            right_rank = sv_dimension(udims, right_mults)
            flat_rank = min(left_rank, right_rank)

            all_rows.append(current.copy())
            all_ranks.append(flat_rank)
            return

        for v in range(nsyms[idx] + 1):
            current.append(v)
            rec(idx + 1, current)
            current.pop()

    rec(0, [])

    if len(all_rows) == 0:
        raise ValueError("No valid flattenings found")

    ranks = torch.tensor(all_ranks, dtype=torch.int64, device=device)
    uflats = torch.tensor(all_rows, dtype=torch.long, device=device)

    idx = torch.argsort(ranks, descending=True)
    idx = idx[: min(m, len(idx))]
    return uflats[idx], ranks[idx]

def find_best_flatpair(dims, n_flats=3):
    """
    find_best_flatpair - Select a good pair of tensor flattenings.

    Chooses two non-complementary flattenings for asymmetric MSPM.
    """
    left_flats, ranks = find_biggest_flattenings(dims, n_flats)

    for j in range(1, n_flats):
        for i in range(j):
            if not torch.all(torch.logical_xor(left_flats[i], left_flats[j])):
                left_flats = left_flats[[i, j], :]
                max_rank = ranks[j]
                return left_flats, max_rank

    # If no good pair found, return the top two
    left_flats = left_flats[:2]
    max_rank = ranks[1]
    return left_flats, max_rank

def search_sorted(v, val):
    """
    search_sorted - Find insertion index in a sorted list.

    Returns the position where val should be inserted.
    """
    n = len(v)
    low = 0
    high = n

    while high > low:
        m = (low + high) // 2
        if val < v[m]:
            low = m + 1
        else:
            high = m

    return low

def find_biggest_flattenings(dims, m=1):
    """
    find_biggest_flattenings - Find top flattenings for an asymmetric tensor.

    Returns the flattenings with the largest admissible rank bounds.
    """
    if isinstance(dims, torch.Tensor):
        dims_list = dims.tolist()
    else:
        dims_list = list(dims)

    log_dims = [math.log(d) for d in dims_list]

    n_dims = len(dims_list)
    subset = [False] * (n_dims - 1) + [True]

    left_flats = []
    k = 0
    ranks = []

    keep_running = True
    while keep_running:
        if sum(subset) >= 1:
            left_dims = [dims_list[i] for i in range(n_dims) if subset[i]]
            right_dims = [dims_list[i] for i in range(n_dims) if not subset[i]]
            left_rank = math.prod(left_dims) - sum(left_dims)
            right_rank = math.prod(right_dims)
            flat_rank = min(left_rank, right_rank)

            if k < m or flat_rank > ranks[m - 1]:
                if k < m:
                    k += 1
                i = search_sorted(ranks, flat_rank)
                ranks = ranks[:i] + [flat_rank] + ranks[i:k - 1]
                left_flats = left_flats[:i] + [subset.copy()] + left_flats[i:k - 1]

        for i in range(n_dims - 1, -1, -1):
            if not subset[i]:
                subset[i] = True
                for j in range(i + 1, n_dims):
                    subset[j] = False
                break
            elif i == 0:
                keep_running = False

    ranks = torch.tensor(ranks, dtype=torch.int64)
    left_flats = torch.tensor(left_flats, dtype=torch.bool)

    return left_flats, ranks

def norm(x):
    """
    norm - Compute Euclidean norm.

    Returns the norm of a tensor or array-like input.
    """
    if isinstance(x, torch.Tensor):
        return torch.linalg.norm(x)
    return torch.linalg.norm(torch.tensor(x, dtype=torch.float64))

def dot(x, y):
    """
    dot - Compute inner product.

    Returns the dot product of two flattened tensors.
    """
    return torch.dot(x.reshape(-1), y.reshape(-1))

def apply_Q_from_QR(A_factor, C, side: str, trans: str):
    """
    apply_Q_from_QR - Apply a Householder orthogonal factor.

    Applies the Q factor defined by A_factor to C from the left or right.
    """
    x = A_factor
    if not isinstance(x, torch.Tensor):
        x = torch.tensor(x, dtype=torch.float64)
    if not isinstance(C, torch.Tensor):
        C = torch.tensor(C, dtype=torch.float64)
    x = x.reshape(-1).to(dtype=C.dtype, device=C.device)
    m = x.numel()
    if m == 0:
        return C
    normx = torch.linalg.norm(x)
    if normx == 0:
        return C
    e1 = torch.zeros_like(x)
    e1[0] = 1.0
    s = 1.0 if x[0] >= 0 else -1.0
    v = x + s * normx * e1
    vnorm = torch.linalg.norm(v)
    if vnorm == 0:
        return C
    v = v / vnorm
    H = torch.eye(m, dtype=C.dtype, device=C.device) - 2.0 * torch.outer(v, v)
    if side == 'L':
        return (H.T @ C) if (trans == 'T') else (H @ C)
    return (C @ H.T) if (trans == 'T') else (C @ H)

def khatri_rao_power(A, n):
    """
    khatri_rao_power - Repeat Khatri-Rao product of a matrix.

    Forms the n-fold Khatri-Rao product of A with itself.
    """
    A = A if isinstance(A, torch.Tensor) else torch.tensor(A, dtype=torch.float64)
    if A.ndim == 1:
        A = A.reshape(-1, 1)
    mats = [A] * n
    return khatri_rao(mats)

def khatri_rao_product(A, B):
    """
    khatri_rao_product - Compute Khatri-Rao product of two matrices.

    Returns the columnwise Kronecker product of A and B.
    """
    A = A if isinstance(A, torch.Tensor) else torch.tensor(A, dtype=torch.float64)
    B = B if isinstance(B, torch.Tensor) else torch.tensor(B, dtype=torch.float64)
    return khatri_rao([A, B])

def symmetric_indices(d, n2):
    """
    symmetric_indices - Build index data for symmetric flattenings.

    Returns representative indices, inverse mapping, and scaling factors.
    """
    grids = torch.meshgrid(*[torch.arange(d) for _ in range(n2)], indexing='ij')
    full = torch.stack(grids, dim=-1).reshape(-1, n2)
    sorted_idx = torch.sort(full, dim=1).values
    uniques, inverse, counts = torch.unique(sorted_idx, dim=0, return_inverse=True, return_counts=True)
    rep = uniques
    findsym = inverse
    symindscale = torch.sqrt(counts.to(dtype=torch.float64))
    powers = (d ** torch.arange(n2 - 1, -1, -1)).to(dtype=torch.long)
    symind_flat = (rep.to(dtype=torch.long) * powers).sum(dim=1)
    return symind_flat, findsym, symindscale

def generate_lowrank_tensor(A, n, w=None):
    """
    generate_lowrank_tensor - Generate a symmetric low-rank tensor.

    Builds a tensor from rank-1 symmetric outer products of the columns of A.
    """
    A = A if isinstance(A, torch.Tensor) else torch.tensor(A, dtype=torch.float64)
    d, r = A.shape
    if w is None:
        w = torch.ones(r, dtype=A.dtype, device=A.device)
    else:
        w = w if isinstance(w, torch.Tensor) else torch.tensor(w, dtype=A.dtype, device=A.device)
        w = w.reshape(-1)
    T = torch.zeros([d] * n, dtype=A.dtype, device=A.device)
    for i in range(r):
        out = A[:, i]
        for _ in range(n - 1):
            out = out[..., None] * A[:, i]
        T = T + w[i] * out
    return T

def tensor_from_columns(A, n):
    """
    tensor_from_columns - Form a symmetric tensor from matrix columns.

    Sums symmetric rank-1 outer products of the columns of A.
    """
    A = A if isinstance(A, torch.Tensor) else torch.tensor(A, dtype=torch.float64)
    d, r = A.shape
    T = torch.zeros([d] * n, dtype=A.dtype, device=A.device)
    for i in range(r):
        out = A[:, i]
        for _ in range(n - 1):
            out = out[..., None] * A[:, i]
        T = T + out
    return T

def eig2(a):
    """
    eig2 - Compute sorted eigendecomposition.

    Returns eigenvalues and eigenvectors in descending order.
    """
    D, V = torch.linalg.eigh(a)
    idx = torch.argsort(D, descending=True)
    return D[idx], V[:, idx]

def pos(x):
    return x > 0

def isbool(x):
    return isinstance(x, bool)

def option_parser(kwargs, *opts):
    """
    option_parser - Parse keyword options.

    Returns a dictionary of validated option values with defaults.
    """
    out = {}
    for name, default, check in opts:
        if name in kwargs:
            if check is not None and not check(kwargs[name]):
                raise ValueError(f"Invalid value for {name}")
            out[name] = kwargs[name]
        else:
            out[name] = default
    return out

