import torch
from time import time
from utils_torch import generate_lowrank_tensor, norm, khatri_rao_power, option_parser, pos, isbool, apply_Q_from_QR, symmetric_indices, dot, eig2

torch.set_default_dtype(torch.float64)

def power_method_iteration(d, n2, Vt, ntries, maxiter, gradtol, ftol):
    """
    power_method_iteration - Estimate a rank-1 factor via power iterations.

    Runs a shifted power method with random restarts to extract a dominant
    component from a projected tensor subspace.
    """
    if n2 <= 4:
        cn = torch.sqrt(torch.tensor(2 * (n2 - 1) / n2))
    else:
        cn = (2 - torch.sqrt(torch.tensor(2.0))) * torch.sqrt(torch.tensor(float(n2)))

    Ak_best = None
    f_best = None

    for tries in range(ntries):
        Ak = torch.randn(d, dtype=Vt.dtype, device=Vt.device)
        Ak = Ak / (norm(Ak) + 1e-12)

        f = torch.tensor(0.0, dtype=Vt.dtype, device=Vt.device)

        for _ in range(maxiter):
            Apow = khatri_rao_power(Ak.reshape((-1, 1)), n2 - 1).reshape(-1)
            VAk = (Vt @ Apow).reshape((-1, d))
            Ak_new = (VAk @ Ak) @ VAk

            f = dot(Ak_new, Ak)

            fcl = torch.clamp(f, min=0.5, max=1.0)
            clambda = torch.sqrt(fcl * (1 - fcl))
            shift = cn * clambda
            Ak_new = Ak_new + shift * Ak

            nrm = norm(Ak_new)
            if nrm < 1e-12:
                Ak_new = torch.randn_like(Ak_new)
                Ak_new = Ak_new / (norm(Ak_new) + 1e-12)
            else:
                Ak_new = Ak_new / nrm

            err = norm(Ak - Ak_new)
            Ak = Ak_new
            if err < gradtol:
                break

        if (Ak_best is None) or (f > f_best):
            Ak_best = Ak
            f_best = f

        if 1 - f < ftol:
            break

    return Ak_best

def subspace_power_method(T, d=None, n=None, r=None, **kwargs):
    """
    subspace_power_method - Decompose a symmetric tensor via SPM.

    Computes a CP decomposition of an even-order symmetric tensor using
    the Subspace Power Method (SPM). The method constructs a symmetric
    flattening, extracts a subspace via eigendecomposition, and recovers
    components through iterative power updates and deflation.
    """
    T = T if isinstance(T, torch.Tensor) else torch.tensor(T, dtype=torch.float64)

    if d is None:
        d = T.shape[0]
    if n is None:
        n = int(torch.log(torch.tensor(T.numel(), dtype=torch.float64)) / torch.log(torch.tensor(float(d))))

    assert n % 2 == 0 and n > 0

    opts = option_parser(
        kwargs,
        ('maxiter', 5000, pos),
        ('ntries', 3, pos),
        ('gradtol', 1e-14, pos),
        ('eigtol', 1e-8, pos),
        ('ftol', 1e-2, pos),
        ('w_out', True, isbool),
    )

    n2 = n // 2
    dn2 = d ** n2

    matT = T.reshape(dn2, dn2)

    symind, findsym, symindscale = symmetric_indices(d, n2)
    findsym = findsym.flatten()

    s = symindscale.to(dtype=T.dtype, device=T.device)
    sym_matT = s.reshape(1, -1) * matT[symind][:, symind] * s.reshape(-1, 1)

    D, symV = eig2(sym_matT)

    if r is None:
        r = int((D > opts['eigtol']).sum().item())

    D = D[:r]
    D_safe = torch.clamp(D, min=1e-12)

    V = (symV[:, :r] / s.reshape(-1, 1))[findsym, :]

    D1 = torch.diag(1.0 / D_safe)

    A = torch.empty((d, r), dtype=T.dtype, device=T.device)
    w = torch.empty((r,), dtype=T.dtype, device=T.device)

    for k in range(r):
        Ak = power_method_iteration(
            d,
            n2,
            V.T.reshape(-1, d ** (n2 - 1)),
            opts['ntries'],
            opts['maxiter'],
            opts['gradtol'],
            opts['ftol'],
        )

        Apow = khatri_rao_power(Ak.reshape(-1, 1), n2)
        alpha = (Apow.T @ V).T
        D1alpha = D1 @ alpha

        A[:, k] = Ak
        w[k] = 1.0 / ((alpha.T @ D1alpha)[0, 0] + 1e-12)

        if k < r - 1:
            if norm(D1alpha) < 1e-12:
                continue
            D1alpha = D1alpha / (norm(D1alpha) + 1e-12)

            D1 = apply_Q_from_QR(D1alpha, D1, side='R', trans='T')
            D1 = apply_Q_from_QR(D1alpha, D1, side='L', trans='N')
            D1 = D1[1:, 1:]

            V = apply_Q_from_QR(D1alpha, V, side='R', trans='T')
            V = V[:, 1:]

    if opts['w_out']:
        return A, w
    return A * w.reshape(1, -1) ** (1.0 / n)

if __name__ == "__main__":
    torch.manual_seed(0)

    d = 20
    r = 120
    n = 4

    A = torch.randn(d, r)
    T = generate_lowrank_tensor(A, n=n)

    start = time()
    A_ = subspace_power_method(T, r=r, w_out=False)
    print(time() - start)

    T_ = torch.zeros_like(T)
    for i in range(r):
        out = A_[:, i]
        for _ in range(n - 1):
            out = out[..., None] * A_[:, i]
        T_ = T_ + out

    print(norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))

    d = 20
    r = 120
    n = 4

    A = torch.randn(d, r)
    T = generate_lowrank_tensor(A, n=n)

    start = time()
    A_ = subspace_power_method(T, r=r, w_out=False)
    print(time() - start)

    T_ = torch.zeros_like(T)
    for i in range(r):
        out = A_[:, i]
        for _ in range(n - 1):
            out = out[..., None] * A_[:, i]
        T_ = T_ + out

    print(norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))

    d = 45
    r = 800
    n = 4

    A = torch.randn(d, r)
    T = generate_lowrank_tensor(A, n=n)

    start = time()
    A_ = subspace_power_method(T, r=r, w_out=False)
    print(time() - start)

    T_ = torch.zeros_like(T)
    for i in range(r):
        out = A_[:, i]
        for _ in range(n - 1):
            out = out[..., None] * A_[:, i]
        T_ = T_ + out

    print(norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))

    print("last example on MPS")

    device = torch.device("mps" if torch.backends.mps.is_available else "cpu")
    T = T.clone().detach()
    start = time()
    A_ = subspace_power_method(T, r=r, w_out=False)
    print(time() - start)
    T_mps = torch.zeros_like(T)
    for i in range(r):
        out = A_[:, i]
        for _ in range(n - 1):
            out = out[..., None] * A_[:, i]
        T_mps = T_mps + out

    print(norm(T.reshape(-1) - T_mps.reshape(-1)) / norm(T.reshape(-1)))



