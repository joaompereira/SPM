import torch
torch.set_default_dtype(torch.float64)
torch.manual_seed(0)
from utils_torch import norm, khatri_rao_power, khatri_rao_product, apply_Q_from_QR, option_parser, pos, isbool
from time import time
from math import sqrt


# @compiler_decorator
def power_method_iteration(Vt, maxiter, gradtol, Ak, Bk):

    r = Vt.shape[0]
    m = Ak.shape[0]
    n = Bk.shape[0]
    V_B = Vt.reshape(r * m, n)
    V_C = Vt.reshape(r, m * n)
    
    for iter in range(maxiter):
        
        VBk = (V_B @ Bk).reshape(r, m)
        Ck = VBk @ Ak
        Ak_new = Ck @ VBk
        
        f = Ak @ Ak_new
        
        Ak = Ak_new / norm(Ak_new)
        
        VCk = (Ck @ V_C).reshape(m, n)
        Bk_new = Ak @ VCk

        Bk_new /= norm(Bk_new)
        err = norm(Bk - Bk_new)
        Bk = Bk_new

        if err < gradtol:
            # Algorithm converged
            break
        
    return Ak, Bk, iter, err, f

# @compiler_decorator
def pm_refinement_iteration(Vt, A_, rho, maxiter, gradtol, Ak, Bk):

    k = A_.shape[1] - 1
    r = Vt.shape[0]
    m = Ak.shape[0]
    n = Bk.shape[0]
    V_B = Vt.reshape(r * m, n)
    V_C = Vt.reshape(r, m * n)
    
    gamma = sqrt(1 - rho ** 2)
    
    for iter in range(maxiter):
        
        VBk = (V_B @ Bk).reshape(r, m)
        Ck = VBk @ Ak
        Ak_new = Ck @ VBk
        
        f = Ak @ Ak_new

        Ak = Ak_new / norm(Ak_new)
        
        if k > 0:
            corr = Ak @ A_
            ucorr = torch.abs(corr)
            ind = int(ucorr.argmax().item())
            if float(ucorr[ind]) > rho:
                Ak = Ak - corr[ind] * A_[:, ind]
                Ak = Ak * (gamma / norm(Ak))
                Ak = Ak + (rho * corr[ind].sign()) * A_[:, ind]

        VCk = (Ck @ V_C).reshape(m, n)
        Bk_new = Ak @ VCk
        Bk_new /= norm(Bk_new)
        
        # corr = np.dot(Bk_new, B_)
        # ucorr = np.abs(corr)
        # ind = np.argmax(ucorr)
        # if ucorr[ind] > rho:
        #     Bk_new -= corr * B_[:, ind]
        #     Bk_new /= norm(Bk_new)
        #     s = np.sign(corr[ind])
        #     Bk_new = np.sqrt(1 - rho ** 2) * Bk_new + (rho * s) * B_[:, ind]
        
        err = norm(Bk - Bk_new)
        Bk = Bk_new

        if err < gradtol:
            # Algorithm converged
            break
        
    return Ak, Bk, iter, err, f

def spm_21sym_robust(T, r=None, **kwargs):
    """
    Decompose symmetric even order tensor using subspace power method.

    Parameters:
        T (ndarray): Tensor of dimension L^n.
        R (int, optional): Tensor rank. If not provided, it will be estimated.
        kwargs: Various SPM options as key-value pairs.
            maxiter (int): Maximum number of iterations of power method (default: 5000).
            ntries (int): Maximum number of tries for initialization (default: 5).
            gradtol (float): Gradient tolerance (default: 1e-15).
            ranksel (float): Tolerance for selecting the rank of T (default: 1e-4).
            ftol (float): Function value tolerance for restarting (default: 1e-2).

    Returns:
        A (ndarray): L x R matrix where the columns are the rank decomposition of T.
        B (ndarray): Scaling factors.
        stat (dict): Various statistics of SPM.
    """

    opts = option_parser(kwargs,
                         ('maxiter', 5000, pos),
                         ('ntries', 3, pos),
                         ('gradtol', 1e-14, pos),
                         ('eigtol', 1e-8, pos),
                         ('ftol', 1e-2, pos),
                         ('rho', .7, pos),
                         ('w_out', True, isbool),
                         ('return_stats', False, isbool),
                         ('O_version', False, isbool),
                         ('lsq_refine', None, None))

    T = T if isinstance(T, torch.Tensor) else torch.tensor(T, dtype=torch.float64)
    m_, m, n = T.shape
    assert m_ == m, "Tensor T must be symmetric."

    # Flatten T
    T = T.reshape(m, -1)

    # Perform SVD
    U, D, Vt = torch.linalg.svd(T, full_matrices=False)

    # Determine tensor rank by the eigenvalues of mat(T)
    if r is None:
        r = int((D > opts['eigtol']).sum().item())

    D1 = torch.diag(1.0 / D[:r])
    D1_deflated = D1.clone()
    Vt = Vt[:r, :]
    Vt_deflated = Vt.clone()
    U = U[:, :r]
    O = torch.eye(r, dtype=T.dtype, device=T.device)

    A = torch.zeros((m, r), dtype=T.dtype, device=T.device)
    B = torch.zeros((n, r), dtype=T.dtype, device=T.device)
    l = torch.zeros((r,), dtype=T.dtype, device=T.device)
    
    stats = []

    for k in range(r):
        
        statsk = []
        
        for tries in range(opts['ntries']):
            # Initialize Ak and Bk
            Ak = torch.randn(m, dtype=T.dtype, device=T.device)
            Ak /= norm(Ak)
            Bk = torch.randn(n, dtype=T.dtype, device=T.device)
            Bk /= norm(Bk)

            Ak, Bk, iter, err, f = power_method_iteration(
                Vt_deflated, opts['maxiter'], opts['gradtol'], Ak, Bk
            )
            
            statsk.append(dict(niter=iter, err=float(err), f=float(f)))

            if 1 - f < opts['ftol']:
                break
            elif tries == 0 or f > f_:
                f_ = f
                Ak_ = Ak
                Bk_ = Bk
            elif tries == opts['ntries'] - 1:
                f = f_
                Ak = Ak_
                Bk = Bk_
                
        Ak, Bk, iter, err, f = pm_refinement_iteration(
            Vt, A[:, :k], opts['rho'], opts['maxiter'], opts['gradtol'], Ak, Bk
        )
        
        statsk.sort(reverse=True, key=lambda stat: stat['f'])
        stats.append({'niter': iter, 'err': float(err), 'f': float(f), 'deflate_stats': statsk})
        
        alphaU = Ak @ U
        alphaV = Vt @ (Ak.reshape(-1, 1) * Bk.reshape(1, -1)).reshape(-1)

        # Solve for lambda
        D1alphaU = alphaU @ D1
        l[k] = norm(alphaU) * norm(alphaV) / (alphaV @ D1alphaU)

        if k < r-1:
                        
            # Update V and D using Householder reflection
            # Calculate the new matrix D and the new subspace
            # Use Householder reflection to update V and D
            if opts['O_version']:
                D1alphaU = (alphaU @ O) @ D1_deflated
                D1alphaV = D1_deflated @ (Vt_deflated @ (Ak.reshape(-1, 1) * Bk.reshape(1, -1)).reshape(-1))
                
                D1_deflated = apply_Q_from_QR(D1alphaU, D1_deflated, side='R', trans='T')
                V_deflated = apply_Q_from_QR(D1alphaU, Vt_deflated.T, side='R', trans='T')
                
                D1_deflated = apply_Q_from_QR(D1alphaV, D1_deflated, side='L', trans='N')
                O = apply_Q_from_QR(D1alphaV, O, side='R', trans='T')

                Vt_deflated = V_deflated[:, 1:].T
                D1_deflated = D1_deflated[1:, 1:]
                O = O[:, 1:]
            else:
                D1alphaU = alphaU @ D1_deflated
                
                D1_deflated = apply_Q_from_QR(D1alphaU, D1_deflated, side='R', trans='T')
                V_deflated = apply_Q_from_QR(D1alphaU, Vt_deflated.T, side='R', trans='T')
                
                Vt_deflated = V_deflated[:, 1:].T
                D1_deflated = D1_deflated[:, 1:]
                
        A[:, k] = Ak
        B[:, k] = Bk
    
    if opts['lsq_refine'] == "B":    
        G = (A.T @ A).pow(2)
        B = torch.linalg.solve(G, (khatri_rao_power(A, 2).T @ T.reshape(-1, n))).T
    elif opts['lsq_refine'] == "lambda":
        G = (A.T @ A).pow(2) * (B.T @ B)
        Z = T.reshape(m, -1) @ khatri_rao_product(A, B)
        u = torch.sum(A * Z, dim=0)
        B *= torch.linalg.solve(G, u).reshape(1, -1)
    else:
        B *= l.reshape(1, -1)
    
    if opts['return_stats']:
        return A, B, stats
    else:
        return A, B

    
def spm_21sym(T, r=None, **kwargs):
    """
    Decompose symmetric even order tensor using subspace power method.

    Parameters:
        T (ndarray): Tensor of dimension L^n.
        R (int, optional): Tensor rank. If not provided, it will be estimated.
        kwargs: Various SPM options as key-value pairs.
            maxiter (int): Maximum number of iterations of power method (default: 5000).
            ntries (int): Maximum number of tries for initialization (default: 5).
            gradtol (float): Gradient tolerance (default: 1e-15).
            ranksel (float): Tolerance for selecting the rank of T (default: 1e-4).
            ftol (float): Function value tolerance for restarting (default: 1e-2).

    Returns:
        A (ndarray): L x R matrix where the columns are the rank decomposition of T.
        B (ndarray): Scaling factors.
        stat (dict): Various statistics of SPM.
    """

    opts = option_parser(kwargs,
                         ('maxiter', 5000, pos),
                         ('ntries', 3, pos),
                         ('gradtol', 1e-14, pos),
                         ('eigtol', 1e-8, pos),
                         ('ftol', 1e-2, pos),
                         ('w_out', True, isbool),
                         ('return_stats', False, isbool),
                         ('lsq_refine', None, None))

    T = T if isinstance(T, torch.Tensor) else torch.tensor(T, dtype=torch.float64)
    m_, m, n = T.shape
    assert m_ == m, "Tensor T must be symmetric."

    dtype = T.dtype
    device = T.device

    # Flatten T
    T = T.reshape(m, -1)

    # Perform SVD
    U, D, Vt = torch.linalg.svd(T, full_matrices=False)

    # Determine tensor rank by the eigenvalues of mat(T)
    if r is None:
        r = int((D > opts['eigtol']).sum().item())

    D1 = torch.diag(1.0 / D[:r])
    V = Vt[:r, :].T
    U = U[:, :r]

    A = torch.zeros((m, r), dtype=dtype, device=device)
    B = torch.zeros((n, r), dtype=dtype, device=device)
    l = torch.zeros((r,), dtype=dtype, device=device)
    
    stats = []

    for k in range(r):
        
        statsk = []
        
        for tries in range(opts['ntries']):
            # Initialize Ak and Bk
            Ak = torch.randn(m, dtype=dtype, device=device)
            Ak /= norm(Ak)
            Bk = torch.randn(n, dtype=dtype, device=device)
            Bk /= norm(Bk)

            Ak, Bk, iter, err, f = power_method_iteration(
                V.T.reshape(r-k, m, n),
                opts['maxiter'], opts['gradtol'], Ak, Bk
            )

            statsk.append(dict(niter=iter, err=float(err), f=float(f)))

            if 1 - f < opts['ftol']:
                break
            elif tries == 0 or f > f_:
                f_ = f
                Ak_ = Ak
                Bk_ = Bk
            elif tries == opts['ntries'] - 1:
                f = f_
                Ak = Ak_
                Bk = Bk_

        statsk.sort(reverse=True, key=lambda stat: stat['f'])
        stats.append(statsk)

        alphaU = Ak @ U
        alphaV = (Ak.reshape(-1, 1) * Bk.reshape(1, -1)).reshape(-1) @ V

        # Solve for lambda
        D1alphaU = alphaU @ D1
        D1alphaV = D1 @ alphaV
        l[k] = norm(alphaU) * norm(alphaV) / (alphaV @ D1alphaU)

        if k < r-1:
            # Update V and D using Householder reflection
            # Calculate the new matrix D and the new subspace
            # Use Householder reflection to update V and D
            D1 = apply_Q_from_QR(D1alphaU, D1, side='R', trans='T')
            V  = apply_Q_from_QR(D1alphaU, V,  side='R', trans='T')

            V = V[:, 1:]

            D1 = apply_Q_from_QR(D1alphaV, D1, side='L', trans='N')
            U  = apply_Q_from_QR(D1alphaV, U,  side='R', trans='T')

            D1 = D1[1:, 1:]
            U  = U[:, 1:]

        A[:, k] = Ak
        B[:, k] = Bk
    
    if opts['lsq_refine'] == "B":    
        G = (A.T @ A).pow(2)
        B = torch.linalg.solve(G, (khatri_rao_power(A, 2).T @ T.reshape(-1, n))).T
    elif opts['lsq_refine'] == "lambda":
        G = (A.T @ A).pow(2) * (B.T @ B)
        Z = T.reshape(m, -1) @ khatri_rao_product(A, B)
        u = torch.sum(A * Z, dim=0)
        B *= torch.linalg.solve(G, u).reshape(1, -1)
    else:
        B *= l.reshape(1, -1)
    
    if opts['return_stats']:
        return A, B, stats
    else:
        return A, B


if __name__ == '__main__':
    # Example usage

    m = 100
    n = 75
    r = 50

    A = torch.randn(m, r, dtype=torch.float64)
    B = torch.randn(n, r, dtype=torch.float64)

    T = (khatri_rao_power(A, 2) @ B.T).reshape(m, m, n)
    
    print("## Example 1 - Noiseless ##")
    print(":: SPM ::")
    start = time()
    A_, B_, stats = spm_21sym(T, r=r, maxiter=1000, ntries=3,
                       gradtol=1e-10, ftol=1e-5, return_stats=True)
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))
    
    print(":: Robust SPM ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T, r=r, maxiter=1000, ntries=3,
                       gradtol=1e-10, ftol=1e-5, return_stats=True)
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))

    print("\n## Example 2 - Noisy ##")
    
    A = torch.randn(m, r, dtype=torch.float64)
    B = torch.randn(n, r, dtype=torch.float64)

    T = (khatri_rao_power(A, 2) @ B.T).reshape(m, m, n)
    T_noisy = T + 2 * torch.randn(m, m, n, dtype=torch.float64) # Add noise
    
    print(":: SPM ::")
    start = time()
    A_, B_, stats = spm_21sym(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True)
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))
    
    print(":: SPM (lambda lsq refine) ::")
    start = time()
    A_, B_, stats = spm_21sym(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True, lsq_refine='lambda')
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))
    
    print(":: SPM (B lsq refine) ::")
    start = time()
    A_, B_, stats = spm_21sym(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True, lsq_refine='B')
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))
    
    print(":: Robust SPM (no refinement) ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True)
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))
    
    print(":: Robust SPM (lambda lsq refine) ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True, lsq_refine='lambda')
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))
    
    print(":: Robust SPM (B refinement) ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True, lsq_refine='B')
    print("Time taken:", time() - start)
    T_ = (khatri_rao_power(A_, 2) @ B_.T).reshape(m, m, n)
    print("Error:", norm(T.reshape(-1) - T_.reshape(-1)) / norm(T.reshape(-1)))
