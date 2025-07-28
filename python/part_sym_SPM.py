import numpy as np
from utils import norm, khatri_rao_power, option_parser, \
    compiler_decorator, prange, pos, isbool, dormqr, lapack
from scipy.linalg import svd
from time import time
from math import sqrt


@compiler_decorator
def power_method_iteration(Vt, maxiter, gradtol, Ak, Bk):

    r = Vt.shape[0]
    m = Ak.shape[0]
    n = Bk.shape[0]
    V_B = Vt.reshape(r * m, n)
    V_C = Vt.reshape(r, m * n)
    
    for iter in range(maxiter):
        
        VBk = np.dot(V_B, Bk).reshape(r, m)
        Ck = np.dot(VBk, Ak)
        Ak_new = np.dot(Ck, VBk)
        
        f = np.dot(Ak, Ak_new)
        
        Ak = Ak_new / norm(Ak_new)
        
        VCk = np.dot(Ck, V_C).reshape(m, n)
        Bk_new = np.dot(Ak, VCk)

        Bk_new /= norm(Bk_new)
        err = norm(Bk - Bk_new)
        Bk = Bk_new

        if err < gradtol:
            # Algorithm converged
            break
        
    return Ak, Bk, iter, err, f

@compiler_decorator
def pm_refinement_iteration(Vt, A_, rho, maxiter, gradtol, Ak, Bk):

    k = A_.shape[1] - 1
    r = Vt.shape[0]
    m = Ak.shape[0]
    n = Bk.shape[0]
    V_B = Vt.reshape(r * m, n)
    V_C = Vt.reshape(r, m * n)
    
    gamma = sqrt(1 - rho ** 2)
    
    for iter in range(maxiter):
        
        VBk = np.dot(V_B, Bk).reshape(r, m)
        Ck = np.dot(VBk, Ak)
        Ak_new = np.dot(Ck, VBk)
        
        f = np.dot(Ak, Ak_new)
        
        Ak = Ak_new / norm(Ak_new)
        
        if k > 0:
            corr = np.dot(Ak, A_)
            ucorr = np.abs(corr)
            ind = np.argmax(ucorr)
            if ucorr[ind] > rho:
                Ak -= corr[ind] * A_[:, ind]
                Ak /= norm(Ak)
                s = np.sign(corr[ind])
                Ak = gamma * Ak + (rho * s) * A_[:, ind] 

        VCk = np.dot(Ck, V_C).reshape(m, n)
        Bk_new = np.dot(Ak, VCk)
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
                         ('O_version', False, isbool))

    m_, m, n = T.shape
    assert m_ == m, "Tensor T must be symmetric."

    # Flatten T
    T = T.reshape(m, -1)

    # Perform SVD
    U, D, Vt = svd(T, full_matrices=False)

    # Determine tensor rank by the eigenvalues of mat(T)
    if r is None:
        r = D.shape[0] - np.searchsorted(D[::-1], opts.eigtol)

    D1 = np.diag(1.0 / D[:r])
    D1_deflated = D1.copy()
    Vt = np.ascontiguousarray(Vt[:r, :])
    Vt_deflated = Vt.copy()
    U = U[:, :r]
    O = np.eye(r)

    A = np.zeros((m, r))
    B = np.zeros((n, r))
    
    stats = []

    for k in range(r):
        
        statsk = []
        
        for tries in range(opts.ntries):
            # Initialize Ak and Bk
            Ak = np.random.randn(m)
            Ak /= norm(Ak)
            Bk = np.random.randn(n)
            Bk /= norm(Bk)

            Ak, Bk, iter, err, f = power_method_iteration(Vt_deflated, opts.maxiter, opts.gradtol, Ak, Bk)
            
            # Ak, Bk, iter, err, f = power_method_iteration(Vt, opts.maxiter, opts.gradtol, Ak, Bk)
            
            statsk.append(dict(niter=iter, err=err, f=f))

            if 1 - f < opts.ftol:
                break
            elif tries == 0 or f > f_:
                f_ = f
                Ak_ = Ak
                Bk_ = Bk
            elif tries == opts.ntries - 1:
                f = f_
                Ak = Ak_
                Bk = Bk_
                
        Ak, Bk, iter, err, f = pm_refinement_iteration(Vt, A[:, :k], opts.rho, opts.maxiter, opts.gradtol, Ak, Bk)
        
        statsk.sort(reverse=True, key=lambda stat: stat['f'])
        #stats.append(statsk)
        stats.append({'niter': iter, 'err': err, 'f': f, 'deflate_stats': statsk})

        if k < r-1:
            alphaU = np.dot(Ak, U)
            
            # Update V and D using Householder reflection
            # Calculate the new matrix D and the new subspace
            # Use Householder reflection to update V and D
            if opts.O_version:
                D1alphaU = np.dot(np.dot(alphaU, O), D1_deflated)
                D1alphaV = np.dot(D1_deflated, np.dot(Vt_deflated, (Ak.reshape(-1, 1) * Bk.reshape(1, -1)).reshape(-1)))
                
                qr, tau, work, info = lapack.dgeqrf(D1alphaU, overwrite_a=1)
                D1_deflated, work, info = dormqr('R', 'T', qr, tau, D1_deflated, overwrite_c=1)
                V_deflated, work, info = dormqr('R', 'T', qr, tau, Vt_deflated.T, overwrite_c=1)
                
                qr, tau, work, info = lapack.dgeqrf(D1alphaV, overwrite_a=1)
                D1_deflated, work, info = dormqr('L', 'N', qr, tau, D1_deflated, overwrite_c=1)
                O, work, info = dormqr('R', 'T', qr, tau, O, overwrite_c=1)

                Vt_deflated = V_deflated[:, 1:].T
                D1_deflated = D1_deflated[1:, 1:]
                O = O[:, 1:]
            else:
                D1alphaU = np.dot(alphaU, D1_deflated)
                
                qr, tau, work, info = lapack.dgeqrf(D1alphaU, overwrite_a=1)
                D1_deflated, work, info = dormqr('R', 'T', qr, tau, D1_deflated, overwrite_c=1)
                V_deflated, work, info = dormqr('R', 'T', qr, tau, Vt_deflated.T, overwrite_c=1)
                
                Vt_deflated = V_deflated[:, 1:].T
                D1_deflated = D1_deflated[:, 1:]
                
        A[:, k] = Ak
        
    G = np.square(np.dot(A.T, A))
    B = np.linalg.solve(G, np.dot(khatri_rao_power(A, 2).T, np.reshape(T, (-1, n)))).T
    
    if opts.return_stats:
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
                         ('return_stats', False, isbool))

    m_, m, n = T.shape
    assert m_ == m, "Tensor T must be symmetric."

    # Flatten T
    T = T.reshape(m, -1)

    # Perform SVD
    U, D, Vt = svd(T, full_matrices=False)

    # Determine tensor rank by the eigenvalues of mat(T)
    if r is None:
        r = D.shape[0] - np.searchsorted(D[::-1], opts.eigtol)

    D1 = np.diag(1.0 / D[:r])
    V = np.ascontiguousarray(Vt[:r, :]).T
    U = U[:, :r]

    A = np.zeros((m, r))
    B = np.zeros((n, r))
    
    stats = []

    for k in range(r):
        
        statsk = []
        
        for tries in range(opts.ntries):
            # Initialize Ak and Bk
            Ak = np.random.randn(m)
            Ak /= norm(Ak)
            Bk = np.random.randn(n)
            Bk /= norm(Bk)

            Ak, Bk, iter, err, f = power_method_iteration(V.T.reshape(r-k, m, n), opts.maxiter, opts.gradtol, Ak, Bk)

            statsk.append(dict(niter=iter, err=err, f=f))

            if 1 - f < opts.ftol:
                break
            elif tries == 0 or f > f_:
                f_ = f
                Ak_ = Ak
                Bk_ = Bk
            elif tries == opts.ntries - 1:
                f = f_
                Ak = Ak_
                Bk = Bk_

        statsk.sort(reverse=True, key=lambda stat: stat['f'])
        stats.append(statsk)

        alphaU = np.dot(Ak, U)
        alphaV = np.dot((Ak.reshape(-1, 1) * Bk.reshape(1, -1)).reshape(-1), V)

        # Solve for lambda
        D1alphaU = np.dot(alphaU, D1)
        D1alphaV = np.dot(D1, alphaV)
        lambdak = norm(alphaU) * norm(alphaV) / np.dot(alphaV, D1alphaU)

        if k < r-1:
            # Update V and D using Householder reflection
            # Calculate the new matrix D and the new subspace
            # Use Householder reflection to update V and D
            qr, tau, work, info = lapack.dgeqrf(D1alphaU, overwrite_a=1)
            D1, work, info = dormqr('R', 'T', qr, tau, D1, overwrite_c=1)
            V, work, info = dormqr('R', 'T', qr, tau, V, overwrite_c=1)

            V = V[:, 1:]

            qr, tau, work, info = lapack.dgeqrf(D1alphaV, overwrite_a=1)
            D1, work, info = dormqr('L', 'N', qr, tau, D1, overwrite_c=1)
            U, work, info = dormqr('R', 'T', qr, tau, U, overwrite_c=1)

            D1 = D1[1:, 1:]
            U = U[:, 1:]

        A[:, k] = Ak
        B[:, k] = lambdak * Bk
        
    G = np.square(np.dot(A.T, A))
    B = np.linalg.solve(G, np.dot(khatri_rao_power(A, 2).T, np.reshape(T, (-1, n)))).T

    if opts.return_stats:
        return A, B, stats
    else:
        return A, B

if __name__ == '__main__':
    # Example usage

    m = 200
    n = 200
    r = 100

    A = np.random.randn(m, r)
    B = np.random.randn(n, r)

    T = np.dot(khatri_rao_power(A, 2), B.T).reshape(m, m, n)
    
    print(":: SPM ::")
    start = time()
    A_, B_, stats = spm_21sym(T, r=r, maxiter=1000, ntries=3,
                       gradtol=1e-10, ftol=1e-5, return_stats=True)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))
    
    print(":: Robust SPM ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T, r=r, maxiter=1000, ntries=3,
                       gradtol=1e-10, ftol=1e-5, return_stats=True)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))
    
    print(":: Robust SPM (O-version) ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T, r=r, maxiter=1000, ntries=3,
                       gradtol=1e-10, ftol=1e-5, return_stats=True, O_version=True)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))

    A = np.random.randn(m, r)
    B = np.random.randn(n, r)

    T = np.dot(khatri_rao_power(A, 2), B.T).reshape(m, m, n)
    T_noisy = T + 2 * np.random.randn(m, m, n) # Add noise
    
    print(":: SPM ::")
    start = time()
    A_, B_, stats = spm_21sym(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))
    
    print(":: Robust SPM ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))
    
    print(":: Robust SPM (O-version) ::")
    start = time()
    A_, B_, stats = spm_21sym_robust(T_noisy, r=r, maxiter=5000, ntries=3,
                       gradtol=1e-10, ftol=1e-1, return_stats=True, O_version=True)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))


