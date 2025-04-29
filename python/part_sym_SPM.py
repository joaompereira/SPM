import numpy as np
from utils import norm, khatri_rao_power, option_parser, \
    compiler_decorator, prange, pos, isbool, dormqr, lapack
from scipy.linalg import svd
from time import time


@compiler_decorator
def power_method_iteration(Vt, ntries, maxiter, gradtol, ftol):

    f_ = 0
    r, m, n = Vt.shape

    for tries in range(ntries):
        # Initialize Ak and Bk
        Ak = np.random.randn(m)
        Ak /= norm(Ak)
        Bk = np.random.randn(n)
        Bk /= norm(Bk)
        V_B = Vt.reshape(r * m, n)
        VAk = np.empty((n, r))

        for iter in range(maxiter):
            for k in prange(r):
                VAk[:, k] = np.dot(Ak, Vt[k])
            Bk = np.dot(VAk, np.dot(Bk, VAk))
            Bk /= norm(Bk)

            VBk = np.dot(V_B, Bk).reshape(r, m)

            Ak_new = np.dot(np.dot(VBk, Ak), VBk)

            f = np.dot(Ak, Ak_new)

            Ak_new /= norm(Ak_new)
            err = norm(Ak - Ak_new)
            Ak = Ak_new

            if err < gradtol:
                # Algorithm converged
                break

        if 1 - f < ftol:
            break
        elif tries == 0 or f > f_:
            f_ = f
            Ak_ = Ak
            Bk_ = Bk
        elif tries == ntries - 1:
            f = f_
            Ak = Ak_
            Bk = Bk_

    return Ak, Bk, f


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
                         ('w_out', True, isbool))

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

    for k in range(r):

        Ak, Bk, f = power_method_iteration(V.T.reshape(r-k, m, n), opts.ntries,
                                           opts.maxiter, opts.gradtol, opts.ftol)

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

    return A, B


if __name__ == '__main__':
    # Example usage

    m = 200
    n = 200
    r = 100

    A = np.random.randn(m, r)
    B = np.random.randn(n, r)

    T = np.dot(khatri_rao_power(A, 2), B.T).reshape(m, m, n)
    start = time()
    A_, B_ = spm_21sym(T, r=r, maxiter=1000, ntries=3,
                       gradtol=1e-10, ftol=1e-5)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))

    A = np.random.randn(m, r)
    B = np.random.randn(n, r)

    T = np.dot(khatri_rao_power(A, 2), B.T).reshape(m, m, n)
    start = time()
    A_, B_ = spm_21sym(T, r=r, maxiter=1000, ntries=3,
                       gradtol=1e-10, ftol=1e-5)
    print("Time taken:", time() - start)
    T_ = np.dot(khatri_rao_power(A_, 2), B_.T).reshape(m, m, n)
    print("Error:", np.linalg.norm(T.reshape(-1) -
          T_.reshape(-1)) / np.linalg.norm(T.reshape(-1)))
