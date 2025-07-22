from warnings import warn

try:
    from numba import njit, prange
    compiler_decorator = njit(cache=True)
    
    NUMBA_COMPILER = True

except ModuleNotFoundError:

    warn("The numba package was not found.\nConsider installing it for improved performance.")

    def compiler_decorator(fun):
        return fun
    
    prange = range

    NUMBA_COMPILER = False

BLAS_DOT = False
if not NUMBA_COMPILER:
    try:
        from scipy.linalg.blas import ddot as dot
        from scipy.linalg.blas import dnrm2 as norm

        BLAS_DOT = True
    except ModuleNotFoundError:
        pass


if not BLAS_DOT:
    from numpy import dot
    from numpy.linalg import norm