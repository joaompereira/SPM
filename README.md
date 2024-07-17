# Subspace power method

An implementation of the subspace power method, in MATLAB and Python,
for decomposing symmetric tensors, as described in 

J. Kileel, J. M. Pereira,
[**Subspace power method for symmetric tensor
decomposition**](
https://arxiv.org/abs/1912.04007), arXiv:1912.04007

## Matlab

### Required External Packages

To compare performance with other tensor decomposition packages, 
and reproduce the results obtained in the paper, the users will have to install
and download external packages themselves. An exception to this is the implementation
of the FOOBI algorithm, from [**Fourth-Order Cumulant-Based Blind Identification
of Underdetermined Mixtures**](https://ieeexplore.ieee.org/document/4203062),
for which we did not find a MATLAB implementation,
and implemented ourselves. If you need assistance setting this up, you can
send an e-mail to [**jpereira@impa.br**](mailto:jpereira@impa.br).
A README for the installation of other packages may also be added in the future.  
The following are required packages:

- [**TensorLAB**](http://www.tensorlab.net/)

- [**Low rank Symmetric Tensor Approximations code**](https://mathweb.ucsd.edu/~njw/CODES/gpstd/symtensor_decm_aprx.html)

- [**GPCA-PDA**](http://www.vision.jhu.edu/gpca.htm) (GPCA only)

- [**GPCA-Voting**](http://people.eecs.berkeley.edu/~yang/software/softwarepage.html) (GPCA only)

### Installation

To install you just need to add the folders
`MATLAB\` and `MATLAB\helper_functions\` to the MATLAB path. Alternatively, you may run the file `MATLAB\setup.m`.

### Reproducibility

- To reproduce the results in [*Subspace power method for symmetric tensor
decomposition*](
https://arxiv.org/abs/1912.04007), run the file `MATLAB\tests\run_SPM_paper.m`.

## Python

### Required Packages

  - `numpy`
  - `scipy`

### Optional Package
  - `numba`: This package provides a just-in-time pre-compiler that can considerably speed up SPM performance. It is used only if the numba package is installed.            
  
### Installation

To install you just need to copy the python files, and from SPM.py, import the method `subspace_power_method`.
Additional methods (such as `generate_lowrank_tensor`) are also available which can be useful for testing SPM.
