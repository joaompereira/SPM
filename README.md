# Subspace power method

Implementations of the subspace power method and multi SPM, in MATLAB and Python, as described in the papers:

- J. Kileel, J. M. Pereira,
[**Subspace power method for symmetric tensor
decomposition**](
https://doi.org/10.1007/s11075-025-02165-y), Numerical Algorithms, 2025. [[preprint]](https://arxiv.org/abs/1912.04007)
- K. Wang, J. M. Pereira, J. Kileel, A. Seigal, [**Multi-subspace power method for decomposing all tensors**](https://www.arxiv.org/abs/2510.18627), ArXiv preprint, arXiv:2510.18627, 2025.

## Matlab

### Installation

To install you just need to add the folders
`MATLAB\` and `MATLAB\helper_functions\` to the MATLAB path. Alternatively, you may run the file `MATLAB\setup.m`.

### Reproducing results of *Subspace power method for symmetric tensor decomposition*

To reproduce all the results in [*Subspace power method for symmetric tensor decomposition*](
https://doi.org/10.1007/s11075-025-02165-y), you must first

- download and install external packages;
- download the ICA dataset.

After these steps, run the file `MATLAB\tests_spm_paper\run_SPM_paper_tests.m`. This will generate the figures in the paper in the `results` folder. If you need assistance setting this up, feel free to open an issue or send an e-mail to [**jpereira@uga.edu**](mailto:jpereira@uga.edu).

### Reproducing results of *Multi-subspace power method for decomposing all tensors*

To reproduce all the results in [*Multi-subspace power method for decomposing all tensors*](https://www.arxiv.org/abs/2510.18627), you must first download and install Tensorlab (see below). After that, run the file `MATLAB\tests_mspm_paper\run_MSPM_paper_tests.m`. This will generate the figures in the paper in the `results` folder. If you need assistance setting this up, feel free to open an issue or send an e-mail to [**jpereira@uga.edu**](mailto:jpereira@uga.edu).

### Required External Packages

The following are required packages:

- [**TensorLAB**](http://www.tensorlab.net/) (SPM and MSPM papers)

- [**Low rank Symmetric Tensor Approximations code**](https://mathweb.ucsd.edu/~njw/CODES/gpstd/symtensor_decm_aprx.html) (SPM paper)

- [**GPCA-PDA**](http://www.vision.jhu.edu/gpca.htm) (GPCA only)

- [**GPCA-Voting**](http://people.eecs.berkeley.edu/~yang/software/softwarepage.html) (GPCA only)

### Datasets

To reproduce the ICA experiment in the SPM paper, you need to download the ICA dataset at

- [https://bnci-horizon-2020.eu/database/data-sets/013-2015/Subject01_s1.mat](https://bnci-horizon-2020.eu/database/data-sets/013-2015/Subject01_s1.mat)

## Python

### Required Packages

- `numpy`
- `scipy`

### Optional Package

- `numba`: This package provides a just-in-time pre-compiler that can considerably speed up SPM performance. It is used only if the numba package is installed.
  
### Installation

To install you just need to copy the python files, and from
`SPM.py`, import the method `subspace_power_method`.
Additional methods (such as `generate_lowrank_tensor`)
are also available which can be useful for testing SPM.
