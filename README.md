# Subspace power method

An implementation of the subspace power method, in MATLAB and Python,
for decomposing symmetric tensors, as described in 

J. Kileel, J. M. Pereira,
[**Subspace power method for symmetric tensor
decomposition and generalized PCA**](
https://arxiv.org/abs/1912.04007), arXiv:1912.04007

<details> <summary>Matlab</font></summary>

### Installation

To install you just need to add
`current_path\` and `current_path\helper_functions\` to the MATLAB path

### Other packages

To compare performance with other tensor decomposition packages, 
and reproduce the results obtained in the paper, the users will have to install
and download external packages themselves. An exception to this is the implementation
of the FOOBI algorithm, from [**Fourth-Order Cumulant-Based Blind Identification
of Underdetermined Mixtures**](https://ieeexplore.ieee.org/document/4203062),
for which we did not find a MATLAB implementation,
and implemented ourselves. If you need assistance setting this up, you can
send an e-mail to [**jpereira@impa.br**](mailto:jpereira@impa.br).
A README for the installation of other packages may also be added in the future.

#### External tensor packages

- [**Tensor Toolbox**](http://www.tensortoolbox.org/)

- [**TensorLAB**](http://www.tensorlab.net/)

- [**GPCA-PDA**](http://www.vision.jhu.edu/gpca.htm)

- [**GPCA-Voting**](http://people.eecs.berkeley.edu/~yang/software/softwarepage.html)

</details>  
  
<details> <summary>Python</summary>

### Required Packages

  - `numpy`
  - `scipy`
  - `numba` (Optional, compiles with numba if available)
  
### Installation

To install you just need to copy the python files, and from SPM.py, import the method `subspace_power_method`.
Additional methods (such as `generate_lowrank_tensor`) are also available which can be useful for testing SPM.
  
</details>  
