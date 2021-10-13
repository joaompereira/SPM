# Subspace power method

An implementation of the subspace power method in MATLAB,
for decomposing symmetric even-order tensors, as described in 

J. Kileel, J. M. Pereira,
[**Subspace power method for symmetric tensor
decomposition and generalized PCA**](
https://arxiv.org/abs/1912.04007), arXiv:1912.04007

## Installation - Matlab

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
send an e-mail to [**joao.pereira@utexas.edu**](mailto:joao.pereira@utexas.edu).
A README for the installation of other packages may also be added in the future.

##### External tensor packages

- [**Tensor Toolbox**](http://www.tensortoolbox.org/)

- [**TensorLAB**](http://www.tensorlab.net/)

- [**GPCA-PDA**](http://www.vision.jhu.edu/gpca.htm)

- [**GPCA-Voting**](http://people.eecs.berkeley.edu/~yang/software/softwarepage.html)
