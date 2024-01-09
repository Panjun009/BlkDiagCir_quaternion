# BlkDiagCir_quaternion
Block Diagonalization of Quaternion Circulant Matrices

Project: BlkDiagCir_quaternion

If you use the code, please cite: 

Reference: J. Pan and M. Ng, "Block Diagonalization of Quaternion Circulant Matrices with Applications "  submitted (will update later)
https://arxiv.org/abs/2302.04086  


Description: a quaternion circulant matrix can be block-diagonalized into 1-by-1 block and 2-by-2 block matrices by permuted discrete quaternion Fourier transform matrix. 
With such a block-diagonalized form, the inverse of a quaternion circulant matrix can be determined efficiently similar to the inverse of a complex circulant matrix. 

The project's code consists of three main applications:

1. computing the inverse of a quaternion circulant matrix

2. solving quaternion Toeplitz system arising from linear prediction of quaternion signals  

3. quaternion tensor singular value decomposition. 
 

Kind reminder: our code is written based on two MATLAB toolboxes: 


MATLAB toolbox QTFM for quaternions: https://qtfm.sourceforge.io/

TensorLab for tensors:      https://www.tensorlab.net/
