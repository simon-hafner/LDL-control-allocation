## LDL decomposition for control allocation
This repository provides a library free C implementation of the LDLT decomposition modified to be used with control allocation algorithms. 
The implementation 
- does not require any math libraries,
- can handle rank-deficient matrices,
- is optimized for a low computational footprint for matrices used in control allocation.

The code specifically solves the (under-) determined equation system 

$$\boldsymbol{\nu} = \mathbf{B}\mathbf{u}$$

### Related Publication 
If you use this code in one of your publications, please cite the related paper

S. Hafner, S. Myschik, J. Bachler, F. Holzapfel, _Fast Pseudoinverse-Based Control Allocation using a Variant of the Cholesky Decomposition_, **Journal of Guidance, Control, and Dynamics** (2025). https://doi.org/10.2514/1.G008767
