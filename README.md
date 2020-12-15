# ewald-summation
Efficient and easy to use fortran implementation of the [Ewald summation](https://en.wikipedia.org/wiki/Ewald_summation) method to calculate electrostatic energies and forces of periodic systems.
The code supports Gaussian shaped atomic charge distributions as well as point charges and is parallelized using OpenMP. 

## Usage
This repository contains a small test program (`src/testEwaldSummation.f90`) that checks the correctness of the Ewald summation implementation using finite differences that can be used as a reference on how to use the code. 

The main subroutine for the Ewald summation is: `ewaldEnergy(nat, ats, lat, q, sigma, e, f [, stress])`.

It takes the following arguments:
* **nat**: Number of atoms in one unit cell.
* **ats**: Array of dimension (3, nat) containing the atomic coordinates.
* **lat**: Array of dimension (3, 3) containing the lattice vectors. lat(:,i) is the i-th lattice vector.
* **q**: Array of dimension (nat) containing the atomic charges. The sum of all charges must be zero.
* **sigma**: Array of dimension (nat) containing the standard deviations of the Gaussian shaped atomic charge distributions. If sigma(1) is < 0 point charges are used. 
* **e**: Will contain the total electrostatic energy.
* **f**: Array of dimension (3, nat) that will contain the atomic forces.
* **stress**: Optional argument of dimension (3, 3) that will contain the stress tensor. A subroutine (`dEdlatFromStress`) is included that allows to convert the stress tensor to the derivatives of the energy with respect to the lattice vectors. 

In the file `src/configuration/parameters.f90` a constant called **ewaldSummationPrecision** is defined. 
This number will be used to determine the Ewald splitting parameter (the size of the auxilary Gaussian charges used in the Ewald summation algorithm) as well as the real- and reciprocal-space cutoff in such a way that the total error of the energy calculated is roughly of that magnitude and optimal performance / scaling is achieved. 
This way the user is never bothered to manually set the splitting parameter or cutoff ranges.

## Compiling the code

To compile the example program you only need a fortran compiler and CMake. 
Using CMake the code can be compiled with the following commands.
Three CMake flags are provided, that allow to compile with or without OpenMP parallelization, with or without debugging flags and with the intel or gnu compiler.

```bash
mkdir build
cd build
cmake -DOPENMP=ON -DDEBUG=OFF -DINTEL=OFF .. # compile with OpenMP parallelization, without debugging flags and gfortran
make
```

## Reference
Most of this code was written during the development of our fourth-generation high-dimensional neural network (4GHDNNP).
I would therefore appreciate if you cite our paper in case you use this code in any academic work.

Ko, Tsz Wai, et al.  
"A fourth-generation high-dimensional neural network potential with accurate electrostatics including non-local charge transfer."   
arXiv preprint [arXiv:2009.06484](https://arxiv.org/abs/2009.06484) (2020).  

