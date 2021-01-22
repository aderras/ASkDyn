# Skyrmion Micromagnetics

Author: Amel Derras-Chouk

The code in this repository computes the evolution of magnetic moments in a two dimensional lattice by evaluating the Landau-Lifshitz-Gilbert equation at every lattice site. Code is written in Julia 1.0.4.

Required libraries:
 - HDF5 (Used to save files in .h5 format.)
 - FFTW (Fourier transform package used in computing dipole-dipole interaction.)
 - PaddedViews (Also used in computing dipole-dipole interaction.)
 - LinearAlgebra (Used for vector operations.)
 - Distributed (To parallelize.)

## Background

This code was written to study the behavior of skyrmions, topologically protected quasiparticles made up of atomic spins in magnetic materials. These spins obey the Landau-Lifshitz-Gilbert equation, which is implemented in this project.

## Usage

To run the program, execute ''julia runScript.jl'' in terminal. A series of prompts requesting input for material parameters, system size, and numerical options will follow.

## Data

When the script is executed, results are saved to the directory "/data/" by default. All data is stored in HDF5 format.

## Works that use this package

% Skyrmions near defects %
A. Derras-Chouk and E. M. Chudnovsky, (2020); https://arxiv.org/abs/2010.14683

% Thermal collapse of a skyrmion %
A. Derras-Chouk, E. M. Chudnovsky, and D. A. Garanin, J. Appl. Phys. 126, 083901 (2019); https://doi.org/10.1063/1.5109728

## References

% Random field XY model in three dimensions %
<a id="1">[1]</a>
D. A. Garanin, E. M. Chudnovsky, and T. Proctor, Phys. Rev. B **88**, 224418 (2013); https://doi.org/10.1103/PhysRevB.88.224418

% Pulse noise approach for classical spin systems %
<a id="1">[2]</a>
D. A. Garanin, Phys. Rev. E **95**, 013306 (2016); https://doi.org/10.1103/PhysRevE.95.013306

% Two-dimensional periodic boundary conditions for demagnetization interactions in micromagnetics %
<a id="1">[3]</a>
W. Wang, C. Mu, B. Zhang, Q. Liu, J. Wang, and D. Xue, Comput. Mater. Sci. **49**, 84
(2010); https://doi.org/10.1016/j.commatsci.2010.04.024
