# Atomistic Skyrmion Dynamics

Author: Amel Derras-Chouk

The code in this repository computes the evolution of magnetic moments in a two dimensional lattice by evaluating the Landau-Lifshitz-Gilbert equation at every lattice site. Code is written in Julia 1.0.4.

Required libraries:
 - HDF5 (Used to save files in .h5 format.)
 - FFTW (Fourier transform package used in computing dipolar interaction.)
 - PaddedViews (Also used in computing dipole-dipole interaction.)
 - LinearAlgebra (Used for vector operations.)
 - Distributed (To parallelize.)
 - BenchmarkTools (For testing.)
 - LoopVectorization (To speed up some loops.)

## Background

This code was written to study the behavior of skyrmions, topologically protected quasiparticles made up of atomic spins in magnetic materials. These spins obey the Landau-Lifshitz-Gilbert equation, which is implemented in this project. An example of the time evolution of spins is shown below.

<p align="center">
  <img src="https://github.com/aderras/ASkDyn/blob/main/spin-precession.gif" />
</p>

## Usage

To run the program, modify values in `01-scripts/UserInputs.jl`, navigate to the `01-scripts/` folder in terminal, and execute `julia Run.jl`

## Data

When the script is executed, results are saved to the directory `/02-data/.` All data is stored in HDF5 format.

## Works using this

A. Derras-Chouk and E. M. Chudnovsky, "Skyrmions near defects" (2020); arXiv: 2010.14683

A. Derras-Chouk, E. M. Chudnovsky, and D. A. Garanin, "Thermal collapse of a skyrmion," *J. Appl. Phys.* **126**, 083901, (2019); doi: 10.1063/1.5109728

## References

<a id="1">[1]</a>
D. A. Garanin, E. M. Chudnovsky, and T. Proctor, "Random field XY model in three dimensions," *Phys. Rev. B* **88**, 224418, (2013); doi: 10.1103/PhysRevB.88.224418

<a id="1">[2]</a>
D. A. Garanin, "Pulse noise approach for classical spin systems," *Phys. Rev. E* **95**, 013306 (2016); doi: 10.1103/PhysRevE.95.013306

<a id="1">[3]</a>
W. Wang, C. Mu, B. Zhang, Q. Liu, J. Wang, and D. Xue, "Two-dimensional periodic boundary conditions for demagnetization interactions in micromagnetics," *Comput. Mater. Sci.* **49**, 84 (2010); doi: 10.1016/j.commatsci.2010.04.024
