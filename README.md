# Fortran Finite-volume Fluid Dynamics Solver (F3DS) Flamework & Collection

F3DS is the modern Fortran (such as Fortran 2003, 2008, ...) library for development fluid dynamics solvers by Finite-volume method (FVM).  
This library support for Object Oriented Desing (OOD) and Structure of Arrays (SoA) layout. 
We provide a fast, maintainable FVM library and solvers.  

Status: **DEVELOPMENT**  

Solvers are specialized. We need to recompile by changing parameters and methods. 
We are currently working on improving to modify methods/parameters from external files.

## Feature

### Equation Collection

- [x] 5 equation model (in collection/five_equation_model. Binary name is "f5eq")
    - Additional terms
        - [x] with -K * div(u) term [Kapila 2001]
        - [ ] with surface tension term [Garrick 2017]
        - [ ] with viscosity term [Coralic 2014]
    - Approximate Riemann solvers
        - [x] HLLC [Toro 1997]
- [ ] Euler equation

### Flamework

#### EoS

- EoS
    - [ ] Ideal gas EoS
- Mixture EoS
    - [x] Stiffened-gas EoS
    - [ ] Noble-Abel Stiffened-gas EoS

#### Reconstruction Methods

- For structure grid (supported orthogonal grid)
    - [x] MUSCL3 (3rd order)
        - [x] with minimod flux limmiter
        - [ ] with monotonized central flux limmiter
    - [ ] MUSCL3 + THINC
    - [x] WENO5 (5th order)
        - [x] with original smoothing indicator (WENO5) [Liu 1994]
        - [x] with smoothing indicator proposed by Jiang & Shu (WENO5-JS) [Jiang 1996]
        - [x] with monotonicity-preserving schemes (MP-WENO5-JS) [Balsara 2000]
- For unstructure grid
    - [ ] Unstructured-MUSCL (3rd order)
    - [ ] Unstructured-MUSCL + Unstructure-THINC (3rd order)
    - [ ] k-expect WENO5 (5th order)

#### Time Stepping

- Explicit Runge-Kutta
    - [x] 2nd order TVD(SPP) Runge-Kutta
    - [x] 3rd order TVD(SPP) Runge-Kutta
    - [ ] 4th order TVD(SPP) Runge-Kutta
    - [ ] Jameson-Baker 4 stage Runge-Kutta [Jameson 1929]

#### Grid System

- [x] Structure grid
- [x] Unstructure grid
- [ ] Addaptive Mesh Refinement (AMR)

#### Grid I/O

- Structure grid format
    - [x] Nishida Lab. legacy format
    - [ ] xyz format
- Unstructure grid format
    - [ ] gmsh
    - [ ] CGNS
    - [ ] OpenFOAM

#### Result I/O

- [ ] F3DS original
- [x] VTK ([VTKFortran](https://github.com/szaghi/VTKFortran) backend)

#### Initial Condition

- [x] Input from Nishida Lab. legacy format

#### Parallelization Backend

- [x] OpenMP
- [ ] OpenMPI (single-node/multi-node)
- [ ] CUDA Fortran
- [ ] OpenACC

#### Grid Decomposition for Parallelization

- [x] Structure (Simple decomposition by OpemMP)
- [ ] [Metis](https://github.com/KarypisLab/METIS)
- [ ] [Scotch](https://gitlab.inria.fr/scotch/scotch)
- [ ] [Zoltan](https://github.com/sandialabs/Zoltan)

#### Measurement Tools

- [x] Sensor (Probe)
- [x] Measurement surface
- [ ] Virtual measurement surface
- [x] Line plotting tool

#### Other tools

- [ ] Nearest Neighbor Search
    - [ ] [FLANN](https://github.com/flann-lib/flann)
    - [ ] [Annoy](https://github.com/spotify/annoy)
- [x] Data Description Language I/O
    - [x] JSON ([JSON-Fortran](https://github.com/jacobwilliams/json-fortran) backend)
    - [x] XML ([FoXy](https://github.com/Fortran-FOSS-Programmers/FoXy) backend)

## How to compile

F3DS only support for Linux system now. We can compile F3DS by Makefile like this.

```:shell
make
```

Default compiler is set "gfortran". And default option is "-O3 -march=native -ffree-line-length-none -fopenmp -cpp".  
If you want to use "ifort" and "-fast -openmp" options, you shoud type following command in your terminal.

```:shell
make FC=ifort FFLAGS=-fast -openmp
```

## How to use

All binaries provided by F3DS is stored in "bins" directory. If you want use "f5eq", you type "./bins/f5eq".  
If you want more infomations, please read README.md in each collection directories.
