# Fortran Fine-volume Fluid Dynamics Solver (F3DS) Flamework & Collection

F3DS is the modern Fortran (such as Fortran 2003, 2008, ...) library for development fluid dynamics solvers by Fine-volume method (FVM).  
This library support for Object Oriented Desing (OOD) and Structure of Arrays (SoA) layout. 
We provide a fast, maintainable FVM library and solvers.  

Status: **DEVELOPMENT**  

Solvers are generalized. We can change parameters by a configuration file witch writen by JSON or XML file format. Of course we can also change a grid and an initial condition.  
We can easily build a solver for any non-viscous fluid using this flamework.  
We are currently working on developing a function for computing viscous flux.

## Feature

### F3DS Collection

- [x] 5 equation model (in collection/five_equation_model. Binary name is "f5eq")
    - Additional terms
        - [x] with -K * div(u) term [Kapila 2001]
        - [x] with surface tension term [Garrick 2017]
        - [ ] with viscosity term [Coralic 2014]
        - [ ] with cavitation model

- [ ] Euler equation

### F3DS Flamework

#### EoS

- [ ] Ideal gas EoS
- [x] Stiffened-gas EoS
- [ ] Noble-Abel Stiffened-gas EoS

#### Approximate Riemann solvers

- HLL family
    - [x] HLLC [Toro 1997]
    - [ ] HLLE (also called HLL)
    - [ ] Rusanov
- AUSM family
    - [ ] AUSM
    - [ ] SLAW
- Rotate Riemann solvers
    - [ ] Rotate
    - [ ] Rotate-hybrid

#### Reconstruction Methods

- For structure grid (supported orthogonal grid)
    - [x] MUSCL3 (3rd order)
        - [x] with minimod flux limmiter
        - [ ] with monotonized central flux limmiter
    - [ ] MUSCL3 + THINC
    - [x] MUSCL3 + rho-THINC (5-equation model only)
    - [x] WENO5 (5th order)
        - [x] with original smoothing indicator (WENO5) [Liu 1994]
        - [x] with smoothing indicator proposed by Jiang & Shu (WENO5-JS) [Jiang 1996]
        - [x] with monotonicity-preserving schemes (MP-WENO5-JS) [Balsara 2000]
- For unstructure grid
    - [ ] UMUSCL3 (3rd order)
    - [ ] UMUSCL3 + THINC (3rd order)

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
    - [x] Nishida Lab. legacy format ([NL Grid Toolbox](https://165.93.124.207/gitlab/tishikawa/nl-grid-toolbox) is helpfull for you.)
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
- [ ] Unifoam flow generation

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

- [ ] Sensor (Probe)
- [ ] Measurement surface
- [ ] Virtual measurement surface
- [ ] Line plotting tool

#### Other tools

- [ ] Nearest Neighbor Search
    - [ ] [FLANN](https://github.com/flann-lib/flann)
    - [ ] [Annoy](https://github.com/spotify/annoy)
- [x] Data Description Language I/O
    - [x] JSON ([JSON-Fortran](https://github.com/jacobwilliams/json-fortran) backend)
    - [X] XML ([FoXy](https://github.com/Fortran-FOSS-Programmers/FoXy) backend)

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