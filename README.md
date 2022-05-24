# Fortran Fine-volume Fluid Dynamics Solver (F3DS) Flamework & Collection

F3DS is the Modern Fortran (such as Fortran 2003, 2008, ...) library for development fluid dynamics solvers by Fine-volume method (FVM).  
This library support for Object Oriented Desing (OOD) and Structure of Arrays (SoA) layout. 
We provide a fast, maintainable FVM library and solvers.  

Status: **DEVELOPMENT**  

Many methods are specialized. We need to recompile by changing parameters and methods. 
We are currently working on improving to modify methods/parameters from external files.

## Feature

### Equation Collection

- [x] 5 equation model (in collection/five_equation_model. Binary name is "f3ds5eq")
    - Additional terms
        - [x] with -K * div(u) term [Kapila 2001]
        - [ ] with surface tension term [Garrick 2017]
        - [ ] with viscosity term [Coralic 2014]
    - Nonviscosity flux approximation
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

- [x] MUSCL3 (3rd order, Support for Cartesian coordinate system only)
    - [x] with minimod flux limmiter
    - [ ] with monotonized central flux limmiter
- [x] WENO5 (5th order, Support for Cartesian coordinate system only)
    - [x] with original smoothing indicator (WENO5)
    - [x] with JS-type smoothing indicator (WENO5-JS)
- [ ] MUSCL3+THINC (3rd order, Support for Cartesian coordinate system only)
- [ ] Unstructure MUSCL (3rd order)
- [ ] Unstructure MUSCL + Unstructure THINC (3rd order)
- [ ] k-expect WENO5 (5th order)

#### Time Stepping

- Explicit Runge-Kutta
    - [x] 2nd order TVD Runge-Kutta
    - [x] 3rd order TVD Runge-Kutta
    - [ ] 4th order TVD Runge-Kutta
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
- [ ] Metis
- [ ] Scotch

#### Measurement Tools

- [x] Sensor (Probe)
- [ ] Measurement surface
- [ ] Virtual measurement surface

#### Other tools

- [ ] Nearest Neighbor Search
    - [ ] [FLANN](https://github.com/flann-lib/flann)
    - [ ] [Annoy](https://github.com/spotify/annoy)
- [x] Data Description Language I/O
    - [x] JSON ([JSON-Fortran](https://github.com/jacobwilliams/json-fortran) backend)

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

All binaries provided by F3DS is stored in "bins" directory. If you want use "f3ds5eq", you type "./bins/f3ds5eq".  
If you want more infomations, please read README.md in each collection directories.