# Fortran Finite volume Fluid Dynamics Solver (F3DS) Package

F3DS package is a modern Fortran (such as Fortran 2003, 2008, ...) software of finite volume method (FVM) for fulid dynamics solvers. F3DS Package is composed of below:

- F3DS Flamework: Flamework for developing fluid dynamics solvers.
- F3DS Resource: Models and schemes for specific solvers.
- F3DS Collection: Solvers built by F3DS Framework & Resource.

Status: **DEVELOPMENT**  

Solvers are generalized. We can change parameters by a configuration file witch writen by JSON or XML file format. Of course we can also change a grid and an initial condition.  
We can easily build a solver for any non-viscous and viscous fluid using this flamework.  
We are working on improving parallel computing and developing features for unstructured mesh.

## Feature

### F3DS Flamework

The source code is located in "framework" directory. The static link library is located in libs/f3ds_framework.a.  
This flamework is desined by Object Oriented Desing (OOD) and supports Structure of Arrays (SoA) layout.
We provide a fast, maintainable FVM library and solvers.

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

- For structure grid
    - [x] MUSCL3 (3rd order)
        - [x] with minimod flux limmiter
        - [ ] with monotonized central flux limmiter
    - [x] WENO5 (5th order)
        - [x] with original smoothing indicator (WENO5) [Liu 1994]
        - [x] with smoothing indicator proposed by Jiang & Shu (WENO5-JS) [Jiang 1996]
        - [x] with monotonicity-preserving schemes (MP-WENO5-JS) [Balsara 2000]
- For unstructure grid
    - [ ] UMUSCL3 (3rd order)

#### Gradient scheme

- [x] Green-Gauss
- [x] Weighted Green-Gauss

#### Divergence scheme

- [x] Gauss divergence
- [x] Weighted Gauss divergence

#### Face gradient interpolation

- [x] Midpoint rule
- [x] Corrected midpoint [Nishikawa 2010]

#### Time Stepping

- Explicit
    - [x] Forward Euler
- Explicit Runge-Kutta
    - [x] 2nd order TVD (SPP) Runge-Kutta
    - [x] 3rd order TVD (SPP) Runge-Kutta
    - [ ] 4th order TVD (SPP) Runge-Kutta
    - [ ] Jameson-Baker 4 stage Runge-Kutta [Jameson 1929]

#### Grid System

- [x] Structure grid
- [x] Unstructure grid
- [ ] Adaptive Mesh Refinement (AMR)

#### Grid I/O

- Structure grid format
    - [x] Nishida Lab. legacy format ([NL Grid Toolbox](https://165.93.124.207/gitlab/tishikawa/nl-grid-toolbox) is helpfull for you.)
    - [ ] xyz
- Unstructure grid format
    - [ ] gmsh

#### Result I/O

- [ ] F3DS original
- [x] VTK ([VTKFortran](https://github.com/szaghi/VTKFortran) backend)

#### Initial Condition

- [x] Input from Nishida Lab. legacy format
- [ ] Unifoam flow generation

#### Parallelization

- [x] single-node multi-thread simd (backends none/openMP/openMP)
- [ ] multi-node multi-thread simd (backends MPI/MPI/openMP)
- [ ] multi-node multi-thread multi-gpu (backends MPI/MPI/{cuda, openACC})

#### Grid Decomposition for Parallelization

- [x] Structure (Simple decomposition)
- [ ] [Metis](https://github.com/KarypisLab/METIS)
- [ ] [Scotch](https://gitlab.inria.fr/scotch/scotch)
- [ ] [Zoltan](https://github.com/sandialabs/Zoltan)

#### Measurement Tools

- [ ] Sensor (Probe)
- [ ] Measurement surface profiler
- [x] Line plotter
- [x] Control volume profiler

#### Other tools

- [ ] Nearest Neighbor Search
    - [ ] [FLANN](https://github.com/flann-lib/flann)
    - [ ] [Annoy](https://github.com/spotify/annoy)
- [x] Data Description Language I/O
    - [x] JSON ([JSON-Fortran](https://github.com/jacobwilliams/json-fortran) backend)
    - [X] XML ([FoXy](https://github.com/Fortran-FOSS-Programmers/FoXy) backend)

### F3DS Resource

#### Five-equation model common tools

The source code is located in resource/five_equation_model_common. The static link library is located in libs/f5eq_common.a.  
This resource provides following schemes:

- [x] rho-THINC reconstruction
### F3DS Collection

### Five-equation model [Kapila 2001] [Allaire 2002]

- Source code location: collection/five_equation_model
- binary location: bins/f5eq.

The solver include the following terms.

- [x] with Kdiv(u) term [Kapila 2001]
- [ ] with cavitation model

### Viscosity five-equation model [Perigaud 2005] [Coralic 2014]

- Source code location: collection/five_equation_model
- binary location: bins/f5eq.

The solver include the following terms.

- [x] with Kdiv(u) term [Kapila 2001]
- [x] with surface tension term [Perigaud 2005] [Garrick 2017]
- [ ] with cavitation model
- [x] with gravity

### Euler equation

Currently under development
## How to compile

At first, we need to clone the repository and download submodules.

```:shell
git clone https://165.93.124.207/gitlab/tishikawa/f3ds.git
git submodule init
git submodule update
```

F3DS only support for Linux system now. We can compile F3DS by Makefile like this.

```:shell
make
```

Default compiler is set "gfortran" with a release build options.  
If you want to use "ifort" and debug options, you can use the following command.

```:shell
make COMPILER=ifort DEBUG=yes
```

More details can be found in 'make help'.

## How to use

### Use F3DS Framework and F3DS Resource

Please link static link libraries and mod files.

```
gfortran your_solver.f90 -o your_solver.exe -L/{your_f3ds_path}/f3ds/libs -I/{your_f3ds_path}/f3ds/mods f3ds_framework.a
```

### Use solvers

All binaries provided by F3DS is stored in "bins" directory. If you want use "f5eq", you type "./bins/f5eq".  
If you want more infomations, please read README.md in each collection directories.
