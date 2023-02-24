# Fortran Finite volume Fluid Dynamics Solver (F3DS) Package

[![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)

F3DS package is a modern Fortran software of finite volume method (FVM) for fulid dynamics solvers. F3DS Package is composed of below:

- F3DS Flamework: Flamework for developing fluid dynamics solvers.
- F3DS Resource: Models and schemes for specific solvers.
- F3DS Collection: Solvers built by F3DS Framework & Resource.

Status: **DEVELOPMENT**  

Solvers are generalized. We can change parameters by a configuration file witch writen by JSON file format. Of course we can also change a grid and an initial condition.  
We can easily build a solver for any non-viscous and viscous fluid using this flamework.  
We are working on improving parallel computing and developing features for unstructured mesh.

## Feature

### F3DS Flamework

The source code is located in "framework" directory. The static link library is located in libs/f3ds_framework.a.  
This flamework is desined by Object Oriented Desing (OOD) and supports Structure of Arrays (SoA) layout.
We provide a fast, maintainable FVM library and solvers.

#### EoS

- [x] Ideal gas EoS
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
        - [x] Minimod flux limmiter (Minmod MUSCL3)
        - [ ] Monotonized central flux limmiter
    - [x] WENO5 (5th order)
        - [ ] Original smoothing indicator (WENO5) [Liu 1994]
        - [x] Smoothing indicator proposed by Jiang & Shu (WENO5-JS) [Jiang 1996]
        - [x] WENO-Z [Borges 2008]
    - [x] Monotonicity preserving scheme (MP)
        - We can conbine WENO5-JS (MP-WENO5-JS) [Balsara 2000] and other reconstruction metheds.
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

#### Face gradient calculation

- [x] Central difference

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
- [x] Surface profiler
- [x] Line plotter
- [x] Control volume profiler

#### Other tools

- [ ] Nearest Neighbor Search
    - [ ] [FLANN](https://github.com/flann-lib/flann)
    - [ ] [Annoy](https://github.com/spotify/annoy)
- [x] Data Description Language I/O
    - [x] JSON ([JSON-Fortran](https://github.com/jacobwilliams/json-fortran) backend)

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

- [x] Kdiv(u) term [Kapila 2001]
- [ ] Cavitation model

### Viscosity five-equation model [Perigaud 2005] [Coralic 2014]

- Source code location: collection/five_equation_model
- binary location: bins/f5eq.

The solver include the following terms.

- [x] Kdiv(u) term [Kapila 2001]
- [x] Surface tension term [Perigaud 2005] [Garrick 2017]
- [ ] Cavitation model
- [x] Gravity term

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

## Lisence

### F3DS Package

[![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)

F3DS Package is an open sorce software, it is distributed under the MIT license. More details of the MIT license available at the following file: [LICENSE](LICENSE).  
Contributors names are listed below:  

- Tatsumasa Ishikawa

### Third party libraries

| Libraries                                                    | Lisence                                                                                                       | Copyright                               |
|--------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------------------------|
|[VTK Fortran](https://github.com/szaghi/VTKFortran)           | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[FACE](https://github.com/szaghi/FACE)                        | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[FoXy](https://github.com/Fortran-FOSS-Programmers/FoXy)      | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[PENF](https://github.com/szaghi/PENF)                        | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[StringiFor](https://github.com/szaghi/StringiFor)            | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[JSON Fortran](https://github.com/jacobwilliams/json-fortran) | [![License](https://img.shields.io/badge/license-BSD3-red.svg)](https://opensource.org/license/bsd-3-clause/) | Copyright (c) 2014-2021, Jacob Williams |

