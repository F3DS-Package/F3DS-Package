# Fortran Finite volume Fluid Dynamics Solver (F3DS) Package

[![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)

Status: **DEVELOPMENT**

## Introduction

F3DS package is a modern Fortran software of finite volume method (FVM) for fluid dynamics solvers. F3DS Package is composed of below:

- F3DS Framework: Framework for developing fluid dynamics solvers.
- F3DS Resource: Models and schemes for specific solvers.
- F3DS Collection: Solvers built by F3DS Framework & Resource.

The framework is designed to discrete any governing equations.
We hope to contribute to many researchers and engineers through this flexible framework-centered software package.
And we also hope that many researchers and engineers can collaborate through this framework.

## F3DS Framework

- Source code location: framework/src
- Static link library location: libs/f3ds_framework.a

This code uses object-oriented design (OOD) and dependency injection strategy to design the generic FVM framework. But, Users do not need to understand OOD.
Variables are expressed as a Structure of Arrays layout, which helps speed up the solver.
You can easily build a useful solver for any non-viscous and viscous fluid using this framework.
We are working on improving parallel computing and developing features for unstructured mesh.

### Features

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

- For structured grid
    - [x] MUSCL3 (3rd order)
        - [x] Minimod flux limmiter (Minmod MUSCL3)
        - [ ] Monotonized central flux limmiter
    - [x] WENO5 (5th order)
        - [ ] Original smoothing indicator (WENO5) [Liu 1994]
        - [x] Smoothing indicator proposed by Jiang & Shu (WENO5-JS) [Jiang 1996]
        - [x] WENO-Z [Borges 2008]
    - [x] Monotonicity preserving scheme (MP)
        - We can conbine WENO5-JS (MP-WENO5-JS) [Balsara 2000] and other reconstruction methods.
- For unstructured grid
    - [ ] UMUSCL3 (3rd order)

#### Gradient scheme

- [x] Green-Gauss
- [x] Weighted Green-Gauss

#### Divergence scheme

- [x] Gauss divergence
- [x] Weighted Gauss divergence

#### Face valiables interpolation

- [x] Midpoint interpolator
- [x] Weighted linear interpolator
- [x] Upwind interpolator
- [x] Linear upwind interpolator

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

- [x] Structured grid
- [x] Unstructured grid
- [ ] Adaptive Mesh Refinement (AMR)

#### Grid I/O

- Structured grid format
    - [x] F3DS legacy format
    - [ ] xyz
- Unstructured grid format
    - [ ] gmsh
    - [ ] CGNS
    - [ ] OpenFOAM

#### Result I/O

- [ ] F3DS original
- [x] VTK ([VTKFortran](https://github.com/szaghi/VTKFortran) backend)

#### Initial Condition

- [x] Input from Nishida Lab. legacy format
- [ ] Uniform flow generation

#### Parallelization

- [x] OpenMP
- [ ] MPI
- [ ] GPU

#### Grid Decomposition for Parallelization

- [x] Simple decomposition according to cell number
- [ ] Recursive Coordinate Bisection
- [ ] [Zoltan](https://github.com/sandialabs/Zoltan) interface

#### Measurement Tools

- [ ] Sensor (Probe)
- [x] Surface profiler
- [x] Line plotter
- [x] Control volume profiler

#### Other tools

- [ ] Nearest Neighbor Search
- [x] Data Description Language I/O
    - [x] JSON ([JSON-Fortran](https://github.com/jacobwilliams/json-fortran) backend)

## F3DS Resource

### Five-equation model common tools

- Source code location: resource/five_equation_model_common
- Static link library location: libs/f5eq_common.a

This resource provides following schemes:

- [x] rho-THINC reconstruction
## F3DS Collection

### Advection Equation

- Source code location: collection/advection_equation
- Binary location: bins/fadvection

This code solves the equation that transports the scaler variable Phi. It will give you a hint as to the taste of the framework.

### Diffusion Equation

- Source code location: collection/diffusion_equation
- Binary location: bins/fdiffusion

This code solves the equation that diffuse the scaler variable Phi.

### Five-equation model [Kapila 2001] [Allaire 2002]

- Source code location: collection/five_equation_model
- Binary location: bins/f5eq

The solver include the following terms.

- [x] Kdiv(u) term [Kapila 2001]
- [ ] Cavitation model

### Viscosity five-equation model [Perigaud 2005] [Coralic 2014]

- Source code location: collection/five_equation_model
- Binary location: bins/f5eq

The solver include the following terms.

- [x] Kdiv(u) term [Kapila 2001]
- [x] Surface tension term [Perigaud 2005] [Garrick 2017]
- [ ] Cavitation model
- [x] Gravity term

### Euler equation

Currently under development
## How to compile & setup

First, you need to clone the repository and download submodules.

```:shell
git clone https://github.com/F3DS-Package/F3DS-Package.git
git submodule init
git submodule update
```

F3DS only support for Linux system now. You can compile F3DS by Makefile like this.

```:shell
make
```

Default compiler is set "gfortran" with a release build options.  
If you want to use "ifort" and debug options, you can use the following command.

```:shell
make COMPILER=ifort DEBUG=yes
```

More details can be found in the script shown by the 'make help' command.  
Finally, to install F3DS Package for your computer, run the bellow:

```:shell
make install
```
The default installation path is '/opt/f3ds-package'. You can change the installation path by the 'PLEFIX={path}' option.   
Executing the script 'setenv.sh', you can easily set the environment variables of F3DS package.
We recommend writing the following script to your bashrc or profile.

```
source /{your_f3ds_package_path}/setenv.sh
```

### Use F3DS Framework and F3DS Resource

Please link static link libraries and mod files.

```
gfortran your_solver.f90 -o your_solver.exe -LF3DS_LIBS -IF3DS_MODS f3ds_framework.a
```

### Use solvers

All binaries provided by F3DS collection are stored in "bins" directory. If you set the environment variables using "setenv.sh", binaries are already set to your environment.
More information can be found in README.md in each collection directory.

## Lisence

### F3DS Package

[![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)

F3DS Package is an open source software, it is distributed under the MIT license. More details of the MIT license available at the following file: [LICENSE](LICENSE).  
Contributors names are listed below:  

- Tatsumasa Ishikawa

### Third party libraries

We use the following software to develop this software.
We would like to thank the authors for the development of useful modern Fortran software.

| Libraries                                                    | License                                                                                                       | Copyright                               |
|--------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------------------------|
|[VTK Fortran](https://github.com/szaghi/VTKFortran)           | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[FACE](https://github.com/szaghi/FACE)                        | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[FoXy](https://github.com/Fortran-FOSS-Programmers/FoXy)      | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[PENF](https://github.com/szaghi/PENF)                        | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[StringiFor](https://github.com/szaghi/StringiFor)            | [![License](https://img.shields.io/badge/license-MIT-red.svg)](https://opensource.org/license/mit/)           | Copyright (c) 2022 Stefano Zaghi        |
|[JSON Fortran](https://github.com/jacobwilliams/json-fortran) | [![License](https://img.shields.io/badge/license-BSD3-red.svg)](https://opensource.org/license/bsd-3-clause/) | Copyright (c) 2014-2021, Jacob Williams |

