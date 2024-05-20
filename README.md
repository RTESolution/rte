# Table of contents
<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [About the project](#about-the-project)
- [Getting Started](#getting-started)
   * [Prerequisites](#prerequisites)
   * [Installation](#installation)
      + [1. Use development version from github](#1-use-development-version-from-github)
      + [2. Install tagged package from PyPI](#2-install-tagged-package-from-pypi)
   * [Running the calculation:](#running-the-calculation)
- [Usage](#usage)
   * [Parameter examples](#parameter-examples)
      + [Time](#time)
      + [Position](#position)
      + [Direction](#direction)
      + [Position again: sphere](#position-again-sphere)
- [License](#license)

<!-- TOC end -->

# About the project

This package implements the calculation of the photon propagation, using the infinite series solution of the Radiative Transfer Equation [[V.Allakhverdian et al.]()].

*TODO*

# Getting Started

## Prerequisites

RTE uses:
* [vegas integrator](https://vegas.readthedocs.io) for numerical calculation of the multidimensional integrals
* [vegas_params](https://github.com/RTESolution/vegas_params) for the definition of the integration space
* (optional) [pyqtgraph](https://www.pyqtgraph.org) for the visualization of the photon trajectories

All the dependencies are installed automatically when you install RTE package.

## Installation

Install using pip:

### 1. Use development version from github
```shell
pip install git+https://https://github.com/RTESolution/rte.git
```
Optionally you can install the dependencies for the 3D visualisation and for running tests
```shell
pip install "git+https://https://github.com/RTESolution/rte.git[git, test]"
```
### 2. Install tagged package from PyPI
*TODO*



## Running the calculation:

```python
import rte #main package module with Target, Source and Process classes
import vegas_params as vp #definition of the integration spaces

#define the photon source
src = rte.Source(R = vp.Vector([0,0,0]), #fixed point
                 T = 0, #photons emitted at fixed time 
                 s = vp.Direction() #uniform direction vector
                )
#define the target
tgt = rte.Target(R = vp.Vector([1,0,1]), #fixed detector point
                 T = 1e-8 #fixed detection time
                )

#define the process calculator
p = rte.Process(src, tgt, 
                Nsteps=1, #number of photon scatterings in the trajectory
                medium=rte.medium.water #the medium properties
               )

#run the calculation
result = p.calculate(nitn=10, neval=10000) #run the numerical calculation using the `vegas` integrator
print(result)
```
The position `R`, time `T` and direction `s` arguments in `Source` and `Target` can either be fixed, or distributed (in which case the final calculation will integrate over it). 
See the [parameter examples](#parameter-examples) section to see how to define them.

# Usage

In order to do the calculation, the user needs to define the following objects:

* `rte.Source` - describing the emission of the photons - their initial time, position and direction
* `rte.Target` - describing the detection of the photons - their final time and position
* `rte.Medium` - contains constants of the uniform medium, in which the photon will be propagating
* `rte.Process` - main object which performs the calculation of the probability of a photon, emitted by the `Source` to reach the `Target` after scattering `Nsteps` times.

In other words, the `Source` and `Target` define the integration spaces of the initial and final state of the photon trajectory, while the `Process` defines also the intermediate steps.

## Parameter examples
RTE uses [vegas_params](https://github.com/RTESolution/vegas_params) for the definition of the integration space:
```python
import vegas_params as vp
```
### Time
Time is a scalar variable, it can either be fixed or uniformely distributed
```python
#exactly at the moment of T=1s
T_fixed = vp.Fixed(1)
#Will be integrated between T=1s and T=10s
T_fixed = vp.Uniform(1,10)
```
### Position
Position should be defined as vector: vp.Vector ensures that output has a correct shape
```python
#fixed at the coordinates (x=1m, y=2m, z=3m)
R_fixed = vp.Vector(xyz = vp.Fixed([1,2,3])) 
#same as above: Vector assumes the xyz is fixed
R_fixed = vp.Vector(xyz = [1,2,3])  
#XY fixed at 0, but Z is uniform in [-10m, 20m]
R_line  = vp.Vector(xyz = vp.Fixed([0,0])|vp.Uniform([-10,20])) 
#X is uniform in [-1,1]m, Y is uniform in [-2,2]m, but Z=0 
R_plane = vp.Vector(xyz = vp.Uniform([[-1,1],[-2,2]])|0) 
#Uniformely distributed in a box at (0,0,0) with size=2m
R_box = vp.Vector(xyz = vp.Uniform([(-1,1),(-1,1),(-1,1)])) 
```
### Direction
Direction is a unit 3D vector, so you can define it as `vp.Vector` directly
```python
S_up = vp.Vector([0,0,1]) 
```
Or use `vp.Direction` to define the unit vector in the spherical coordinates:
```python
#same as previous example: fixed vector along Z axis
S_up = vp.Direction(cos_theta=1, phi=0)
#fully isotropic case: all directions are equal
S_uniform = vp.Direction(cos_theta=vp.Uniform([-1,1]), phi=vp.Uniform([0,2*np.pi]))
#same as above: default arguments are to make uniform direction
S_uniform = vp.Direction()
#looking up at fixed zenith angle. It will still be integrated over phi=[0, 2*pi]
S_fixed_zenith = vp.Direction(cos_theta=0.9) 
```
### Position again: sphere
One can combine the `vp.Vector`, `vp.Fixed`, `vp.Direction` etc. to make a uniform position on a given sphere - for example, an optical detector
```python
detector_radius = 0.5
detector_center = vp.Vector([1,0,1])
R_sphere = detector_radius * vp.Direction() + detector_position
```

# License

*TODO*





