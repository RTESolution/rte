![PyPI - Version](https://img.shields.io/pypi/v/rte)
[![Test](https://github.com/RTESolution/rte/actions/workflows/test.yml/badge.svg)](https://github.com/RTESolution/rte/actions/workflows/test.yml)

# Table of contents

- [About the project](#about-the-project)
- [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
    + [1. Install tagged package from PyPI](#1-install-tagged-package-from-pypi)
    + [2. Use development version from github](#2-use-development-version-from-github)
  * [Running the calculation](#running-the-calculation)
  * [Running the calculation in parallel](#running-the-calculation-in-parallel)
- [How to define parameter expressions](#how-to-define-parameter-expressions)
    + [Time](#time)
    + [Position](#position)
    + [Direction](#direction)
    + [Position again: point on sphere](#position-again--point-on-sphere)
  * [Working with parameter expressions](#working-with-parameter-expressions)
    + [Sampling](#sampling)
    + [Tree structure of expressions](#tree-structure-of-expressions)
    + [Inspecting](#inspecting)
    + [Modifying](#modifying)
- [Inspecting trajectories](#inspecting-trajectories)
    + [3D viewer](#3d-viewer)
- [RTE Method details](#rte-method-details)
- [License](#license)

# About the project

This package implements the calculation of the photon propagation, using the infinite series solution of the Radiative Transfer Equation (RTE) [[V. Allakhverdian, D. Naumov](https://arxiv.org/abs/2401.15698)].

RTE is a method developed for calculating light fluxes at an arbitrary point in space with an arbitrary source and an arbitrary detector at the registration point. See [this section](#rte-method-details) for details 

# Getting Started

## Prerequisites

RTE uses:
* [vegas integrator](https://vegas.readthedocs.io) for numerical calculation of the multidimensional integrals
* [vegas_params](https://github.com/RTESolution/vegas_params) for the definition of the integration space
* (optional) [pyqtgraph](https://www.pyqtgraph.org) for the visualization of the photon trajectories

All the dependencies are installed automatically when you install RTE package.

## Installation

Install using pip:

### 1. Install tagged package from PyPI
```shell
pip install rte
```
Optionally you can install the dependencies for the 3D visualisation and for running tests
```shell
pip install rte[gui, test]
```

### 2. Use development version from github
```shell
pip install git+https://github.com/RTESolution/rte.git
```

## Running the calculation

In order to do the calculation, the user needs to define the following objects:

* `rte.Source` - describing the emission of the photons - their initial time, position and direction. See [Sources section](#sources) for more information.
* `rte.Target` - describing the detection of the photons - their final time and position
* `rte.Medium` - contains constants of the uniform medium, in which the photon will be propagating
* `rte.Process` - main object which performs the calculation of the probability of a photon, emitted by the `Source` to reach the `Target` after scattering `Nsteps` times.

In other words, the `Source` and `Target` define the integration spaces of the initial and final state of the photon trajectory, while the `Process` defines also the intermediate steps.

Here is a full example script of such calculation:
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

#optionally: set the integrator parameters
p.vegas_kwargs = {'nitn':10, 'neval':10000} 
#run the numerical calculation using the `vegas` integrator
result = p.calculate() 
print(result) 
#65505.6(8.4) - this number may vary because of random sampling
```
The position `R`, time `T` and direction `s` arguments in `Source` and `Target` can either be fixed, or distributed (in which case the final calculation will integrate over it). 
See the [How to define parameter expressions](#how-to-define-parameter-expressions) section to see the examples.

The `Process.calculate` function evaluates the integral for the given setup of the source and target distributions.

### Running the calculation in parallel

The `Process.caluclate_map` function allows to perform parallel calculation of several processes, varying some of the parameters.

For example if we want to calculate for many fixed final state times for the example above:
```python
import numpy as np
T = np.geomspace(1e-8, 1e-7,21)
#calculate the results for give time
results = p.calculate_map(override={'tgt.T':T}, output='reshape')
#get only the mean values
import gvar #to get mean for each value in array
print(gvar.mean(results))
#[65515.26584496 42490.06720284 29069.30460075 20556.71848331
# 14857.42366258 10894.38730677  8063.57276645  6002.31205735
#  4481.53863444  3348.59624871  2499.77893137  1861.40452433
#  1380.65335953  1018.65576917   746.67916569   543.01496697
#   391.23942716   278.88939372   196.39585017   136.38307807
#    93.22038259]
```

### Sources
Currently we provide two types of sources, all of them can be used for constructing `rte.Process` objects:
#### Basic source
`rte.sources.Source`. position, time and direction are sampled ___independently___ of each other.
Initilaized as in the example above:
```python
src = rte.Source(R = vp.Vector([0,0,0]), #fixed point
                 T = 0, #photons emitted at fixed time 
                 s = vp.Direction() #uniform direction vector
                )
```
#### `rte.sources.TrackSource`:  creation of photons along a given track.
It samples track position, direction, speed, time and photon direction relative to track. And in the output it provides the photon starting points, similartly to the basic source.
```python
src = rte.source.TrackSource(
    T = vp.Uniform([0,1e-8]),
    track_pos = Vector([0,0,0]), #starting point 
    track_dir = vp.Direction(0.5,0), #45 degrees between X and Z axes
    track_speed = 3e8, #m/s
    photon_dir = Direction(cos_theta=0.9)#photon emission angles
)
```
# How to define parameter expressions
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
#same as above: ** operator used to repeat one value
R_box = vp.Vector(xyz = vp.Uniform([-1,1])**3)
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
### Position again: point on sphere
One can combine the `vp.Vector`, `vp.Fixed`, `vp.Direction` etc. to make a uniform position on a given sphere - for example, an optical detector
```python
detector_radius = 0.5
detector_center = vp.Vector([1,0,1])
R_sphere = detector_radius * vp.Direction() + detector_center
```
## Working with parameter expressions
In this section we explain parameter expressions in more detail.

The parameter expression is a recipee of how to construct the required values based on a set of random numbers from vegas.

### Sampling
You can request a random sample of any parameter expression:
```python
R = vp.Vector([0,0,1])
#get 3 "random" saples of a fixed value
xyz = R.sample(3)
#vector([[[0., 0., 1.]],
#        [[0., 0., 1.]],
#        [[0., 0., 1.]]])

#get 100 values in range [-1,1]
x = vp.Uniform([-1,1]).sample(100)
#np.array with 100 values

#get 1 random direction
s = vp.Direction().sample(1)
#vector([[[ 0.00657508, -0.06280371,  0.99800424]]])
```
This is useful for inspecting and debugging.

### Tree structure of expressions
`vp.Source`, `vp.Target` and even `vp.Process` are parameter expressions, which depend on other expressions: `R`,`T`,`s` etc.

This way the expression can be seen as a dependency tree of the parameters. You can see this tree just by printing the expression.

### Inspecting
The `print(src)` in the very first example shows:
```python
Source[2](
     --> R=Vector[0](
         --> xyz=Fixed[0]([[0 0 0]])
    )
     --> T=Fixed[0]([[0]])
     --> s=Direction[2](
         --> cos_theta=Uniform[1]([-1, 1])
         --> phi=Uniform[1]([0.0, 6.283185307179586])
    )
)
```

The number after the class name in `Source[2]` shows the number of degrees of freedom of this expression.
This is the number of dimensions our integral will have when using this source.
### Modifying

You can access and change the inner parameters with the following syntax:
```python
#get the parameter
source_time = src['T'] #vp.Fixed(0)
#set the parameter
src['R'] = R_plane
#access the nested parameter with dot separated path
src['s.phi'] = vp.Fixed(np.pi)
```
## Inspecting trajectories
`vp.Process` calculates the weights of many photon trajectories.
They are not stored by default, but enabling the `save_trajectory` flag will allow you to access them after the calculation
```python
p.save_trajectory=True #enable saving trajectories
p.sample(100) #generate 100 trajectories
trajectories = p.trajectory #get the results
```
Trajectories contain the position `R`, time `T` and direction `s` of each photon trajectory segment.

### 3D viewer 
This package provides a simple 3D viewer functions to inspect trajectories and photon sources. 
In order to use it, make sure to install this module with the "gui" dependencies:
```shell
pip install -U "rte[gui]"
```

and then you can run the viewer:
```python
from rte.viewer import plot_trajectory
#get the trajectories from the process
plot_trajectory(p.trajectory, 
                p.factor**0.95)
```

additionally one can inspect the photon source samples:
```python
from rte.viewer import plot_photons
plot_photons(src.sample(100), #generate 100 photons from source
             photon_size=1) #show each of the photons as 1m line
```


## RTE Method details

We solve the radiative transfer equation (RTE) in anisotropically scattering media as an infinite series. Each series term represents a distinct number of scattering events, with analytical solutions derived for zero and single scattering. Higher-order corrections are addressed through numerical calculations or approximations. The RTE solution corresponds to Monte Carlo sampling of photon trajectories with fixed start and end points. 
The medium in which light propagates is specified by the following parameters: $n(\lambda)$ - refraction index, $\mu_s(\lambda)$ - inverse scattering length, $\mu_a(\lambda)$ - inverse absorption length, $p_g(\hat{\bf{s}},\hat{\bf{s}}_1)$ - scattering indicatrix. For Mi scattering we are using Heney-Greenstein function as scattering indicatrix:

$p_g(\hat{\bf{s}},\hat{\bf{s}}_1) = \frac{1-g^2}{4\pi}(1+g^2-2g(\hat{\bf{s}},\hat{\bf{s}}_1))^{-3/2}$,

but the method in the article can work with an arbirtary indicatrix. 

In the article Green function for radiative transfer equation was found with iteration method. For each iteration was found exact formula for fluxes: 
    
$L=L_0+\sum\limits_{n=1}^\infty\delta L^{(n)}$, where

$L_0  = e^{-\mu_t ct}c \delta^3(\bf{r}-c\hat{\bf{s}}_0 t)\delta^2(\hat{\bf{s}}-\hat{\bf{s}}_0)$

and

$\delta L^{(n)} = e^{-\mu_a ct} P_n(\mu_s ct) \delta L_s^{(n)}$.

$P_n(\mu_s ct) = e^{-\mu_s ct}\frac{(\mu_s ct)^n}{n!}$ - Poisson probability for $n$ scattering events, $\delta L_s^{(n)}$ - factor containing all scattering information. Full expression one can find in [[V. Allakhverdian, D. Naumov](https://arxiv.org/abs/2401.15698)], formula (34). 

Knowing the green function, you can convolve it with an arbitrary source and detector, obtaining the complete signal recorded in the detector. This procedure is carried out in this program.


# License

`rte` is distributed under the terms of MIT License
