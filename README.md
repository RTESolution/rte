# About the project

This package implements the calculation of the photon propagation, using the infinite series solution of the Radiative Transfer Equation [[V.Allakhverdian et al.]()].

*TODO*

# Installation

Install using pip:

## 1. Use development version from github
```shell
pip install git+https://https://github.com/RTESolution/rte.git
```
## 2. Install tagged package from PyPI
*TODO*

# Usage

### Basic usage example:
```python
import rte #main package module with Target, Source and Process classes
import vegas_params as vp #

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




