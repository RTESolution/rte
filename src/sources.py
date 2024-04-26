import numpy as np
from vegas_params import expression, Uniform
from vegas_params.vector import vector,Vector, Scalar, Direction, expression

from dataclasses import dataclass

@dataclass 
class Point:
    """A container for the photon trajectory points"""
    R: vector
    T: np.array
    s: vector
    def __getitem__(self, n):
        return Point(self.R[n],self.T[n],self.s[n])
    def __len__(self):
        return len(np.atleast_1d(self.R))
    def iter(self):
        for n in range(len(self)):
            yield self[n]

@expression
class Spherical:
    r:Scalar = Scalar(Uniform([0,1]))
    s:Direction = Direction()
    def __call__(self, r,s):
        self.factor=r**2
        return r*s
        
@expression
class Target:
    R: Vector = Spherical()
    T: Scalar = Uniform([0,1])
    s: Direction = Direction(1,0)#fixed: not used
    
    def __call__(self, R, T, s):
        self.factor=1
        return Point(np.atleast_3d(R),
                     np.atleast_3d(T),
                     np.atleast_3d(s))

class Source(Target):
    pass