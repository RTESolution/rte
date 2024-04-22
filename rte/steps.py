import numpy as np
from scipy.spatial.transform import Rotation

import vegas_params as vp
from vegas_params.vector import vector
from vegas_params.utils import swapaxes
from .medium import Medium

from loguru import logger

#intermediate steps

#times
@vp.expression
@swapaxes
def Times_from_xi(self, xi):
    r"""Construct normalized times (t/T) from given xi values:

    Parameters:
    -----------
        xi: np.ndarray
            Array of shape (Nsamples,N,...) of values :math:`\xi\in[0,1]`
    Factor:
    -------
    .. math ::
        f = \prod_{i=1}^{N} (1-\xi_i)^{i-1}

    Returns: 
    --------
    array of normalized time of shape (Nsamples,N,...), t[0]=0, t[N]=1
    .. math ::
        \tau_i = \prod_{k=i}^{N}(1-\xi_k)
    """
    #calculate factor
    powers = np.arange(len(xi))-1
    powers = powers.reshape(-1,1)
    self.factor = np.prod(np.power((1.-xi)[1:],powers[1:]), axis=0)
    #calculate T
    xi_rev = np.flip(xi, axis=0)
    t_rev = np.cumprod(1-xi_rev, axis=0)
    t = np.flip(t_rev, axis=0)
    t = np.atleast_3d(t)
    #append last ending point
    t = np.append(t, np.ones(shape=(1,*t.shape[1:])), axis=0)
    #return result
    return t

#directions
class Directions_from_scattering_distr(vp.Expression):
    def __init__(self, medium, Nsteps:int, s0=vp.Direction(1,0)):
        self.medium = medium
        self.Nsteps = Nsteps
        super().__init__(s0=s0, 
                         #cos_theta = Fixed([0.5]*(Nsteps-1)),
                         cos_theta = vp.FromDistribution(medium.h_g.icdf)**(Nsteps-1),
                         phi = vp.Uniform([0,2*np.pi])**(Nsteps-1)
                        )
    @swapaxes
    def __call__(self, s0, cos_theta, phi):
        theta = np.arccos(cos_theta)
        step_dir = s0.reshape(-1,3)
        s = [step_dir]
        for theta_step,phi_step in zip(theta, phi):
            #transversal vector along X(if possible)
            s_T = np.cross(step_dir, [1,0,0]).view(vector)
            #recalculate for Y axis, if the cross with X axis was zero
            is_zero = np.isclose(s_T.mag(),0).squeeze()
            s_T[is_zero,:] = np.cross(step_dir[is_zero], [0,1,0])
            #apply the rotation
            rot_vec0 = (phi_step[:,np.newaxis]*step_dir)
            rot_axis = Rotation.from_rotvec(rot_vec0).apply(s_T).view(vector)
            #additionally normalize rot_axis
            rot_axis = rot_axis/rot_axis.mag()
            rot_vec1 = (theta_step*rot_axis)
            #
            step_new = Rotation.from_rotvec(rot_vec1).apply(step_dir) #apply rotation to the previous step
            step_dir = step_new
            s.append(step_new)
        s = np.asarray(s).view(vector)
        return s

#steps variants
class StepsUniform(vp.Expression):
    """Expression that describes the intermediate photon steps"""
    def __init__(self, Nsteps:int, medium:Medium, s0:vp.Direction=vp.Direction(1,0)):
        """
        Parameters
        ----------
        Nsteps:int
            number of intermediate photon steps (Number of scatterings)
        medium:Medium
            scattering medium, to define the probability pistribution
        s0:Expression
            Direction of the first step (default is fixed along positive Z)
        Returns
        -------
            Dictionary:
                't': time of each step, shape=(Nsamples,Nsteps,1),
                's': direction at each step, shape=(Nsamples,Nsteps,3)}
        """
        self.Nsteps = Nsteps
        self.medium = medium
        super().__init__(t=Times_from_xi(xi=1|vp.Uniform([0,1])**(Nsteps-1)),
                         s=s0|vp.Direction()**(Nsteps-1)
                        )
        
    @swapaxes
    def __call__(self, t, s):
        s = s.view(vector)
        s_dot_s = s[1:].dot(s[:-1])
        self.factor = np.prod(self.medium.scatter(s_dot_s), axis=0).swapaxes(0,1)
        return {'t':t, 's':s}

class StepsDistribution(vp.Expression):
    """Expression that describes the intermediate photon steps"""
    def __init__(self, Nsteps:int, medium:Medium, s0:vp.Direction=vp.Direction(1,0)):
  
        self.Nsteps = Nsteps
        super().__init__(t=Times_from_xi(xi=1|vp.Uniform([0,1])**(Nsteps-1)),
                         s=Directions_from_scattering_distr(medium, Nsteps, s0=s0)
                        )
        
    def __call__(self, t, s):
        self.factor = 1/(2*np.pi)**(self.Nsteps-1)
        return {'t':t, 's':s}
            
