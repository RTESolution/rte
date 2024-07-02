import numpy as np
import vegas_params as vp
#from vegas_params import expression, Expression, Uniform, Fixed, FromDistribution
from vegas_params.vector import vector
from vegas_params.utils import swapaxes, swap_axes
from vegas_params.utils import mask_arrays, unmask_arrays, save_input
#from vegas_params import integral
from scipy.spatial.transform import Rotation

from itertools import product
from functools import partial

from .sources import Point #dataclass which holds trajectory points
from .steps import StepsUniform, StepsDistribution
from .medium import Medium

from loguru import logger

def align_rotation_to_vector(v:vector)->Rotation:
    v = np.asarray(v).view(vector)
    n = v/v.mag()
    theta = np.arccos(n.z)
    phi = np.arctan2(n.y,n.x)
    return Rotation.from_euler('yz',np.array([theta,phi]).T.squeeze())

# --- Expressions are defined here
class Process(vp.Expression):
    r"""The calculator of the RTE term :math:`\delta L^{(n)}` of the orrder `Nsteps`:

    .. math ::
        \delta L^{(n)} = 
        A_{src}(\vec{R}_{src}, T_{src}, \hat{s}_{src}) \cdot A_{det}(\vec{R}_{det}, T_{det}, \hat{s}_{det}) \cdot e^{-\mu_t\,c\,(T_{det}-T_{src})}\cdot \left[\mu_s\,c\,(T_{det}-T_{src})\right]^n \cdot \prod_{i=1}^{n} f(\hat{s}_i\cdot\hat{s}_{i-1})\cdot \delta^3\left(\vec{R}_{tgt}-\vec{R}_{src}-\sum_{i=1}^{n} \hat{s}_{i} c (t_{i+1}-t_{i})\right)
    """
    def __init__(self, 
                 src:vp.Expression, 
                 tgt:vp.Expression, 
                 medium:Medium, 
                 Nsteps:int,
                 save_trajectory:bool=False, 
                 use_masking:bool=True,
                 use_uniform_sampling=False
                ):
        self.medium = medium
        self.Nsteps = Nsteps
        self.save_trajectory = save_trajectory
        self.use_masking = use_masking

        if use_uniform_sampling:
            super().__init__(src=src,tgt=tgt,
                             steps=StepsUniform(Nsteps,medium)
                            )
        else:
            super().__init__(src=src,tgt=tgt,
                             steps=StepsDistribution(Nsteps,medium)
                            )

    def calculate_last_step(self, R, T, t_i, s_i):
        """calculate the length and direction of the last step, needed to get to the given point R,T"""
        dt_i = np.diff(t_i, axis=0)
        N = np.sum(s_i*dt_i, axis=0, keepdims=True).view(vector)
        R1 = R-N*T*self.medium.c
        R1_2 = R1.mag2()
        gamma = 1-N.mag2()
        R1N = R1.dot(N)
        #calculate the missing step length
        #the equation is 
        # gamma*(l**2) - 2*R1N*l - R1**2 = 0
        is_linear = np.isclose(gamma, 0)#in these cases the solution is linear
        #determinant
        D = R1N**2+R1_2*gamma
        
        step_length = np.where(is_linear, -R1_2/(2*R1N), (R1N+np.sqrt(D))/(gamma))
        step_direction = R1/step_length+N
        #calculate factor
        factor = 1/(step_length*np.abs(R1N - step_length*gamma))
        factor *= 1/T
        factor *= 1/self.medium.c
        return step_length, step_direction, factor
        
    @swapaxes
    def __call__(self, src, tgt, steps):
        #target points
        R = tgt.R-src.R
        T = tgt.T-src.T
        #convenience stuff
        s_i = steps['s']
        t_i = steps['t']
        c = self.medium.c
        Nsamples = T.shape[1]
        #set the initial values
        #s_i[0] = src.s
        #rotate the trajectories
        r = align_rotation_to_vector(src.s)
        s_rot = np.asarray([r.apply(s) for s in s_i]).view(vector)
        s_i = s_rot
        #apply the delta-function - calculate the last step
        l_last,s_last, delta_factor = self.calculate_last_step(R,T, t_i, s_i)
        T1 = T-l_last/c
        # --- mask the unphysical values
        mask = (T*c< R.mag()).flatten()
        mask |= (l_last<0).flatten()
        mask |= ((T1<0)|(T1>T)).flatten()
        
        if self.use_masking:
            if np.all(mask):
                return np.zeros((1,Nsamples))
            l_last, s_last, delta_factor = mask_arrays([l_last, s_last, delta_factor], mask=mask, axis=1)
            s_i, t_i, T1, T = mask_arrays([s_i, t_i, T1, T], mask=mask, axis=1)
            
        #gather the resulting values
        s = np.concatenate([s_i, s_last], axis=0).view(vector)
        t = np.concatenate([t_i*T1, T], axis=0)
        #calculate factors
        scat_factor = self.medium.scatter(s[-2].dot(s[-1]))[np.newaxis,:]
        att_factor = self.medium.attenuation(T, n_scattering=self.Nsteps)
        step_factor = (T1/T)**(self.Nsteps-1)
        self.factor = scat_factor*att_factor*step_factor*delta_factor
        self.factor = self.factor.reshape(s.shape[1])
        self.factor*= self.medium.c
        #optional stuff
        if self.use_masking:
            # --- unmask all arrays --- 
            t,s = unmask_arrays([t,s], mask=mask, axis=1)
            self.factor = unmask_arrays(self.factor, mask=mask, axis=0)
        else:
            self.factor[mask]=0
        if self.save_trajectory:
            #calculate the trajectory (for debugging, can be disabled)
            dt = np.diff(t, axis=0)
            dr = c*np.cumsum(dt*s, axis=0)
            r = np.concatenate([np.zeros((1,*dr.shape[1:])),dr],axis=0)+src.R
            self.trajectory = swap_axes(Point(r,t,s))#store the trajectory
        
        #return ones - they will be multiplied by the factor automatically
        return np.ones(shape=(1,Nsamples))

    def __getitem__(self, key):
        expr = super()
        for token in key.split('.'):
            expr = expr.__getitem__(token)
        return expr
    
    def __setitem__(self, key, value):
        try:
            path, key = key.rsplit('.',1)
            expr = self[path]
        except ValueError:
            expr = super()
        expr.__setitem__(key, value)
        
    def calculate(self, **vegas_kwargs):
        return vp.integral(self)(**vegas_kwargs)

    def calculate_with(self, override:dict, **vegas_kwargs):
        for key, value in override.items():
            if not key.startswith('_'):
                self[key]=value
        return self.calculate(**vegas_kwargs)

    def calculate_map(self, override, output='dict', map_function=map, **vegas_kwargs):
        keys,values = override.keys(), override.values()
        kwargs = [dict(zip(keys,v)) for v in product(*values)]
        #run the actual calculation
        result =  list(map_function(partial(self.calculate_with, **vegas_kwargs), kwargs))
        #output the result
        if output=='reshape':
            return np.asarray(result).reshape([len(v) for v in values])
        elif output=='dict':
            return [dict(**d,result=r) for d,r in zip(kwargs,result)]
        else:
            return result