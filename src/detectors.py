from .sources import Point
import vegas_params as vp
import numpy as np

class DetectorSpherical(vp.Expression):
    """A spherical detector with a given radius and center position"""
    def __init__(self,
                 center = vp.Vector([0,0,0]),
                 radius = vp.Scalar(0.216),
                 T = vp.Uniform([0,1])
                ):
        super().__init__(center = vp.Vector(center),
                         radius = vp.Scalar(radius),
                         T = vp.Scalar(T),
                         _R_local = vp.Direction(),
                         _s_local = vp.Vector([0,0,1]) #not used
                        )
        self.soften_parameter = 0.01 #a parameter for making a soft aiming function
        
    def __call__(self, center, radius, T, _R_local, _s_local):
        R = vp.Vector.__call__(_R_local) * vp.Scalar.__call__(radius) + vp.Vector.__call__(center)
        s = vp.Vector.__call__(_s_local)
        self.factor = (radius**2).squeeze()
        return Point(R, T, s)

    def efficiency(self, p: Point)->np.array:
        #convert to local RF
        R_local = p.R-self['center'].sample()
        N_local = R_local/R_local.mag()
        #discard rays coming from inside
        cosTheta_to_normal = N_local.dot(p.s).squeeze()
        eff_valid = 1.*(cosTheta_to_normal<0)
        #multiply by -cosTheta to project rays on the surface
        eff_valid *= -cosTheta_to_normal
        return eff_valid

    def get_intersection(self,
                         p: Point,
                         speed_of_light:np.array,
                         soft=False
                        )->(Point, np.array):
        r0, s = p.R, p.s
        R = self['radius'].sample()
        r = self['center'].sample()-r0
        sr = s.dot(r)
        #solving equation for the intersection point
        #equation is:
        # l^2 - 2*sr*l + r^2 - R^2 = 0
        #calculate discriminant
        D = sr**2 - (r.mag2() - R**2)
        #calculate l - and we need the intersection with the lowest possible l
        is_hit = (sr>0) & (D>0) & (D<=sr**2)
        l = np.where(is_hit,  sr-np.sqrt(D), r.mag())
        #get the final point
        r1 = r0+s*l
        t1 = p.T + l/speed_of_light
        p1 = Point(R=r1, T=t1, s=s)
        if(soft):
            is_hit = np.exp(self.soften_parameter*D)
        return p1, is_hit.squeeze()