import rte
import vegas_params as vp
import numpy as np

class DetectorSpherical(vp.Expression):
    def __init__(self,
                 center = vp.Vector([0,0,0]),
                 radius = vp.Scalar(0.216),
                 T = vp.Uniform([0,1])
                ):
        super().__init__(center = vp.Vector(center),
                         radius = vp.Scalar(radius),
                         T = vp.Scalar(T),
                         _R_local = vp.Direction(),
                         _s_local = np.nan #vp.Direction()
                        )
    def __call__(self, center, radius, T, _R_local, _s_local):
        R = vp.Vector.__call__(_R_local) * vp.Scalar.__call__(radius) + vp.Vector.__call__(center)
        s = vp.Vector(_s_local)
        return rte.Point(R, T, s)

    def efficiency(self, p: rte.Point)->np.array:
        return np.ones(len(p))

    def get_intersection(self,
                         p: rte.Point,
                         speed_of_light:np.array,
                         return_is_hit_array:bool=True
                        )->rte.Point:
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
        is_hit = (D>0) & (D<=sr**2)
        l = np.where(is_hit,  sr-np.sqrt(D), r.mag())
        #get the final point
        r1 = r0+s*l
        t1 = p.T + l/speed_of_light
        p1 = rte.Point(R=r1, T=t1, s=s)
        if return_is_hit_array:
            return p1, is_hit.squeeze()
        else:
            #mark the invalid points as NaN
            t1[is_hit==False] = np.nan
            return p1