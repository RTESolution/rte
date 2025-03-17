import rte
import vegas_params as vp
import numpy as np

class DetectorSpherical(vp.Expression):
    def __init__(self, position, radius, T):
        self.R = radius
        super().__init__(position=position, time=T)

    def efficiency(p: rte.Point)->np.array:
        return np.ones(len(p))

    def get_intersection(self,
                         p: rte.Point, 
                         speed_of_light: np.array)->rte.Point:

        r0, s = p.R, p.s
        R = self.R
        r = self['position'].value-r0
        sr = s.dot(r)
        #solving equation for the intersection point
        #equation is:
        # l^2 - 2*sr*l + r^2 - R^2 = 0
        #calculate discriminant
        D = sr**2 - (r.mag2() - R**2)
        #calculate l - and we need the intersection with the lowest possible l
        valid = (D>0)& (D<=sr**2)
        l = np.where(valid,  sr-np.sqrt(D), r.mag())
        #get the final point
        r1 = r0+s*l
        t1 = p.T + l/speed_of_light
        #mark the invalid points as NaN
        t1[valid==False] = np.nan
        return rte.Point(R=r1, T=t1, s=s)