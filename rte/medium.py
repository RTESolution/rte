from dataclasses import dataclass
import numpy as np

class HeneyGreenshtein:
    def __init__(self, g=0.9):
        self.g = g
    
    def pdf(self, x:np.ndarray)->np.ndarray:
        "probability density function"
        g = self.g
        return 0.5*(1 - g**2) / ((1 + g**2 - 2*g*x)**1.5)
    
    def cdf(self, x:np.ndarray)->np.ndarray:
        "Cumulative distirbution function"
        g = self.g
        return 0.5*(1-g**2)/g*(1/np.sqrt(1 + g**2 - 2*g*x) - 1/(1+g))

    def icdf(self, y:np.ndarray):
        "Inverse of cdf"
        g = self.g
        return (1+g**2 - ((1-g**2)/(2*y*g+1-g))**2)/(2*g)

@dataclass
class Medium:
    mu_a:float #inverse absorption length
    mu_s:float #inverse scattering length
    c:float  #speed of light in the medium
    g:float # average cosine of scattering (asymmetry parameter)
    
    def __post_init__(self):
        self.mu_t=self.mu_a+self.mu_s
        self.h_g = HeneyGreenshtein(self.g)

    def attenuation_factor(self, time:np.ndarray)->np.ndarray:
        return np.exp(-self.mu_t*self.c*time)

    def attenuation(self, time:np.ndarray, n_scattering:int)->np.ndarray:
        return ((self.c*self.mu_s*time)**n_scattering) * self.attenuation_factor(time)

    def scatter(self, x:np.ndarray)->np.ndarray:
        return self.h_g.pdf(x)/(2*np.pi)

#define the medium parameters
water = Medium(mu_a = 0.047846, mu_s = 0.014439, c = 219662673.8, g = 0.9)