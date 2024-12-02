import vegas_params as vp
import rte
import numpy as np
import pytest


v_expected = {1e-6: 3.34332961e-07, 1e-8:1.10669159e+04}

vegas_kwargs = dict(nitn=10, neval=10000)

@pytest.mark.parametrize(
    'R, T, Nsteps, expected',
    [
        ([1,0,1], 1e-6, 1, 3.34332961e-07),
        ([1,0,1], 1e-8, 1, 1.10669159e+04)
    ])
def test_fixed_time(R, T, Nsteps, expected):
    #initial setup for source and target
    src = rte.Source(R=vp.Vector([0,0,0]),
                     T=0,
                     s=vp.Vector([0,0,1]))
    
    tgt = rte.Target(R=vp.Vector(R),
                     T=T,
                     s=vp.Vector([0,0,1]))

    p = rte.Process(src, tgt, medium=rte.medium.water, Nsteps=Nsteps)
    p.vegas_kwargs = vegas_kwargs
    result = p.calculate()
    assert np.isclose(result.mean, expected)
