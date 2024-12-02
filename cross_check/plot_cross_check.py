#!/bin/env python
import numpy as np
import rte
import vegas_params as vp
from tqdm.auto import tqdm
import pylab as plt
import sys
import warnings
import gvar

#read the argument name
if len(sys.argv)<2:
    print(f"Usage:`python {sys.argv[0]} <filename>.npz`")
    exit()
fname = sys.argv[1] #"Light_1_0_1_0.9.npz" 

#read file
data = np.load(fname)
r_tgt = data['r']
#select some points, we don't want to process all of them
n = 5
#take every n-th point
times = data['time'][1::n]
photons = data['photons'][:,:,1::n]
#skip 0-th order
photons = photons[1:]

#initial setup for source and target
src = rte.Source(R=vp.Vector([0,0,0]),
                 T=0,
                 s=vp.Vector([0,0,1]))

tgt = rte.Target(R=vp.Vector(data['r']),
                 T=1e-6,
                 s=vp.Vector([0,0,1]))#will be ignored

#prepare the orders of calculation
Nsteps_total = photons.shape[0]
Ns = np.arange(1,Nsteps_total)

#do the calculation
result = []
for Nsteps in tqdm(Ns):
    p = rte.Process(src, tgt, medium=rte.medium.water, Nsteps=Nsteps)
    p.vegas_kwargs = dict(nitn=10, neval=10000)
    res = p.calculate_map(override={'tgt.T':times}, output='reshape')
    result.append(res)

#make it the same shape as photons
result = np.array(result)
result = np.stack([gvar.mean(result), gvar.sdev(result)], axis=1)

#import pylab as plt
#plot the comparison
fig, axes = plt.subplots(2,int(np.ceil(len(Ns)/2)), figsize=(8,6))#, sharey=True)
axes = axes.flatten()
ax = iter(axes)

for N in Ns:
    plt.sca(next(ax))
    d_mean = photons[N-1,0,:]
    d_sdev = photons[N-1,1,:]
    plt.errorbar(x=times, y=d_mean, yerr=d_sdev, fmt='.k', label='old')
    c_mean = result[N-1,0,:]
    c_sdev = result[N-1,1,:]
    plt.errorbar(x=times, y=c_mean, yerr=c_sdev, fmt='r-', label='new')
    #plt.errorbar(x=times, y=c_mean*S_factor, yerr=c_sdev*S_factor,fmt='r:', label='new*S')
    plt.yscale('log')
    plt.grid(ls=':')
    plt.ylabel(r'$\Phi^{'+f'({N})'+'}$')
    
axes[0].legend()
plt.xlabel('Time, s')
plt.suptitle(fname.rsplit('/',1)[-1])
plt.tight_layout()
#save the figure
out_fname = f'{fname}.check.png'
plt.savefig(out_fname)
print(f"Generated figure {out_fname}. Bye!")
plt.show()