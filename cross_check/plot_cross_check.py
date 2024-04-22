#!/bin/env python
import numpy as np
import rte
import vegas_params as vp
from tqdm.auto import tqdm
import pylab as plt
import sys
import warnings

#read the argument name
if len(sys.argv)<2:
    print(f"Usage:`python {sys.argv[0]} <filename>.npz`")
    exit()
fname = sys.argv[1]

#calculate the factor
R = 0.21
S_factor = np.pi * R**2

#read file
data = np.load(fname)
r_tgt = data['r']
#select some points, we don't want to process all of them
n = 5 
#take every n-th point
times = data['time'][::n]
photons = data['photons'][:,:,::n]
#skip 0-th order, stop at 4th 
photons = photons[1:5]

#initial setup for source and target
src = rte.Source(R=vp.Vector([0,0,0]),
                 T=0,
                 s=vp.Vector([0,0,1]))

tgt = rte.Target(R=vp.Vector(data['r']),
                 T=1e-6,
                 s=vp.Vector([0,0,1]))#will be ignored

#prepare the orders of calculation
Nsteps_total = photons.shape[0]
Ns = np.arange(Nsteps_total)+1

#do the calculation
result = []
for Nsteps in Ns:
    result_n = []
    p = rte.Process(src, tgt, medium=rte.medium.water, Nsteps=Nsteps)
    for t in tqdm(times, desc=f'{Nsteps} order'):
        p['tgt']['T']=vp.Fixed(t)
        res = p.calculate(nitn=10, neval=10000)
        result_n.append([res.mean, res.sdev])
    result.append(result_n)
#make it the same shape as photons
result = np.array(result)
result = np.swapaxes(result, 1,2)

#plot the comparison
fig, axes = plt.subplots(2,len(Ns)//2, figsize=(8,6))#, sharey=True)
axes = axes.flatten()
ax = iter(axes)
import pylab as plt
for N in Ns:
    plt.sca(next(ax))
    d_mean = photons[N,0,:]
    d_sdev = photons[N,1,:]
    plt.errorbar(x=times, y=d_mean, yerr=d_sdev, fmt='.k', label='old')
    c_mean = np.array([v.mean for v in result[N-1]])#*S_factor
    c_sdev = np.array([v.sdev for v in result[N-1]])#*S_factor
    plt.errorbar(x=times, y=c_mean, yerr=c_sdev, fmt='r-', label='new')
    plt.errorbar(x=times, y=c_mean*S_factor, yerr=c_sdev*S_factor, 
                 fmt='r:', label='new*S')
    plt.yscale('log')
    plt.grid(ls=':')
    plt.ylabel('$\Phi^{'+f'({N})'+'}$')
    
axes[0].legend()
plt.xlabel('Time, s')
plt.suptitle(fname.rsplit('/',1)[-1])
plt.tight_layout()
#save the figure
out_fname = f'{fname}.check.png'
plt.savefig(out_fname)
print(f"Generated figure {out_fname}. Bye!")
plt.show()