import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np

def normalize(a, axis=None, v_min=None, v_max=None):
    v_min = v_min or a.min(axis=axis)
    v_max = v_max or a.max(axis=axis)
    if v_min==v_max:
        return np.ones_like(a)
    return (a - v_min)/(v_max-v_min)

def plot_trajectory(points, factor):
    view = gl.GLViewWidget()
    #grids
    #view.addItem(grid:=gl.GLGridItem())
    #axes
    view.addItem(zaxis:=gl.GLGridItem())
    zaxis.setSize(0,0,10)
    #colormap
    points.T = normalize(points.T, axis=None).squeeze()
    cm = pg.colormap.get('CET-D8')
    for pt, f in zip(points.iter(),normalize(factor)):
        width = 1+f*3
        alpha = f**0.3
        #logger.debug(p.R.shape)
        if width<0.75 or alpha<0.1:
            continue
        color = cm.map(1-pt.T, mode='Float')
        color[:,3]=alpha#set the alpha
        
        line=gl.GLLinePlotItem(pos=pt.R, 
                               color=color,
                               width=width,
                               antialias=False)
        view.addItem(line)
    view.show()

def plot_photons(points, photon_size=1):
    view = gl.GLViewWidget()
    #grids
    #view.addItem(grid:=gl.GLGridItem())
    #axes
    view.addItem(zaxis:=gl.GLGridItem())
    zaxis.setSize(0,0,10)
    #colormap
    alpha=0.7
    width=2
    factor=1
    points.T = normalize(points.T, axis=None).squeeze()
    cm = pg.colormap.get('CET-D8')
    for pt in points.iter():
        color = cm.map([1-pt.T]*2, mode='Float')
        #color[:,3]=alpha#set the alpha
        line=gl.GLLinePlotItem(pos=[pt.R, pt.R+pt.s*photon_size], 
                               color=color,
                               width=width,
                               antialias=False)
        #plot all the photon directions
        view.addItem(line)
    view.show()