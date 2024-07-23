import numpy as np
import vegas_params as vp
from vegas_params.vector import vector
from scipy.spatial.transform import Rotation

def combine_rotations(dir_base:vector, dir_relative:vector)->vector:
    rot = align_rotation_to_vector(dir_base)
    return rot.apply(dir_relative.reshape(-1,3)).view(vector)

def align_rotation_to_vector(v:vector)->Rotation:
    v = np.asarray(v).view(vector)
    n = v/v.mag()
    theta = np.arccos(n.z)
    phi = np.arctan2(n.y,n.x)
    return Rotation.from_euler('yz',np.array([theta,phi]).T.squeeze())
