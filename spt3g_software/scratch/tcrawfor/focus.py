import __future__
import numpy as np

def bench2optical(x, y, z):
    '''
    Convert bench offsets to offsets relative to the optical axis.
    
    INPUTS
    ------
    x, y, z: float
        The x, y and z coordinates (or offets) of the optical bench.
    RETURNS
    -------
    vp: array-like
        An array giving the rotated position. 
        The coordinate system (as viewed from the primary):
        vp[0]    left (+) or right (-)
        vp[1]    up (+)   or down (-)
        vp[2]    away (+) or towards (-)
    '''
    # angle between the bench coordinate system and optical axis
    theta = np.deg2rad(41.47)
    m = np.matrix([[1., 0., 0.],
                  [0., np.sin(theta), -np.cos(theta)],
                  [0., np.cos(theta), np.sin(theta)]])
    v = np.matrix([z, x, y]).transpose()
    return m * v
    
def optical2bench(v_optical):
    '''
    Convert from an offset measured from the optical axis to bench offsets
    
    INPUTS
    ------
    v_optical: array-like
        Optical axis coordinates, as defined in bench2optical
    RETURNS
    -------
    v_bench: array-like
        Bench coordinates:
        v_bench[0] = y1/y2/y3
        v_bench[1] = x4/x5
        v_bench[2] = z6
    '''
    # angle between the bench coordinate system and optical axis
    theta = np.deg2rad(41.47)
    m = np.matrix([[1., 0., 0.],
                  [0., np.sin(theta), -np.cos(theta)],
                  [0., np.cos(theta), np.sin(theta)]]).I # not inverse
    v_optical = np.matrix(v_optical)
    if np.shape(v_optical) == (1, 3):
        v_optical = v_optical.transpose()
    v_bench = m * v_optical
    return np.array([v_bench[2], v_bench[1], v_bench[0]]).flatten()

def bench_command(v_bench):
    print('benchOffset {y:0.2f}, {y:0.2f}, {y:0.2f}, {x:0.2f}, {x:0.2f}, {z:0.2f}'.format(y = v_bench[0], x = v_bench[1], z = v_bench[2]))
