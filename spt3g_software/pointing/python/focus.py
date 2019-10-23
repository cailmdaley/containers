import numpy as np

bench_zero = {'y1': 24, 'y2': 25, 'y3': 25,
              'x4': 24, 'x5': 25,
              'z6': 33.5}

def g3bench2xyz(benchpos, gcp_offsets=None):
    '''
    Convert from G3MapDoubles containing bench position to mean x, y
    and z offset.

    NB: For most of 2018 (and onwards)  `gcp_offsets` are stored in the
    Observation frame.  For pre-2018 data, the offset is hardcoded at
    the top of this file.  The offets are unlikely to change, so it is
    probably safe to use the hard-coded values.

    Arguments
    ---------
    benchpos : dictionary
        A dictionary of optical bench positions in units of mm. Must
        contain the keys:
        'y1', 'y2', 'y3', 'x4', 'x5', 'z6'
    gcp_offsets : dictionary, optional
        A dictionary of offsets (in units of mm) that GCP applies to
        the commanded optical bench positions. Must contain the keys:
        'y1', 'y2', 'y3', 'x4', 'x5', 'z6'

    Returns
    -------
    x, y, z : tuple, 3 elements
        The mean x, y, and z offsets of the optical bench in units of
        mm.
    '''
    if gcp_offsets is None:
        gcp_offsets = bench_zero
    x = np.mean([benchpos['x4'] - gcp_offsets['x4'], 
                 benchpos['x5'] - gcp_offsets['x5']])
    y = np.mean([benchpos['y1'] - gcp_offsets['y1'], 
                 benchpos['y2'] - gcp_offsets['y2'], 
                 benchpos['y3'] - gcp_offsets['y3']])
    z = benchpos['z6'] - gcp_offsets['z6']
    return x, y, z


def bench2optical(x, y, z, offset=True):
    '''
    Convert bench offsets to offsets relative to the optical axis.
    
    Arguments
    ---------
    x, y, z: float
        The x, y and z coordinates (or offets) of the optical bench.
    offset: boolean, optional
        If `True` `x`, `y`, and `z` are offsets from nominal bench
        center.

    Returns
    -------
    vp: array-like
        An array giving the rotated position. The coordinate system
        (as viewed from the primary):
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
    Convert from an offset measured from the optical axis to optical
    bench offsets.
    
    Arguments
    ------
    v_optical: array-like
        Optical axis coordinates, as defined in `bench2optical()`

    Returns
    -------
    v_bench: array-like
        Bench coordinates:
        v_bench[0] = y1 or y2 or y3
        v_bench[1] = x4 or x5
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
    """
    Print out the optical bench offsets.

    Arguments
    ---------
    v_bench : array-like
        An array of optical bench coordinates in units of mm. The array
        should be of the same form as that returned by `optical2bench()`.
    """
    print('benchOffset {y:0.2f}, {y:0.2f}, {y:0.2f}, {x:0.2f}, {x:0.2f}, {z:0.2f}'.format(y = v_bench[0], x = v_bench[1], z = v_bench[2]))
