import numpy as np

def flat_map_to_png(flat_map, file_name ):
    import scipy.misc
    '''
    Convert a FlatSkyMap  into a png with a 1 to 1 mapping of the pixels
    '''
    scipy.misc.imsave(file_name, flat_map)

def mask_array_from_png(file_name):
    import scipy.misc
    '''
    Convert a png into a 2d array
    '''
    arr = 1.0 - np.asarray((np.asarray(scipy.misc.imread(file_name), dtype = 'float')/255.0)[:,:,3] == 1, dtype = 'float')
    return arr
    
