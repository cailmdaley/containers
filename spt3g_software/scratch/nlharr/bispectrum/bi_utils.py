from __future__ import division 
import numpy as np
import os, sys, itertools, copy
import bi_anal_utils as anal_utils
from spt3g import mapspectra, core, util
from spt3g.util.genericutils import uniquify_list
from functools import reduce

'''
Current status:  testing for auto spectra
                 cross spectra is not well tested

'''


DEBUG_SIMS = True
GLOBAL_SIM_DEBUG_NUM = 0



DO_MEMORY_CHECK = False
if DO_MEMORY_CHECK:
    import psutil
    def print_mem_usage(info):
        if DO_MEMORY_CHECK:
            print( info, ' Memory Usage: ', psutil.Process().memory_info().rss/1024.0**3, 'GB')
else:
    def print_mem_usage(info):
        pass

'''
debug_args = {'plot': bool,
              'verbose': bool,
              'dir':str,
              'tag':str,
              }

'''
def _debug_plot_grid(mp, dargs, spec_tag, extent = None, vmin = None, vmax  = None, use_std = False, std_scale = 1.):
    if not dargs['plot']:
        return

    import pylab as pl

    out_name = dargs['dir']+'/'+dargs['tag']+'_'+spec_tag+'.png'
    pl.clf()
    if use_std:
        vmin = -1*np.std(mp)*std_scale
        vmax = -1 * vmin
    pl.imshow(mp,  vmin = vmin, vmax = vmax, extent = extent, cmap = 'gray')
    pl.savefig(out_name)

def debug_plot_lin(y, dargs, spec_tag, x=None):
    if not dargs['plot']:
        return
    import pylab as pl
    out_name = dargs['dir']+'/'+dargs['tag']+'_'+spec_tag+'.png'
    pl.clf()
    if (x == None):
        pl.plot(y)
    else:
        pl.plot(x,y)
    pl.savefig(out_name)

def validate_debug_args(darg):
    if ((not 'plot' in darg) or
        (not 'dir' in darg) or
        (not 'verbose' in darg) or
        (not 'tag' in darg)):
        raise RuntimeError("needed key not in debug_args")
    if not type(darg['plot']) is bool:
        raise RuntimeError("in debug args plot needs to be bool")
    if not type(darg['verbose']) is bool:
        raise RuntimeError("in debug args verbose needs to be bool")
    if not type(darg['tag']) is str:
        raise RuntimeError("in debug args tag needs to be str")
    if not type(darg['dir']) is str:
        raise RuntimeError("in debug args dir needs to be str")



def get_filtered_ffts(fft_map, els, delta_el, reso_rad):
    fft_dims = np.shape(fft_map)
    el_grid = anal_utils.makeEllGrid( fft_dims, reso_rad)# + .01  uncomment to make it line up
    npts = len(els) * reduce(lambda x,y: x*y, fft_dims,1)

    out_maps = np.zeros( (len(els),)+tuple(fft_dims), dtype = 'complex' )

    for i, el in enumerate(els):
        valid_inds = np.where( np.logical_and(el_grid > el - delta_el/2., 
                                              el_grid < el + delta_el/2. ))
        out_maps[i][valid_inds] = fft_map[valid_inds]
    return out_maps


def fft_maps_numpy(mps, is_forward = True):
    if (mps.ndim == 2):
        mps = np.array([mps])
    if is_forward:
        fft_func = np.fft.fft2
    else:
        fft_func = np.fft.ifft2
    return np.array(map(fft_func, mps))

fft_maps = fft_maps_numpy

def cpp_inpaint_operation(inpaint_mask, mp, n_iters=2000):
    assert(inpaint_mask.shape == mp.shape)
    ny, nx = np.shape(inpaint_mask)
    inp_mask = core.G3VectorInt(inpaint_mask.flatten().astype('int32'))
    mp = core.G3VectorDouble(mp.flatten().astype('float64'))
    mapspectra.inpaint_map_laplace(inp_mask, nx, ny, n_iters, mp)
    mp = np.array(mp)
    return mp.reshape( (ny,nx))


def cpp_filter_map(f_mps, mp, map_scaling, els, delta_el, reso_rad, do_forward):
    n_maps = f_mps.shape[0]
    n0 = f_mps.shape[1]
    n1 = f_mps.shape[2]
    f_mps = core.G3DoubleVector(f_mps.flatten())
    mp = core.G3DoubleVector(mp.flatten())
    map_scaling = core.G3DoubleVector(map_scaling.flatten())
    els = core.G3DoubleVector(els.flatten())
    ellgrid = core.G3DoubleVector(anal_utils.makeEllGrid( np.shape(mp), reso_rad).astype('float64').flatten())
    mapspectra.bispec_filter_maps(n_maps, n0, n1, do_forward, 
                                  mp, map_scaling, ell_grid, els, delta_el,
                                  f_mps)
    return np.array(f_mps).reshape( (n_maps, n0, n1) )
    

def cpp_bispec_call(bi_in_maps, bands, els):
    
    tshape = bi_in_maps[0].shape
    n_maps = tshape[0]
    map_size = tshape[1] * tshape[2]

    n_bands = len(bi_in_maps)
    assert(n_bands > 0 and n_bands <= 3)
    assert( len(bi_in_maps) == len(bands))
    assert( len(uniquify_list(bands)) == len(bands))

    input_bands = core.G3VectorInt([0,0,0])
    m0 = core.G3VectorDouble(bi_in_maps[0].flatten())
    input_bands[0] = bands[0]
    if len(bi_in_maps) > 1):
        m1 = core.G3VectorDouble(bi_in_maps[1].flatten())
        input_bands[1] = bands[1]
    else:
        m1 = m0
    if len(bi_in_maps) > 2):
        m2 = core.G3VectorDouble(bi_in_maps[2].flatten())
        input_bands[2] = bands[2]
    else:
        m2 = m0
    els = core.G3VectorDouble(els)
    bispec_vals = core.G3VectorDouble()
    band_0 = core.G3VectorInt()
    band_1 = core.G3VectorInt()
    band_2 = core.G3VectorInt()
    el_0 = core.G3VectorInt()
    el_1 = core.G3VectorInt()
    el_2 = core.G3VectorInt()
    mapspectra.calc_full_bispec_est(m0,m1,m2,n_maps,map_size,bands,els,
                                    bispec_vals,band_0,band_1,band_2,
                                    el_0,el_1,el_2)
    bispec_vals = np.array(bispec_vals, dtype='int32')
    band_0 = np.array(band_0, dtype='int32')
    band_1 = np.array(band_1, dtype='int32')
    band_2 = np.array(band_2, dtype='int32')
    el_0 = np.array( el_0, dtype='int32')
    el_1 = np.array(el_1, dtype='int32')
    el_2 = np.array(el_2, dtype='int32')
    return bispec_vals,band_0,band_1,band_2, el_0, el_1, el_2



    
#################################
# Actual bispectrum code begins #
#################################
def con_auto_bi_list_to_grid(n_el_inds, el_inds_0, el_inds_1, el_inds_2, bispec_vals):
    out_grid = np.zeros( (n_el_inds, n_el_inds, n_el_inds))
    for i in range(len(bispec_vals)):
        for perm in itertools.permutations( (el_inds_0[i], el_inds_1[i], el_inds_2[i])):
            out_grid[ perm[0], perm[1], perm[2] ] = bispec_vals[i]
    return out_grid

#@VPU.timeit
def con_multi_bi_list_to_grids(n_el_inds, el_inds_0, el_inds_1, el_inds_2, 
                               n_bands, bands_0, bands_1, bands_2, bispec_vals,
                               scale_repeated_indices = False):
    #vs = [  ( (b0, el0), (b1,el1), (b2,el2)),
    #            ...

    vs = map(sorted, zip(zip(bands_0, el_inds_0),
                         zip(bands_1, el_inds_1),
                         zip(bands_2, el_inds_2)))

    bands = map(sorted, zip(bands_0, bands_1, bands_2))

    if scale_repeated_indices:
        for i, v in enumerate(vs):
            if (v[0]==v[1] and v[1]==v[2]):
                bispec_vals[i] *= (1.0/6.0)
            elif (v[0]==v[1] or v[1]==v[2]):
                bispec_vals[i] *= 1.0/2.0
        
    unique_bands = sorted(uniquify_list(bands))

    out_grid = np.zeros((len(unique_bands), n_el_inds, n_el_inds, n_el_inds))
    grid_indices = map(lambda x: unique_bands.index(x), bands)

    AUTOBI = 0
    FIRSTTWOSHARED = 1
    SECTWOSHARED = 2
    ALLDIFF = 3

    iter_type = []
    for ub in unique_bands:
        if ub[0] == ub[1] and ub[1] == ub[2]:
            iter_type.append(AUTOBI)
        elif ub[0] == ub[1]:
            iter_type.append(FIRSTTWOSHARED)
        elif ub[1] == ub[2]:
            iter_type.append(SECTWOSHARED)
        else:
            iter_type.append(ALLDIFF)

    for i,v  in enumerate(vs):
        grid_ind = grid_indices[i]
        it = iter_type[grid_ind]
        if it == AUTOBI:
            for perm in itertools.permutations( (v[0][1], v[1][1], v[2][1])):
                out_grid[grid_ind, perm[0], perm[1], perm[2] ] = bispec_vals[i]
        elif it == FIRSTTWOSHARED:
            for perm in itertools.permutations( (v[0][1], v[1][1])):
                out_grid[grid_ind, perm[0], perm[1], v[2][1] ] = bispec_vals[i]
        elif it == SECTWOSHARED:
            for perm in itertools.permutations( (v[2][1], v[1][1])):
                out_grid[grid_ind, v[0][1], perm[0], perm[1] ] = bispec_vals[i]
        elif it == ALLDIFF:
            out_grid[grid_ind, v[0][1], v[1][1], v[2][1]] = bispec_vals[i]
    return out_grid, unique_bands



def expandify_map(m_in, new_shape):
    in_shape = np.shape(m_in)
    m_out = np.zeros(new_shape)
    xlim = min(in_shape[1], new_shape[1])
    ylim = min(in_shape[0], new_shape[0])
    m_out[ :ylim, :xlim] = m_in[ :ylim, :xlim  ]
    return m_out








def process_map_for_bs_nofft(mp, reso_rad, apod_mask, tf_grid, tf_thresh, 
                             cmb_ells, cmb_cls, psd_grid, 
                             inpaint_mask, dargs, kmask, inpaint_n_iters,
                             is_sim = False):
    '''
    Should not modify any of the input values now.  Frack.
    '''

    mp = copy.copy(mp)
    kmask = copy.copy(kmask)
    #sets the kmask if we don't supply it
    if kmask == None:
        kmask = np.zeros(np.shape(mp)) + 1.0

    #first let's inpaint the map
    if dargs['verbose']:
        print( 'inpainting')
    #_debug_plot_grid(mp, dargs, 'pre_inpaint', vmin=-100, vmax = 100)
    mp = cpp_inpaint_call(inpaint_mask, mp, n_iters=inpaint_n_iters)
    #_debug_plot_grid(mp, dargs, 'post_inpaint', vmin = -100, vmax = 100)
    mp *= apod_mask


    if dargs['verbose']:
        print( 'done inpainting')

    #assigns a bunch of constants we will later use
    ngrid = float(np.shape(mp)[0])

    dngrid = float(ngrid)
    maskfac_bl = np.sum(apod_mask**3.0)/ngrid**2.0
    maskfac_cl = np.sum(apod_mask**2.0)/ngrid**2.0
    nsr = (reso_rad * ngrid)**2.0


    if dargs['verbose']:
        print( 'cor fac', maskfac_bl*nsr )

    #totvar, kmask, itfgrid = get_totvar_and_kmask(reso_rad, cmb_ells, cmb_cls, psd_grid, tf_grid, tf_thresh)

    ellgrid = anal_utils.makeEllGrid( np.shape(psd_grid), reso_rad)
    #gets the inverse of the transfer function
    bad_thresh_inds = np.where( tf_grid < tf_thresh)
    itfgrid = tf_grid.copy()
    itfgrid[bad_thresh_inds] = 1e6 #don't worry about this being big we filter it later
    itfgrid = 1. / itfgrid
    itfgrid[bad_thresh_inds] = 0


    #itfgrid[bad_thresh_inds] = 0
    kmask[bad_thresh_inds] = 0

    if not np.all( (cmb_ells[1:]-cmb_ells[:-1]) == 1):
        raise RuntimeError("This code assumes that the cmb_ells have a unit step between them")
    
    cmb_offset_ind = cmb_ells[0]    
    #grid values of the cmb_grid
    cmb_inds = (ellgrid+cmb_offset_ind).astype('int')
    bad_inds = np.where(cmb_inds >= len(cmb_cls))
    cmb_inds[bad_inds] = -1
    cmb_cl_grid = cmb_cls[cmb_inds]
    
    #includes the transfer function with the noise
    totvar = (psd_grid*itfgrid)**2.0 + cmb_cl_grid
    # correct totvar for the fact that input clcmb and noise psd are
    # assumed to be referred to full sky
    totvar *= maskfac_cl * nsr
    
    #invert the values in the variance grid to get the weight grid
    weight_grid = anal_utils.safe_inverse(totvar)
    weight_grid *= kmask #filter out the funky portions from the transfer function inversion
    ###################################################
    #warning, mp and fft_mp should point to same object
    ###################################################
    
    if dargs['verbose']:
        print( 'not ffting map')
    if is_sim:
        scale_fac = weight_grid * reso_rad**2.0
    else:
        scale_fac = weight_grid * itfgrid * reso_rad**2.0
    return mp, scale_fac, weight_grid, maskfac_bl, nsr, ngrid


def float_arr_eq(a,b,eps=1e-4):
    return np.all(np.abs(a-b)<eps)

def get_bispectrum_multi(els, delta_el, reso_arcmin, 
                         mps, bands, apod_masks, tf_grids, tf_threshs, cmb_ells, cmb_cls, 
                         psd_grids, inpaint_masks,  #map dependent variables
                         dargs, kmask = None, no_area_corr = False, inpaint_n_iters = 2000,
                         is_sim = False, skip_inpaint = False, cached_weight = None):

    #check for basic validity and that the apodization mask is the same for all maps
    validate_debug_args(dargs)

    if skip_inpaint:
        inpaint_masks = map(lambda x: 0*x, inpaint_masks)

    len_x = len(mps)
    shp  = np.shape(mps[0])
    #check that shape of all the maps are same and all the lists have the same length
    for maps in [mps, apod_masks, tf_grids, psd_grids, inpaint_masks]:
        assert(len(maps) == len_x)
        for m in maps:
            assert(np.shape(m) == shp)

    assert( len(tf_threshs)==len(cmb_ells)==len(cmb_cls))

    #for now check that the maps and apod_masks are same size
    if len(apod_masks)>=2:
        assert(float_arr_eq(apod_masks[0], apod_masks[1]))
    if len(apod_masks)==3:
        assert(float_arr_eq(apod_masks[1], apod_masks[2]))

    ######
    ###### And now some code
    ######

    reso_rad = reso_arcmin/60. * np.pi/180.

    apod_mps = []
    scale_facs = []
    weight_grids = []
    #maskfac_bl, nsr, ngrid should be the same.

    #print_mem_usage('Get SF')



    #_debug_plot_grid(mps[0], dargs, 'pre_inpaint', vmin = -100, vmax = 100)
    #_debug_plot_grid(inpaint_masks[0], dargs, 'inpaint_mask', vmin = 0, vmax = 1)

    for i in range(len_x):
        apdmp, scale_fac, weight_grid, maskfac_bl, nsr, ngrid = process_map_for_bs_nofft(
            mps[i], reso_rad, apod_masks[i], tf_grids[i], tf_threshs[i],
            cmb_ells[i], cmb_cls[i], psd_grids[i], inpaint_masks[i], 
            dargs, kmask, inpaint_n_iters,
            is_sim = is_sim)

        apod_mps.append(apdmp)
        scale_facs.append(scale_fac)
        weight_grids.append(weight_grid)

    if DEBUG_SIMS:
        import pickle
        global GLOBAL_SIM_DEBUG_NUM
        for i in range(len_x):
            print( np.shape(apod_mps[i]) )
            pickle.dump(apod_mps[i], open('sim_debug_map_%d_%d.pkl'%(i, GLOBAL_SIM_DEBUG_NUM), 'w'), protocol=2)
            #pickle.dump(scale_facs[i], open('sim_scalefacs_map_%d_%d.pkl'%(i, GLOBAL_SIM_DEBUG_NUM), 'w'), protocol=2)
        GLOBAL_SIM_DEBUG_NUM += 1


    #fmp, els, delta_el, reso_rad, fmps, els, ngrid, 

    #_debug_plot_grid(apod_mps[0], dargs, 'post_inpaint', vmin = -100, vmax = 100)

    if dargs['verbose']:
        print( 'removing kmask' )
    out_vals = {}

    fmps_shape = (len(els), np.shape(mps[0])[0], np.shape(mps[0])[1])

    for amp_lst, tag in zip([apod_mps, weight_grids], ['bi', 'wt']):
        if tag == 'wt' and cached_weight != None:
            if dargs['verbose']: print( 'using cached weight' )
            out_vals['wt'] = cached_weight
            continue    
        if tag == 'bi':
            do_forward = True
        else:
            do_forward = False

        bi_in_maps = []

        print_mem_usage('Making fmps')
        if dargs['verbose']: print( 'filtering maps' )
        for amp, sf in zip(amp_lst, scale_facs):
            fmps = np.zeros(fmps_shape, dtype = 'float64')
            print_mem_usage('one fmps')
            fmps = cpp_filter_map(fmps, amp, sf, els, delta_el, reso_rad, do_forward)
            print_mem_usage('filtering')
            bi_in_maps.append(fmps)

        print_mem_usage('Done making fmps')
        if dargs['verbose']: print( 'bispecing maps' )
        bispec_vals, band_0, band_1, band_2, el_0, el_1, el_2 = cpp_bispec_call(
            ms = bi_in_maps, bands = bands,  els = els) 
        
        bispec_vals /=  ngrid**2.0 
        if dargs['verbose']: print( 'converting to grid format' )
        bi_grid, grid_bands = con_multi_bi_list_to_grids(len(els), el_0, el_1, el_2, 
                                                         len(bands), band_0, band_1, band_2,
                                                         bispec_vals,
                                                         scale_repeated_indices = True)
        out_vals[tag] = bi_grid
        if tag == 'bi':
            out_vals[tag+'_grid_bands'] = grid_bands
    returned_cached_weight = copy.copy(out_vals['wt'])

    out_vals['bi'] *= anal_utils.safe_inverse(out_vals['wt'])
    out_vals['wt'] *= maskfac_bl**1.5 
    if no_area_corr:
        out_vals['bi'] *= (maskfac_bl*nsr)
        out_vals['wt'] /= ((maskfac_bl*nsr))**2.0

    print_mem_usage('Done with BS')
    return out_vals, returned_cached_weight


def validate_bispec_grid(grid, nels):
    if grid.ndim != 3:
        raise RuntimeError("Improperly formatted bispec grid, must be np array with 3 dims")
    if (np.shape(grid)[0] != nels or
        np.shape(grid)[1] != nels or
        np.shape(grid)[2] != nels):
        raise RuntimeError("Improperly formatted bispec grid, 3 dims must be same size")


def generate_noise_map_pair(noise_psd, reso_arcmin, tf_grid, tf_thresh):
    d2r = 0.01745329251
    pixel_size = d2r*reso_arcmin/60.0
    shape = np.shape(noise_psd)
    n1,n2 = shape
    assert n1 == n2
    npix = n1*n2

    bad_thresh_inds = np.where( tf_grid < tf_thresh)
    itfgrid = tf_grid.copy()
    itfgrid[bad_thresh_inds] = 1e6 #don't worry about this being big we filter it later
    itfgrid = 1. / itfgrid
    itfgrid[bad_thresh_inds] = 0


    if DEBUG_SIMS:
        import pickle
        pickle.dump(noise_psd, open('sim_noise_powspec.pkl', 'w'), protocol=2)

    uf = np.random.normal(size=npix).reshape(shape)+1.0j*np.random.normal(size=npix).reshape(shape)
    tmap = noise_psd * uf  * itfgrid

    #import pdb; pdb.set_trace()
    tmap  = np.fft.ifft2( tmap )*npix**.5 / pixel_size
    return np.real(tmap), np.imag(tmap)


'''
returns fast,slow
'''
def get_cr_sim_in_map_dims(fn):
    sn = os.path.basename(fn).split('_')
    return int(sn[5]),int(sn[6])

def load_sim_t_in_map(fn, sim_num, slow_n,fast_n):
    offset = 4 * fast_n * slow_n * sim_num
    f = open(fn)
    f.seek(offset)
    a = np.fromfile(f, dtype='float32', count = fast_n * slow_n)
    return a.reshape(slow_n, fast_n) * 1e6

def load_cr_sim_t_map(fn, sim_num):
    nf, ns = get_cr_sim_in_map_dims(fn)
    return load_sim_t_in_map(fn, sim_num, ns,nf)


    

#noise file
#bispec file


def safe_inverse(arr):
    arr = arr.copy()
    bad_inds = np.where(arr == 0)
    arr[bad_inds] = 1
    arr = 1./arr
    arr[bad_inds] = 0
    return arr


def makeFFTGrid(shape, resolution=1.):
    """
    Creates a grid in which each element is proportional to the
    distance from the point to the origin (with wrapping).
    The distance between neighboring grid squares is 1/(n*resolution).
    
    >>> utils.math.makeFFTGrid(10, 2.)
    array([ 0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25, -0.2 , -0.15, -0.1 , -0.05])
    
    >>> utils.math.makeFFTGrid([5,5], 1.)
    array([[ 0.        ,  0.2       ,  0.4       ,  0.4       ,  0.2       ],
           [ 0.2       ,  0.28284271,  0.4472136 ,  0.4472136 ,  0.28284271],
           [ 0.4       ,  0.4472136 ,  0.56568542,  0.56568542,  0.4472136 ],
           [ 0.4       ,  0.4472136 ,  0.56568542,  0.56568542,  0.4472136 ],
           [ 0.2       ,  0.28284271,  0.4472136 ,  0.4472136 ,  0.28284271]])
    """
    # Make sure the input "shape" is an iterable.
    try:
        len(shape)
    except TypeError:
        shape = [shape]

    # Create a set of coordinates. We can then find the "distance" of
    # each point from the origin.
    slices = [slice( -np.floor((dim-1)/2.), np.ceil((dim-1)/2.)+1 ) for dim in shape]
    coordinates = np.mgrid[slices]

    # Get the correct normalization - divide each dimension by
    # its size and by the overall resolution.
    for index, dim_size in enumerate(shape):
        coordinates[index] /= (dim_size*resolution)

    # If we want a 1D array, then we already have what we want.
    # Don't apply the square and square-root, so that we can
    # keep the negative signs.
    if len(shape)!=1:
        grid = np.sqrt( np.sum( coordinates**2, axis=0 ) )
    else:
        grid = coordinates[0]
    
    for index, dim_size in enumerate(shape):
        grid = np.roll(grid, int(np.ceil((dim_size-1)/2.)+1), axis=index)

    return grid

def makeEllGrid(shape, resolution=1.):
    """
    OUTPUT 
       2.*np.pi*makeFFTGrid(shape=shape, resolution=resolution)
    """
    return 2.*np.pi*makeFFTGrid(shape=shape, resolution=resolution)


def combine_bispecs_and_variances(bispec_vals, variances):
    wgts = map(get_bi_field_weight, variances)
    wsum = sum(wgts)
    wgts = map(lambda w: w/wsum, wgts)
    #weighted mean
    wmean = sum(map(lambda i: wgts[i] * bispec_vals[i], range(len(bispec_vals))))

    #weighted variance
    wvar = sum(map(lambda i: wgts[i]**2 * variances[i], range(len(bispec_vals))))
    return wmean, wvar



def get_bispec_chisq_diff(b0, b_var0, b1, b_var1):
    bad_inds = np.where(b_var0 == 0)
    n_vals = np.size(b0) - np.size(bad_inds)/3

    sq_diff = (b0-b1)**2/(b_var0+b_var1)
    sq_diff[bad_inds] = 0

    #import pdb; pdb.set_trace()

    chisq = np.sum(sq_diff)
    red_chisq = chisq/n_vals
    return chisq, red_chisq



def get_bi_field_weight(var_matrix):
    good_inds = np.where(var_matrix != 0)
    return np.median( 1.0/(var_matrix[good_inds]))


