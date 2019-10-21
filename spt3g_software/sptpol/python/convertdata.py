from spt3g import core, coordinateutils
from sptpol_software.observation import sky

import numpy as np

def convert_map_frame_to_sptpol_object(fr):
    '''
    Takes a single map frame and converts it to an sptpol-style Map or PolarizedMap object.
    '''

#    map_shape_npix = fr['T'].shape
#
# actually, just in case, don't assume any attributes are there. but
# collect some warnings to throw if they're not
    warnstr = '\n'
    nwarn = 0
    try:
        map_shape_npix = np.shape(fr['T'])
    except: 
        try:
            map_shape_npix = np.shape(fr['Q'])
        except:
            try:
                map_shape_npix = np.shape(fr['U'])
            except:
                core.log_error('Bad Map frame, has no T, Q, or U information.')
                map_shape_npix = (1,1)
                out_object= sky.Map(shape=map_shape_npix)
                return out_object
    try:
        alpha_center = fr['T'].alpha_center
    except:
        alpha_center = 0.
        warnstr += ' No value found for alpha_center, setting to 0.\n'
        nwarn += 1
    try:
        delta_center = fr['T'].delta_center
    except:
        delta_center = 0.
        warnstr += ' No value found for delta_center, setting to 0.\n'
        nwarn += 1
    try:
        proj = fr['T'].proj
    except:
        proj = coordinateutils.MapProjection.ProjNone
        warnstr += ' No value found for proj, setting to None (integer value 42).\n'
        nwarn += 1
    projint = [key for key in (coordinateutils.MapProjection.values).items() if key[1] == proj][0][0]
    try:
        res = fr['T'].res
    except:
        res = 0.
        warnstr += ' No value found for res, setting to 0.\n'
        nwarn += 1

# is this a polarized or unpolarized map?
    if (('T' in fr and 'Q' in fr) or ('T' in fr and 'U' in fr) or 
        ('Q' in fr and 'U' in fr)):
        try: 
            is_weighted = fr['T'].is_weighted
        except:
            is_weighted = False
            warnstr += ' No value found for is_weighted, setting to False.\n'
            nwarn += 1
        out_object= sky.PolarizedMap(shape=map_shape_npix,
                                     data_type=np.float64,
                                     projection=projint, 
                                     center=[alpha_center,delta_center],
                                     reso_arcmin=res/core.G3Units.arcmin,
                                     weighted_map=is_weighted)
        try: 
            out_object['T'] = np.asarray(fr['T'])
        except:
            out_object['T'] = np.zeros(map_shape_npix)
            warnstr += ' No values found for T map, setting to 0.\n'
            nwarn += 1
        try: 
            out_object['Q'] = np.asarray(fr['Q'])
        except:
            out_object['Q'] = np.zeros(map_shape_npix)
            warnstr += ' No values found for Q map, setting to 0.\n'
            nwarn += 1
        try: 
            out_object['U'] = np.asarray(fr['U'])
        except:
            out_object['U'] = np.zeros(map_shape_npix)
            warnstr += ' No values found for U map, setting to 0.\n'
            nwarn += 1
        wtemp2 = np.zeros(map_shape_npix[0],map_shape_npix[1],3,3)
        try:
            wtemp = fr['Wpol']
            wtemp2[:,:,0,0] = wtemp.TT
            wtemp2[:,:,0,1] = wtemp.TQ
            wtemp2[:,:,0,2] = wtemp.TU
            wtemp2[:,:,1,1] = wtemp.QQ
            wtemp2[:,:,1,2] = wtemp.QU
            wtemp2[:,:,2,2] = wtemp.UU
        except:
            try: 
                wtemp = fr['Wunpol']
                wtemp2[:,:,0,0] = wtemp.TT
                wtemp2[:,:,0,1] = wtemp.TQ
                wtemp2[:,:,0,2] = wtemp.TU
                wtemp2[:,:,1,1] = wtemp.QQ
                wtemp2[:,:,1,2] = wtemp.QU
                wtemp2[:,:,2,2] = wtemp.UU
                warnstr += ' Map frame has multiple polarizations, but weights are stored as Wunpol.\n'
                nwarn += 1
            except:
                warnstr += ' No weight values found, setting to 0.\n'
                nwarn += 1
        wtemp2[:,:,1,0] = wtemp2[:,:,0,1]
        wtemp2[:,:,2,0] = wtemp2[:,:,0,2]
        wtemp2[:,:,2,1] = wtemp2[:,:,1,2]
    else:
        if 'T' in fr:
            try: 
                is_weighted = fr['T'].is_weighted
            except:
                is_weighted = False
                warnstr += ' No value found for is_weighted, setting to False.\n'
                nwarn += 1
            out_object= sky.Map(shape=map_shape_npix,
                                data_type=np.float64,
                                projection=projint, 
                                center=[alpha_center,delta_center],
                                reso_arcmin=res/core.G3Units.arcmin,
                                polarization='T',
                                weighted_map=is_weighted)
            out_object.map = np.asarray(fr['T'])
            try:
                out_object.weight = np.asarray((fr['Wunpol']).TT)
            except:
                try:
                    out_object.weight = np.asarray((fr['Wpol']).TT)
                    warnstr += ' Map frame only contains T information but has polarized weights. Using TT component.\n'
                    nwarn += 1
                except:
                    out_object.weight = np.zeros(map_shape_npix)
                    warnstr += ' No weight values found, setting to 0.\n'
                    nwarn += 1
        if 'Q' in fr:
            try: 
                is_weighted = fr['Q'].is_weighted
            except:
                is_weighted = False
                warnstr += ' No value found for is_weighted, setting to False.\n'
                nwarn += 1
            out_object= sky.Map(shape=map_shape_npix,
                                data_type=np.float64,
                                projection=projint, 
                                center=[alpha_center,delta_center],
                                reso_arcmin=res/core.G3Units.arcmin,
                                polarization='Q',
                                weighted_map=is_weighted)
            out_object.map = np.asarray(fr['Q'])
            try:
                wtemp = fr['Wunpol']
                if np.max(wtemp.QQ) == 0 and np.max(wtemp.TT) > 0:
                    out_object.weight = np.asarray(wtemp.TT)
                    warnstr += ' Map frame only contains Q information but weights are stored as TT.\n'
                    nwarn += 1
                else:
                    out_object.weight = np.asarray(wtemp.QQ)
            except:
                try:
                    out_object.weight = np.asarray((fr['Wpol']).QQ)
                    warnstr += ' Map frame only contains Q information but has polarized weights. Using QQ component.\n'
                    nwarn += 1
                except:
                    out_object.weight = np.zeros(map_shape_npix)
                    warnstr += ' No weight values found, setting to 0.\n'
                    nwarn += 1
        if 'U' in fr:
            try: 
                is_weighted = fr['U'].is_weighted
            except:
                is_weighted = False
                warnstr += ' No value found for is_weighted, setting to False.\n'
                nwarn += 1
            out_object= sky.Map(shape=map_shape_npix,
                                data_type=np.float64,
                                projection=projint, 
                                center=[alpha_center,delta_center],
                                reso_arcmin=res/core.G3Units.arcmin,
                                polarization='U',
                                weighted_map=is_weighted)
            out_object.map = np.asarray(fr['U'])
            try:
                wtemp = fr['Wunpol']
                if np.max(wtemp.QQ) == 0 and np.max(wtemp.TT) > 0:
                    out_object.weight = np.asarray(wtemp.TT)
                    warnstr += ' Map frame only contains Q information but weights are stored as TT.\n'
                    nwarn += 1
                else:
                    out_object.weight = np.asarray(wtemp.QQ)
            except:
                try:
                    out_object.weight = np.asarray((fr['Wpol']).QQ)
                    warnstr += ' Map frame only contains Q information but has polarized weights. Using QQ component.\n'
                    nwarn += 1
                except:
                    out_object.weight = np.zeros(map_shape_npix)
                    warnstr += ' No weight values found, setting to 0.\n'
                    nwarn += 1

    if nwarn > 0:
        core.log_warn(warnstr)

    return out_object


