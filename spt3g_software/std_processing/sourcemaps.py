from spt3g import core, mapmaker
from spt3g.coordinateutils import FlatSkyMap, MapProjection, MapCoordReference, MapPolType
from spt3g.core import G3Units as U
import time

# Make a globally accessible list of sources to mask
# Format is {source_name: [ra1, dec1, radius1], <[ra2, dec2, radius2]>}
# By default, puts a 4 arcmin mask around the main source
def get_src_mask(at_time=None):
    rcw38 = [( 8.0*U.rahour + 59.0*U.raminute + 4.8*U.rasecond,
               -(47.0*U.degrees + 30.0*U.arcmin + 36.0*U.arcsec),
               7 * core.G3Units.arcmin),
             ( 8.0*U.rahour + 59.0*U.raminute + 40.0*U.rasecond,
               -(47.0*U.degrees + 29.0*U.arcmin + 36.0*U.arcsec),
               7 * core.G3Units.arcmin)]
    mat5a = [(11 * U.rahour + 11 * U.raminute + 52.8 * U.rasecond,
              -(61 * U.degrees + 18 * U.arcmin + 47 * U.arcsec),
              7 * U.arcmin),
             (11 * U.rahour + 14 * U.raminute + 48.0 * U.rasecond,
              -(61 * U.degrees + 18 * U.arcmin + 0.0 * U.arcsec),
              10 * core.G3Units.arcmin)]
    mask_map = {
        'rcw38':rcw38,
        'mat5a':mat5a
    }
    src_map = get_src_map(at_time)
    for src, pos in src_map.items():
        #XXX: 4 arcmin is a pretty arbitrary choice
        if src not in mask_map:
            mask_map[src] = [(pos[0], pos[1], 4 * U.arcmin)]
    return mask_map
    

# Make a globally accessible source list
def get_src_map(at_time=None):
    def get_planet_at(planet, at):
        import ephem
        if at is None:
            core.log_fatal('Cannot find planets without knowing the time')
        p = getattr(ephem, planet.capitalize())()
        t = time.gmtime(float(at)/core.G3Units.s)
        t = time.strftime('%Y/%m/%d %H:%M:%S', t)
        p.compute(t)
        return ( p.ra*U.rad, p.dec*U.rad )

    rcw38 = ( 8.0*U.rahour + 59.0*U.raminute + 4.8*U.rasecond,
              -(47.0*U.degrees + 30.0*U.arcmin + 36.0*U.arcsec))
    cena =  ( 13*U.rahour + 25*U.raminute + 27.8*U.rasecond,
              -(43*U.degrees +  1*U.arcmin + 5.7*U.arcsec))
    pica =  ( 5*U.rahour +19*U.raminute + 49.7*U.rasecond,
              -(45*U.degrees +  46*U.arcmin + 43.8*U.arcsec))
    focusquasar =  ( 84.710*U.degrees, -44.086*U.degrees)
    newfocusquasar =  ( 32.6925*U.degrees, -51.0172*U.degrees)
    ra23h30 = (352.5 * U.degrees, -55 * U.degrees)
    ra0hdec575 = (0 * U.degrees, -57.5 * U.degrees)
    ra0hdec4475 = (0 * U.degrees, -44.75 * U.degrees)
    ra0hdec5225 = (0 * U.degrees, -52.25 * U.degrees)
    ra0hdec5975 = (0 * U.degrees, -59.75 * U.degrees)
    ra0hdec6725 = (0 * U.degrees, -67.25 * U.degrees)

    mat5a = (11 * U.rahour + 11 * U.raminute + 52.8 * U.rasecond,
             -(61 * U.degrees + 18 * U.arcmin + 47 * U.arcsec))

# radio day sources
    j00060623 = (  1.5579*U.degrees, -6.3931*U.degrees)
    j02105101 = ( 32.6925*U.degrees, -51.0172*U.degrees)
    j04508101 = ( 72.5227*U.degrees, -81.0173*U.degrees)
    j04532807 = ( 73.3110*U.degrees, -28.1270*U.degrees)
    j05223627 = ( 80.7416*U.degrees, -36.4586*U.degrees)
    j05384405 =  ( 84.7098*U.degrees, -44.0858*U.degrees)
    j06091542 = ( 92.4206*U.degrees, -15.71129*U.degrees)
    j07301141 = (112.5796*U.degrees, -11.6868*U.degrees)
    j09045735 = ( 136.2216*U.degrees, -57.5849*U.degrees)
    j10372934 = ( 159.3170*U.degrees, -29.5674*U.degrees)
    j10588003 = ( 164.6805*U.degrees, -80.0650*U.degrees)
    j11074449 = ( 166.7862*U.degrees, -44.8188*U.degrees)
    j11476753 = ( 176.8892*U.degrees, -67.8949*U.degrees)
    j12560547 = ( 194.0465*U.degrees, -5.7893*U.degrees)
    j13265256 = ( 201.7051*U.degrees, -52.9399*U.degrees)
    j13295608 = ( 202.2548*U.degrees, -56.1341*U.degrees)
    j13371257 = ( 204.4158*U.degrees, -12.9569*U.degrees)
    j14274206 = ( 216.9846*U.degrees, -42.1054*U.degrees)
    j15120905 = ( 228.2106*U.degrees, -9.1000*U.degrees)
    j16175848 = ( 244.3245*U.degrees, -58.8022*U.degrees)
    j16262951 = ( 246.5251*U.degrees, -29.8575*U.degrees)
    j17331304 = (263.2612*U.degrees, -13.0804*U.degrees)
    j22354835 = ( 338.8052*U.degrees, -48.5997*U.degrees)
    j22461206 = (341.5759*U.degrees, -12.1142*U.degrees)
    j22582758 = ( 344.5248*U.degrees, -27.9726*U.degrees)


    src_map = {
        'rcw38':rcw38,
        'cena':cena,
        'centaurusa':cena,
        'pica': pica,
        'pictora':pica,
        'ra23h30': ra23h30,
        'ra0hdec575': ra0hdec575,
        'ra0hdec4475': ra0hdec4475,
        'ra0hdec5225': ra0hdec5225,
        'ra0hdec5975': ra0hdec5975,
        'ra0hdec6725': ra0hdec6725,
        '0537441': focusquasar,
        'pmnj02105101': newfocusquasar,
        'mat5a': mat5a,
        'ngc 3576':mat5a,
        'ngc3576':mat5a,
        'j00060623': j00060623,
        'j02105101': j02105101,
        'j04508101': j04508101,
        'j04532807': j04532807,
        'j05223627': j05223627,
        'j05384405': j05384405,
        'j06091542': j06091542,
        'j07301141': j07301141,
        'j09045735': j09045735,
        'j10372934': j10372934,
        'j10588003': j10588003,
        'j11074449': j11074449,
        'j11476753': j11476753,
        'j12560547': j12560547,
        'j13265256': j13265256,
        'j13295608': j13295608,
        'j13371257': j13371257,
        'j14274206': j14274206,
        'j15120905': j15120905,
        'j16175848': j16175848,
        'j16262951': j16262951,
        'j17331304': j17331304,
        'j22354835': j22354835,
        'j22461206': j22461206,
        'j22582758': j22582758
   }
    if at_time is not None:
        # Add planets if we can
        try:
            for src in ['mercury', 'venus', 'mars', 'jupiter', 'saturn',
              'uranus', 'neptune', 'pluto', 'moon', 'sun']:
                src_map[src] = get_planet_at(src, at_time)
        except ImportError:
            pass
    return src_map

def get_source_ra_dec(src, at_time=None):
    src_map = get_src_map(at_time)
    src = src.lower().strip().replace(' ','').replace('-','').replace('.','')
    if src in src_map:
        return src_map[src]
    else:
        raise RuntimeError("Source %s not recognized"%src)

def get_source_mask_ra_dec_radius(src, at_time=None):
    src_map = get_src_mask(at_time)
    src = src.lower().strip().replace(' ','').replace('-','').replace('.','')
    if src in src_map:
        masks = src_map[src]
        ra = [mask[0] for mask in masks]
        dec = [mask[1] for mask in masks]
        radius = [mask[2] for mask in masks]
        return ra, dec, radius
    else:
        raise RuntimeError("Source %s not recognized"%src)

def CreateSourceMapStub(src, x_len, y_len, res, proj, at_time=None):
    ra, dec = get_source_ra_dec(src, at_time)
    m = FlatSkyMap(x_len = int(x_len), y_len = int(y_len), 
                   res = res, proj = proj,
                   alpha_center=ra,
                   delta_center=dec,
                   pol_type = MapPolType.T,
                   coord_ref = MapCoordReference.Equatorial)
    return m

def CreateSourceMap(src, data, res, proj):
    ra, dec = get_source_ra_dec(src)
    m = FlatSkyMap(data, 
                   res = res, proj = proj,
                   alpha_center=ra,
                   delta_center=dec,
                   coord_ref = MapCoordReference.Equatorial,
                   pol_type = MapPolType.T   )
    return m


def CreateGroundMap(res, proj = MapProjection.ProjCAR, 
                    delta_center = 45.0 * core.G3Units.deg, 
                    delta_angular_extent = 90 * core.G3Units.deg):
    m = FlatSkyMap(x_len = int(360 * core.G3Units.deg / res), 
                   y_len = int(delta_angular_extent/res), 
                   res = res, proj = proj,
                   alpha_center = 0.0,
                   delta_center = delta_center,
                   pol_type = MapPolType.T,
                   coord_ref = MapCoordReference.Local)
    return m


def CreateFieldMapStub(res=2 * core.G3Units.arcmin,
                       proj=MapProjection.ProjLambertAzimuthalEqualArea,
                       width=75 * core.G3Units.deg,
                       height=50 * core.G3Units.deg,
                       ra_center=0.0,
                       dec_center=-57.5 * core.G3Units.deg):
    m = FlatSkyMap(x_len=int(width / res),
                   y_len=int(height / res),
                   res=res,
                   alpha_center=ra_center,
                   delta_center=dec_center,
                   proj=proj,
                   pol_type=MapPolType.T,
                   coord_ref=MapCoordReference.Equatorial)
    return m

def CreateSourceMask(src, *args, **kwargs):
    '''
    Create a mask around the source with radius `radius`.  
    All other arguments are passed to `CreateSourceMap`.
    Returns a FlatSkyMap with ones in the pixels to be masked.
    '''
    maskmap = CreateSourceMapStub(src, *args, **kwargs)
    if 'at_time' in kwargs:
        at_time = kwargs['at_time']
    else:
        at_time = None
    tomask = get_source_mask_ra_dec_radius(src, at_time = at_time)
    mapmaker.make_point_source_mask_cpp(*tomask, False, maskmap)
    return maskmap
    
