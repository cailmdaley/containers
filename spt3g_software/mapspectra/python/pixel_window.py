import numpy as np
'''
Has functions for getting the pixel window function for flat sky

The healpy function healpy.sphtfunc.pixwin returns the pixel window function for
healpix skies.
'''
def calc_pl(ell, reso_arcmin ):
    '''
    returns the flat sky pixel window function.  Was stolen from idl code which
    was stolen from a paper which was stolen from the math gods.  They punished
    the titan that brought us this knowledge by having birds eat his liver every
    day.  We should remember to include the pixel window function so his
    sacrifice was worth it.
    '''

    pix = reso_arcmin/60. * np.pi/180.
    return np.sqrt( ( np.exp(-1./18.1*(ell*pix)**2.04 ) * ( 1. - .0272*(ell*pix)**2 ) )  )

def pixwin_flatsky(ells, reso_arcmin, npts=100):
    '''
    This is a port of pixwin_flatsky.pro.  Note, however, this outputs the pixel
    window, NOT the pixel window squared as the IDL version does.  I believe
    calc_pl above is an approximation at low ell for this function, while this
    function should be exact.
    '''

    nl = len(ells)

    res_factor = reso_arcmin/60. * np.pi/180.

    phi = (np.arange(npts, dtype=np.float64)/npts*2. - 1. + 1./npts)*np.pi

    pixwin = np.zeros(nl)

    if int(ells[0]) == 0:
        pixwin[0] = res_factor**4.
        start_index = 1
    else:
        start_index = 0

    for i in range(start_index, nl):
        this_ell = ells[i]

        u = this_ell/2./np.pi*np.cos(phi)
        v = this_ell/2./np.pi*np.sin(phi)

        wtemp = np.sin(np.pi*res_factor*u) * np.sin(np.pi*res_factor*v) / (np.pi**2.*u*v)

        pixwin[i] = np.mean(wtemp**2.)

    pixwin /= res_factor**4.

    return np.sqrt(pixwin)

