from spt3g import core
import numpy as np
from astropy.io import fits

def bl2fits(bl_array, fitsfile):
    '''
    Python implementation of the IDL routine bl2fits.

    bl_array
      real or double array of Bl coefficients to be written to file. 
      This has dimension (n, lmax+1) with n == 1 or n == 3 given in the sequence T E B.
    fitsfile
      String containing the name of the file to be written.
    '''
    
    ndims = len(np.shape(bl_array))
    assert(ndims == 1 or ndims == 2)
    if ndims == 2:
        assert(np.shape(bl_array)[0] == 1 or 
               np.shape(bl_array)[0] == 3)
    is_pol = (ndims == 2 and np.shape(bl_array)[0] == 3)
    
    hdu_lst = fits.HDUList()
    beam_header = fits.Header()
    hdu_lst.append(fits.PrimaryHDU(header = beam_header))

    if is_pol:
        tcol = fits.Column(name='TEMPERATURE', array=bl_array[0], format='K')
        ecol = fits.Column(name='GRADIENT', array=bl_array[1], format='K')
        bcol = fits.Column(name='CURL', array=bl_array[2], format='K')
        thdu = fits.TableHDU.from_columns([tcol, ecol, bcol])
    else:
        t_values = bl_array if ndims == 1 else bl_array[0]
        tcol = fits.Column(name='TEMPERATURE', array=np.array(t_values), format='K')
        thdu = fits.TableHDU.from_columns([tcol])

    thdu.header['EXTNAME'] = 'WINDOW FUNCTION'
    thdu.header['MAX-LPOL'] = np.shape(bl_array)[-1] - 1
    thdu.header['POLAR'] = is_pol
    thdu.header['BCROSS'] = False
    thdu.header['ASYMCL'] = False

    hdu_lst.append(thdu)
    hdu_lst.writeto(fitsfile)

def bl_for_lenspix(bl_array, outfile):
    '''
    Generates a beamfile for lenspix simulations
    bl_array: a 1d array of length lmax+1
    outfile: filename that this is stored in
    '''

    assert(len(np.shape(bl_array == 1)))
    f = open(outfile, 'w')
    for i, l in enumerate(bl_array):
        f.write( '%d %f\n'%(i,l))
    f.close()

