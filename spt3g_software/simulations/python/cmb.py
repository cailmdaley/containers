'''
simulations/python/cmb.py

Read in CMB power spectra from CAMB outputs
'''
import os
import numpy as np


def read_camb(filename, as_cls=False, lmax=np.inf, lmin=None):
    """
    Read CAMB .dat file into a dictionary.

    Parameters
    ----------
    filename : str
        The .dat file to be read
    as_cls : bool
        If True, converts data to Cls
    lmax : int, optional
        Omit any values for higher l than this.
    lmin : optional
        The minimum l to include. If this is lower than
        found in the file, pad spectra with zeros.

    Returns
    -------
    camb_cls : dict
        Dictionary of labelled el and Dls (or Cls)

    Notes
    -----
    By default, CAMB outputs are in uK^2 as Dl = l*(l+1)/(2pi) * Cl,
    except for the lensing potential column which is l^4 * Cl and 
    the cltp column which is l^3 * Cl. Tcmb also needs to be 
    corrected for. If do_lensing=T and lens_potential_output_file
    will be a separate file. In this file, Cdd=[l(l+1)]^2*Cl/2pi, 
    CdT=[l(l+1)]^(3/2)Cl/2pi, and CdE=[l(l+1)]^(3/2)Cl/2pi.
    If `scalar_amp` and `CMB_outputscale` are both set to 1 in the
    CAMB params.ini file, the outputs will be dimensionless.
    Also note that usually the monopole (l=0) and dipole (l=1) are
    omitted, but e.g. healpy.synfast expects values for these.
    """
    camb_cls = {}
    # look in the camb data directory if the input file is a relative path that is
    # otherwise not found
    if not os.path.exists(filename):
        if not os.path.isabs(filename):
            filename = os.path.join(os.path.dirname(__file__), 'data/camb', filename)
    with open(filename) as f:
        header = f.readline().strip('\n').strip('#').split(' ')
        header = list(filter(('').__ne__, header))
    fmt = [float] * len(header)
    fmt[0] = float
    dat = np.loadtxt(filename, unpack=True, dtype={'names': header, 'formats': fmt})
    tcmb = 2.726 * 1e6  # uK
    for i, col in enumerate(header):
        if as_cls:
            if len(header) == 5:
                # *lensedCls.dat
                if i != 0:
                    conv = 2.0 * np.pi / (dat[0] ** 2 + dat[0])
                else:
                    conv = 1.0
            elif len(header) == 6:
                # *lensing_scalCls.dat
                if 'pp' in col.lower():
                    conv = 1.0 / (dat[0] ** 4) / tcmb ** 2
                elif 'p' in col.lower():
                    conv = 1.0 / (dat[0] ** 3) / tcmb
                elif i != 0:
                    conv = 2.0 * np.pi / (dat[0] ** 2 + dat[0])
                else:
                    conv = 1.0
            elif len(header) == 8:
                # *lenspotentialCls.dat
                if 'pp' in col.lower():
                    conv = 2.0 * np.pi / (dat[0] ** 2 + dat[0]) ** 2
                elif 'p' in col.lower():
                    conv = 2.0 * np.pi / (dat[0] ** 2 + dat[0]) ** 1.5
                elif i != 0:
                    conv = 2.0 * np.pi / (dat[0] ** 2 + dat[0])
                else:
                    conv = 1.0
        else:
            conv = 1.0
        camb_cls[col] = (conv * dat[i])[dat[0] < lmax + 1]

    if lmin is not None:
        orig = camb_cls['L']
        diff = int(min(camb_cls['L']) - lmin)
        for spec, ps in camb_cls.items():
            if diff < 0:
                camb_cls[spec] = ps[orig > lmin - 1]
            else:
                if spec == 'L':
                    camb_cls[spec] = np.concatenate((np.arange(diff), ps))
                else:
                    camb_cls[spec] = np.concatenate((np.zeros(diff), ps))

    return camb_cls
