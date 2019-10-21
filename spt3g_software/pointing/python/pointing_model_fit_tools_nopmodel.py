from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
import re
import os, pickle
from scipy.optimize import leastsq

def sex2deg(val):
    """
    Convert a sexagesimal number to degrees.

    Arguments
    ---------
    val : string? or float?
        A sexagesimal value to be converted to degrees

    Returns
    -------
    thing : float
        A value in degrees.
    """
    if ':' not in str(val):
        return float(val)
    try:
        val = val.decode()
    except AttributeError:
        pass
    sign = -1 if val.strip()[0] == '-' else 1
    try:
        d, m, s = [np.abs(float(x)) for x in val.split(':')]
    except ValueError:
        raise ValueError("Cannot parse sexadecimal value {}".format(val))
    return sign * (d + m / 60. + s / 3600.)

# Pointing model functions


def calculate_pointing_model_offsets(az, el, a_vector):
    """
    Takes perfect (topocentric) az/el, applies pointing model with 
    selected parameters, returns difference between result and original
    values. Actually, return negative of that, because input offsets are 
    defined as actual measured position minus topocentric position.

    Arguments
    ---------
    az : number
        An azimuth position in units of degrees.
    el : number
        An elevation position in units of degrees.
    a_vector : vector
        Vector of pointing model parameters, ordered as in official documentation (all in units of degrees):
      a[0] : high-el flexure
      a[1] : low-el flexure
      a[2] : azimuth bearing tilt in az=0 direction
      a[3] : azimuth bearing tilt in az=90 direction
      a[4] : elevation bearing tilt relative to plane of az bearing
      a[5] : cross-elevation collimation
      a[6] : elevation collimation
      a[7] : azimuth encoder offset
      a[8] : elevation encoder offset

    Returns
    -------
    d_az, d_el : tuple of floats
    """
    el_rad = el*np.pi/180.
    az_rad = az*np.pi/180.
    a_vector_rad = np.asarray(a_vector)*np.pi/180.
    a = a_vector_rad

    el1 = el_rad - a[0]*np.sin(el_rad) - a[1]*np.cos(el_rad)
    theta_tilt = np.arccos(np.cos(a[3])*np.cos(a[2]*np.cos(a[3])))
    phi_tilt = -np.pi - np.arctan2(np.sin(a[2]*np.cos(a[3])),(np.cos(a[2]*np.cos(a[3]))*np.sin(a[3])))
    w = phi_tilt - np.pi/2. - az_rad
    el2 = np.arcsin(np.sin(el1)*np.cos(theta_tilt) - np.cos(el1)*np.sin(theta_tilt)*np.sin(w))
    el3 = np.arcsin(np.sin(el2)/np.cos(a[4]))
    el4 = el3 + a[6]
    el5 = np.arcsin(np.sin(el4)*np.cos(a[5]))
    el6 = el5 + a[8]

    az1 = phi_tilt - np.arctan2(-(np.cos(el2)*np.cos(w)),(np.cos(theta_tilt)*np.sin(w)*np.cos(el2) + np.sin(theta_tilt)*np.sin(el2))) + np.pi
    az2 = az1 - np.arcsin(np.tan(a[4])*np.tan(el3))
    az3 = az2 + np.arctan2(np.tan(a[5]),np.cos(el4))
    az4 = az3 + a[7]

    az_new = az4*180./np.pi
    if np.size(az_new) == 1:
        if az_new < 0.:
            az_new += 360.
    else:
        az_new[np.where(az_new < 0.)] += 360.
    el_new = el6*180./np.pi

    return az_new-az, el_new-el


def model_daz(az, el, a2, a3, a4, a5, a7, a8):
    """
    Calculate the first-order approximation to the pointing model
    applied in GCP to go from topocentric azimuth to raw encoder
    azimuth. (I.e., returns az_actual minus az_topo for these pointing
    parameters.)

    Arguments
    ---------
    az : number
        An azimuth position in units of degrees.
    el : number
        An elevation position in units of degrees. 
    a2 : first component of az tilt [degrees]
    a3 : second component of az tilt [degrees]
    a4 : el tilt [degrees]
    a5 : cross-elevation collimation [degrees]
    a7 : az encoder offset [degrees]
    a8 : el encoder offset [degrees]

    Returns
    -------
    d_az : float
        actual minus topocentric az in degrees

    """
    el_rad = np.radians(el)
    az_rad = np.radians(az)
    a8_rad = np.radians(a8)

    d_az = -(a2 * np.cos(az_rad) + a3 * np.sin(az_rad)) * np.tan(el_rad) - a4 * np.tan(el_rad) + a5 / np.cos(el_rad) + a7

    return d_az

def model_del(az, el, a0, a1, a2, a3, a6, a8):
    """
    Calculate the first-order approximation to the pointing model
    applied in GCP to go from topocentric elevation to raw encoder
    elevation. (I.e., returns el_actual minus el_topo for these pointing
    parameters.)

    Arguments
    ---------
    az : number
        An azimuth position in units of degrees.
    el : number
        An elevation position in units of degrees.
    a0 : first component of flexure [degrees]
    a1 : second component of flexure [degrees]
    a2 : first component of az tilt [degrees]
    a3 : second component of az tilt [degrees]
    a6 : elevation collimation [degrees]
    a8 : el encoder offset [degrees]

    Returns
    -------
    d_el : float
        actual minus topocentric el in degrees
    """
    el_rad = np.radians(el)
    az_rad = np.radians(az)

    d_el = -a0 * np.sin(el_rad) - a1 * np.cos(el_rad) + a2 * np.sin(az_rad) - a3 * np.cos(az_rad) + a6 + a8

    return d_el

# SPT model

def model_factory(az, el, parlist):
    """
    Generate a model function for SPT pointing

    Arguments
    ---------
    az : array-like
        Azimuth points at which observations were taken
    el : array-like
        Elevation points at which observatins were taken
    parlist : list of strings
        A list of string names for the parameters in the order they will
        be supplied to the model function

    Returns
    -------
    fun : callable
        A function that takes as input a list of parameters and returns
        the resulting model. Fixed parameters can be supplied as optional
        keyword arguments or an input dictionary.
    """
    def f(p, *args, **kwargs):
        split = kwargs.pop('split', False)

        # get parameter order from parlist
        pp = {}
        for idx, par in enumerate(parlist):
            pp[par] = p[idx]

        # fixed parameters
        if len(args) and isinstance(args[0], dict):
            kwargs.update(**args[0])
        for k, v in kwargs.items():
            pp[k] = v

        thisdaz = np.zeros_like(az)
        thisdel = np.zeros_like(az)
        
        newmethod = True
#        newmethod = False
        if newmethod:
            thisdaz, thisdel = calculate_pointing_model_offsets(az, el, [pp['a0'], pp['a1'], pp['a2'], pp['a3'], pp['a4'], pp['a5'], pp['a6'], pp['a7'], pp['a8']])
        else:
            thisdaz = model_daz(az, el, pp['a2'], pp['a3'], pp['a4'], pp['a5'], pp['a7'], pp['a8'])
            thisdel = model_del(az, el, pp['a0'], pp['a1'], pp['a2'], pp['a3'], pp['a6'], pp['a8'])
        if split:
            return thisdaz, thisdel
        return np.concatenate((thisdaz, thisdel))

    return f

def lmfitfun(model, data, err):
    """
    Construct a fit function for use by scipy.optimize.least_squares.

    Arguments
    ---------
    model : callable
        A model function that takes a list of parameters
    data : array-like
        A list of data points to input to the model
    err : array-like
        A list of errors on each data point

    Returns
    -------
    fitfun : callable
        A fit function that can be passed to scipy.optimize.least_squares
    """
    def f(p, *args, **kwargs):
        return (model(p, *args, **kwargs) - data) / err
    return f

# Define parameter order
parlist = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8']

# ================================================================================
# Read data

def read_data(spt_files=[]):

    az = []
    el = []
    az_offset = []
    el_offset = []
    is_spt = []

    if not isinstance(spt_files, list):
        spt_files = [spt_files]

    # SPT positions and offsets
    #   - change magic file format numbers for 2019
    for spt_file in spt_files:
        # Read in the data
        data = np.loadtxt(spt_file, delimiter=', ', usecols=[4, 5, 6, 7, 8, 9, 10, 11, 12], unpack=True)
        this_az = data[4]
        this_el = data[5]
        # this is making the approximation that az_off=ra_off and el_off=-dec_off
        this_az_off = data[0] / 60.
        this_el_off = -data[1] / 60.

        # store
        az.append(this_az)
        el.append(this_el)
        az_offset.append(this_az_off)
        el_offset.append(this_el_off)

    # Cleanup and sort by elevation
    az = np.concatenate(az)
    el = np.concatenate(el)
    az_offset = np.concatenate(az_offset)
    el_offset = np.concatenate(el_offset)

    sel = np.argsort(el)
    az = az[sel]
    el = el[sel]
    az_offset = az_offset[sel]
    el_offset = el_offset[sel]

    # Fake errors
    errors = np.zeros(len(az)) + 0.05

    return dict(az=az, el=el, az_offset=az_offset, el_offset=el_offset,
                errors=errors)

# ================================================================================
# Fit the data

def fit_data(in_data, fit_params=None, fix_boom_terms=False, fix_az_tilts=True,
             fix_el_tilt=False, fixed_params=None):

    this_parlist = list(parlist) # local changes

    online_params = {}
    online_params['a0'] = 0.05471194
    online_params['a1'] = 0.05054306
    online_params['a2'] = 0.00693388
    online_params['a3'] = -0.00793388
    online_params['a4'] = -0.02206277
    online_params['a5'] = 0.04253917
    online_params['a6'] = -0.13203667
    online_params['a7'] = -0.304527
    online_params['a8'] = 0.

    if fixed_params is None:
        input_params = {}
        for parname in this_parlist:
            input_params[parname] = online_params[parname]
    else:
        input_params = fixed_params.copy()
        for parname in this_parlist:
            if parname not in input_params:
                input_params[parname] = online_params[parname]

    # Construct fit parameters dictionary
    # Use fixed params as the starting point
    pardict = {
        'a0' : {'value': input_params['a0'],
                'fixed': int(fix_boom_terms)},
        'a1' : {'value': input_params['a1'],
                'fixed': int(fix_boom_terms)},
        'a2' : {'value': input_params['a2'],
                'fixed': int(fix_az_tilts)},
        'a3' : {'value': input_params['a3'],
                'fixed': int(fix_az_tilts)},
        'a4' : {'value': input_params['a4'],
                'fixed': int(fix_el_tilt)},
        'a5' : {'value': input_params['a5'],
                'fixed': 0},
        'a6': {'value': input_params['a6'],
               'fixed': 0},
#        'a7': {'value': input_params['a7'],
#               'fixed': 1},
        'a7': {'value': input_params['a7'],
               'fixed': 0},
        'a8': {'value': input_params['a8'],
               'fixed': 1}
    }

    # collect fixed parameters
    kwargs = {}
    for k, d in pardict.items():
        if d['fixed']:
            this_parlist.pop(this_parlist.index(k))
            kwargs[k] = d['value']

    # Generate the full pointing model
    model = model_factory(in_data['az'], in_data['el'], this_parlist)

    # Combine data for the fit function
    data = np.concatenate((in_data['az_offset'], in_data['el_offset']))
    err = np.concatenate((in_data['errors'], in_data['errors']))

    # Compute the best fit
    par0 = [pardict[x]['value'] for x in this_parlist]
    fitfun = lmfitfun(model, data, err)
    # result = least_squares(fitfun, par0, method='lm', kwargs=kwargs)
    # par_fit = result.x
    result = leastsq(fitfun, par0, args=(kwargs,))
    par_fit = result[0]

    # Extract fit parameters
    popt = {k: v for k, v in zip(this_parlist, par_fit)}
    popt.update(kwargs)
    res = fitfun(par_fit, **kwargs) * err

    N = len(in_data['az_offset'])
    az_residuals = res[:N]
    el_residuals = res[N:]

    fit = res + data
    az_fit = fit[:N]
    el_fit = fit[N:]

    return dict(az_residuals=az_residuals, el_residuals=el_residuals,
                az_fit=az_fit, el_fit=el_fit, fit_params=popt, fitres=result)

# ================================================================================
# Process results

def process_fit(name, in_data, result, output_dir=".", verbose=True, plot=True):

    az = in_data['az']
    el = in_data['el']
    az_offset = in_data['az_offset']
    el_offset = in_data['el_offset']

    az_residuals = result['az_residuals']
    el_residuals = result['el_residuals']
    az_fit = result['az_fit']
    el_fit = result['el_fit']
    fit_params = result['fit_params']

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    fileroot = os.path.join(output_dir, name)

    with open("{}.txt".format(fileroot), 'w') as f:
        f.write("""# Pointing model {name} parameters:
az encoder offset {a7}
tilts {a2}, {a3}, {a4}
flexure radio, {a0}, {a1}
""".format(name=name, **fit_params))
        f.write("""# for SPT
collimate_fixed radio, {a5}, {a6}
""".format(**fit_params))

        f.write('\n# Fitting {} points:\n'.format(len(az)))
        f.write('# Az residuals rms (arcsec): {}\n'.format(np.std(az_residuals) * 3600))
        f.write('# El residuals rms (arcsec): {}\n'.format(np.std(el_residuals) * 3600))

        bad_idx, = np.where( (np.abs(az_residuals) * 3600 > 30) |
                             (np.abs(el_residuals) * 3600 > 30) )
        f.write('\n# Found {} sources with large outliers:\n'.format(len(bad_idx)))
        for idx in bad_idx:
            f.write("# {} {} {} {}\n".format(az[idx], el[idx], az_residuals[idx] * 3600, el_residuals[idx] * 3600))

    if verbose:
        os.system('cat {}.txt'.format(fileroot))

    if plot:

        def plot_stuff(fig, x, y, y_fit, xlabel, ylabel, title, filename):
            ix = np.argsort(x)
            x = x[ix]
            y = y[ix]
            if y_fit is not None:
                y_fit = y_fit[ix]

            plt.figure(fig)
            plt.clf()
            plt.plot(x, y, 'ko')
            if y_fit is not None:
                plt.plot(x, y_fit, 'b-')
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
            plt.savefig('{}_{}'.format(fileroot, filename), bbox_inches='tight')

        for x, xname, xlabel in [(az, 'az', 'Azimuth [deg]'),
                                 (el, 'el', 'Elevation [deg]')]:
            plot_stuff(1, x, az_offset, az_fit, xlabel, 'Az Offset [deg]',
                       'Az Offsets', '{}_vs_az_off'.format(xname))
            plot_stuff(1, x, az_residuals * 3600, None, xlabel, 'Az Residuals [arcsec]',
                       'Az Residuals', '{}_vs_az_res'.format(xname))
            plot_stuff(1, x, el_offset, el_fit, xlabel, 'El Offset [deg]',
                       'El Offsets', '{}_vs_el_off'.format(xname))
            plot_stuff(1, x, el_residuals * 3600, None, xlabel, 'El Residuals [arcsec]',
                       'El Residuals', '{}_vs_el_res'.format(xname))


# ================================================================================
# Main run function
def run(name, spt_files=[], fix_boom_terms=False, fix_az_tilts=True, fix_el_tilt=False,
        output_dir=".", process=True, verbose=True, plot=True):
    """
    Perform a least squares fit of the pointing model to various SPT 
    pointing observations.

    Arguments
    ---------
    name : string
        Indentifier for the output data.
    spt_files : list of strings
        Data files stored in the appropriate format, containing measured
        pointing offsets for each instrument.
        TODO: details on data format.
    fix_boom_terms : bool
        Fix the flexure parameters to the online parameter values.
    fix_az_tilts, fix_el_tilt : bool
        Fix az or el tilt terms in the fit to the online parameter values.
    output_dir : string
        Output directory where processed data should be stored.
    process : bool
        If True, post-process the fit data to print a summary text file
        and make several summary plots.
    verbose : bool
        If True, print the fit summary to screen.
    plot : bool
        If True, make summary plots.
    """

    in_data = read_data(spt_files)

    # !!!
    fixed_params = {}
    fixed_params['a2'] = 0.006729
    fixed_params['a3'] = -0.009464
    # !!!

    fit_result = fit_data(in_data, fix_boom_terms=fix_boom_terms,
                          fix_az_tilts=fix_az_tilts, fix_el_tilt=fix_el_tilt, fixed_params=fixed_params)
    if process:
        process_fit(name, in_data, fit_result, output_dir, verbose=verbose, plot=plot)

    return in_data, fit_result
