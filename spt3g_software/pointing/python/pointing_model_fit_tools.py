from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
import re
import os
# from scipy.optimize import least_squares
from scipy.optimize import leastsq
from spt3g import core

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

def model_daz(az, el, a4, a5, a7, az0, a2, a3):
    """
    Calculate the azimuth pointing model.

    Arguments
    ---------
    az : number
        An azimuth position in units of degrees.
    el : number
        An elevation position in units of degrees.
    a4 : type?
        Description
    a5 : type?
        Description
    a7 : type?
        Description
    az0 : type?
        Description
    a2 : type?
        Description
    a3 : type?
        Description

    Returns
    -------
    d_az : float
        The azimuth offset caculated from the model in units of radians.
    """
    el = np.radians(el)
    az = np.radians(az)
    a7 = np.radians(a7)

    d_az = a4 * np.tan(el) - a5 / np.cos(el - a7) - az0
    d_az += (a2 * np.cos(az) + a3 * np.sin(az)) * np.tan(el) # check the sign of this term

    return d_az

def model_del(az, el, a0, a1, a6, a2, a3):
    """
    Calculate the elevation pointing model.

    Arguments
    ---------
    az : number
        An azimuth position in units of degrees.
    el : number
        An elevation position in units of degrees.
    a0 : type?
        Description
    a6 : type?
        Description
    a2 : type?
        Description
    a3 : type?
        Description

    Returns
    -------
    d_el : float
        The elevation offset caculated from the model in units of radians.
    """
    el = np.radians(el)
    az = np.radians(az)

    d_el = a0 * np.sin(el) + a1 * np.cos(el) - a6
    d_el += -(a2 * np.sin(az) - a3 * np.cos(az)) # check the sign of this term

    return d_el

# All-purpose combined SPT/EHT model

def model_factory(az, el, is_spt, parlist):
    """
    Generate a model function for combined SPT and EHT pointing.

    Arguments
    ---------
    az : array-like
        Azimuth points at which observations were taken
    el : array-like
        Elevation points at which observatins were taken
    is_spt : boolean, array-like
        Whether each point in the az/el arrays corresponds to SPT data
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

        model_az = np.zeros_like(az)
        model_el = np.zeros_like(az)

        if np.any(is_spt):
            model_az[is_spt] = model_daz(az[is_spt], el[is_spt], pp['a4'], pp['a5_spt'],
                                         pp['a7'], pp['az0'], pp['a2'], pp['a3'])
            model_el[is_spt] = model_del(az[is_spt], el[is_spt], pp['a0'], pp['a1'],
                                         pp['a6_spt'], pp['a2'], pp['a3'])

        if np.any(~is_spt):
            model_az[~is_spt] = model_daz(az[~is_spt], el[~is_spt], pp['a4'], pp['a5_eht'],
                                          pp['a7'], pp['az0'], pp['a2'], pp['a3'])
            model_el[~is_spt] = model_del(az[~is_spt], el[~is_spt], pp['a0'], pp['a1'],
                                          pp['a6_eht'], pp['a2'], pp['a3'])

        if split:
            return model_az, model_el
        return np.concatenate((model_az, model_el))

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

def read_model(input, use_offline_model=False):
    """
    Read the online pointing model from the file header or a frame 
    in a g3 file or a frame extracted from a g3 file.

    If header, lines should look like:

        # tilts -0.0114374830026, 0.00372502737393, -0.0224453462368
        # flexure radio, 0.0631999760291, 0.0592164524428
        # collimate_fixed radio, 0.0434865395172, -0.125024298806

    Arguments
    ---------
    input : string
        Either the absolute path to a file that contains either az/el offsets 
        and online pointing model (as strings) or a g3 file that 
        has a field 'OnlinePointingModel' in one of the frames OR a 
        g3 frame with 'OnlinePointingModel' as a key.

    use_offline_model: use the key OfflinePointingModel (rather than
        OnlinePointingModel) in frame. Default is False.

    Returns
    -------
    model : dictionary
        A dictionary containing keys corresponding to model parameters:
        'a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'az0'
    """
    #    model = {'a7': 0., 'az0': -0.304527}
    # Ugh, replace old magic number with new 2019 magic number
    model = {'a7': 0., 'az0': -0.337629}
    
    if isinstance(input,str):
        fn2use = input
        input_is_file = True
    else:
        fn2use = ''
        input_is_file = False

    model_from_frame = None
    if isinstance(input,core.G3Frame):
        if use_offline_model:
            if 'OfflinePointingModel' in input.keys():
                model_from_frame = input['OfflinePointingModel']
            else:            
                core.log_warn("No offline pointing model found in frame.")
                return model
        else:
            if 'OnlinePointingModel' in input.keys():
                model_from_frame = input['OnlinePointingModel']
            else:            
                core.log_warn("No online pointing model found in frame.")
                return model

    if model_from_frame is None and not input_is_file:
        core.log_warn("Input must be string filename or g3 frame with OnlinePointingModel field.")
        return model

    if model_from_frame is None and (fn2use.split('.')[-1] == 'g3' or fn2use.split('.')[-2] == 'g3'):
        for frame in core.G3File(fn2use):
            if 'OnlinePointingModel' in frame.keys():
                model_from_frame = frame['OnlinePointingModel']
        if model_from_frame is None:
            core.log_warn("No online pointing model found in file.")
            return model

    if model_from_frame is not None:
        model['a0'] = model_from_frame['flexure'][0]/core.G3Units.deg
        model['a1'] = model_from_frame['flexure'][1]/core.G3Units.deg
        model['a2'] = model_from_frame['tilts'][0]/core.G3Units.deg
        model['a3'] = model_from_frame['tilts'][1]/core.G3Units.deg
        model['a4'] = model_from_frame['tilts'][2]/core.G3Units.deg
        model['a5'] = model_from_frame['fixedCollimation'][0]/core.G3Units.deg
        model['a6'] = model_from_frame['fixedCollimation'][1]/core.G3Units.deg
        return model

    with open(fn2use,'r') as f1:
        for line in f1.readlines():
            if not line.startswith('#'):
                continue

            m = re.search('tilts ([0-9\.\-\:]+), ([0-9\.\-\:]+), ([0-9\.\-\:]+)', line.strip())
            if m:
                model['a2'], model['a3'], model['a4'] = [sex2deg(x) for x in m.groups()]

            m = re.search('flexure radio, ([0-9\.\-\:]+), ([0-9\.\-\:]+)', line.strip())
            if m:
                model['a0'], model['a1'] = [sex2deg(x) for x in m.groups()]

            m = re.search('collimate_fixed radio, ([0-9\.\-\:]+), ([0-9\.\-\:]+)', line.strip())
            if m:
                model['a5'], model['a6'] = [sex2deg(x) for x in m.groups()]

    return model

def remove_model(az_offset, el_offset, model, params, parlist):
    """
    Remove the pointing model from ???

    Arguments
    ---------
    az_offset : type?
        Description
    el_offset : type?
        Description
    model : type?
        Description
    params : type?
        Description
    parlist : type?
        Description

    Returns
    -------
    -az_offset, -el_offset : tuple, 2 elements
        Description
    """

    if isinstance(params, dict):
        params = [params[k.split('_')[0]] for k in parlist]
    model_az, model_el = model(params, split=True)
    # check signs here
    az_offset -= model_az
    el_offset -= model_el
    return -az_offset, -el_offset

# Define parameter order
parlist = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5_spt', 'a5_eht', 'a6_spt', 'a6_eht', 'a7', 'az0']

# ================================================================================
# Read data

def read_data(spt_files=[], eht_files=[], cut_hi_el=4, el_range=None, online_params=None):

    az = []
    el = []
    az_offset = []
    el_offset = []
    is_spt = []

    if not isinstance(spt_files, list):
        spt_files = [spt_files]
    if not isinstance(eht_files, list):
        eht_files = [eht_files]

    # SPT positions and offsets
    #   - change magic file format numbers for 2019
    for spt_file in spt_files:
        params = read_model(spt_file)
        if not online_params:
            online_params = params

#        # Read in the data and remove the online pointing model
#        data = np.loadtxt(spt_file, delimiter=', ', usecols=[4, 5, 6, 8, 10, 12, 14, 15, 16], unpack=True)
#        this_az = data[0]
#        this_el = data[1]
#        this_is_spt = np.ones_like(this_az, dtype=bool)
#        this_az_off = data[4] / 60.
#        this_el_off = -data[5] / 60.
#        tilts0 = data[6]
#        tilts1 = data[7]

        # Read in the data and remove the online pointing model
        data = np.loadtxt(spt_file, delimiter=', ', usecols=[4, 5, 6, 7, 8, 9, 10, 11], unpack=True)
        this_az = data[0]
        this_el = data[1]
        this_meas_az = data[2]
        this_az_off = (this_meas_az - this_az)
        this_meas_daz = data[3]
        this_meas_el = data[4]
        this_el_off = (this_meas_el - this_el)
        this_meas_del = data[5]
        this_is_spt = np.ones_like(this_az, dtype=bool)
        tilts0 = data[6]
        tilts1 = data[7]

        this_az_off_orig = this_az_off.copy()
        this_el_off_orig = this_el_off.copy()

        # Actually, we now measure tilts for every observation independently, so remove model for every source individually.
        this_az_off = np.zeros(len(this_az_off_orig))
        this_el_off = np.zeros(len(this_az_off_orig))
        for q in np.arange(len(this_az_off_orig)):
            this_model = model_factory(this_az[q], this_el[q], this_is_spt[q], parlist)
            params_temp = params
            params_temp['a2'] = tilts0[q]
            params_temp['a3'] = tilts1[q]
            this_az_off_temp, this_el_off_temp = remove_model(this_az_off_orig[q], this_el_off_orig[q], this_model, params_temp, parlist)
            this_az_off[q] = this_az_off_temp
            this_el_off[q] = this_el_off_temp
                
        # store
        az.append(this_az)
        el.append(this_el)
        az_offset.append(this_az_off)
        el_offset.append(this_el_off)
        is_spt.append(this_is_spt)

    # EHT positions and offsets
    for eht_file in eht_files:
        # Read in the online pointing model data
        params = read_model(eht_file)
        if not len(spt_files) and not online_params:
            online_params = params

        # Read in the data and remove the online pointing model
        this_az, this_el, this_az_off, this_el_off = np.loadtxt(
            eht_file, converters={x: sex2deg for x in [1, 2, 3, 4]},
            usecols=[1,2,3,4], unpack=True)

        this_az_off /= np.cos(np.radians(this_el))
        this_is_spt = np.zeros_like(this_az, dtype=bool)
        this_model = model_factory(this_az, this_el, this_is_spt, parlist)
        this_az_off, this_el_off = remove_model(this_az_off, this_el_off, this_model,
                                                params, parlist)

        # store
        az.append(this_az)
        el.append(this_el)
        az_offset.append(this_az_off)
        el_offset.append(this_el_off)
        is_spt.append(this_is_spt)

    # Cleanup and sort by elevation
    az = np.concatenate(az)
    el = np.concatenate(el)
    az_offset = np.concatenate(az_offset)
    el_offset = np.concatenate(el_offset)
    is_spt = np.concatenate(is_spt)

    sel = np.argsort(el)
    az = az[sel]
    el = el[sel]
    az_offset = az_offset[sel]
    el_offset = el_offset[sel]
    is_spt = is_spt[sel]

    # !!! cut low elevation points
    if cut_hi_el:
        az = az[:-cut_hi_el]
        el = el[:-cut_hi_el]
        az_offset = az_offset[:-cut_hi_el]
        el_offset = el_offset[:-cut_hi_el]
        is_spt = is_spt[:-cut_hi_el]

    if el_range is not None:
        idx = np.logical_and(el >= el_range[0], el <= el_range[1])
        az = az[idx]
        el = el[idx]
        az_offset = az_offset[idx]
        el_offset = el_offset[idx]
        is_spt = is_spt[idx]

    ## wrap az to -180 < az < 180
    # az_orig = az.copy()
    # az[np.where(az > 180.)] -= 360.

    # Fake errors
    errors = np.zeros(len(az)) + 0.05

    return dict(az=az, el=el, az_offset=az_offset, el_offset=el_offset,
                is_spt=is_spt, errors=errors, online_params=online_params)

# ================================================================================
# Fit the data

def fit_data(in_data, fit_params=None, fixed_params=['boom_terms', 'az_tilts', 'x_collimation', 'az_encoder_offset']):

    online_params = in_data['online_params']
    use_spt = np.any(in_data['is_spt'])
    use_eht = not np.all(in_data['is_spt'])
    this_parlist = list(parlist) # local changes

    input_params = online_params.copy()
    if fit_params:
        input_params.update(**fit_params)

    # Construct fit parameters dictionary
    # Use the online model as the starting point
    pardict = {
        'a0' : {'value': input_params['a0'],
                'fixed': int('boom_terms' in fixed_params)},
        'a1' : {'value': input_params['a1'],
                'fixed': int('boom_terms' in fixed_params)},
        'a2' : {'value': input_params['a2'],
                'fixed': int('az_tilts' in fixed_params)},
        'a3' : {'value': input_params['a3'],
                'fixed': int('az_tilts' in fixed_params)},
        'a4' : {'value': input_params['a4'],
                'fixed': int('el_tilt' in fixed_params)},
        'a5_spt' : {'value': input_params['a5'],
                    'fixed': int(not use_spt)},
        'a5_eht': {'value': input_params['a5'],
                   'fixed': int(not use_eht)},
        'a6_spt': {'value': input_params['a6'],
                   'fixed': int(not use_spt)},
        'a6_eht': {'value': input_params['a6'],
                   'fixed': int(not use_eht)},
        'a7': {'value': input_params['a7'],
               'fixed': 1},
        'az0': {'value': input_params['az0'],
                'fixed': int('az_encoder_offset' in fixed_params)},
    }
    if 'x_collimation' in fixed_params:
        pardict['a5_spt']['fixed'] = 1
    if 'y_collimation' in fixed_params:
        pardict['a6_spt']['fixed'] = 1

    # collect fixed parameters
    kwargs = {}
    for k, d in pardict.items():
        if d['fixed']:
            this_parlist.pop(this_parlist.index(k))
            kwargs[k] = d['value']

    # Generate the full pointing model
    model = model_factory(in_data['az'], in_data['el'], in_data['is_spt'], this_parlist)

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
    is_spt = in_data['is_spt']

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
az encoder offset {az0}
tilts {a2}, {a3}, {a4}
flexure radio, {a0}, {a1}
""".format(name=name, **fit_params))
        if np.any(is_spt):
            f.write("""# for SPT
collimate_fixed radio, {a5_spt}, {a6_spt}
""".format(**fit_params))
        if not np.all(is_spt):
            f.write("""# for EHT
collimate_fixed radio, {a5_eht}, {a6_eht}
""".format(**fit_params))

        f.write('\n# Fitting {} points:\n'.format(len(az)))
        f.write('# Az residuals rms (arcsec): {}\n'.format(np.std(az_residuals) * 3600))
        f.write('# El residuals rms (arcsec): {}\n'.format(np.std(el_residuals) * 3600))

        bad_idx, = np.where( (np.abs(az_residuals) * 3600 > 30) |
                             (np.abs(el_residuals) * 3600 > 30) )
        f.write('\n# Found {} sources with large outliers:\n'.format(len(bad_idx)))
        for idx in bad_idx:
            f.write("# {} {} {} {} {}\n".format(az[idx], el[idx], az_residuals[idx] * 3600,
                                                el_residuals[idx] * 3600, is_spt[idx]))

    if verbose:
        os.system('cat {}.txt'.format(fileroot))

    if plot:

        def plot_stuff(fig, x, y, y_fit, xlabel, ylabel, title, filename):
            ix = np.argsort(x)
            x = x[ix]
            y = y[ix]
            if y_fit is not None:
                y_fit = y_fit[ix]
            ispt = is_spt
            ispt = ispt[ix]

            plt.figure(fig)
            plt.clf()
            plt.plot(x[ispt], y[ispt], 'ko')
            plt.plot(x[~ispt], y[~ispt], 'ro')
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
            plot_stuff(2, x, az_residuals * 3600, None, xlabel, 'Az Residuals [arcsec]',
                       'Az Residuals', '{}_vs_az_res'.format(xname))
            plot_stuff(3, x, el_offset, el_fit, xlabel, 'El Offset [deg]',
                       'El Offsets', '{}_vs_el_off'.format(xname))
            plot_stuff(4, x, el_residuals * 3600, None, xlabel, 'El Residuals [arcsec]',
                       'El Residuals', '{}_vs_el_res'.format(xname))


# ================================================================================
# Main run function
def run(name, spt_files=[], eht_files=[], cut_hi_el=0, el_range=None,
        fit_params=None, 
        fixed_params=['boom_terms', 'az_tilts', 'x_collimation', 'az_encoder_offset'],
        output_dir=".", process=True, verbose=True, plot=True):

    """Perform a least squares fit of the pointing model to various SPT and EHT
    pointing observations.

    Arguments
    ---------
    name : string
        Indentifier for the output data.
    spt_files, eht_files : list of strings
        Data files stored in the appropriate format, containing measured
        pointing offsets for each instrument.
        TODO: details on data format.
    cut_hi_el : integer
        Cut this number of points at the lowest elevation from the combined
        SPT+EHT dataset
    fit_params : dict
        Dictionary of pointing parameters to input to the least squares fitter,
        either as values to be fixed (if any of the options below are True) or
        as new starting points for the fitter.  Values not supplied here
        are taken from the set of parameters loaded from the first input data
        file.
    fixed_params : List of parameters to fix in fit. Choices are 'boom_terms',
        'az_tilts', 'el_tilt', 'x_collimation', 'y_collimation', and
        'az_encoder_offset'. Default=['boom_terms', 'az_tilts',
        'x_collimation', 'az_encoder_offset']. We assume the user will never
        want to fix just one of the boom or az tilt terms.
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

    in_data = read_data(spt_files, eht_files, cut_hi_el=cut_hi_el, el_range=el_range)
    fit_result = fit_data(in_data, fit_params=fit_params, fixed_params=fixed_params)
    if process:
        process_fit(name, in_data, fit_result, output_dir, verbose=verbose, plot=plot)

    return in_data, fit_result
