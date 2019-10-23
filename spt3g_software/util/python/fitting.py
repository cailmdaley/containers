"""
This contains different fitting and minimization functions that might be
useful for maps or timestreams
"""
import numpy as np
import itertools
import scipy.optimize as opt
import scipy.stats
from spt3g import core
from . import maths
from . import stats

@core.usefulfunc
def fit_gaussian2d(data, fit_offset=False, verbose=False, **leastsq_kwargs):
    '''
    Returns (amp, x_0, y_0, sigma_x, sigma_y, theta, height),
    the parameters of a 2D gaussian found by a fit.
    '''
    params = stats.gaussian2d_moments(data)

    # Remove the vertical offset from the parameters if we want it fixed.
    if fit_offset is False:
        params=params[:-1]

    errorfunction = lambda p: np.ravel(
        maths.gaussian2d(*p)(*np.indices(data.shape)) - data)
    p, success = opt.leastsq(errorfunction, params, **leastsq_kwargs)

    if verbose:
        print("Finished fit: %s" % str(success))

    # Add a zero offset onto the end of the parameter list, if we kept it fixed.
    if fit_offset is False:
        p = np.append(p, 0.)

    return p

def fit_gaussian2d_fixedparam(data, verbose=False, **fixed_args):
    '''

    Returns (height, amplitude, center_x, center_y, width_x, width_y, rota),
    the parameters of a 2D gaussian found by a fit.
    
    This version allows you to fix some parameters of the 2D Gaussian by
    adding extra keyword arguments.
    
    Valid argument names for "fixed_args":
        ['amp', 'y_0', 'x_0', 'sigma_x', 'sigma_y', 'theta', 'height']
    
    AUTHOR
        Stephen Hoover, 10 October 2013
    '''
    # Check the input fixed arguments, and figure out what's left to fit. 
    arg_names = ['height', 'amplitude', 'center_x', 'center_y', 'width_x', 
                ' width_y', 'rota']
    for name in fixed_args: assert name in arg_names
    free_arg_names = filter(lambda x: x not in fixed_args, arg_names)
    
    params = stats.gaussian2d_moments(data) # Calculate initial parameters.
    
    # Remove fixed parameters from the list of initial parameters.
    params = [par for i_par, par in enumerate(params) if
              arg_names[i_par] not in fixed_args]
    
    # Fit the data.
    def errorfunction(*params):
        return np.ravel(
            maths.gauss2d( **dict(zip(free_arg_names, *params),
                               **fixed_args))(*np.indices(data.shape)) - data)
    params, success = opt.leastsq(errorfunction, params)

    if verbose:
        print("Finished fit: %s" % str(success))
        
    # Now rearrange the output into the expected order.
    param_dict = dict(zip(free_arg_names, params), **fixed_args)
    ordered_params = [param_dict[name] for name in arg_names]

    return ordered_params

@core.usefulfunc
def curve_fit_constrained(func, x, y, p0, sigma=None, fprime=None, bounds=[],
                          debug_level=0):
    """
    A wrapper to perform least-squares fitting using scipy.optimize.fmin_slsqp. 
    Fit the data points (x, y), with error on y of sigma, to the function func. 
    p0 is our initial guess at the values of the fit parameters.

    Note Bene: The default scipy implmentation of slsqp is bugged!
    You'll get segmentation faults if you don't fix it manually after
    installing scipy. See the sptpol_software README for fix info.

    Parameters:
    -----------
    func: The function which we wish to fit. Must take arguments `(x, *params)`, 
        where x is possibly (probably) a NumPy array.

    x: Independent variable of data.

    y: Dependent variable of data.

    p0: Initial guess of values for the function parameters.

    sigma [None]: Errors on y. If None (not recommended!),
        each point will be assigned error=1.

    fprime [None]: A function which returns a tuple of partial derivates of
        func. If None, the derivatives will be estimated inside
        scipy.optimize.fmin_slsqp.

    bounds [ [] ]: list
       A list of tuples specifying the lower and upper bound
       for each independent variable [(xl0, xu0),(xl1, xu1),...]

    debug_level [0]: Sets the amount of feedback printed to the screen by
        scipy.optimize.fmin_slsqp.
    """
    def sum_residual_sq(p, y, x, sigma):
        # Sum of residuals, squared. Compare the data to the fit function.
        err = y - func(x, *p)
        return np.sum((err/sigma)**2)

    if fprime is not None:
        # Assume that fprime gives the partial derivatives of func with
        # respect to each parameter. Now take the partial derivative of
        # sum_residual_sq with respect to each parameter.
        def sum_residual_sq_deriv(p, y, x, sigma):
            err = func(x, *p) - y
            return np.sum(2*err/(sigma**2)*fprime(x,*p), axis=1)
    else:
        sum_residual_sq_deriv = None

    if bounds is None: bounds=[]
        
    if len(x)!=len(y):
        raise ValueError("Input x and input y must have the same length.")

    if sigma is None: sigma = np.ones(len(x))
    if len(sigma)!=len(x):
        raise ValueError("Input errors must have the same length as the data.")
    
    fit = opt.fmin_slsqp(sum_residual_sq, np.ascontiguousarray(p0),
                         args=(y, x, sigma), iprint=debug_level,
                         fprime=sum_residual_sq_deriv, bounds=bounds,
                         full_output=True, iter=200)

    return fit

@core.usefulfunc
def grid_search_min(func,x,y,sigma,grid_ranges):
    """
    Performs a grid search minimization of y-func(x,params)

    Parameters:
    -----------
    func: objective function to which y values are being fit.
        Should take x and a list of parameters as its input

    x: independent variable array

    y: array of data points

    sigma: array of errors on y

    func: function of x and params that is modeling y(x)

    grid_ranges: ranges of values for parameters for func to be iterated over. 
        Should be a nested list of lists/arrays.  The ranges should be in the
        same list as the parameters sent to "func", so
        [[a1_0,...,a1_n],[a2_0,...,a2_n],...,[an_0,...,an_n]]
        if the params for func are f(x,a1,a2,...,an)
    """
    param_grid = itertools.product(*grid_ranges)
    
    chi_sq_min=sum(y**2)
    cnt=0
    best_fit_params= None
    stat = 0
    while True:
        try:
            params = param_grid.next()
            chi_sq = sum(((y - func(x,*params))/sigma)**2)
            if chi_sq < chi_sq_min:
                chi_sq_min = chi_sq
                best_fit_params = list(params)
                stat_message = 'success'
            cnt+=1
        except StopIteration:
            break
    if best_fit_params is None:
        best_fit_params = params
        stat = -1
        stat_message = 'failure'

    return best_fit_params, chi_sq, cnt, stat, stat_message

@core.usefulfunc
def gaussfit_hist(values, nbins, minval, maxval, do_plot=False):
    """
    histogram some values, fit the histogram to a gaussian, and return
    the best-fit parameters.
    """

    binsize = (maxval-minval)/(np.float(nbins-1))
    ht = np.histogram(values*1.0,bins=nbins,range=(minval,maxval))
    dx = ht[1][1] - ht[1][0]
    x_out = ht[1][0:nbins] + dx/2.
    hist_out = ht[0]

    # throw out some outliers before fittin
    (mutemp,sigmatemp) = scipy.stats.norm.fit(values)
    whgood = np.where(np.abs(values-mutemp)/sigmatemp < 5.)
    (mu,sigma) = scipy.stats.norm.fit(values[whgood])
    yfit_unnorm = np.exp(-(x_out-mu)**2/2./sigma**2)
    ampl = np.sum(yfit_unnorm*hist_out)/np.sum(yfit_unnorm**2)
    yfit = ampl*yfit_unnorm
    hfit_out = yfit

    params = [ampl,mu,sigma]

    if do_plot:
        import matplotlib.pyplot as plt
        plt.hist(values, bins=nbins, range=(minval,maxval), color='g')
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = ampl*np.exp(-(x-mu)**2/2./sigma**2)
        plt.plot(x, p, 'k', linewidth=2)
        title = "Fit results: mu = %.4g,  std = %.4g" % (mu, sigma)
        plt.title(title)

    dict_out = {'params':params, 'x_out':x_out, 'hist_out':hist_out,
                'hfit_out':hfit_out}

    return dict_out
