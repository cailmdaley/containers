import numpy as np
from scipy import optimize

def twoDgaussian(height, amp, y_0, x_0, sigma_y, sigma_x, theta):
    '''Returns a 2D gaussian with the given parameters.'''
    sigma_x = np.float(sigma_x)
    sigma_y = np.float(sigma_y)
    
    a = np.cos(theta)**2./(2.*sigma_x**2.) + np.sin(theta)**2./(2.*sigma_y**2.)
    b = np.sin(2.*theta)/(4.*sigma_x**2.) - np.sin(2.*theta)/(4.*sigma_y**2.)
    c = np.sin(theta)**2./(2.*sigma_x**2.) + np.cos(theta)**2./(2.*sigma_y**2.)

    return lambda x,y: height + amp*np.exp(-(a*(x-x_0)**2. + 2.*b*(x-x_0)*(y-y_0) + c*(y-y_0)**2.))

def moments(data):
    '''Returns (height, amp, x_0, y_0, sigma_x, sigma_y, theta), the gaussian
       parameters of a 2D distribution by calculating its moments.
    '''
    total = np.abs(data).sum()
    Y,X = np.indices(data.shape)
    y_0 = np.argmax((X*np.abs(data)).sum(axis=1)/total)
    x_0 = np.argmax((Y*np.abs(data)).sum(axis=0)/total)

    col = data[int(y_0),:]
    sigma_x = np.sqrt(np.abs((np.arange(col.size)-y_0)*col).sum()/np.abs(col).sum())

    row = data[:, int(x_0)]
    sigma_y = np.sqrt(np.abs((np.arange(row.size)-x_0)*row).sum()/np.abs(row).sum())

    height = np.median(data.ravel())
    amp = data.max() - height
       
    if np.isnan(sigma_x) or np.isnan(sigma_y) or np.isnan(height) or np.isnan(amp):
        #raise ValueError('Something is nan...')
        print 'Something is nan... no good...'
        mylist = [0.,0.,0.,0.,0.,0.,0.]

    mylist = [height,amp,x_0,y_0,sigma_x,sigma_y,0.]

    return mylist

def fit2Dgaussian(data):
    '''
    Returns (height, amp, x_0, y_0, sigma_x, sigma_y, theta), the parameters of a 2D
    gaussian found by a fit.
    '''
    params = moments(data)
    
    errorfunction = lambda p: np.ravel(twoDgaussian(*p)(*np.indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)

    return p
