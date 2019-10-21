'''
Calculate the mode-mode coupling kernel as in the MASTER paper Hivon et al.
astro-ph/0105302v1 Appendix A1. This def returns the coupling kernel M_k1_k2,
which is used to account for mode mixing as follows:
  Cl_meas[k1] = M_k1_k2 * Cl_true[k2]
  IDL> Cl_meas = M_k1_k2 ## Cl_true
  IDL> Cl_true = inv_M ## Cl_meas

Note that this is not exactly equation A14 from Hivon et al.  It has
been scaled by the determinant factor k2 in equation A13 in order
for the matrix-formalism above to work.  See Martin's notes
for more details:
  https://spt.uchicago.edu/trac/attachment/wiki/AnalysisNotesFeb2007/Pseudo-Cl-081021.pdf 

For details on the polarized mode-mixing term, see:
  https://pole.uchicago.edu/sptpol/index.php/Polarized_Coupling_Kernel_for_EE_and_TE
  https://pole.uchicago.edu/sptpol/images/Mode_mixing_qu_201308.pdf
  Appendix to Crites el al. SPTpol 100d EE/TE astro-ph/1411.1042

MODIFICATION HISTORY:
    2008/09    - authored by Martin Lueker
    2010/06/11 - added lowmem keyword, much slower but less ram usage.
    2010/06/15 - replaced lowmem keyword with "fast" keyword which does the
            opposite, now the default is bitwise identical to Lueker 09 version
    2013/07/31 - (KTS) Re-write the "grid-shift" definition of v in
             coupling_kernel. The output is identical to the previous version.
    2013/08/15 - (KTS) Add the determinant scaling factor into this function.
            Note, this no longer returns the same thing as the original
            spt_analysis code.
    2013/08/20 - (KTS) Add option to calculate the EE polarized coupling kernel
    2013/09/11 - (KTS) Add option to calculate the TE polarized coupling kernel
    2019/05/18 - (DPD) Convert from IDL to python. The output is not identical to
            the IDL version due to a difference in integration methods.
    2019/10/01 - (SR and DPD) Significant speed improvements. Got rid of
            unused options cheby=False, changevar=False, fast=True, angle_interp=False
'''
import time
import numpy as np
from scipy.interpolate import interp1d
import scipy.integrate
from spt3g import core


def coupling_kernel(
    apodization_mask,
    res=None,
    spec_type='TT',
    curlyw=None,
    lmax=10000,
    oversampfact=8,
    interpfact=1000,
    verbose=False,
):
    '''
    Calculate the mode-mode coupling kernel of `apodization_mask`.

    This function performs the analytical calculation of the
    TT mode-mode coupling kernel as found in the
    MASTER paper Hivon et al. astro-ph/0105302v1 Appendix A1,
    and of the polarized kernels TE, EE as in the appendix to
    the Crites et al. 100d paper, astro-ph/1411.1042.

    Parameters
    ----------
    apodization_mask : array-like
        Must be square
    res : float, optional
        The resolution of `apodization_mask` in G3 units [radians].
        Must be specifed if `apodization_mask` is an ndarray.
    spec_type : ['TT', 'TE', 'EE', 'BB']
        The spectrum for which to compute the coupling kernel.
        Note: 'BB' uses the same calculation as 'EE', and only considers
        B-to-B coupling, and not e.g. E-to-B coupling.
    curlyw : array-like, optional
        Integral over theta of the fourier transform of `apodization_mask`.
        Calculated by calc_curly_w()
    lmax : int
        The maximum ell up to which to compute the coupling kernel
    oversampfact : int
        How much to oversample the curlyw integral
    interpfact : int
        interpolation factor
    verbose : bool
        Reduce log level to 'INFO'

    Returns
    -------
    Dictionary with the following keys :
        'kernel' : 2-d array
            The coupling kernel M_k1_k2, which is used
            to account for mode mixing as follows:
                Cl_meas[k1] = M_k1_k2 * Cl_true[k2]
        'ellkern' : 1-d array
            The central values of the ell bins corresponding
            to the rows and columns of M_k1_k2
        'curlyw' : 1-d array
            Integral over theta of the fourier transform of
            `apodization_mask`.
    '''
    if verbose:
        core.set_log_level(core.G3LogLevel.LOG_INFO, unit='CouplingKernel')
    core.log_info('Spec_type is: ' + spec_type, unit='CouplingKernel')

    if hasattr(apodization_mask, 'res'):
        res = apodization_mask.res
    else:
        if res is None:
            raise ValueError("Input has no 'res' attribute and kwarg `res` not set.")

    kmax = lmax / 2.0 / np.pi

    s = np.shape(apodization_mask)
    # check the mask
    if len(s) != 2:
        raise ValueError('Apodization mask must be a 2-d array')
    if s[0] != s[1]:
        raise ValueError('Apodization mask must be square')

    ngrid = s[0]
    delta_k = 1.0 / (ngrid * res)

    if curlyw is None:
        curlyw = calc_curly_w(
            apodization_mask, delta_k=delta_k, oversampfact=oversampfact
        )

    wincnt = len(curlyw)
    k_wfact = np.arange(len(curlyw)) * (delta_k / oversampfact) ** 2 * 2.0 * np.pi
    k_wfact[0] = (delta_k / oversampfact / 2.0) ** 2 * np.pi
    k_1d = np.arange(int(ngrid / 2))
    kidx = np.where(k_1d < (kmax / delta_k) * 2)[0]
    nk = len(kidx)
    mkk = np.zeros((nk, nk))
    wweighted = curlyw * k_wfact

    ellkern = (np.arange(nk) + 0.5) * delta_k * 2.0 * np.pi

    start = time.time()
    lastreport = start
    reportinterval = 60
    core.log_info('Calculating coupling kernel', unit='CouplingKernel')

    ncurlyw = len(curlyw)
    k3 = np.arange(ncurlyw)
    npts = oversampfact * interpfact
    # points chosen for interpolation by gauss chebyshev quadrature
    # http://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature
    tmp = (2 * (np.arange(npts) + 1.0) - 1) / 2.0 / npts * np.pi
    vdesired = np.cos(tmp)
    # added 0 to iteration 2/20/09
    for i in np.arange(nk):
        now = time.time()
        if (now - lastreport) > reportinterval:
            # this algorithm should be roughly quadratic in npix
            progress = (float(i) / nk) ** 2
            remaining = (1.0 - progress) * (now - start) / progress
            remainingmin = remaining / 60
            core.log_info(
                'Estimated time remaining: %.2f min' % remainingmin,
                unit='CouplingKernel',
            )
            lastreport = now
            # report every 5minutes from now on
            reportinterval = 300

        # added 0 to iteration 2/20/09
        k1 = (i + 0.5) * oversampfact
        for j in np.arange(i + 1):
            k2 = (j + 0.5) * oversampfact
            # v here is the same as "u" in Eqn A11 from Hivon et al.
            # It is derrived by completing the square of (k3**2)**2
            # in the definition of J()
            # k3 defined outside loop
            v = (k3 ** 2 - k1 ** 2 - k2 ** 2) / (2.0 * k1 * k2)
            vmax = v[ncurlyw - 1]
            # If k1+k2 > max(k3) : some points of Vdes
            # will be outside of the interpolation range of
            # curly w.  We identify these points to avoid
            # extrapolation.
            idx = np.where(vdesired < vmax)[0]
            icnt = len(idx)
            if icnt > 0:
                if spec_type == 'TT':
                    integrand = curlyw
                elif spec_type == 'EE' or spec_type == 'BB':
                    pol_term = (
                        (k1 ** 2 + k2 ** 2 - k3 ** 2) ** 2 / (2 * k1 ** 2 * k2 ** 2) - 1
                    ) ** 2.0
                    integrand = curlyw * pol_term
                elif spec_type == 'TE':
                    pol_term = (k1 ** 2 + k2 ** 2 - k3 ** 2) ** 2 / (
                        2 * k1 ** 2 * k2 ** 2
                    ) - 1
                    integrand = curlyw * pol_term
                else:
                    raise ValueError("spec_type not recognized: ", spec_type)
                # Note: this interpolation onto vdesired enforces the
                # correct limits on the integral: see Hivon eqn A12.
                int_func = interp1d(v, integrand, kind='linear')
                tmp = (2.0 * np.pi / npts) * np.sum(int_func(vdesired[idx]))
                mkk[i, j] = tmp
                mkk[j, i] = tmp
            else:
                mkk[i, j] = 0
                mkk[j, i] = 0

    core.log_info('Done', unit='CouplingKernel')

    # scale the kernel by the Jacobian determinant conversion:
    kernsize = nk
    winsize = ngrid
    u = np.outer(
        (np.zeros(kernsize) + (1.0 / (res * winsize) ** 4)), (np.arange(kernsize) + 0.5)
    )
    mkk *= u

    return {'kernel': mkk, 'ellkern': ellkern, 'curlyw': curlyw}


def calc_curly_w(apod_mask, delta_k, siglim=1e-9, oversampfact=1):
    '''
    Integrate the fourier transform of the `apod_mask` along theta.
    
    This is W(k) defined in Hivon et al. just after eqn A12.

    Parameters
    ----------
    apod_mask : array-like
    delta_k : float
        The resolution of `apod_mask` in k-space
    siglim : float
        significance limit. Affects what elements of curlyw are kept.
    oversampfact : int
        How much to oversample the integral

    See Also
    --------
    angle_integral
    '''
    ngrid = np.shape(apod_mask)[0]
    delta_k = delta_k / oversampfact

    if oversampfact <= 1:
        w_k = abs(np.fft.fft2(apod_mask)) / apod_mask.size / delta_k ** 2
    else:
        core.log_info(
            'Padded mask edge size: %s' % (ngrid * oversampfact), unit='CouplingKernel'
        )
        starti = int(np.floor(0.5 * (oversampfact - 1) * ngrid))
        padded = np.zeros((ngrid * oversampfact, ngrid * oversampfact))
        padded[starti : starti + ngrid, starti : starti + ngrid] = apod_mask
        w_k = (abs(np.fft.fft2(padded)) / padded.size) / delta_k ** 2
        del padded

    curlyw = angle_integral(abs(w_k) ** 2)

    curlywtot = np.cumsum(curlyw * np.arange(ngrid * oversampfact / 2))
    curlywtot /= curlywtot[int(ngrid * oversampfact / 2 - 1)]
    idx = np.where((1.0 - curlywtot) >= siglim)[0]
    wincnt = len(idx)
    if not wincnt:
        raise ValueError('Found 0 elements of curly w within threashold.  Empty mask?')
    nidx = np.setdiff1d(np.arange(len(curlywtot)), idx)

    if wincnt != ngrid / 2:
        curlyw[nidx] = 0
    core.log_info('Keeping %s elements of curly w' % wincnt, unit='CouplingKernel')
    return curlyw[0:wincnt]


def angle_integral(griddedfunc):
    '''
    Do the integral over theta used in coupling kernel calculations.
    
    Given an NxN square array, this will perform the radial integral
    of the array along theta for each of the N/2 values of radii.
    '''
    # Since each circle only intersects a small number of grid points
    # in addition to the four on the kx, ky axes, the integrals are
    # very sparse.  To get more values in the integral, this function also
    # uses off-grid points calculated from weighted averages of on-grid points.
    s = np.shape(griddedfunc)
    npix = s[0]
    result = np.zeros(int(npix / 2))
    result[0] = griddedfunc[0, 0]

    start = time.time()
    lastreport = start
    reportinterval = 60
    core.log_info("Entering angle integral", unit='CouplingKernel')
    for i in np.arange(1, int(npix / 2)):
        now = time.time()
        if (now - lastreport) > reportinterval:
            # this algorithm should be roughly quadratic in npix
            progress = (float(i) / (npix / 2)) ** 2
            remaining = (1.0 - progress) * (now - start) / progress
            remainingmin = remaining / 60
            core.log_info(
                'Estimated time remaining: %.2f min' % remainingmin,
                unit='CouplingKernel',
            )
            lastreport = now
            # report every 5 minutes from now on
            reportinterval = 300
        x, y, a = gridpts(i)
        # Where x or y are negative,
        # shift x and y to the upper half of the array, where all the
        # negative values are stored (Remember that we are assuming
        # that the gridded function is gridded like an fft, negative
        # values in the upper half of the array)
        x[x < 0] += npix
        y[y < 0] += npix
        ycnt = 0
        xcnt = 0
        xintidx = np.where((np.round(x) - x) == 0)[0]
        xcnt = len(xintidx)
        nidx = np.setdiff1d(np.arange(len(x)), xintidx)
        if len(nidx) > 0:
            yintidx = nidx[(np.round(y[nidx]) - y[nidx]) == 0]
            ycnt = len(yintidx)
        values = np.zeros(len(x))
        if (ycnt + xcnt) != len(x):
            # all points on the grid should either have an integer value
            # of X or an integer value of Y. If some points do not fall
            # in either of these categories make an error announcement
            raise ValueError('X and Y integer indices do not add up?')
        # Get the value of the gridded function at each set of coordinates.
        # If one ordinate is not an integer, use the weighted average of
        # the function at the neighboring ordinates.
        if xcnt:
            y0 = np.floor(y[xintidx]).astype(int)
            y1 = (y0 + 1) % npix
            alpha = y[xintidx] - y0
            ix = (x[xintidx]).astype(int)
            values[xintidx] = (1.0 - alpha) * griddedfunc[ix, y0] + alpha * griddedfunc[
                ix, y1
            ]
        if ycnt:
            x0 = np.floor(x[yintidx]).astype(int)
            x1 = (x0 + 1) % npix
            alpha = x[yintidx] - x0
            iy = (y[yintidx]).astype(int)
            values[yintidx] = (1.0 - alpha) * griddedfunc[x0, iy] + alpha * griddedfunc[
                x1, iy
            ]

        result[i] = scipy.integrate.simps(values, a) / 2 / np.pi

    core.log_info('Done', unit='CouplingKernel')

    return result


def gridpts(radius):
    '''
    Given an integer radius from the origin, this returns the
    x,y coordinates of points on the circle with either integer x
    or integer y, along with the angles of these points in radians.
    '''
    xord = np.arange(2 * radius + 1) - radius
    xabs = np.sqrt(radius ** 2 - xord ** 2)
    xord = np.concatenate((xord, xord[1 : 2 * radius][::-1]))
    xabs = np.concatenate((xabs, -1 * xabs[1 : 2 * radius][::-1]))

    x = np.concatenate((xabs, xord))
    y = np.concatenate((xord, xabs))
    ang = np.arctan2(y, x)

    # output is sorted
    uidx = np.unique(ang, return_index=True)[1]
    ang = ang[uidx]
    x = x[uidx]
    y = y[uidx]

    return x, y, ang
