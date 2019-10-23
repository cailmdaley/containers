def get_foreground_power_george_2015(component, freq1=150, freq2=None, units='uk'):
    """
    Foreground powers from George et al. 2015 results.

    Uses .sav file generated by Christain Reichardt.

    Parameters
    ----------
    component : str
        The foreground component to use. Must be one of
        'all', 'tSZ', 'kSZ', 'DG-Cl', 'DG-Po', 'RG', 'tSZ-CIB', 'Total', 'CMB'
    freq1 : int
        Frequency band. If `freq2` is specified, the cross-spectrum between
        the two frequencies will be returned. Otherwise autospectrum of freq1.
    freq2 : int, optional
        Frequency band for cross-spectrum with `freq1`
    units : str
        'k' or 'uk'. Note: default savfile is Dls in uK

    Returns
    -------
    fgnd_cls : array
        Power spectrum of `component` at specified frequency band.
    """
    components = [
        'all',
        'tSZ',
        'kSZ',
        'DG-Cl',
        'DG-Po',
        'RG',
        'tSZ-CIB',
        'Total',
        'CMB',
    ]
    if component not in components:
        raise ValueError(
            '{} not in list of possible foregrounds, must be one of {}'.format(
                component, components
            )
        )

    filename = os.path.join(
        os.path.dirname(__file__), 'data/foregrounds/george_plot_bestfit_line.sav'
    )
    from scipy.io import readsav

    data = readsav(filename)

    if freq2 is None:
        freq2 = freq1
    if freq1 == 90:
        freq1 = 95
    if freq2 == 90:
        freq2 = 95

    freqs = np.asarray(
        [(95, 95), (95, 150), (95, 220), (150, 150), (150, 220), (220, 220)]
    )
    dls_all = data['ml_dls'][(freqs[:, 0] == freq1) & (freqs[:, 1] == freq2)][0]
    labels = data['ml_dl_labels'].astype('str')
    ells = np.asarray(data['ml_l'], dtype=int)

    if component == 'all':
        spec = ells * 0.0
        for fg in components:
            if fg in ['all', 'tSZ-CIB', 'Total', 'CMB']:
                continue
            spec += dls_all[labels == fg][0]
    else:
        spec = dls_all[labels == component][0]

    # Changing Dls to Cls
    spec /= ells * (ells + 1.0) / 2.0 / np.pi
    if units.lower() == 'k':
        spec /= 1e12

    # Pad to l=0
    spec = np.concatenate((np.zeros(min(ells)), spec))

    return spec

def get_bl(beamval, els):

    fwhm_radians = np.radians(beamval/60.)
    sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
    sigma2 = sigma ** 2
    bl = np.exp(els * (els+1) * sigma2)

    return bl

def get_nl(noiseval, els, beamval, use_beam_window = 1, units='uk'): #here we will smooth noise by the output required beam

    if units.lower() == 'k': noiseval /= 1e6

    if use_beam_window:
        fwhm_radians = np.radians(beamval/60.)
        sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
        sigma2 = sigma ** 2
        bl = np.exp(els * (els+1) * sigma2)

    delta_T_radians = noiseval * np.radians(1./60.)
    nl = np.tile(delta_T_radians**2., int(max(els)) + 1 )

    nl = np.asarray( [nl[int(el)] for el in els] )

    if use_beam_window: nl *= bl

    return nl

def apply_ILC_weights(maparr, weight_arr, dx):
    ny, nx = maparr[0].shape
    mapparams = [nx, ny, dx, dx]
    maparr_weighted = []
    for mm in range(len(maparr)):
        currW = Cls2CLS(weight_arr[mm], mapparams)
        currmap = np.fft.ifft2( np.fft.fft2(maparr[mm]) * currW ).real
        maparr_weighted.append(currmap)
    MAP = np.sum(maparr_weighted, axis = 0)[0]

    return MAP

def get_Cls_mat(nuarr, beamarr, noisearr, units = 'uk'):

    Cls_dic = {}
    for mcnt1, (nuval1, beamval1, noiseval1) in enumerate(zip(nuarr, beamarr, noisearr)):
        for mcnt2, (nuval2, beamval2, noiseval2) in enumerate(zip(nuarr, beamarr, noisearr)):

            #print nuval1, nuval2
            if (nuval2, nuval1) in Cls_dic:
                Cls_dic[(nuval1, nuval2)] = Cls_dic[(nuval2, nuval1)]
                continue

            els, Cls_CMB = get_foreground_power_george_2015('CMB', freq1 = nuval1, freq2 = nuval2)
            els, Cls_tSZ = get_foreground_power_george_2015('tSZ', 1, just_Cls = 1, freq1 = nuval1, freq2 = nuval2)
            els, Cls_kSZ = get_foreground_power_george_2015('kSZ', 1, just_Cls = 1, freq1 = nuval1, freq2 = nuval2)
            els, Cls_RG = get_foreground_power_george_2015('RG', 1, just_Cls = 1, freq1 = nuval1, freq2 = nuval2)
            els, Cls_DGPo = get_foreground_power_george_2015('DG-Po', 1, just_Cls = 1, freq1 = nuval1, freq2 = nuval2)
            els, Cls_DGclus = get_foreground_power_george_2015('DG-Cl', 1, just_Cls = 1, freq1 = nuval1, freq2 = nuval2)

            Cls = Cls_CMB + Cls_tSZ + Cls_kSZ + Cls_RG + Cls_DGPo + Cls_DGclus 
            if units.lower() == 'k': Cls /= 1e12
            
            bl = get_bl(beamval2, els)

            #if autospectrum: add noise convolved with the final output beam
            if nuval1 == nuval2: 
                nl = get_nl(noiseval1, els, beamval1)
            else:
                nl = np.zeros(len(Cls))

            Cls = Cls + nl

            Cls_dic[(nuval1, nuval2)] = Cls

    return Cls_dic

def create_Clmat(elcnt):
    Clmat = np.zeros( (nc, nc) )
    for ncnt1, nuval1 in enumerate(nuarr):
        for ncnt2, nuval2 in enumerate(nuarr):
            Clmat[ncnt2, ncnt1] = Cls_dic[(nuval1, nuval2)][elcnt]
    return Clmat

def make_ilc_map(maparr, dx, nuarr, beamarr, noisearr, final_comp = 'CMB', weight_arr = None, nu_calib_fac = None):

    if weight_arr is not None:
        ilc_map = apply_ILC_weights(maparr, weight_arr)
        return ilc_map

    nchannels = len(nuarr)
    
    if final_comp.lower() == 'cmb':
        nu_scale_fac = np.ones(nc)
    if final_comp.lower() == 'comptony':
        nu_scale_fac = [ysz_Tsz_conv_fac_90, ysz_Tsz_conv_fac_150, ysz_Tsz_conv_fac_220]

    if nu_calib_fac is None:
        nu_calib_fac = np.ones( nchannels )

    acap = np.zeros(nc) + (nu_scale_fac * nu_calib_fac) #assuming CMB is the same and calibrations factors are same for all channels
    acap = np.mat(acap).T #should be nc x 1

    #get weights in harmonic space
    weight_arr = np.zeros( (nc, len(els)) )
    for elcnt, el in enumerate(els):
        Clmat = np.mat( create_Clmat(elcnt) )
        Clinv = sc.linalg.pinv2(Clmat)
        
        nr = np.dot(Clinv, acap)
        dr = np.dot( acap.T, np.dot(Clinv, acap) )
        weight_mat = np.asarray(nr/dr).squeeze()

        weight_arr[:, elcnt] = weight_mat

    if debug: #plot weights
        plot(weight_arr[0], 'r', label = r'\textbf{90}'); plot(weight_arr[1], 'k' , label = r'\textbf{150}'); plot(weight_arr[2], 'g', label = r'\textbf{220}'); plot(weight_arr[3], 'm', label = r'\textbf{270}')
        ylim(-5.,5);legend(loc=1);show()

    ilc_map = apply_ILC_weights(maparr, weight_arr, dx)

    return ilc_map

