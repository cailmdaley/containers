from spt3g import core, mapmaker, coordinateutils
from spt3g.util import fitting, framecombiner, files
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from y1_ee.null_tests import calculate_cl_chi_square

def plot_ps_with_err(data, spectrum='EE', ell_range=[300, 3000], plot_dl = False,
                     use_sigma_errors = True, n_spectra = 1, label_pte = False,
                     title = '', save_plots=False, filename='', 
                     return_pte = False, **kwargs):
    if isinstance(data, str):
        data = files.load_pickle(data)
    elif not isinstance(data, dict):
        raise TypeError("Input must be dict or path to a dict, "
                        "not %s"%type(data))
    # Figure out if these are Dls, mostly for labelling purposes
    dls = data[spectrum].spec_type == 'dl'
    # Get number of spectra used
    if 'n_spectra' in data.keys():
        n_spectra = data['n_spectra']
    # Get ell range index
    ell = data[spectrum].bin_centers
    idx = np.where((ell >= ell_range[0]) & (ell <= ell_range[1]))
    if dls and plot_dl:
        conv = 1
    if not dls and not plot_dl:
        conv = 1
    if not dls and plot_dl:
        conv = ell[idx] * (ell[idx] + 1)/(2 * np.pi)
    if dls and not plot_dl:
        conv = (2 * np.pi) / (ell[idx] * (ell[idx] + 1))

    fig, ax = plt.subplots()
    if 'error_'+spectrum in data and not use_sigma_errors:
        yerr = data['error_'+spectrum][idx]
    elif 'sigma_'+spectrum in data:
        if not use_sigma_errors:
            print('%s cov errors not found, using sigma errors'%spectrum)
        yerr = conv*data['sigma_'+spectrum][idx] * (1 / core.G3Units.uK**2)
    else:
        raise ValueError('Cannot find uncertainties on input Cls')
    ax.errorbar(ell[idx], conv*data[spectrum][idx] * (1 / core.G3Units.uK**2),
                 yerr = yerr, **kwargs)
    ax.set_xlabel('Ell')
    if plot_dl:
        ax.set_ylabel('$D_\ell$ $[\mu K ^2]$',fontsize=12)
    else:
        ax.set_ylabel('$C_\ell$ $[\mu K ^2]$',fontsize=12)
    ax.set_title(title)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(1e-2,1e4))

    if label_pte:
        chi2, dof, pte = calculate_cl_chi_square(
            data, spectrum, use_sigma_errors = use_sigma_errors)
        ax.axhline(0,linestyle='--', color='gray')
        ax.text(0.55, 0.9, '$\chi^2$ = %.1f, PTE = %.4f'%(chi2, pte),
                transform=ax.transAxes, fontsize=12)
    if save_plots:
        plt.savefig(filename)
    if label_pte and return_pte:
        return pte
    
    
def plot_spt3g_map(input_map, map_id = None, pixel_mask = None,
                   combine_map_frames = False, fig = None,
                   title = '', suptitle = ''):
    if isinstance(input_map, str):
        g3file = core.G3File(input_map)
        if combine_map_frames == True:
            fc = framecombiner.FrameCombiner(
                type = core.G3FrameType.Map, fr_id = map_id, key_op_mapping = {
                    'T' : framecombiner.fc_add,
                    'Q' : framecombiner.fc_add,
                    'U' : framecombiner.fc_add,
                    'Wpol' : framecombiner.fc_add,
                    'Wunpol' : framecombiner.fc_add})
            for frame in g3file:
                fc(frame)
            input_map = fc.output_frame
        elif map_id is not None:
            found = False
            for frame in g3file:
                if 'Id' in frame and frame['Id']==map_id:
                    input_map = frame
                    found = True
                    break
            if not found:
                print("No Map frame with Id matching %s found."%map_id)
                return
        else:
            print("Using first valid Map frame found in .g3 file")
            for frame in g3file:
                if frame.type == core.G3FrameType.Map:
                    if 'Wpol' in frame or 'Wunpol' in frame:
                        input_map = frame
                        break

    if isinstance(input_map, core.G3Frame):
        wgt = None
        if 'Wunpol' in input_map:
            wgt = input_map['Wunpol']
        elif 'Wpol' in input_map:
            wgt = input_map['Wpol']
        if input_map['T'].is_weighted:
            tmap = mapmaker.mapmakerutils.remove_weight_t(
                input_map['T'], wgt)
        else:
            tmap = input_map['T']

    if isinstance(input_map, coordinateutils.FlatSkyMap):
        if input_map.is_weighted:
            print('Input map is marked as weighted, '
                  'but plotting will proceed')
        tmap = input_map
    
    if pixel_mask is None:
        pixel_mask = np.isfinite(np.asarray(tmap))
    
    #plot in uK
    tmap /= core.G3Units.uK
    
    # Figure out an appropriate colorscale
    sdtemp = np.nanstd(np.array(tmap)[np.where(pixel_mask)])
    if sdtemp == 0 or not np.isfinite(sdtemp):
        raise ValueError(
            "Map has zero or non-finite RMS.")
    atemp = fitting.gaussfit_hist(
        np.array(tmap)[np.where(pixel_mask)], 1000,
        -5.*sdtemp, 5.*sdtemp, do_plot = False)
    
    # Now that we're done using it as an indexing bool,
    # use pixel_mask as plotting mask
    pixel_mask = np.array(pixel_mask, dtype = float)
    pixel_mask[np.where(pixel_mask == 0.)] = np.nan

    # Make the plot
    plt.figure(fig)
    ax = plt.gca()
    im = ax.imshow(tmap*pixel_mask,
                   origin='lower', cmap=plt.cm.gray,
                   interpolation = None,
                   vmin = -5.*atemp['params'][2],
                   vmax = 5.*atemp['params'][2])
    plt.xticks([])
    plt.yticks([])
    plt.title(title)
    plt.suptitle(suptitle)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=0)

    cbar = plt.colorbar(im, cax = cax)
    cbar.set_label("$\mu K_{CMB}$", labelpad = -10)
    cbar.ax.tick_params(axis='y',width=0,length=0)
    plt.tight_layout()
