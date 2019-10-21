# calvel_tools.py
# Tools useful in calibrator response vs. elevation analaysis.
#
# Angelina Harke-Hosemann <aharkehosemann@gmail.com>
# November 2017
#
#
# Analysis workflow ususally goes : 
# - make list of observation IDs for cal vs el scan
# - get data dictionary from get_data_dict(obs, norm_ob)   [norm_ob usually chosen to be the middle index]
# - maybe make histogram of delta response with hist_delres()
# - maybe plot delta response on the focal plane with delres_fullplane()




import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import os, glob
import cPickle as pkl
from spt3g import core, mapmaker, std_processing
import matplotlib.gridspec as gridspec
import pdb


C = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']   # matplotlib 2.0 default colors


def plot_calvel(data_dict, legend=False, title='Calibrator Response vs Elevation', xlim=(43, 67), ylim=(0,0.0043)):
    
    '''
    Plots Cal v El Data, probably pulled from get_cal_v_el.
    
    INPUT
    bolos: array of bolometer names of interest
    cals: array of bolometer responses, indexed cals[bolometer #, obs_ID]
    els: array of elevation for response measurements, indexed els[bolometer #, obs_ID]
    
    RETURNS
    None, inline plots the dang plot.
    '''
    
    bolos = data_dict['data'].keys()

    calvel_fig = plt.figure(figsize=(6,4))
    for bb, bolo in enumerate(bolos):
        plt.plot(data_dict['scan']['els'], data_dict['data'][bolo]['cal_res'], label=bolo, color=C[bb%9], alpha=.3)
    if legend == True:
        plt.legend(loc='upper left', fontsize='x-small', framealpha=.8)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(title)
    plt.xlabel('Elevation [deg]')
    plt.ylabel('Bolometer Response [pW]')

    plt.tight_layout()


def get_normd_cals(cals, norm_ind, plot=False):
    
    '''
    Gives normalized calibration responses.
    
    INPUT
    cals: array of unnormalized calibration responses. format: cals[observation, bolometer]
    ob: observation response to normalize to.
    
    RETURNS
    normd_cals: normalized calibration responses
    if plot=True, plots the normalized responses of first 5 bolos.
    '''
    
    norm = np.array(cals[norm_ind])
    norm = np.tile(norm, [len(cals), 1]).T
    normd_cals = cals/norm
    normd_cals = normd_cals[0]   # dumb hack
    
    ''' Normalized Plot '''
    if plot == True:
        for bb, bolo in enumerate(bolos[0:4]):
            plt.plot(els[bb], normd_cals[bb], label=bolo)
        plt.title('Calibrator Response vs Elevation')
        plt.xlabel('Elevation [deg]')
        plt.ylabel("Bolometer Response [Norm'd at El = 55 deg]")
        plt.legend(loc='lower right', fontsize='x-small')
        plt.xlim(min(els[0,:])-2., max(els[0,:])+2.)
    
    return normd_cals


def sn_cut(data_dict, cut, ob):
    
    '''
    Makes SN cut, returns indices of uncut bolos.
    
    INPUTS:
    bolos: all bolometers in scan
    cut: cut threshold in sigma
    ob: observation ID used to find calibration file
    
    RETURNS:
    abcut: array of indices of bolometers that survived the cut  
    '''
    
    bolos = data_dict['data'].keys()
    abcut = np.array([bolo for bolo in bolos if data_dict['data'][bolo]['res_sn'] > cut])
    sn_ab = [data_dict['data'][bolo]['res_sn'] for abolo in abcut]

    print 'Mean SN for bolos that survived cut: ', np.mean(sn_ab)
    return abcut


def get_delres(normd_cals, bolos, els, plot=False):

    # delta res = (norm cal response at HIGHEST el - norm cal response at LOWEST el)
    del_res = np.empty(len(bolos))
    for bb, bolo in enumerate(bolos):
        del_res[bb] = np.abs(normd_cals[bb,-1] - normd_cals[bb,0])
        
    return del_res


def hist_delres(data_dict, plot_type='all_bolos', **kwargs):
    

    '''
    Plot histograms of delta response. Shows total population and break down by band. 
    Optional breakdowns : all_bolos, by_wafer

    INPUT
    data_dict : output of get_data_dict()
    plot_type [optional] : breakdown population, default all_bolos
    kwargs : range, bins, ylim, xlabel, title

    OUTPUT
    returns histogram figure
    '''
    
    # check plot type
    plot_types = ['all_bolos', 'by_wafer']
    
    if plot_type not in plot_types:
        print "plot_type not supported. Available plot_type's: ", plot_types
        return

    
    # useful things to know
    bolos = data_dict['data'].keys()
    els = data_dict['scan']['els']

    # 3G has the following band centers:
    bands = ['90.0', '150.0', '220.0']

    if plot_type == 'all_bolos':
        
        nrows, ncols = 1, 1
        
        delres = {}
        delres['all'] = {}

        dres = [data_dict['data'][bolo]['del_res'] for bolo in bolos]
        delres['all']['all'] = dres

        # sort by band
        for band in bands:
            bandbolos = [bolo for bolo in bolos if 'band' in data_dict['data'][bolo] and str(data_dict['data'][bolo]['band']) == band]
            dres_band = [data_dict['data'][bbolo]['del_res'] for bbolo in bandbolos]
            delres['all'][band] = dres_band

        delres = delres   
        #titles = [str(band)+' GHz' for band in bands]
    
        empty = []
        titles = ['']
        fsize = (4,3)

        
    if plot_type == 'by_wafer':
        
        # find wafers 
        wafers=[]
        for bolo in bolos:
            wafers.append(data_dict['data'][bolo]['wafer'])
        wafers = np.array(list(set(wafers)))

        nrows, ncols = int(np.ceil(len(wafers)/3.)), 3

        # sort delres by wafer                                                        
        delres_bywafer = {wafer: dict() for wafer in wafers}
        for ww, wafer in enumerate(wafers):
    
            # find all in wafer
            delres_bywafer[wafer] = {band:[] for band in bands}
            waferbolos = [bolo for bolo in bolos if data_dict['data'][bolo]['wafer'] == wafer]
            delres_bywafer[wafer]['all'] = [data_dict['data'][wbolo]['del_res'] for wbolo in waferbolos]

            # separate by band
            for band in bands:
                bandbolos = [wbolo for wbolo in waferbolos if 'band' in data_dict['data'][wbolo] and str(data_dict['data'][wbolo]['band']) == band]
                dres_band = [data_dict['data'][bbolo]['del_res'] for bbolo in bandbolos]
                delres_bywafer[wafer][band] = dres_band

        delres = delres_bywafer                
        titles = [wafer for wafer in wafers]
        fsize = (3*ncols, 3*nrows)

        if len(wafers) % 3 == 1:   # one subplot in last row, two left over
            empty = [-1, -2]
        elif len(wafers) % 3 == 2:   # two subplots in last row, one left over
            empty = [-1]
        else:
            empty = []


    # initialize figure
    hgrid = gridspec.GridSpec(nrows,ncols)
    hfig = plt.figure(figsize=fsize)
    #hfig=plt.figure(figsize=(4,3))
    
    axes = []
    for rr in np.arange(nrows):
        for cc in np.arange(ncols): 
            axes = np.append(axes, plt.subplot(hgrid[rr,cc]))
            
    for emp in empty:
        axes[emp].set_axis_off()    
    
    # read in user-defined plot parameters
    range = kwargs.pop('range', None)
    bins = kwargs.pop('bins', 30)
    ylim = kwargs.pop('ylim', None)
    xlabel = kwargs.pop('xlabel', "Delta Response [Norm'd at 55 deg]")

    stitle = 'Observation on ' + str(data_dict['scan']['date'])
    stitle = kwargs.pop('title', stitle)
    
    # plot histogram
    
    for kk, key in enumerate(delres.keys()):
        ax = axes[kk]
        dres_all = delres[key]['all']
        dres_all = [dres for dres in dres_all if not math.isnan(dres)]   # handle nans
        ax.hist(dres_all, bins=bins, color=C[3], range=range, alpha=0.8, histtype='step')
        
        
        # by band
        for bb, band in enumerate(bands):
            dres_band = delres[key][band]
            dres_band = [dres for dres in dres_band if not math.isnan(dres)]   # handle nans
            ax.hist(dres_band, bins=bins,color=C[bb], range=range, alpha=0.6, histtype='stepfilled', label=band, edgecolor='none')
        
        ax.legend(fontsize='small')
        ax.set_xlabel(xlabel)
        ax.set_title(titles[kk])
        ax.set_ylim(ylim)
        ax.set_xlim(range)

    hfig.suptitle(stitle, fontsize=14, y=1.01)
    hfig.tight_layout()
    
    return hfig


def delres_fullplane(data_dict, dres_cutoff, sn_cut=-1, **kwargs):

    """
  
    A function for plotting delta response as a colormap on the SPT-3G focal plane. 
    Outputs a figure with six subplots corresponding each SPT-3G detector band and polarization. 
    Takes delta response caculations for every bolometer, employs an optional SN cut (user-defined theshold),
    then makes the plots.

    Assumes 3G bands and polarizations, S/N cut is desired only at one elevation, that the user thinks
    plasma colormap is cool, probably other things. 

    INPUT
    data_dict : output of get_data_dict()
    dres_cutoff : at what delta response to cap the colorbar (outliers can easily dominate colormap)
    sn_cut [optional] : at what sigma S/N to cut bolometers (exclusive)
    kwargs : title
    
    RETURNS
    figure of colormapped focal plane for each band, polarization

    """

    # SPT-3G bands and polarizations
    bands = ['90.0', '150.0', '220.0']
    pols = ['x', 'y']

    bolos = [bolo for bolo in data_dict['data'].keys() if 'res_sn' in data_dict['data'][bolo] and data_dict['data'][bolo]['res_sn'] > sn_cut]

    ### initialize figure
    nrows, ncols = 3, 2
    planefig = plt.figure()
    planefig.set_size_inches(10,10)
    planegrid = gridspec.GridSpec(nrows, ncols+1, width_ratios=[25]*ncols+[1])
    
    ### get colormap
    cmap = plt.get_cmap('plasma_r') 
    cbar_norm = matplotlib.colors.Normalize(vmin=0, vmax=dres_cutoff)       

    ### make plots
    for bnd, band in enumerate(bands):
        bbolos = [bolo for bolo in bolos if 'band' in data_dict['data'][bolo].keys() and str(data_dict['data'][bolo]['band']) == band]
        for pl, pol in enumerate(pols): 
            pbolos = [bolo for bolo in bbolos if 'pol' in data_dict['data'][bolo].keys() and data_dict['data'][bolo]['pol'] == pol]

            ax = planefig.add_subplot(planegrid[bnd, pl])
            ax.set_title(str(band) + ' GHz, ' + str(pol) + ' pol')

            for bolo in pbolos:
                x_pos = data_dict['data'][bolo]['x_offset']
                y_pos = data_dict['data'][bolo]['y_offset']

                c = data_dict['data'][bolo]['del_res'] / dres_cutoff
                color = cmap(c)

                ax.plot(x_pos, y_pos, '.', color=color)
                ax.axis('off')

    ### colorbar
    for row in np.arange(nrows):           
        cax = planefig.add_subplot(planegrid[row, ncols])
        cbase = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=cbar_norm)

    title = kwargs.pop('title', data_dict['scan']['date'])
    planefig.suptitle(title, fontsize=13)
    planefig.tight_layout()
    plt.subplots_adjust(left=0.12, bottom=0.08, right=0.85, top=0.93, wspace=0.01, hspace=0.2)
    
    return planefig


def get_data_dict(obs, norm_ind):
    
    """
    Function designed to do all file I/O at once. Creates the data dictionary used by other functions 
    in this analysis. 

    INPUT
    obs : list of observation ID's in cal vs el scan. 
    norm_ind : index of observation in obs with which to normalize responses. 

    OUTPUT
    dictionary of data and scan information
    data_dict['data'] is indexed by bolometer, with properties : cal_res, normd_res, del_res, timestream, res_sn, band, pol, x_offset, y_offset, wafer
    data_dict['scan'] has properties : els, date
    """


    data_dict = dict()
    data_dict['data'] = dict()
    data_dict['scan'] = dict()
    data_dict['scan']['els'] = np.array([])
    
    els = []

    for ob in obs:
        
        obsdir = '/spt/data/bolodata/downsampled/calibrator/' + str(ob) + '/'
        calfile = obsdir + '/offline_calibration.g3'
        scanfile = obsdir + '0000.g3'
        onlinefile = obsdir + '/nominal_online_cal.g3'

        # Cal has Calibrator Data, Scan has El Data                                      
        calibrator = list(core.G3File(calfile))
        scan = list(core.G3File(scanfile))
        online = list(core.G3File(onlinefile))

        # get bolos
        if len(data_dict['data']) == 0:
            bolos = scan[-1]['RawTimestreams_I'].keys()
            for bolo in bolos:
                data_dict['data'][bolo] = dict()
                data_dict['data'][bolo]['cal_res'] = np.array([])
                data_dict['data'][bolo]['timestream'] = np.array([])

        # get scan info
        data_dict['scan']['els'] = np.append(data_dict['scan']['els'], scan[-1]['OnlineBoresightEl'][2000]/core.G3Units.deg)
        data_dict['scan']['date'] = str(scan[0]['ObservationStart']).split(':')[0]   # gets overwritten every ob

        # get cals
        for bolo in bolos:
            cal_res = calibrator[0]['CalibratorResponse'][bolo]/core.G3Units.pW
            data_dict['data'][bolo]['cal_res'] = np.append(data_dict['data'][bolo]['cal_res'], cal_res)
            timestream = scan[-1]['RawTimestreams_I'][bolo][10:-10]
            data_dict['data'][bolo]['timestream'] = np.append(data_dict['data'][bolo]['timestream'], timestream)


    for bolo in bolos:
        cals = data_dict['data'][bolo]['cal_res']
        data_dict['data'][bolo]['normd_res'] = np.array(get_normd_cals(cals, norm_ind))
        data_dict['data'][bolo]['del_res'] = data_dict['data'][bolo]['normd_res'][-1] - data_dict['data'][bolo]['normd_res'][0]
        
        # get bolo info
        if bolo in calibrator[0]['BolometerProperties'].keys():
            data_dict['data'][bolo]['res_sn'] = calibrator[0]['CalibratorResponseSN'][bolo]
            data_dict['data'][bolo]['band'] = calibrator[0]['BolometerProperties'][bolo].band/core.G3Units.GHz
        data_dict['data'][bolo]['pol'] = online[0]['NominalBolometerProperties'][bolo].physical_name.split('.')[-1]
        data_dict['data'][bolo]['x_offset'] = online[0]['NominalBolometerProperties'][bolo].x_offset
        data_dict['data'][bolo]['y_offset'] = online[0]['NominalBolometerProperties'][bolo].y_offset
        data_dict['data'][bolo]['wafer'] = online[0]['NominalBolometerProperties'][bolo].wafer_id

        
    return data_dict
