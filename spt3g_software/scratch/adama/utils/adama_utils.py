import numpy as np
import scipy.constants as consts
import scipy.signal as ss
from spt3g import core, dfmux, todfilter, calibration
import inspect
import types
import matplotlib.pyplot as plt
from glob import glob
import os.path
import fnmatch
from datetime import datetime
from spt3g.std_processing import obsid_to_g3time, time_to_obsid

# typical values taken from Z. Pan's FTS summary page:
# https://spt-trac.grid.uchicago.edu/trac_south_pole/wiki/Spt3gFTS
bandwidths_3g = {90.0*core.G3Units.GHz: 25.5e9,
                 150.0*core.G3Units.GHz: 31.6e9,
                 220.0*core.G3Units.GHz: 47.9e9}
bandcenters_3g = {90.0*core.G3Units.GHz: 94e9,
                  150.0*core.G3Units.GHz: 146e9,
                  220.0*core.G3Units.GHz: 218e9}

def wattsPerKcmb(band):
    x = consts.h * bandcenters_3g[band] / consts.k / 2.7
    return consts.k * bandwidths_3g[band] * x**2 * np.exp(x) / ((np.exp(x) - 1)**2) * (core.G3Units.watt / core.G3Units.K)

def wattsPerKrj(band):
    return consts.k * bandwidths_3g[band] * (core.G3Units.watt / core.G3Units.K)

def get_obsids_from_times(start_time, end_time, sources, downsampled=True):
    '''
    Get obsids between start and end times for a particular source or set of
    sources.

    Parameters
    ----------
    start_time : G3Time object
        Start of the desired time range for which to return obsids
    end_time : G3Time object
        End of the time range
    source : string or list of strings
        Sources for which to get obsids
    '''
    data_path = '/spt/data/bolodata/'
    if downsampled:
        data_path = os.path.join(data_path, 'downsampled')
    else:
        data_path = os.path.join(data_path, 'fullrate')
        
    if isinstance(sources, list) == False:
        sources = [sources]

    obsids = {}
    for source in sources:
        source_path = os.path.join(data_path, source)
        obsids_all = [int(os.path.basename(path)) for path in glob(os.path.join(source_path, '*'))]
        obsids[source] = [obsid for obsid in obsids_all \
                          if obsid_to_g3time(obsid) > start_time and \
                          obsid_to_g3time(obsid) < end_time]
        obsids[source] = np.sort(obsids[source])
    return obsids
        
def extract_array(value, source):
    '''
    Extra

    Parameters
    ----------
    value: str or func
        Two possibilities here:
            1.) A key to an entry in a frame, which you wish to extract
            2.) A function of the form f(key1, key2, ..., keyN), where each
            'keyj' are the names of fields in a G3Frame you wish to analyze,
            and are assumed to be floats by the function f. In other words,
            f is used to calculate a number, per bolometer, using each of the
            keys in the argument list.
    source: str or G3Frame
        Two possibilities here:
            1.) Name of a G3File to load
            2.) Specific G3Frame object
            
    Returns
    -------
    data: numpy array
        Array of data extracted from the g3 file or g3 frame
    '''
    if type(value) == str:
        argnames = [value]
    elif type(value) == types.FunctionType:
        argnames = inspect.getargspec(value).args
        fvalue = value

    data = {}
    if type(source) == str:
        # If supplied a filename instead of a frame, then open the file and
        # loop through the frames until we encounter a frame with all the
        # required keys.
        ff = core.G3File(source)
        fr = ff.next()
        while fr.type != core.G3FrameType.EndProcessing:
            if np.all([k in fr.keys() for k in argnames]):
                source = fr
                break
            fr = ff.next()
        if fr.type == core.G3FrameType.EndProcessing:
            print('Could not find frame with all keys!')
            return None

    if type(source) == core.G3Frame:
        bololist_intersect = []
        for arg in argnames:
            bololist = [bolo for bolo in source[arg].keys()]
            if bololist_intersect:
                bololist_intersect = bololist
            else:
                bololist_intersect = np.intersect1d(bololist, bololist_intersect)
        for bolo in bololist:
            arglist = {arg: source[arg][bolo] for arg in argnames}
            if type(value) == str:
                data[bolo] = source[value][bolo]
            else:
                data[bolo] = fvalue(**arglist)

    return data

def focal_plane_plot(data_dict, band, pol, vrange=None,
                     nominal_bolo_props=None):
    '''
    Parameters
    ----------
    data_dict : dict-like
        A dictionary-like object keyed by bolometer name, which has a single
        float for each key.
    band : float
        90, 150, or 220
    pol : str
        Nominal polarization state ('x' or 'y')
    vrange : 2-element list
        [vmin, vmax] where vmin and vmax are the min and max limits of the color
        scale
    nominal_bolo_props : map of BolometerProperties objects
        The thing called "NominalBolometerProperties" stored in a
        nominal_online_cal.g3 file. If None, then we use properties from an
        arbitrary observation.

    Returns
    -------
    None
    '''
    if nominal_bolo_props:
        bp = nominal_bolo_props
    else:
        bp_filename = '/spt/data/bolodata/fullrate/ra0hdec-59.75/36844428/nominal_online_cal.g3'
        bp = [fr for fr in core.G3File(bp_filename)][0]["NominalBolometerProperties"]

    values   = np.array([data_dict[bolo] for bolo in data_dict.keys()
                         if bp[bolo].band / core.G3Units.GHz == band and 
                            bp[bolo].physical_name.split('.')[-1] == pol])
    x_offset = np.array([bp[bolo].x_offset / core.G3Units.degree for bolo in data_dict.keys()
                         if bp[bolo].band / core.G3Units.GHz == band and 
                            bp[bolo].physical_name.split('.')[-1] == pol])
    y_offset = np.array([bp[bolo].y_offset / core.G3Units.degree for bolo in data_dict.keys()
                         if bp[bolo].band / core.G3Units.GHz == band and 
                            bp[bolo].physical_name.split('.')[-1] == pol])

    plt.figure(figsize=(12,10))
    scatter_args = {'x': x_offset, 'y': y_offset, 'c': values}
    if vrange:
        scatter_args['vmin'] = vrange[0]
        scatter_args['vmax'] = vrange[1]
    plt.scatter(**scatter_args)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.xlabel('x offset [deg]')
    plt.ylabel('y offset [deg]')


def plot_noise(frame, boloprops, obsid,
               noise_type='NET',
               band='10.0Hz_to_15.0Hz', wafer=None, bywafer=False,
               filestub=None):
    units   = {'NEI': core.G3Units.amp*1e-12 / np.sqrt(core.G3Units.Hz),
               'NET': core.G3Units.microkelvin * np.sqrt(core.G3Units.sec),
               'NEP': core.G3Units.attowatt / np.sqrt(core.G3Units.Hz)}
    bins    = {'NEI': np.linspace(0, 100, 60),
               'NET': np.linspace(0, 3000, 60),
               'NEP': np.linspace(0, 100, 60)}
    labels  = {'NEI': 'NEI [pA/sqrt(Hz)]',
               'NET': 'NET [uK rtsec]',
               'NEP': 'NEP [aW/sqrt(Hz)]'}
    wafer_list = np.sort(np.unique([boloprops[bolo].wafer_id for bolo in boloprops.keys()]))

    noise_key = '{}_{}'.format(noise_type, band)

    data = {}
    if bywafer:
        plt.figure(figsize=(20,8))
        for jwafer, wafer in enumerate(wafer_list):
            plt.subplot(2,5,jwafer+1)
            data[wafer] = {}
            for obs_band in [90, 150, 220]:
                noise = np.array([frame[noise_key][bolo]
                                  for bolo in frame[noise_key].keys()
                                  if boloprops[bolo].band / core.G3Units.GHz == obs_band and \
                                  boloprops[bolo].wafer_id == wafer]) / units[noise_type]
                plt.hist(noise[np.isfinite(noise)],
                         bins=bins[noise_type],
                         histtype='step',
                         label='{} GHz'.format(obs_band)) 
                data[wafer][obs_band] = noise[np.isfinite(noise)]
            plt.legend()
            plt.title('observation {}: {}'.format(obsid, wafer))
            plt.xlabel(labels[noise_type])
            plt.xlim([np.min(bins[noise_type]), np.max(bins[noise_type])])
    
    else:
        plt.figure(figsize=(8,6))
        for obs_band in [90, 150, 220]:
            if wafer is not None:
                noise = np.array([frame[noise_key][bolo]
                                  for bolo in frame[noise_key].keys()
                                  if boloprops[bolo].band / core.G3Units.GHz == obs_band and
                                  boloprops[bolo].wafer_id == wafer]) / units[noise_type]
                data[obs_band] = noise[np.isfinite(noise)]
            else:
                noise = np.array([frame[noise_key][bolo]
                                 for bolo in frame[noise_key].keys()
                                 if boloprops[bolo].band / core.G3Units.GHz == obs_band]) / units[noise_type]
                data[obs_band] = noise[np.isfinite(noise)]
            plt.hist(noise[np.isfinite(noise)],
                     bins=bins[noise_type],
                     histtype='step',
                     label='{} GHz'.format(obs_band)) 
        plt.legend()
        band_label = band.replace('_', ' ')
        plt.title('observation {}: {}'.format(obsid, band_label))
        plt.xlabel(labels[noise_type])
        plt.xlim([np.min(bins[noise_type]), np.max(bins[noise_type])])
    plt.tight_layout()
    if filestub != None:
        plt.savefig('{}_{}.png'.format(filestub, obsid), dpi=200)

    return data


def plot_noise_comparison(frames, boloprops, obsids,
                          noise_type='NET',
                          band='10.0Hz_to_15.0Hz', bywafer=False,
                          legend_labels=None, filestub=None):
    units   = {'NEI': core.G3Units.amp*1e-12 / np.sqrt(core.G3Units.Hz),
               'NET': core.G3Units.microkelvin * np.sqrt(core.G3Units.sec),
               'NEP': core.G3Units.attowatt / np.sqrt(core.G3Units.Hz)}
    bins    = {'NEI': np.linspace(0, 100, 60),
               'NET': np.linspace(0, 3000, 60),
               'NEP': np.linspace(0, 100, 60)}
    bins_byband = {'NET': {90*core.G3Units.GHz: np.linspace(0, 1500, 60),
                           150*core.G3Units.GHz: np.linspace(0, 1000, 60),
                           220*core.G3Units.GHz: np.linspace(0, 3000, 60)}}
    labels  = {'NEI': 'NEI [pA/sqrt(Hz)]',
               'NET': 'NET [uK rtsec]',
               'NEP': 'NEP [aW/sqrt(Hz)]'}
    wafer_list = np.sort(np.unique([boloprops[0][bolo].wafer_id for bolo in boloprops[0].keys()]))
    noise_key = '{}_{}'.format(noise_type, band)
    linelist = ['-', '--', ':', '-.']

    if bywafer:
        for jband, obs_band in enumerate([90, 150, 220]):
            plt.figure(figsize=(20,8))
            for jwafer, wafer in enumerate(wafer_list):
                plt.subplot(2,5,jwafer+1)
                for jframe, (frame, bps, obsid) in enumerate(zip(frames, boloprops, obsids)):
                    noise = np.array([frame[noise_key][bolo]
                                      for bolo in frame[noise_key].keys()
                                      if bps[bolo].band / core.G3Units.GHz == obs_band and \
                                      bps[bolo].wafer_id == wafer]) / units[noise_type]
                    if labels:
                        legend_label = legend_labels[jframe]
                    else:
                        legend_label = '{}'.format(obsid)
                    plt.hist(noise[np.isfinite(noise)],
                             bins=bins_byband[noise_type][obs_band*core.G3Units.GHz],
                             histtype='step',
                             label=legend_label,
                             color='C{}'.format(jband),
                             linestyle=linelist[jframe])
                plt.legend()
                plt.title('{}: {} GHz'.format(wafer, obs_band))
                plt.xlabel(labels[noise_type])
                plt.xlim([np.min(bins_byband[noise_type][obs_band*core.G3Units.GHz]),
                          np.max(bins_byband[noise_type][obs_band*core.G3Units.GHz])])
            plt.tight_layout()
            if filestub != None:
                plt.savefig('{}_{}_{}.png'.format(filestub,
                                                  '_'.join([str(oid) for oid in obsids]),
                                                  obs_band), dpi=200)
    #else: # something...
    plt.tight_layout()
    

    

def plot_noise_from_file(fname, cal_fname, noise_type='NET',
                         band='10.0Hz_to_15.0Hz'):
    d = [fr for fr in core.G3File(fname)][0]
    dcal = core.G3File(cal_fname)
    bp = dcal.next()['BolometerProperties']
    obsid = fname.split('/')[-1].split('.')[0]
    
    if len(d) > 0:
        plot_noise(d, bp, obsid,
                   noise_type=noise_type,
                   band=band)


def plot_autoproc_data(obsid, field, xlims, units, plot_name, xlabel,
                       nbins=100, bololist=None):
    '''
    Plots a histogram of a quantity in an autoprocessing output file, over
    a specific list of bolometers, if requested.

    Parameters
    ----------
    obsid : int
        Observation ID
    field : str
        Name of the quantity to plot
    xlims : 2-tuple
        Min and max range of histogram
    units : core.G3Units object
        The units in which you would like the data plotted
    plot_name : str
        File name of plot to save
    xlabel : str
        X axis label
    nbins : int
        Number of bins to include in the histogram
    bolos : array of str
        Names of bolometers to include in histogram

    Returns
    -------
    None
    '''
    # load autoproc data
    file_path = '/spt/user/production/calibration/*/{}.g3'.format(obsid)
    fname = glob(file_path)
    if len(fname)!=1 or os.path.exists(fname[0])==False:
        print('ERROR: Cannot find file.')
        return
    fname_to_read = fname[0]
    d = [fr for fr in core.G3File(fname_to_read)][0]

    # load boloprops
    bp_path = '/spt/data/bolodata/downsampled/*/{}/nominal_online_cal.g3'.format(obsid)
    fname = glob(bp_path)
    if len(fname)!=1 or os.path.exists(fname[0])==False:
        print('ERROR: Cannot find file {}.'.format(fname[0]))
        return
    fname_to_read = fname[0]    
    bp = [fr for fr in core.G3File(fname_to_read)][0]['NominalBolometerProperties']

    plt.figure(figsize=(8,6))
    for band in [90, 150, 220]:
        if bololist is not None:
            plot_data = np.array([d[field][bolo] for bolo in d[field].keys() \
                                 if bolo in bololist and \
                                 bp[bolo].band / core.G3Units.GHz == band])
        else:
            plot_data = np.array([d[field][bolo] for bolo in d[field].keys() \
                                 if bp[bolo].band / core.G3Units.GHz == band])

        plt.hist(plot_data[np.isfinite(plot_data)] / units,
                 bins=np.linspace(xlims[0], xlims[1], nbins),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.xlim(xlims)
    plt.xlabel(xlabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_name, dpi=200)


def get_scan_bounds(obsid, source, filenames=['0*.g3'], fullrate=True):
    '''
    Utilty function for getting times of scan boundaries and whether the scans
    correspond to turnaround frames.

    Parameters:
    -----------
    obsid : str, int
        Observation ID for which to get data.
    source : str
        Name of source of data.
    filenames : arr of str
        File names (or globbable string) in observation ID for which to get
        data.
    fullrate : str
        Use fullrate or downsampled data? Shouldn't matter in this case because
        the times of scan boundaries should not care about whether data is
        downsampled or not.

    Outputs:
    --------
    time_bounds : numpy array of datetime
        Times of the boundaries between scan frames.
    turnarounds : numpy array of bool
        Whether the scan frames corresponding to time_bounds are turnaround
        frames or not.
    '''
    if fullrate:
        data_basepath = '/spt/data/bolodata/fullrate/'
    else:
        data_basepath = '/spt/data/bolodata/downsampled/'

    data_fullpaths = []
    for fname in filenames:
        data_fullpaths.extend(glob(os.path.join(data_basepath, source, obsid, fname)))

    pipe = core.G3Pipeline()
    scan_bounds = ScanBounds()

    pipe.Add(core.G3Reader, filename=data_fullpaths)
    pipe.Add(scan_bounds)
    pipe.Run()

    return scan_bounds.time_bounds, scan_bounds.turnarounds

def get_raw_timestreams(bolos, obsid, file_name='0000.g3', scan_num=None,
                        fullrate=True, n_bolos_per_fig=1,
                        plot=True, data=False, cut_turnarounds=True,
                        psd=False, units=None, phase='I', timerange=None,
                        freqrange=None, plot_name_stub=''):
    if fullrate:
        file_path = '/spt/data/bolodata/fullrate/*/{}/{}'.format(obsid, file_name)
        boloprops_path = '/spt/data/bolodata/fullrate/*/{}/nominal_online_cal.g3'.format(obsid)
    else:
        file_path = '/spt/data/bolodata/downsampled/*/{}/{}'.format(obsid, file_name)
        boloprops_path = '/spt/data/bolodata/downsampled/*/{}/nominal_online_cal.g3'.format(obsid)
    fname = glob(file_path)
    bp_fname = glob(boloprops_path)

    # check for existence of unique data file
    if len(fname)!=1 or os.path.exists(fname[0])==False:
        print('ERROR: Cannot find file {}.'.format(fname[0]))
        return
    fnames_to_read = [fname[0]]

    # check for existence of unique file containing bolometer properties
    if len(bp_fname)!=1 or os.path.exists(bp_fname[0])==False:
        print('ERROR: Cannot find file {}.'.format(boloprops_path[0]))
        return
    bp_fname_to_read = bp_fname[0]
    
    if units == core.G3TimestreamUnits.Tcmb:
        calfile_path = '/spt/user/production/calibration/calframe/*/{}.g3'.format(obsid)
        calfname = glob(calfile_path)

        if len(calfname)!=1 or os.path.exists(calfname[0])==False:
            print('ERROR: Cannot find file {}.'.format(calfname[0]))
            return

        fnames_to_read.insert(0, calfname[0])

    # convert `bolos` to single-element list if it is a string
    if type(bolos) == str:
        bolos = [bolos]

    # do wildcard search with each element of `bolos` argument over all 
    # bolometer names
    bolonames = [fr for fr in core.G3File(bp_fname_to_read)][0]['NominalBolometerProperties'].keys()
    bolonames_to_grab = []
    for bname in bolos:
        bolonames_to_grab.extend(fnmatch.filter(bolonames, bname))

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=fnames_to_read)
    if cut_turnarounds:
        pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])
    data_key = 'RawTimestreams_' + phase
    if scan_num != None:
        pipe.Add(UseScanNums, frame_list=scan_num, frame_type=core.G3FrameType.Scan)
    if units == core.G3TimestreamUnits.Power:
        pipe.Add(dfmux.ConvertTimestreamUnits, Input=data_key,
                 Output='TimestreamsWatts_' + phase, Units=core.G3TimestreamUnits.Power)
        data_key = 'TimestreamsWatts_' + phase
    elif units == core.G3TimestreamUnits.Tcmb:
        pipe.Add(dfmux.ConvertTimestreamUnits, Input=data_key,
                 Output='TimestreamsWatts_' + phase, Units=core.G3TimestreamUnits.Power)
        pipe.Add(calibration.ApplyTCalibration, Input='TimestreamsWatts_' + phase, 
                 Output='TimestreamsKcmb', OpacityCorrection=False)
        data_key = 'TimestreamsKcmb'
    if data:
        accumulator = AccumulateTimestreams(bololist=bolonames_to_grab, Input=data_key)
        pipe.Add(core.Dump)
        pipe.Add(accumulator)
    if plot:
        plotter = PlotTimestreams(bololist=bolonames_to_grab,
                                  n_bolos_per_fig=n_bolos_per_fig,
                                  Input=data_key,
                                  psd=psd,
                                  timerange=timerange,
                                  freqrange=freqrange,
                                  name_stub=plot_name_stub,
                                  units=units)
        pipe.Add(plotter)
    pipe.Run()

    # extract the stored timestreams from the accumulator
    if data:
        return accumulator.saved_ts
    

class PlotTimestreams(object):
    def __init__(self, bololist=[], n_bolos_per_fig=1,
                 Input='RawTimestreams_I', psd=False,
                 timerange=None, freqrange=None, name_stub='',
                 units=None, units_factor=None):
        self.bololist = bololist
        self.obsid = None
        self.jscan = 0
        self.n_bolos_per_fig = n_bolos_per_fig
        self.Input = Input
        self.psd = psd
        self.timerange = timerange
        self.freqrange = freqrange
        self.name_stub = name_stub
        self.units = units
        self.units_factor = units_factor

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Observation:
            self.obsid = frame["ObservationID"]
        if frame.type == core.G3FrameType.Scan:
            if self.psd:
                psd, f = todfilter.dftutils.get_psd_of_ts_map(frame[self.Input])

            jbolo = 0
            for bolo in self.bololist:
                if jbolo == self.n_bolos_per_fig:
                    plt.close()
                    plt.figure()
                    jbolo = 0
                if not self.psd:
                    try:
                        times = np.array([ts.time / core.G3Units.second for ts in frame[self.Input][bolo].times()])
                        times = times - np.min(times)
                        if self.timerange != None:
                            cTime = (times > self.timerange[0]) & (times < self.timerange[1])
                            times = times[cTime]
                            timestream = frame[self.Input][bolo][cTime]
                        else:
                            timestream = frame[self.Input][bolo]

                        # deal with units
                        if self.units_factor is None:
                            if self.units is core.G3TimestreamUnits.Power:
                                timestream_plot = timestream / (core.G3Units.watt * 1e-12)
                                ylabel = 'power [pW]'
                            elif self.units is core.G3TimestreamUnits.Tcmb:
                                timestream_plot = timestream / (core.G3Units.K)
                                ylabel = 'temperature [K]'
                            else:
                                timestream_plot = timestream
                                ylabel = ''
                        else:
                            timestream_plot = timestream / self.units_factor
                            ylabel = ''
                        # finally plot
                        plt.plot(times, timestream, linewidth=0.5)
                    except:
                        pass
                    plt.title('{}: obs {}, scan {}'.format(bolo,
                                                           self.obsid,
                                                           self.jscan))
                    plt.xlabel('time [sec]')
                    domain_stub = 'timestream'
                else:
                    if bolo in frame[self.Input].keys():
                        # deal with units
                        if self.units_factor is None:
                            if self.units is core.G3TimestreamUnits.Power:
                                asd_plot = np.sqrt(psd[bolo]) / (core.G3Units.watt * 1e-18 / np.sqrt(core.G3Units.Hz))
                                ylabel = 'NEP [aW / rtHz]'
                            elif self.units is core.G3TimestreamUnits.Tcmb:
                                asd_plot = np.sqrt(psd[bolo]) / (core.G3Units.K * 1e-6 * np.sqrt(core.G3Units.second))
                                ylabel = 'NET [uK rtsec]'
                            else:
                                timestream_plot = timestream
                                ylabel = ''
                        else:
                            asd_plot = np.sqrt(psd[bolo]) / self.units_factor
                            ylabel = ''

                        # finally plot                        
                        plt.semilogy(f / core.G3Units.Hz, asd_plot, linewidth=0.5)
                        plt.title('{}: obs {}, scan {}, length {}'.format(bolo,
                                                                          self.obsid,
                                                                          self.jscan,
                                                                          len(f)))
                    if self.freqrange != None:
                        plt.xlim(self.freqrange)
                    plt.xlabel('frequency [Hz]')
                    domain_stub = 'psd'
                plt.ylabel(ylabel)

                plt.tight_layout()
                if self.name_stub == '':
                    plt.savefig('{}_{}_{}_{:03d}.png'
                                .format(self.obsid, domain_stub, bolo, self.jscan),
                                dpi=200)
                else:
                    plt.savefig('{}_{}_{}_{}_{:03d}.png'
                                .format(self.name_stub, self.obsid, domain_stub, bolo, self.jscan),
                                dpi=200)
                plt.close()
                jbolo += 1
            self.jscan += 1


class AccumulateTimestreams(object):
    def __init__(self, bololist=[], n_bolos=0, Input='RawTimestreams_I'):
        self.bololist = bololist
        self.n_bolos = n_bolos
        self.saved_ts = {'time': np.array([])}
        self.Input = Input

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration and \
                "NominalBolometerProperties" in frame.keys():
            if len(self.bololist) == 0:
                self.bololist = np.random.choice(list(frame["NominalBolometerProperties"].keys()), self.n_bolos)

        if frame.type == core.G3FrameType.Scan:
            for bolo in self.bololist:
                if bolo in frame[self.Input].keys() and \
                        bolo not in self.saved_ts.keys():
                    self.saved_ts[bolo] = np.array(frame[self.Input][bolo])
                elif bolo in frame[self.Input].keys():
                    self.saved_ts[bolo] = np.append(self.saved_ts[bolo], frame[self.Input][bolo])
            if len(frame[self.Input].keys()) > 0:
                bolo = list(frame[self.Input].keys())[0]
                frame[self.Input][bolo].times()
                times_datetimes = [datetime.utcfromtimestamp(g3tt.time / core.G3Units.second)
                                   for g3tt in frame[self.Input][bolo].times()]
                self.saved_ts['time'] = np.append(self.saved_ts['time'],
                                                  times_datetimes)
        

def Keep(frame, keys=[]):
    ''' 
    Deletes all but the enumerated keys from the frame.
    '''
    for key in frame:
        if key not in keys:
            del frame[key]

class SkipNFrames(object):
    def __init__(self, n_frames, type_to_skip=None):
        self.j_frame = 0
        self.n_frames = n_frames
        self.type_to_skip = type_to_skip

    def __call__(self, frame):
        if self.type_to_skip != None and frame.type == self.type_to_skip:
            self.j_frame += 1
        if self.type_to_skip == None:
            self.j_frame += 1
            
        if self.j_frame < self.n_frames:
            return []

class StopAfterNFrames(object):
    def __init__(self, n_frames, type_to_skip=None):
        self.j_frame = 0
        self.n_frames = n_frames
        self.type_to_skip = type_to_skip

    def __call__(self, frame):
        if self.type_to_skip != None and frame.type == self.type_to_skip:
            self.j_frame += 1
        if self.type_to_skip == None:
            self.j_frame += 1

        if self.j_frame > self.n_frames:
            core.G3Pipeline.halt_processing()


class UseScanNums(object):
    def __init__(self, frame_list, frame_type=None):
        self.jframe = 0
        self.frame_list = frame_list
        self.frame_type = frame_type

    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            return
        if self.frame_type == None or frame.type == self.frame_type:
            if self.jframe not in self.frame_list:
                self.jframe += 1
                return []
            else:
                self.jframe += 1

class ScanBounds(object):
    def __init__(self):
        self.time_bounds = np.array([])
        self.turnarounds = np.array([])

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            if 'DetectorSampleTimes' in frame.keys():
                self.time_bounds = np.append(self.time_bounds,
                                             datetime.utcfromtimestamp(frame['DetectorSampleTimes'][0].time \
                                                                       / core.G3Units.second))
            if 'Turnaround' in frame.keys() and frame['Turnaround'] is True:
                self.turnarounds = np.append(self.turnarounds,
                                             True)
            else:
                self.turnarounds = np.append(self.turnarounds,
                                             False)
