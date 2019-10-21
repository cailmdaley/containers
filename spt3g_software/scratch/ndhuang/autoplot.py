import os
import argparse
import glob
from datetime import datetime, timedelta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
from spt3g import core, calibration

# some configs
# For each source, we need to know which key to check against, and the
# appropriate threshold
caltypes = {'elnod': {'framekey': 'ElnodSNSlopes',
                      'thresh': 5},
            'calibrator': {'framekey': 'CalibratorResponseSN',
                           'thresh': 5}}
class AutoBoloCounter(object):
    '''
    Provides a dict-like interface to access (or automatically read and then access)
    calibration frames, and then calculate their number of live bolos

    TODO: how do we handle the fact that most, but not all, elnods have negative response?
    '''
    def __init__(self, directory, source, 
                 boloprops_dir = '/poleanalysis/sptdaq/calresult/calibration/boloproperties/'):
        '''
        directory: str
            The directory of autoprocessed data files for `source`.
            Usually /spt/user/production/calibration/<source>
        source: str
            The source type
        '''
        self.dir = directory
        self.framekey = caltypes[source]['framekey']
        self.thresh = caltypes[source]['thresh']
        self.livecount = {}
        self.source = source
        self.boloprops_dir = boloprops_dir
        self.bprop_ids = np.array(sorted([int(bp.split('.')[0]) for bp in 
                                          os.listdir(self.boloprops_dir)]))

    def __getitem__(self, key):
        try:
            return self.livecount[key]
        except KeyError:
            self.load_obsid(key)
            return self.livecount[key]

    def load_obsid(self, obsid):
        '''
        Load the data for obsid, and count the number of live bolos.
        Store in self.livecount[obsid].
        '''
        bprop = max(self.bprop_ids[self.bprop_ids <= obsid])
        print(bprop, obsid, self.bprop_ids)
        bpropfile = core.G3File(os.path.join(self.boloprops_dir, 
                                             str(bprop) + '.g3'))
        boloproperties = bpropfile.next()['BolometerProperties']
        calfile = os.path.join(self.dir, str(obsid) + '.g3')
        calframe = core.G3File(calfile).next()
        if self.source == 'calibrator':
            # Some extra checks on the calibrator
            # We might consider checking el, but that requires a different info source
            if 'CalibratorResponseFrequency' not in calframe:
                # Don't know the calibrator frequency?  drop it
                self.livecount[obsid] = None
                return
            if calframe['CalibratorResponseFrequency'] / core.G3Units.Hz > 6.5:
                # Drop high frequency calibrator stares
                self.livecount[obsid] = None
                return
        self.livecount[obsid] = self.live_bolos_by_wafer(calframe[self.framekey],
                                                         boloproperties)

    @staticmethod
    def live_bolos_by_wafer(metric, boloproperties, thresh = 20):
        '''
        Count the number of live bolometers on a wafer.
    
        metric: dict-like
            A dict-like object that maps bolo name to a liveness
            metric (calibrator SN, for example)
        thresh: float
            The threshold above which a bolometer is considered alive
        '''
        live = {}
        for bolo, v in metric.iteritems():
            try:
                wafer = boloproperties[bolo].wafer_id
            except KeyError as err:
                print('Boloproperties and observation don\'t have the same bolometers!')
                return None
            if wafer not in live:
                live[wafer] = 0
            if v > thresh:
                live[wafer] += 1
        return live



def make_plot(time, livecount, total_live = False, linear = True):
    '''
    Plot the number of live bolometers over time, by wafer.
    
    time: array-like
        The times at which each measurement of liveness was made
    livecount: dict-like
        A dict-like mapping wafer name to an array of the same length
        as `time` containing the number of live bolometers on said
        wafer.
    '''
    # Plot configuration
    # This tries to make this plot comprehensible in black and white,
    # as well as color
    wafer_specs = [{'color': 'blue',
                    'linestyle': 'solid',
                    'marker': 'v'},
                   {'color': 'green',
                    'linestyle': 'solid',
                    'marker': '^'},
                   {'color': 'cyan',
                    'linestyle': 'solid',
                    'marker': 'o'},
                   {'color': 'black',
                    'linestyle': 'solid',
                    'marker': '*'},
                   {'color': 'magenta',
                    'linestyle': 'solid',
                    'marker': 'd'},
                   {'color': 'blue',
                    'linestyle': 'dashed',
                    'marker': 'o'},
                   {'color': 'green',
                    'linestyle': 'dashed',
                    'marker': 'v'},
                   {'color': 'cyan',
                    'linestyle': 'dashed',
                    'marker': 'd'},
                   {'color': 'magenta',
                    'linestyle': 'dashed',
                    'marker': '*'},
                   {'color': 'black',
                            'linestyle': 'dashed',
                            'marker': '^'}]
    fig = pl.figure(figsize = (14, 12))
    if total_live:
        all_wafers = np.array([livecount[w] for w in livecount.keys()])
        total = np.sum(all_wafers, axis = 0)
        wafer_specs[0]['linestyle'] = ''
        pl.plot(time, total, markeredgewidth = 0, label = 'All Wafers', **wafer_specs[0])
    else:
        for i, wafer in enumerate(sorted(livecount.keys())):
            wafer_specs[i]['linestyle'] = ''
            pl.plot(time, livecount[wafer], label = wafer, markeredgewidth = 0, **wafer_specs[i])
    pl.yscale('symlog')
    pl.ylabel('# of live bolos')
    pl.legend(loc = 'best')
    # Do some prettyfying
    ax = fig.axes[0]
    # Reasonable date format
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    fig.autofmt_xdate()
    # Make the logscale look nice
    if not linear:
        pl.yscale('symlog')
        # Setup the ticks we want on the y-axis
        # Regular logspacing for major ticks
        # a few minor ticks, labelled by their multiplier, not exponent
        maxtick = np.ceil(np.log10(np.max(ax.yaxis.get_ticklocs())))
        if maxtick > 0:
            maxtick = int(maxtick)
            majorticks = np.zeros(maxtick + 2)
            majorticks[1:] = np.logspace(0, maxtick, num = maxtick + 1, base = 10.)
            minortick_locs = [2, 4, 6]
            minortick_locs = np.round(minortick_locs)
            minorticks = np.concatenate([i * majorticks[1:] for i in minortick_locs])
            ax.yaxis.set_ticks(majorticks)
            ax.yaxis.set_major_formatter((matplotlib.ticker.LogFormatterMathtext(
                        base = 10., labelOnlyBase = True)))
            ax.yaxis.set_ticks(minorticks, minor = True)
            ax.yaxis.set_ticklabels(['{:d}'.format(mt) for mt in 
                                     np.repeat(minortick_locs, len(majorticks) - 1)],
                                    minor = True)
        # reset the y-axis in case we screwed it up while ticking
        yl = pl.ylim()
        pl.ylim(-.5, yl[1])
    # gridlines
    pl.grid(axis = 'y', which = 'both')
    pl.grid(axis = 'x')
    return fig
    

def live_count_for_time(time, obsids, counter, start, stop, do_all = False,
                        linear = False):
    calinds = np.where((time <= stop) & (time > start))[0]
    if len(calinds) == 0:
        if do_all:
            return None, None
    livecount = {}
    for j, cali in enumerate(calinds):
        try:
            live = counter[obsids[cali]]
        except StopIteration:
            # Sometimes there are empty g3 files
            continue
        if live is not None:
            for wafer in live:
                if wafer not in livecount.keys():
                    livecount[wafer] = np.zeros(len(calinds))
                livecount[wafer][j] = live[wafer]
    if livecount == {}:
        if do_all:
            return None, None
    if do_all:
        fig_all = make_plot(time[calinds], livecount, True, linear = linear)
        fig_wafer = make_plot(time[calinds], livecount, linear = linear)
        return fig_all, fig_wafer
    else:
        return make_plot(time[calinds], livecount)

'''
The following 3 functions make a full set of plots for the given time range.
They take the same arguments:
source: str
    The source to calculate liveness from.  The source must
      a) exist in caltypes (defined at the top of this file)
      b) be a valid GCP source that we have observed, and autoprocessed
obsids: (int) array-like
    A list of all the obsids that have been autoprocessed for `source`.
counter: AutoBoloCounter
    An instance of the AutoBoloCounter, initialized for `source`
start: datetime (Jan 1, 2017)
    The start date.  If this is changed, we will probably generate an
    entirely new set of plots, rather than just making the newest
    round.
dir: str
    The directory in which to save the plots.  A subdirectory for
    `source` will be created and the timeframe will be created.
'''

def make_monthlies(source, obsids, counter, start = datetime(2017, 1, 1, 0, 0, 0),
                   dir = '/home/ndhuang/plots/autoplots', linear = False):
    if linear:
        tag = '_linear'
    else:
        tag = ''
    dir = os.path.join(dir, source, 'monthly')
    if not os.path.exists(dir):
        os.makedirs(os.path.join(dir))
    t0 = datetime(2017, 1, 1, 0, 0)
    time = np.array([t0 + timedelta(seconds = int(oid)) for oid in obsids])
    if len(time) < 1:
        return
    if start.month < 12:
        next = datetime(start.year, start.month + 1, start.day, 
                        start.hour, start.minute)
    else:
        next = datetime(start.year + 1, 1, start.day, 
                        start.hour, start.minute)
    while next < time[-1]:
        filename = os.path.join(dir, '{year}-{month}{tag}.png'.format(year = start.year,
                                                                      month = start.month,
                                                                      tag = tag))
        if os.path.exists(filename):
            start = next
            if start.month < 12:
                next = datetime(start.year, start.month + 1, start.day, 
                                start.hour, start.minute)
            else:
                next = datetime(start.year + 1, 1, start.day, 
                                start.hour, start.minute)
            continue
        fig_all, fig_wafer = live_count_for_time(time, obsids, counter, 
                                                 start, next, do_all = True,
                                                 linear = linear)
        if fig_all is None:
            start = next
            if start.month < 12:
                next = datetime(start.year, start.month + 1, start.day, 
                                start.hour, start.minute)
            else:
                next = datetime(start.year + 1, 1, start.day, 
                                start.hour, start.minute)
            continue
        
        fig_all.gca().set_title('# Bolos with {source} SN > {thresh} for {year}-{month}'.format(
                source = source,
                thresh = caltypes[source]['thresh'],
                year = start.year,
                month = start.month))
        fig_all.savefig(filename.replace(str(start.month), str(start.month) + '_all'))
        pl.close(fig_all)
        fig_wafer.gca().set_title('# Bolos with {source} SN > {thresh} for {year}-{month}'.format(
                source = source,
                thresh = caltypes[source]['thresh'],
                year = start.year,
                month = start.month))
        fig_wafer.savefig(filename)
        pl.close(fig_wafer)
        start = next
        if start.month < 12:
            next = datetime(start.year, start.month + 1, start.day, 
                            start.hour, start.minute)
        else:
            next = datetime(start.year + 1, 1, start.day, 
                            start.hour, start.minute)

    now = datetime.utcnow()
    then = now - timedelta(days = 30)
    fig_all, fig_wafer = live_count_for_time(time, obsids, counter, then, now, 
                                             do_all = True, linear = linear)
    if fig_all is None:
        return
    fig_all.gca().set_title('# Bolos with {source} SN > {thresh} for the Last 30 Days'.format(
            source = source,
            thresh = caltypes[source]['thresh']))
    fig_all.savefig(os.path.join(dir, 'last30_all{}.png'.format(tag)))
    pl.close(fig_all)
    fig_wafer.gca().set_title('# Bolos with {source} SN > {thresh} for the Last 30 Days'.format(
            source = source,
            thresh = caltypes[source]['thresh']))
    fig_wafer.savefig(os.path.join(dir, 'last30{}.png'.format(tag)))
    pl.close(fig_all)

def make_weeklies(source, obsids, counter, start = datetime(2017, 1, 1, 0, 0, 0),
                  dir = '/home/ndhuang/plots/autoplots/', linear = False):
    if linear:
        tag = '_linear'
    else:
        tag = ''
    dir = os.path.join(dir, source, 'weekly')
    if not os.path.exists(dir):
        os.makedirs(dir)
    t0 = datetime(2017, 1, 1, 0, 0)
    time = np.array([t0 + timedelta(seconds = int(oid)) for oid in obsids])
    if len(time) < 1:
        return
    week = timedelta(days = 7)
    next = start + week
    while next < time[-1]:
        filename = os.path.join(dir, '{year}-{month}-{day}{tag}.png'.format(
                                   year = start.year, month = start.month,
                                   day = start.day, tag = tag))
        if os.path.exists(filename):
            start = next
            next = start + week
            continue
        fig_all, fig_wafer = live_count_for_time(time, obsids, counter, 
                                                 start, next, do_all = True,
                                                 linear = linear)
        if fig_all is None:
            start = next
            next = start + week
            continue
        fig_all.gca().set_title('# Bolos with {source} SN > {thresh} for week of {year}-{month}-{day})'.format(
                source = source,
                thresh = caltypes[source]['thresh'],
                year = start.year,
                month = start.month,
                day = start.day))
        fig_all.savefig(filename.replace(str(start.day), str(start.day) + '_all'))
        pl.close(fig_all)
        fig_wafer.gca().set_title('# Bolos with {source} SN > {thresh} for week of {year}-{month}-{day})'.format(
                source = source,
                thresh = caltypes[source]['thresh'],
                year = start.year,
                month = start.month,
                day = start.day))
        fig_wafer.savefig(filename)
        pl.close(fig_wafer)
        start = next
        next = start + week

    now = datetime.utcnow()
    then = now - week
    fig_all, fig_wafer = live_count_for_time(time, obsids, counter, then, now, 
                                             do_all = True, linear = linear)
    if fig_all is None:
        return
    fig_all.gca().set_title('# Bolos with {source} SN > {thresh} for the Last Week'.format(
            source = source,
            thresh = caltypes[source]['thresh']))
    fig_all.savefig(os.path.join(dir, 'last7_all{}.png'.format(tag)))
    pl.close(fig_all)
    fig_wafer.gca().set_title('# Bolos with {source} SN > {thresh} for the Last Week'.format(
            source = source,
            thresh = caltypes[source]['thresh']))
    fig_wafer.savefig(os.path.join(dir, 'last7{}.png'.format(tag)))
    pl.close(fig_wafer)

def make_last(source, obsids, counter, dir = '/home/ndhuang/plots/autoplots/',
              linear = False):
    '''
    Make liveness plots for 24 hours, ending at the last observation we have
    '''
    if not os.path.exists(os.path.join(dir, source)):
        os.makedirs(os.path.join(dir, source))
    t0 = datetime(2017, 1, 1, 0, 0)
    time = np.array([t0 + timedelta(seconds = int(oid)) for oid in obsids])
    if len(time) < 1:
        return
    now = time[-1]
    then = now - timedelta(days = 1)
    fig_all, fig_wafer = live_count_for_time(time, obsids, counter, then, now, 
                                             do_all = True, linear = linear)
    if fig_all is None:
        return
    fig_all.gca().set_title('# Bolos with {source} SN > {thresh} for the Last Day'.format(
            source = source,
            thresh = caltypes[source]['thresh']))
    if linear:
        fig_all.savefig(os.path.join(dir, source + '_last_all.png'))
    else:
        fig_all.savefig(os.path.join(dir, source + '_last_all_linear.png'))
    pl.close(fig_all)
    fig_wafer.gca().set_title('# Bolos with {source} SN > {thresh} for the Last Day'.format(
            source = source,
            thresh = caltypes[source]['thresh']))
    if linear:
        fig_wafer.savefig(os.path.join(dir, source + '_last.png'))
    else:
        fig_wafer.savefig(os.path.join(dir, source + '_last_linear.png'))
    pl.close(fig_wafer)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot-dir', default = None)
    args = parser.parse_args()
    if args.plot_dir is None:
        # try to do something sensible depending on where we run
        plot_dir = os.path.join(os.getenv('HOME'), 'plots', 'data_quality', 'liveness')
    else:
        plot_dir = args.plot_dir
    for source in caltypes.keys():
        calfiles = np.array(glob.glob(os.path.join('/spt/user/production/calibration/', source, '*.g3')))
        # calfiles = np.array(glob.glob(os.path.join(
        #       '/poleanalysis/sptdaq/calresult/calibration/', source, '*.g3')))
        obsids = np.array([int(os.path.basename(cf).split('.')[0]) for cf in calfiles])
        inds = np.argsort(obsids)
        calfiles = calfiles[inds]
        obsids = obsids[inds]
        counter = AutoBoloCounter(os.path.join('/spt/user/production/calibration/', source), 
                                  source, 
                                  boloprops_dir = '/spt/user/production/calibration/boloproperties/')
        # counter = AutoBoloCounter(os.path.join('/poleanalysis/sptdaq/calresult/calibration/', source), source)
        for linear in [True, False]:
            make_weeklies(source, obsids, counter, dir = plot_dir, linear = linear)
            make_monthlies(source, obsids, counter, dir = plot_dir, linear = linear)
            make_last(source, obsids, counter, dir = plot_dir, linear = linear)
