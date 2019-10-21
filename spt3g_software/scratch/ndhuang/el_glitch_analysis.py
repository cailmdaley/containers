import os
import glob
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg', False)
from matplotlib import pyplot as pl
from spt3g import core

'''
Hokay, so.
This is a mess of tools I use to analyze el glitches.
There's a lot of copy-pasta to do the cuts, but the only way I thought
of to do that better would have taken effort.  It wouldn't be too
crazy to repack all the data with the cuts (i.e. use repack_dir).
But you figure out your own optimal workflow.
'''

colors = {0: 'blue',
          1: 'red',
          'el': 'orange'}

##########################################################################
# Plotting tools

def plot_glitch(gl, xaxis_time = False, ax = None):
    '''
    Plot a single glitch.  This is mostly useful as a tool
    for overlaying all the glitches at once.  See below.
    '''
    if ax is None:
        ax = pl.gca()
    if xaxis_time:
        time = [datetime.utcfromtimestamp(t.time / core.G3Units.s) 
                for t in gl['time']]
        alpha = .7
    else:
        time = np.arange(len(gl['time']), dtype = int)
        alpha = .1
    for key in [0, 1, 'el']:
        c = colors[key]
        ax.plot(time, gl[key] - np.median(gl[key]),
                color = c, alpha = alpha, lw = 2)
    ax.set_ylim(-.01, .01)
    ax.set_ylabel('Degrees (arbitrary offset)')
    return ax

def lotsaplots(outdir = '/home/ndhuang/plots/el_glitches/lotsa/'):
    '''
    Plot every single glitch.  Note the hard coded path below.
    This was useful for designing cuts, but it probably not useful any more.
    '''
    files = sorted(glob.glob('/poleanalysis/ndhuang/el_glitches/*.pkl'))
    i_fig = 0
    last_t = 0
    for f in files:
        try:
            out = np.load(f)
        except OSError:
            continue
        for gl in out[1]:
            if gl[1].time - last_t <= .02 * core.G3Units.s:
                continue
            gl[0]['el'] /= core.G3Units.deg
            dl = {k: np.diff(gl[0][k]) for k in [0, 1, 'el']}
            single_sample = False
            for k in [0, 1, 'el']:
                if dl[k][600] > 1:
                    single_sample = single_sample or is_single_sample(dl[k])
            if not single_sample:
                continue
            pl.figure()
            plot_glitch(gl[0], True)
            pl.legend(['Encoder 0', 'Encoder 1', 'El'])
            pl.savefig(os.path.join(outdir, '{:04d}.png'.format(i_fig)))
            pl.close()
            i_fig += 1
            last_t = gl[1].time

def plot_classes(outdir = '/home/ndhuang/plots/el_glitches/classified/'):
    '''
    Classify all the glitches, then overplot all glitches of each class.
    This is useful for determining facts such as 
    "All encoder 1 glitches have kicks".
    Classification is defined by `classify`, below.
    '''
    files = sorted(glob.glob('/spt/user/ndhuang/el_glitches/2018/*.pkl'))
    last_t = 0
    keys = [0, 1, 'el']
    figs = {k: pl.figure() for k in keys}
    axes = {k: figs[k].add_axes((.15, .1, .75, .8)) for k in keys}
    for k, ax in axes.items():
        ax.set_title(str(k))
    for f in files:
        try:
            out = np.load(f)
        except OSError:
            continue
        for gl in out[1]:
            if gl[1].time - last_t <= .02 * core.G3Units.s:
                continue
            if len(gl[0]['el']) == 0:
                continue
            gl[0]['el'] /= core.G3Units.deg
            dl = {k: np.diff(gl[0][k]) for k in keys}
            single_sample = False
            for k in [0, 1, 'el']:
                if dl[k][600] > 1:
                    single_sample = single_sample or is_single_sample(dl[k])
            if not single_sample:
                continue
            cl = classify(dl)
            for k in keys:
                if cl[k]:
                    plot_glitch(gl[0], False, axes[k])
            last_t = gl[1].time
    for k in keys:
        axes[k].legend(keys, loc = 'best')
        figs[k].savefig(os.path.join(outdir, str(k) + '_old.png'))

def get_class_and_time(pkldir = '/poleanalysis/ndhuang/el_glitches/sept-oct2018/'):
    '''
    This is an excercise in repacking the glitches.
    Returns a dictionary mapping time (in seconds since the epoch) to 
    glitch classificatoin.  Largely useful for feeding into other functions.
    '''
    files = sorted(glob.glob(os.path.join(pkldir, '*.pkl')))
    last_t = 0
    keys = [0, 1, 'el']
    classtime = {}
    for f in files:
        try:
            out = np.load(f)
        except OSError:
            continue
        for gl in out[1]:
            if gl[1].time - last_t <= .02 * core.G3Units.s:
                continue
            if len(gl[0]['el']) == 0:
                continue
            gl[0]['el'] /= core.G3Units.deg
            dl = {k: np.diff(gl[0][k]) for k in keys}
            single_sample = False
            for k in [0, 1, 'el']:
                if dl[k][600] > 1:
                    single_sample = single_sample or is_single_sample(dl[k])
            if not single_sample:
                continue
            classtime[gl[1].time / core.G3Units.s] = classify(dl)
            last_t = gl[1].time
    return classtime

def to_int(bools):
    '''
    Don't worry about this.
    '''
    out = 0
    for i, b in enumerate(bools):
        out += 2**i * b
    return out

classes = {1: 'Encoder 1',
           2: 'Encoder 2',
           3: 'Encoders 1 and 2',
           4: 'El',
           5: 'Encoder 1 and El',
           6: 'Encoder 2 and El',
           7: 'All'}

def get_frequency(times_by_class, dt, t_min, t_max):
    '''
    Get the number of glitches in each class for every `dt` seconds
    between t = `t_min` and t = `t_max`.  I.e. use this function
    if you want to get the number of glitches per hour as a fucntion
    of time.
    '''
    out = {}
    bins = np.arange(t_min, t_max + dt, dt)
    for key, v in times_by_class.items():
        out[key] = np.histogram(v, bins = bins)[0]
    out['Total'] = np.sum(list(out.values()), axis = 0)
    t = bins[:-1] + (bins[1] - bins[0]) / 2
    out['t'] = np.array([datetime.utcfromtimestamp(_t) for _t in t])
    return out    

def plot_frequency(freq):
    '''
    Plot glitch frequency for all classes.  Becomes indecipherable
    if you have too many bins.
    '''
    for k in ['Total', 'El', 'Encoder 1', 'Encoder 2', 'Encoder 1 and El']:
        pl.step(freq['t'], freq[k], label = k, where = 'mid', alpha = .7)

def analyze_classtime(classtime):
    '''
    Analyzes the output of get_class_and_time.
    Returns the number of glitches per hour, per day, and counts of glitches by
    class.
    '''
    keys = [0, 1, 'el']
    counts = {}
    times_by_class = {}
    for t in sorted(classtime.keys()):
        bools = [classtime[t][k] for k in keys]
        derpy = classes[to_int(bools)]
        if derpy in counts:
            counts[derpy] += 1
        else:
            counts[derpy] = 1
        if derpy in times_by_class:
            times_by_class[derpy].append(t)
        else:
            times_by_class[derpy] = [t]
    t_min = np.min(list(classtime.keys()))
    t_max = np.max(list(classtime.keys()))
    # hourly = get_frequency(times_by_class, 3600, t_min, t_max)
    daily = get_frequency(times_by_class, 3600 * 24, t_min, t_max)
    weekly = get_frequency(times_by_class, 3600 * 24 * 7, t_min, t_max)
    return daily, weekly, counts

def classify(dl):
    '''
    For a glitch, figure out which fields glitched.
    Returns a dictionary with True in each field that glitched.
    '''
    return {k: any(dl[k][599:602] > 1) for k in [0, 1, 'el']}

def filter_and_do(do):
    '''
    This is a trap.  Don't use it.
    '''
    files = sorted(glob.glob('/poleanalysis/ndhuang/el_glitches/*.pkl'))
    i_glitch = 0
    last_t = 0
    for f in files:
        try:
            out = np.load(f)
        except OSError:
            continue
        for gl in out[1]:
            if gl[1].time - last_t <= .02 * core.G3Units.s:
                continue
            if len(gl[0]['el']) == 0:
                continue
            gl[0]['el'] /= core.G3Units.deg
            dl = {k: np.diff(gl[0][k]) for k in [0, 1, 'el']}
            single_sample = False
            for k in [0, 1, 'el']:
                if dl[k][600] > 1:
                    single_sample = single_sample or is_single_sample(dl[k])
            if not single_sample:
                continue
            do(gl, i_glitch)
            i_glitch += 1
            last_t = gl[1].time

def get_jump_size_and_expected(el):
    '''
    For a glitchy timestream, figure out how far off the glitch
    was from the expected value.
    Returns (glitch_size, expected_value).
    '''
    dl = np.diff(el)
    inds = np.where(abs(dl) > 1)[0]
    m = 0
    n = 0
    for i in [598, 599, 600, 601]:
        if i in inds:
            m += abs(dl[i])
            n += 1
    assert(n == 2)
    for i in [598, 599, 600, 601]:
        if i in inds:
            expected = np.mean([el[i], el[i + 2]])
            break
    return m / 2., expected
        
def get_rel_jump(pkldir = '/poleanalysis/ndhuang/el_glitches/sept-oct2018/',
                 tocheck = 0):
    '''
    Prove that encoder 1 glitches look like bit shifts.
    For each glitch, calculate the glitch size and the expected value.
    Store expected - 2 * glitch, which should be identically 0 for a perfect
    bit shift.
    '''
    # aka prove to steve that bit shifts are happening
    files = sorted(glob.glob(os.path.join(pkldir, '*.pkl')))
    last_t = 0
    keys = [0, 1, 'el']
    results = {}
    for f in files:
        try:
            out = np.load(f)
        except OSError:
            continue
        for gl in out[1]:
            if gl[1].time - last_t <= .02 * core.G3Units.s:
                continue
            if len(gl[0]['el']) == 0:
                continue
            gl[0]['el'] /= core.G3Units.deg
            dl = {k: np.diff(gl[0][k]) for k in keys}
            single_sample = False
            for k in [0, 1, 'el']:
                if dl[k][600] > 1:
                    single_sample = single_sample or is_single_sample(dl[k])
            if not single_sample:
                continue
            cls = classify(dl)
            if cls[tocheck]:
                if tocheck == 1:
                    # The second encoder doesn't actually have its raw value 
                    # stored in arcfiles.  Correct for that.
                    thingy = 360 - gl[0][tocheck]
                else:
                    thingy = gl[0][tocheck]
                jump, expected = get_jump_size_and_expected(thingy)
                results[gl[1].time / core.G3Units.s] = expected - 2 * jump
    return results

def repack_dir(pkldir):
    '''
    Unpack all the glitches in pickle files in a directory, and
    return them as a dictionary mapping the glitch time to the glitch data.
    '''
    files = sorted(glob.glob(os.path.join(pkldir, '*.pkl')))
    last_t = 0
    keys = [0, 1, 'el']
    results = {}
    for f in files:
        try:
            out = np.load(f)
        except OSError:
            continue
        for gl in out[1]:
            if gl[1].time - last_t <= .02 * core.G3Units.s:
                continue
            if len(gl[0]['el']) == 0:
                continue
            gl[0]['el'] /= core.G3Units.deg
            dl = {k: np.diff(gl[0][k]) for k in keys}
            single_sample = False
            for k in [0, 1, 'el']:
                if dl[k][600] > 1:
                    single_sample = single_sample or is_single_sample(dl[k])
            if not single_sample:
                continue
            results[gl[1].time / core.G3Units.s] = gl[0]
    return results
##########################################################################
# Event cuts

def is_single_sample(dl):
    '''
    Return True if the glitch is a single sample, False otherwise.
    '''
    inds = np.where(abs(dl) > 1)[0]
    # 600 or 599 must be in
    return (599 in inds and (600 in inds or 598 in inds)) or \
        (600 in inds and (601 in inds or 599 in inds))
    

if __name__ == '__main__':
    '''
    laziness
    '''
    plot_classes()
    # lotsaplots()
            
