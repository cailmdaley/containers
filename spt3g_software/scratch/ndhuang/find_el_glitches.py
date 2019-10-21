import pickle
import numpy as np
from spt3g import core
from spt3g.util.extractdata import extract_keys

keys = {0: ['antenna0', 'tracker', 'raw_encoder', 0],
        1: ['antenna0', 'tracker', 'raw_encoder', 1],
        'el': ['antenna0', 'tracker', 'actual', 1],
        'time': ['antenna0', 'tracker', 'utc']}


def get_glitch(encoders, el, time, ind):
    '''
    Return a dictionary containing glitchy data.
    '''
    s = slice(ind - 600, ind + 600)
    return {0: encoders[0][s],
            1: encoders[1][s],
            'el': el[s],
            'time': time[s]}

def grub_for_glitches(encoders, el, time):
    '''
    Find glitches in pointing timestreams, and return 
    small chunks of data around the glitches.
    
    INPUTS
    ------
    encoders: array-like
        The elevation encoder timestreams, in degrees.
        `encoders[0]` should probably contain data from encoder 1
        and `encoders[1]` should probably contain data from encoder 2
    el: array-like
        Elevation timestream in native G3Units.
    time: array-like
        Sample times

    RETURNS
    -------
    A list of tuples of the form [(glitch, time), ...]
    where `glitch` is a dictionary containing data from the 2 encoders,
    elevation and the sample time for 6 seconds before and after 
    a glitch.  `time` is the time of the glitchy sample.
    '''
    enc_diff = np.diff(encoders, axis = 1)
    el_diff = np.diff(el)
    bad_enc = [np.where(d > 1)[0] for d in enc_diff]
    bad_el = np.where(el_diff > 1)[0]
    allbad = np.array([bad_el, bad_enc[0], bad_enc[1]])
    glitch_samples = np.unique(sorted(np.concatenate(allbad)))
    glitch_times = time[glitch_samples]
    glitches = [(get_glitch(encoders, el, time, i), t) 
                for i, t in zip(glitch_samples, glitch_times)]
    return glitches

description = 'Grub through a list of arcfiles looking for encoder and \
elevation glitches.  The arcfiles should be in (chronological) order.  \
The result is saved in a file specified by the -o flag (which is required).'

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('infiles', type = str, nargs = '+',
                        help = 'Arc files to search')
    parser.add_argument('-o', type = str, required = True,
                        help = 'Output file name')
    args = parser.parse_args()
    glitches = {}
    for af in args.infiles:
        data = extract_keys(af, keys)
        glitches.append(grub_for_glitches([data[0], data[1]], data['el'],
                                     data['time']))
    with open(args.o, 'wb') as f:
        pickle.dump((list(args.infiles), glitches), f)
