#!/usr/bin/env python

import os
import glob
import re
import numpy as np
import subprocess as sp
from spt3g import core, dfmux

def check_obs_files(obs_root, cal_file_name='nominal_online_cal.g3',
                    data_root='/spt_data/bolodata'):
    """
    Verify observation data structure on disk:

    (1) cal file must exist and be verified
    (2) frame file sequence must be continuous starting from 0

    Arguments
    ---------
    obs_root : string
        Observation pathstring
    cal_file_name : string
        Name of nominal calibration file
    data_root : string
        Observation data file root

    Returns
    -------
    files : list of strings
        Ordered sequence of frame files to verify
    """

    files = sorted(glob.glob(os.path.join(data_root, obs_root, '*.g3')))
    file_bases = [os.path.basename(x) for x in files]

    # ensure cal file exists
    if cal_file_name not in file_bases:
        msg = 'Observation {} is corrupt! Missing {}.'
        raise RuntimeError(msg.format(obs_root, cal_file_name))

    # ensure that cal file is not corrupt
    cal_file = os.path.join(data_root, obs_root, cal_file_name)
    sp.check_call(['spt3g-verify', cal_file])

    # keep only the observation frame files
    regs = [re.search('(\d+.g3)', f) for f in file_bases]
    files = [f for f, r in zip(files, regs) if r and f.endswith(r.group(0))]
    file_bases = [os.path.basename(x) for x in files]

    # ensure at least one frame file
    if not len(files):
        msg = 'Observation {} is corrupt! No frame files found.'
        raise RuntimeError(msg.format(obs_root))

    # ensure file IDs are contiguous and start at zero
    file_nums = np.array([int(f.split('.')[0]) for f in file_bases])
    if file_nums[0] != 0 or (len(file_nums) > 1 and (np.diff(file_nums) != 1).any()):
        msg = 'Observation {} is corrupt! Found non-contiguous files {}.'
        raise RuntimeError(msg.format(obs_root, file_bases))

    return files

class OnlineVerifier:
    """
    Class for verifying data frames.
    """

    def __init__(self, obs, check_continuity=True, check_hk=True):
        """
        Arguments
        ---------
        obs : string
            Observation pathstring
        """
        self.obs = obs
        self.check_continuity = check_continuity
        self.check_hk = check_hk
        self.last = None
        self.count = 0

    def __call__(self, frame):
        """
        Verify that:

        (1) Scan frames are contiguous within the data sample rate
        (2) Housekeeping is present in every scan frame
        """
        if frame.type != core.G3FrameType.Scan:
            return

        self.count += 1

        # check continuity
        if self.check_continuity:
            tsmap = frame['RawTimestreams_I']
            if len(tsmap.keys()) == 0:
                raise RuntimeError(
                    'Observation {} is missing all bolometers at frame {} ({})'.format(
                        self.obs, self.count - 1, str(tsmap.start)))
            if self.last is not None:
                if tsmap.sample_rate <= 0:
                    raise RuntimeError(
                        'Observation {}: Scan frame {} at {} has a nonsense sample rate {:.4f} Hz'.format(
                            self.obs, self.count - 1, str(tsmap.start), tsmap.sample_rate / core.G3Units.Hz))
                rdelta = 1. / tsmap.sample_rate / core.G3Units.s
                delta = (tsmap.start.time - self.last.time) / core.G3Units.s
                if delta > 1.5 * rdelta:
                    raise RuntimeError(
                        ('Observation {} is not contiguous at frame {} ({}).\n'
                         'Scan frame separation is {} s, expected {} s.').format(
                             self.obs, self.count - 1, str(tsmap.start), delta, rdelta))
            self.last = tsmap.stop

        # check housekeeping
        if self.check_hk:
            hkmap = frame['DfMuxHousekeeping']
            if hkmap is None or not len(hkmap.keys()):
                raise RuntimeError(
                    'Observation {} is missing housekeeping data at frame {} ({}).'.format(
                        self.obs, self.count - 1, str(tsmap.start)))

        return


if __name__ == "__main__":

    import argparse as ap

    P = ap.ArgumentParser(
        description='Check processed observations',
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    P.add_argument('rate', choices=['fullrate', 'downsampled'],
                   help='Observation data rate')
    P.add_argument('source', help='Observation source')
    P.add_argument('observation', type=int, help='Observation ID')
    P.add_argument('--cal-file-name', default='nominal_online_cal.g3',
                   help='Calibration file to check')
    P.add_argument('--data-root', default='/spt_data/bolodata',
                   help='Observation data file root')
    P.add_argument('--no-check-hk', default=True, dest='check_hk',
                   action='store_false',
                   help='Check for housekeeping data in scan frames.')
    P.add_argument('--no-check-continuity', default=True, dest='check_continuity',
                   action='store_false',
                   help='Check for continuity between scan frames.')
    P.add_argument('-v', '--verbose', default=False, action='store_true',
                   help='More words')

    args = P.parse_args()

    # collect frame file list
    obs_root = os.path.join(args.rate, args.source, str(args.observation))

    # check files structure on disk
    files = check_obs_files(obs_root, cal_file_name=args.cal_file_name,
                            data_root=args.data_root)

    verifier = OnlineVerifier(obs_root, check_hk=args.check_hk,
                              check_continuity=args.check_continuity)

    # verify data frames
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=files)
    pipe.Add(verifier)
    if args.verbose:
        pipe.Add(core.Dump)
    pipe.Run()

    print('Observation {} OK'.format(obs_root))
