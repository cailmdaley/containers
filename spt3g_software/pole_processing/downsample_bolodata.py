#!/usr/bin/env python
import numpy as np
from spt3g import core, dfmux, pointing, std_processing
from spt3g.todfilter.downsampler import Downsampler
from spt3g.calibration.elnod_analysis import RotateIQ
from spt3g.calibration.calibrator_analysis import CalibratorRotateIQ

@core.scan_func_cache_data(wiringmap = 'WiringMap', hkmap = 'DfMuxHousekeeping')
def check_housekeeping(frame, keys = ['RawTimestreams_I'], 
                       wiringmap = None, hkmap = None):
    """
    Check for housekeeping data.
    Drop channels with zero carrier amplitude from the frame
    """
    # elif frame.type != core.G3FrameType.Scan:
    #     return

    # check that hkmap exists and is populated
    if hkmap is None or not len(hkmap.keys()):
        raise RuntimeError(
            'Observation {} is missing housekeeping data'.format(
                args.output_root))

    # drop unbiased channels
    for k in keys:
        data = frame.pop(k, None)
        if data is None:
            raise ValueError('Missing data {}'.format(k))
        new_data = core.G3TimestreamMap()
        for bolo in data.keys():
            status = dfmux.HousekeepingForBolo(hkmap, wiringmap, bolo)
            if status.carrier_amplitude != 0:
                new_data[bolo] = data[bolo]
        frame[k] = new_data

def downsample_nan_times(fr, nantimes_key = 'NanSampleTimes',
                         downsampledtsmap_key = 'RawTimestreams_I'):
    '''
    "Downsample" the times at which we have interpolated over NaNs.
    This is a fast approximation assuming a 33-sample sinc kernel, with a 
    downsampling factor of 2.  It will be *wrong* for any other kernel.

    The algorithm is simple: if the timestamp is exactly that of an output sample
    (i.e. the NaN is in phase with the decimation), keep the timestamp as is.
    If not, mark the two closest samples as interpolated NaNs.
    '''
    if nantimes_key not in fr:
        return
    bolotimes = fr[downsampledtsmap_key].times()
    nantimes = fr.pop(nantimes_key, None)
    for bolo, this_time in nantimes:
        new_time = []
        for t in this_time:
            if t not in bolotimes:
                i = np.searchsorted(bolotimes, t)
                new_time += [bolotimes[i - 1], bolotimes[i]]
            else:
                new_time.append(t)
        nantimes[bolo] = core.G3VectorTime(new_time)
    fr[nantimes_key] = nantimes
        
if __name__ == "__main__":

    import argparse as ap
    import os, glob, re, shutil
    import subprocess as sp

    P = ap.ArgumentParser(
        description='Downsample processed observations',
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    P.add_argument('input_root',
                   help='Input directory where the fullrate frame files are stored')
    P.add_argument('output_root',
                   help='Output directory where the downsampled frame files are stored')
    P.add_argument('-d', '--downsample-factor', default=2, type=int,
                   help='Factor by which to downsample input data')
    P.add_argument('--keep-unbiased', default=False, action='store_true',
                   help='Keep unbiased channels in the frame')
    P.add_argument('--keep-q', default=False, action='store_true',
                   help='Keep Q-phase data')
    P.add_argument('--scratch', default=False, action='store_true',
                   help="Copy input files to the condor scratch directory, and copy "
                   "completed jobs back to the output directory when complete.")
    P.add_argument('--calframe', help="Calframe to use for IQ rotation")
    P.add_argument('-v', '--verbose', default=False, action='store_true',
                   help='More words')

    args = P.parse_args()

    # These are the frame keys we want to downsample
    ds_keys = ['RawTimestreams_I'] + ['RawTimestreams_Q'] * args.keep_q

    pipe = core.G3Pipeline()

    # find all observation files
    files = sorted(glob.glob(os.path.join(args.input_root, '*.g3')))
    regs = [re.search('([0-9]+.g3)', os.path.basename(f)) for f in files]
    files = [f for f, r in zip(files, regs) if r and f.endswith(r.group(0))]

    scratch = None
    if args.scratch:
        scratch = os.getenv('_CONDOR_SCRATCH_DIR')

    if scratch:
        remote_input = args.input_root
        remote_output =args.output_root

        # create temporary directory structure on local disks
        args.input_root = os.path.join(scratch, 'raw')
        try:
            os.makedirs(args.input_root)
        except OSError:
            pass
        args.output_root = os.path.join(scratch, 'output')
        try:
            os.makedirs(args.output_root)
        except OSError:
            pass

        # copy input files over
        sp.check_call(
            'rsync -aviP {} {}'.format(' '.join(files), os.path.join(args.input_root, '.')),
            shell=True,
        )
        files = [os.path.join(args.input_root, os.path.basename(f)) for f in files]

    if args.calframe:
        files = [args.calframe] + files

    # read in data
    pipe.Add(core.G3Reader, filename=files)

    # interpolate over nans
    pipe.Add(std_processing.InterpOverNans)

    # rotate IQ phase
    if args.calframe:
        kwargs = dict(i_rotate_key='RawTimestreams_I',
                      q_rotate_key='RawTimestreams_Q' if args.keep_q else None,
                      destructive=True)
        if 'calibrator' in args.calframe:
            pipe.Add(CalibratorRotateIQ, **kwargs)
        else:
            pipe.Add(RotateIQ, **kwargs)
        # drop calibration frame once it's been used
        pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)

    # drop some data early to reduce memory usage
    pipe.Add(core.Delete,
             keys=['OnlineRaDecRotation'] + ['RawTimestreams_Q'] * (not args.keep_q))

    # filter, decimate, and clean up the outputs
    pipe.Add(Downsampler, ds_factor=args.downsample_factor,
             keys=ds_keys, compress=True, 
             key_suffix='_downsampled')

    # "Downsample" the timestamps of interpolated samples.  
    # This is hardcoded for downsampling by a factor of 2, so just refuse
    # to do this if `downsample_factor` is not 2.
    if args.downsample_factor == 2:
        pipe.Add(downsample_nan_times)
    else:
        core.log_warn("Downsample factor is not 2, so fr['NanSampleTimes'] will not be downsampled")

    # check for HK data and drop channels that have zero carrier amplitude
    if not args.keep_unbiased:
        pipe.Add(check_housekeeping, keys = ds_keys)

    # remove duplicate housekeeping frames
    pipe.Add(core.DeduplicateMetadata)

    # Add Online pointing to scans.
    # Clean up pre-existing timestreams
    pipe.Add(
        core.Delete,
        keys=['OnlineBoresightAz', 'OnlineBoresightEl',
              'OnlineBoresightRa', 'OnlineBoresightDec', 'OnlineRaDecRotation']
    )
    pipe.Add(
        pointing.CalculateCoordTransRotations,
        raw_az_key='RawBoresightAz',
        raw_el_key='RawBoresightEl',
        output='OnlineBoresight',
        transform_store_key='OnlineRaDecRotation',
        model='OnlinePointingModel',
        flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
        #, 'thermolin'] # Thermoline broken as of 4/13/17
    )

    # sanity check
    if args.verbose:
        pipe.Add(core.Dump)

    # store downsampled data in re-sequenced frame files
    output_file = os.path.join(args.output_root, '%04d.g3')
    pipe.Add(core.G3MultiFileWriter, filename=output_file, size_limit=1024**3)

    pipe.Run()

    if scratch:
        # copy output files back to network storage
        sp.check_call(
            'rsync -aviP {} {}'.format(os.path.join(args.output_root, '*'), remote_output),
            shell=True,
        )
