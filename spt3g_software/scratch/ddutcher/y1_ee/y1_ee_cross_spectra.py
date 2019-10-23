"""
This module contains functions for pre-processing the maps
(T-P deprojection, scaling the field halves to match Planck)
and for computing their cross spectra.
"""
import os
import numpy as np
import argparse
from spt3g import core, mapmaker, coordinateutils, std_processing
from spt3g.util import files
from spt3g.mapspectra import map_analysis


def make_planck_calibration_template(
    map_template=None, band=None, high_dec=1, low_dec=1):
    """
    Creates a template to calibrate 1500d field maps to Planck.
    
    Takes input map_template fills all points with
    declination <= -56deg by low_dec number, and all points with
    declination > -56deg with high_dec number.
    Multiply maps by this template to calibrate them.
    """
    if map_template is None:
        map_template = std_processing.CreateFieldMapStub()
        # temporary work-around before sparse_maps is merged
        map_template[0] = 0

    if band is not None:
        if int(band) == 90:
            high_dec = 0.9934
            low_dec = 0.9331
        elif int(band) == 150:
            high_dec = 0.9648
            low_dec = 0.8882
        elif int(band) == 220:
            high_dec = 0.9551
            low_dec = 0.9400
        else:
            raise ValueError(band)
        
    template = coordinateutils.maputils.get_ra_dec_map(map_template)[1]
    del map_template
    np.asarray(template)[
        np.asarray(template) > -56 * core.G3Units.deg
    ] = high_dec
    np.asarray(template)[
        np.asarray(template) <= -56 * core.G3Units.deg
    ] = low_dec
    
    return template


def apply_planck_calibration(mp, cal_template):
    for k in ['T', 'Q', 'U']:
        old = mp.pop(k)
        mp[k] = old * cal_template
        del old
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('maps', nargs=2, action='store',
                        help='Filepaths of input maps')
    parser.add_argument('--apod', action='store', default=None, type=str,
                        help='Filepath to apod mask')
    parser.add_argument('--band', action='store', type=int,
                        choices=[90, 150, 220])
    parser.add_argument('--deproject-tp', action='store_true', default=False,
                        help='Perform monopole T-to-P deprojection.')
    parser.add_argument('--planck-calibration', action='store_true', default=False,
                        help='Multiply field halves by Planck calibration numbers.')
    parser.add_argument('--output', '-o', action='store', default='output.pkl',
                        help='Location to store outputs')
    args = parser.parse_args()

    map1 = list(core.G3File(args.maps[0]))[-1]
    map2 = list(core.G3File(args.maps[1]))[-1]
    
    if args.apod is not None:
        if str(args.apod).endswith('.pkl'):
            apod = files.load_pickle(args.apod)
        elif str(args.apod).endswith('.g3') or str(args.apod).endswith('.gz'):
            apod = list(core.G3File(args.apod))[-1]
        else:
            raise OSError(args.apod)
    else:
        apod = None

    if args.deproject_tp:
        map1 = map_analysis.deproject_tp(map1, band=args.band)
        map2 = map_analysis.deproject_tp(map2, band=args.band)

    if args.planck_calibration:
        template = make_planck_calibration_template(map1['T'], band=args.band)
        apply_planck_calibration(map1, template)
        apply_planck_calibration(map2, template)

    cross_spectrum = map_analysis.calculate_powerspectra(
        map1, input2=map2, apod_mask=apod, b_mode_method='basic', calculate_dls=True)
    
    files.save_pickle(cross_spectrum, args.output)