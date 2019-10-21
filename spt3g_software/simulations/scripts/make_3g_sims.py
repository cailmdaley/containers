#!/usr/bin/env python
import os, sys
import argparse
import glob
import time
import yaml

# Create a parser for first handling the config file
cfg_parser = argparse.ArgumentParser(
    description='Generate sky simulations', add_help=False
)
cfg_parser.add_argument(
    '-c',
    '--config-file',
    action='store',
    help='YAML file containing configuration parameters for sims. '
    'Any of the arguments below can be included in the file, and '
    'will be overridden if supplied at the command line.',
)
args, _ = cfg_parser.parse_known_args()

# https://stackoverflow.com/a/20422915
class ActionNoYes(argparse.Action):
    def __init__(self, option_strings, dest, default=None, required=False, help=None):

        if default is None:
            raise ValueError('You must provide a default with Yes/No action')
        if len(option_strings) != 1:
            raise ValueError('Only single argument is allowed with YesNo action')
        opt = option_strings[0]
        if not opt.startswith('--'):
            raise ValueError('Yes/No arguments must be prefixed with --')

        opt = opt[2:]
        opts = ['--' + opt, '--no-' + opt]
        super(ActionNoYes, self).__init__(
            opts,
            dest,
            nargs=0,
            const=None,
            default=default,
            required=required,
            help=help,
        )

    def __call__(self, parser, namespace, values, option_strings=None):
        if option_strings.startswith('--no-'):
            setattr(namespace, self.dest, False)
        else:
            setattr(namespace, self.dest, True)


parser = argparse.ArgumentParser(
    parents=[cfg_parser], formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
group = parser.add_argument_group('Configuration parameters')
group.add_argument(
    '--output-root',
    '-o',
    action='store',
    required=True,
    help='Output directory for storing sims data files',
)
group.add_argument(
    '--num-threads', action='store', type=int, default=1, help='number of threads'
)
group.add_argument(
    '--sim-index',
    action='store',
    type=int,
    default=0,
    help='Index of first sim to generate. The random seed is uniquely '
    'determined for each foreground component from the sim number, so '
    'use this parameter to add more sims to an existing database.',
)
group.add_argument(
    '--num-sims',
    action='store',
    type=int,
    default=1,
    help='Number of simulations.  Sims will be numbered from '
    'sim_index to sim_index + num_sims',
)
group.add_argument(
    '--freqs',
    action='store',
    type=int,
    nargs='+',
    choices=[90, 150, 220],
    default=[150],
    help='Frequency channel(s) to include. Important for foregrounds '
    'and instrumental beam/noise.',
)
group.add_argument(
    '--nside',
    action='store',
    type=int,
    default=512,
    help='Healpix map resolution parameter',
)
group.add_argument(
    '--lmax',
    action='store',
    type=int,
    default=1000,
    help='Maximum multipole (ell) to simulate',
)
group.add_argument(
    '--pol', action=ActionNoYes, default=True, help='Include polarization (Q + U)'
)
group.add_argument('--cmb', action=ActionNoYes, default=False, help='Make cmb skies')
group.add_argument(
    '--foregrounds', action=ActionNoYes, default=False, help='Make foreground skies'
)
group.add_argument('--noise', action=ActionNoYes, default=False, help='Make noise sims')
group.add_argument(
    '--combine-gaussian-fg',
    action=ActionNoYes,
    default=False,
    help='Combine gaussian foregrounds into one file per frequency per sim.',
)
group.add_argument(
    '--beam',
    action=ActionNoYes,
    default=False,
    help='Smooth the final combined map with the beam',
)
group.add_argument(
    '--mask-file', action='store', help='Path to sky mask to apply to the final map(s)'
)
group.add_argument(
    '--debug',
    default=False,
    action=ActionNoYes,
    help='Show plots of the sim CMB and foregrounds.',
)
group.add_argument(
    '--verbose',
    default=False,
    action=ActionNoYes,
    help='Print healpy outputs, reduce log level',
)

group = parser.add_argument_group('CMB parameters')
group.add_argument(
    '--camb-file',
    default=os.path.join(
        'planck18_TTEEEE_lowl_lowE_lensing',
        'base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat',
    ),
    help='CAMB power spectra to use for CMB simulations',
)

group = parser.add_argument_group('Beam parameters')
group.add_argument(
    '--beam-file',
    action='store',
    help='Path to file containing B_ell for each frequency',
)
group.add_argument(
    '--fwhm-90',
    action='store',
    type=float,
    default=1.7,
    help='Default FWHM in arcmin for 90 GHz band, '
    'used for a gaussian beam if beam file is not supplied',
)
group.add_argument(
    '--fwhm-150',
    action='store',
    type=float,
    default=1.2,
    help='Default FWHM in arcmin for 150 GHz band, '
    'used for a gaussian beam if beam file is not supplied',
)
group.add_argument(
    '--fwhm-220',
    action='store',
    type=float,
    default=1.0,
    help='Default FWHM in arcmin for 220 GHz band, '
    'used for a gaussian beam if beam file is not supplied',
)

group = parser.add_argument_group('Foreground parameters')
group.add_argument(
    '--gaussian-fg-model',
    action='store',
    default='george',
    choices=['george'],
    help='Gaussian foreground model to use',
)
group.add_argument(
    '--gaussian-thermal-sz',
    action=ActionNoYes,
    default=True,
    help='Include a gaussian thermal SZ foreground component, '
    'using values from SPT-SZ (arXiv: 1408.3161)',
)
group.add_argument(
    '--gaussian-kinetic-sz',
    action=ActionNoYes,
    default=True,
    help='Include a gaussian kinetic SZ foreground component, '
    'using values from SPT-SZ (arXiv: 1408.3161)',
)
group.add_argument(
    '--gaussian-radio-galaxies',
    action=ActionNoYes,
    default=True,
    help='Include a gaussian radio galaxy foreground component, '
    'using values from SPT-SZ (arXiv: 1408.3161)',
)
group.add_argument(
    '--gaussian-dusty-galaxies',
    action=ActionNoYes,
    default=True,
    help='Include a gaussian dusty galaxy foreground component, '
    'using values from SPT-SZ (arXiv: 1408.3161)',
)
group.add_argument(
    '--gaussian-dg-clustering',
    action=ActionNoYes,
    default=True,
    help='Include the clustering term in the gaussian dusty galaxy '
    'foreground component, using values from SPT-SZ (arXiv: 1408.3161). '
    'Only used if --gaussian-dusty-galaxies is supplied.',
)
group.add_argument(
    '--poisson-radio-galaxies',
    action=ActionNoYes,
    default=True,
    help='Include poisson radio galaxy foreground component.',
)
group.add_argument(
    '--poisson-rg-model',
    action='store',
    default='dezotti',
    choices=['dezotti'],
    help='Model giving source counts and fluxes for radio galaxies.',
)
group.add_argument(
    '--poisson-dusty-galaxies',
    action=ActionNoYes,
    default=True,
    help='Include poisson dusty galaxy foreground component.',
)
group.add_argument(
    '--poisson-dg-model',
    action='store',
    default='bethermin',
    choices=['bethermin'],
    help='Model giving source counts and fluxes for dusty galaxies.',
)
group.add_argument(
    '--sz-pol-fraction',
    action='store',
    type=float,
    default=0.0,
    help='Polarization fraction of SZ components.',
)
group.add_argument(
    '--dg-pol-fraction',
    action='store',
    type=float,
    default=0.02,
    help='Polarization fraction of dusty galaxy component, '
    'Default from arXiv: astro-ph/0610485',
)
group.add_argument(
    '--rg-pol-fraction',
    action='store',
    type=float,
    default=0.03,
    help='Polarization fraction of radio galaxy component, '
    'Default from ACTPol: arXiv: 1811.01854; '
    'SPTpol: Gupta et al. in prep., etc.',
)
group.add_argument(
    '--min-flux-limit',
    action='store',
    type=float,
    default=6.4e-3,
    help='Poisson minimum source flux in Janskys.  Must agree with '
    'value used to calculate foreground power spectra.',
)
group.add_argument(
    '--max-flux-limit',
    action='store',
    type=float,
    default=5.0e-2,
    help='Poission maximum source flux in Janskys.',
)
group.add_argument(
    '--spec-index-radio-90-150',
    action='store',
    type=float,
    default=-0.7,
    help='Power-law index for scaling radio source flux at 150 GHz'
    'to 90 GHz, default from SPT-SZ W. Everett paper Fig. 3',
)
group.add_argument(
    '--spec-index-radio-220-150',
    action='store',
    type=float,
    default=-0.6,
    help='Power-law index for scaling radio source flux at 150 GHz '
    'to 220 GHz, default from SPT-SZ W. Everett paper Fig. 3',
)
group.add_argument(
    '--spec-index-dust-90-150',
    action='store',
    type=float,
    default=3.4,
    help='Power-law index for scaling dusty source flux at 150 GHz '
    'to 90 GHz, default from SPT-SZ W. Everett paper Fig. 3',
)
group.add_argument(
    '--spec-index-dust-220-150',
    action='store',
    type=float,
    default=3.4,
    help='Power-law index for scaling dusty source flux at 150 GHz '
    'to 220 GHz, default from SPT-SZ W. Everett paper Fig. 3',
)

group = parser.add_argument_group('Noise parameters')
group.add_argument(
    '--delta-t-90',
    action='store',
    type=float,
    default=3.0,
    help='Map depth at 90 GHz in uK-arcmin, default is ' 'SPT-3G 5-year forecast',
)
group.add_argument(
    '--delta-t-150',
    action='store',
    type=float,
    default=2.2,
    help='Map depth at 150 GHz in uK-arcmin, default is ' 'SPT-3G 5-year forecast',
)
group.add_argument(
    '--delta-t-220',
    action='store',
    type=float,
    default=8.8,
    help='Map depth at 220 GHz in uK-arcmin, default is ' 'SPT-3G 5-year forecast',
)
group.add_argument(
    '--lknee-t-90',
    action='store',
    type=int,
    default=1200,
    help='Multipole where 1/f noise flattens out in '
    'temperature at 90 GHz, default based on SPT-SZ',
)
group.add_argument(
    '--lknee-t-150',
    action='store',
    type=int,
    default=2200,
    help='Multipole where 1/f noise flattens out in '
    'temperature at 150 GHz, default based on SPT-SZ',
)
group.add_argument(
    '--lknee-t-220',
    action='store',
    type=int,
    default=2300,
    help='Multipole where 1/f noise flattens out in '
    'temperature at 220 GHz, default based on SPT-SZ',
)
group.add_argument(
    '--lknee-p-90',
    action='store',
    type=int,
    default=300,
    help='Multipole where 1/f noise flattens out in '
    'polarization at 90 GHz, default is a conservative guess',
)
group.add_argument(
    '--lknee-p-150',
    action='store',
    type=int,
    default=300,
    help='Multipole where 1/f noise flattens out in '
    'polarization at 150 GHz, default is a conservative guess',
)
group.add_argument(
    '--lknee-p-220',
    action='store',
    type=int,
    default=300,
    help='Multipole where 1/f noise flattens out in '
    'polarization at 220 GHz, default is a conservative guess',
)

# override default values from the config file
if args.config_file:
    cfg = yaml.safe_load(open(args.config_file, 'r'))
    parser.set_defaults(**cfg)

# parse all the arguments now
args = parser.parse_args()

# set number of threads before importing anything else
os.putenv('OMP_NUM_THREADS', str(args.num_threads))

# import all the big packages
import numpy as np
import healpy as hp
from spt3g import core
from spt3g.simulations import sim_tools, cmb, instrument, foregrounds

# for basic debugging purposes we will make simulations of 90 and 150.
if args.debug:
    args.freqs = [90, 150]

# print timestamps in log messages
core.G3Logger.global_logger.timestamps = True
if args.verbose:
    core.set_log_level(core.G3LogLevel.LOG_INFO)

# standardize file structure relative to output_root
def get_filename(subdir, name, sim_num, freq=None):
    if freq is not None:
        # ensure a unique filename
        name = name.replace('_alms', '_%dghz_alms' % freq)
        name = name.replace('_map', '_%dghz_map' % freq)
        if 'ghz' not in name:
            name = '%s_%dghz' % (name, freq)
    # ordered alphabetically by sim number using zero padding
    filename = os.path.join(args.output_root, subdir, '%s_%04d.fits' % (name, sim_num))
    # ensure output directory exists
    save_dir = os.path.dirname(filename)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return filename


# ===========================================================
# Set up simulation parameters
# ===========================================================
# list of sim numbers to generate
sim_range = list(range(args.sim_index, args.sim_index + args.num_sims))

# Set up dictionary to keep track of alms made for each sim.
# These will be added into a combined map if requested.
coadd_alms = {}
for freq in args.freqs:
    coadd_alms[freq] = {}
    for sim_num in sim_range:
        coadd_alms[freq][sim_num] = []

# Lower lmax to that allowed by nside
args.lmax = min([3 * args.nside - 1, args.lmax])

if args.debug:
    import pylab as plt

    def debug_plot(map_file, name, sim_num, freq=None, ref_cls=None):
        m = hp.read_map(map_file, verbose=False)
        title = '%s, Sim %d' % (name, sim_num)
        if freq:
            title = '%s, %s GHz' % (title, freq)
        hp.mollview(m, title=title)
        plt.savefig(map_file.replace('.fits', '.png'), bbox_inches='tight')

        cls = hp.anafast(m, lmax=args.lmax)

        if ref_cls is None:
            plt.figure()
            plt.plot(cls, 'r', label='Sim %d' % sim_num)
            plt.xscale('log')
            plt.yscale('log')
            plt.title(title)
        else:
            cls_diff = (ref_cls[: args.lmax + 1] - cls) / ref_cls[: args.lmax + 1]
            fig, axs = plt.subplots(2, 1, sharex=True)
            axs[0].plot(ref_cls, 'k', label=name)
            axs[0].plot(cls, 'r', label='Sim %d' % sim_num)
            axs[0].legend()
            axs[0].set_xscale('log')
            axs[0].set_yscale('log')
            axs[0].set_title(title)
            axs[1].plot(cls_diff, 'r', label='Sim %d' % sim_num)
            axs[1].axhline(color='k')
            axs[1].set_ylim(-1, 1)
        plt.savefig(map_file.replace('.fits', '_spec.png'), bbox_inches='tight')

        plt.show()


# ===========================================================
# Set up random seed
# -----------------------------------------------------------
# maximum number of unique components per sim number
max_components = 100

# dictionary of offsets for each component that requires a
# unique realization of random numbers
comp_offsets = {
    'cmb': 0,
    'tsz': 1,
    'ksz': 2,
    'rg': 3,
    'dg-cl': 4,
    'dg-po': 5,
    'poisson_rg': 6,
    'poisson_dg': 7,
    'inst_noise_90': 8,
    'inst_noise_150': 9,
    'inst_noise_220': 10,
    'atm': 11,
}


def set_seed(sim_num, comp):
    # random seed is determined uniquely for each sim index and
    # independent sim component
    np.random.seed(sim_num * max_components + comp_offsets[comp])


# ===========================================================
# Make simulated skies
# ===========================================================

# ===========================================================
# Make CMB realizations
# -----------------------------------------------------------
if args.cmb:
    core.log_notice('Creating CMB skies', unit='CMB')

    # Read CAMB spectrum
    # By default these units will be uK^2
    cls_dict = cmb.read_camb(args.camb_file, as_cls=True, lmax=args.lmax, lmin=0)
    if args.pol:
        camb_cls = np.asarray([cls_dict[x] for x in ['TT', 'EE', 'BB', 'TE']])
    else:
        camb_cls = cls_dict['TT']

    for sim_num in sim_range:
        core.log_info("Generating sim %d" % sim_num, unit='CMB')

        set_seed(sim_num, 'cmb')

        alms = hp.synalm(camb_cls, lmax=args.lmax, new=True, verbose=args.verbose)

        alm_file = get_filename('cmb', 'cmb_alms', sim_num)
        map_file = get_filename('cmb', 'cmb_map', sim_num)

        sim_tools.save_sims(
            alms=alms,
            store_alm=True,
            store_map=args.debug,
            alm_out=alm_file,
            map_out=map_file,
            lmax=args.lmax,
            nside=args.nside,
            pol=args.pol,
        )

        for freq in args.freqs:
            coadd_alms[freq][sim_num].append(alm_file)

        del alms

        if args.debug and sim_num == sim_range[0]:
            debug_plot(map_file, 'CMB', sim_num, ref_cls=cls_dict['TT'])

# ===========================================================
# Make foreground realizations
# -----------------------------------------------------------
if args.foregrounds:

    # -----------------------------------------------------------
    # Gaussian foregrounds
    # -----------------------------------------------------------
    for freq in args.freqs:
        core.log_notice(
            'Creating gaussian foregrounds for %sGHz channel' % freq, unit='Foregrounds'
        )

        fg_cls = foregrounds.get_foreground_sim_spectra(
            freq,
            model=args.gaussian_fg_model,
            thermal_sz=args.gaussian_thermal_sz,
            kinetic_sz=args.gaussian_kinetic_sz,
            radio_galaxies=args.gaussian_radio_galaxies,
            dusty_galaxies=args.gaussian_dusty_galaxies,
            dg_clustering=args.gaussian_dg_clustering,
        )

        for sim_num in sim_range:
            for key, cls in fg_cls.items():
                core.log_info(
                    "Generating %d GHz %s sim %d" % (freq, key, sim_num),
                    unit='Foregrounds',
                )

                # these foregrounds should be correlated between frequencies
                # so the seed is independent of freq
                set_seed(sim_num, key.lower())

                if args.pol:
                    # get polarization fractions
                    if 'sz' in key.lower():
                        pol_fraction = args.sz_pol_fraction
                    elif key.lower() == 'dg-cl':
                        # DG clustering term is not polarized
                        pol_fraction = 0
                    elif key.lower() == 'dg-po':
                        pol_fraction = args.dg_pol_fraction
                    elif 'rg' in key.lower():
                        pol_fraction = args.rg_pol_fraction
                    else:
                        raise KeyError(
                            'Cannot get pol_fraction for unrecognized %s' % key
                        )

                    # Assuming TE, TB, EB are all 0 for foregrounds
                    cls = [
                        cls,
                        cls * (pol_fraction ** 2) / 2,
                        cls * (pol_fraction ** 2) / 2,
                        cls * 0,
                    ]

                alms = hp.synalm(cls, lmax=args.lmax, new=True, verbose=args.verbose)

                if not args.combine_gaussian_fg:
                    alm_file = get_filename(
                        'foregrounds', '%s_alms' % key.lower(), sim_num, freq
                    )
                    map_file = get_filename(
                        'foregrounds', '%s_map' % key.lower(), sim_num, freq
                    )

                    sim_tools.save_sims(
                        alms=alms,
                        store_alm=True,
                        store_map=args.debug,
                        alm_out=alm_file,
                        map_out=map_file,
                        lmax=args.lmax,
                        nside=args.nside,
                        pol=args.pol,
                    )

                    coadd_alms[freq][sim_num].append(alm_file)

                    del alms

                    if args.debug and sim_num == sim_range[0]:
                        debug_plot(
                            map_file,
                            'Foreground: %s' % key,
                            sim_num,
                            freq=freq,
                            ref_cls=cls,
                        )

                else:
                    try:
                        gauss_fg_alms += alms
                    except:
                        gauss_fg_alms = np.copy(alms)
                    del alms

            if args.combine_gaussian_fg:
                alm_file = get_filename(
                    'foregrounds','combined_gaussian_alms', sim_num, freq
                )

                map_file = get_filename(
                    'foregrounds', 'combined_gaussian_map', sim_num, freq
                )

                sim_tools.save_sims(
                    alms=gauss_fg_alms,
                    store_alm=True,
                    store_map=args.debug,
                    alm_out=alm_file,
                    map_out=map_file,
                    lmax=args.lmax,
                    nside=args.nside,
                    pol=args.pol,
                )

                coadd_alms[freq][sim_num].append(alm_file)

                del gauss_fg_alms

                if args.debug and sim_num == sim_range[0]:
                    debug_plot(
                        map_file,
                        'Foreground: %s' % key,
                        sim_num,
                        freq=freq,
                        ref_cls=cls,
                    )
        # cleanup
        fg_cls.clear()

    # -----------------------------------------------------------
    # Poisson foregrounds
    # -----------------------------------------------------------
    # Generate one realization per source population at 150GHz,
    # then scale fluxes appropriately.
    for comp in ['radio_galaxies', 'dusty_galaxies']:

        if not getattr(args, 'poisson_%s' % comp):
            continue

        core.log_notice(
            'Creating Poisson sims of %s' % comp.replace('_', ' '), unit='Foregrounds'
        )

        comp_short = 'rg' if comp == 'radio_galaxies' else 'dg'
        comp_spec = 'radio' if comp == 'radio_galaxies' else 'dust'

        pol_fraction = getattr(args, '%s_pol_fraction' % comp_short)

        spec_index = {
            90: getattr(args, 'spec_index_%s_90_150' % comp_spec),
            150: 1.0,
            220: getattr(args, 'spec_index_%s_220_150' % comp_spec),
        }

        fluxes, counts = foregrounds.get_poisson_source_counts(
            comp,
            freq=150,
            rg_model=args.poisson_rg_model,
            dg_model=args.poisson_dg_model,
            min_flux_limit=args.min_flux_limit,
            max_flux_limit=args.max_flux_limit,
        )

        for sim_num in sim_range:
            core.log_info(
                "Generating Poisson %s sim %d" % (comp_short.upper(), sim_num),
                unit='Foregrounds',
            )

            # foregrounds correlated between frequencies
            set_seed(sim_num, 'poisson_%s' % comp_short)

            maps = foregrounds.make_poisson_source_sim(
                fluxes,
                counts,
                freq=150,
                pol=args.pol,
                pol_fraction=pol_fraction,
                nside=args.nside,
            )

            # Compute alms once here, and scale them for saving
            core.log_info("Converting %s poisson maps to alms..."%comp_spec)
            alms = np.asarray(hp.map2alm(maps, lmax=args.lmax, iter=0))
            core.log_info("...done converting %s poisson maps to alms."%comp_spec)

            for freq in args.freqs:
                scaling = (freq / 150.0) ** spec_index[freq]

                alm_file = get_filename(
                    'foregrounds', 'poisson_%s_alms' % comp_short, sim_num, freq
                )
                map_file = get_filename(
                    'foregrounds', 'poisson_%s_map' % comp_short, sim_num, freq
                )

                sim_tools.save_sims(
                    maps=maps * scaling,
                    alms=alms * scaling,
                    store_alm=True,
                    store_map=args.debug,
                    alm_out=alm_file,
                    map_out=map_file,
                    lmax=args.lmax,
                    nside=args.nside,
                    pol=args.pol,
                )

                coadd_alms[freq][sim_num].append(alm_file)

                if args.debug and sim_num == sim_range[0]:
                    debug_plot(
                        map_file,
                        'Foregrounds: Poisson %s' % comp_short.upper(),
                        sim_num,
                        freq=freq,
                    )

            # cleanup
            del maps, alms

# ===========================================================
# Make noise realizations
# -----------------------------------------------------------
if args.noise:
    # pending: atmosphere T/Q/U are currently uncorrelated. must be fixed.
    for freq in args.freqs:
        for comp in ['instrument', 'atmosphere']:
            core.log_notice(
                'Creating %s noise skies for %sGHz channel' % (comp, freq), unit='Noise'
            )

            comp_short = 'inst' if comp == 'instrument' else 'atm'

            nls = instrument.get_noise_sim_spectra(
                freq,
                component=comp,
                lmax=args.lmax,
                delta_t=getattr(args, 'delta_t_%d' % freq),
                lknee_t=getattr(args, 'lknee_t_%d' % freq),
                lknee_p=getattr(args, 'lknee_p_%d' % freq),
            )

            for sim_num in sim_range:
                core.log_info(
                    'Generating %d GHz %s sim %d' % (freq, comp, sim_num), unit='Noise'
                )

                if comp == 'instrument':
                    # independent noise realizations per frequency
                    set_seed(sim_num, 'inst_noise_%d' % freq)
                else:
                    # noise correlated between frequencies
                    set_seed(sim_num, 'atm')

                # temperature noise
                alms = hp.synalm(nls[0], lmax=args.lmax, new=True, verbose=args.verbose)
                if args.pol:
                    # pol noise
                    alms_q = hp.synalm(
                        nls[1], lmax=args.lmax, new=True, verbose=args.verbose
                    )
                    alms_u = hp.synalm(
                        nls[1], lmax=args.lmax, new=True, verbose=args.verbose
                    )
                    alms = np.asarray([alms, alms_q, alms_u])
                    del alms_q, alms_u

                alm_file = get_filename('noise', '%s_alms' % comp_short, sim_num, freq)
                map_file = get_filename('noise', '%s_map' % comp_short, sim_num, freq)

                sim_tools.save_sims(
                    alms=alms,
                    store_alm=True,
                    store_map=args.debug,
                    alm_out=alm_file,
                    map_out=map_file,
                    lmax=args.lmax,
                    nside=args.nside,
                    pol=args.pol,
                )

                coadd_alms[freq][sim_num].append(alm_file)

                # cleanup
                del alms

                if args.debug and sim_num == sim_range[0]:
                    debug_plot(
                        map_file, '%s Noise' % comp.capitalize(), sim_num, freq=freq
                    )

# ===========================================================
# Combine components and save the final maps
# ===========================================================
core.log_notice('Combining sim components into maps', unit='Combined')

for freq in args.freqs:
    for sim_num in sim_range:

        combined_map = sim_tools.combine_alms_into_map(
            coadd_alms[freq][sim_num],
            freq=freq,
            nside=args.nside,
            lmax=args.lmax,
            pol=args.pol,
            add_beam=args.beam,
            beamfile=args.beam_file,
            verbose=args.verbose,
        )

        filename = get_filename('total', 'total_map_3g', sim_num, freq)
        core.log_info('Saving spt3g map {}'.format(filename), unit='Combined')
        sim_tools.save_healpix_as_spt3g_map(
            combined_map, filename, maskfile=args.mask_file
        )

# storing the params in a timestamped file
args.alms = coadd_alms
tstamp = time.strftime('%Y%m%d_%H%M%S')
if args.config_file:
    prefix = os.path.splitext(os.path.basename(args.config_file))[0]
else:
    prefix = 'default_sims'
filename = os.path.join(args.output_root, '%s_%s.yaml' % (prefix, tstamp))
with open(filename, 'w') as f:
    yaml.dump(vars(args), f, default_flow_style=False)
