-----------
simulations
-----------

The simulations directory contains code for creating simulated skies. The key
components of simulations are:

*CMB*
    The external software package CAMB is used to generate Cls from a specific
    cosmology that are then stored within the ``data/camb`` directory.  These
    Cls can be handed to e.g. healpy.synalm or healpy.synfast to generate a
    random CMB realization.

*Foregrounds*
    Foregrounds are cosmological sources of radiative power that are not the
    CMB. They include gaussian-distributed fields such as the thermal and
    kinetic Sunyaev-Zel'dovitch effects (tSZ and kSZ) and the unresolved flux
    from synchrotron-emitting galaxies, as well as poisson-distributed bright
    point sources.  Data on the power spectrum and/or number counts of these
    sources from previous experiments, such as SPT-SZ, is located in the
    ``data/foregrounds`` directory.

*Noise*
   The noise spectrum of our measurements can be roughly broken down into two
   parts: a flat white-noise component at higher frequencies coming from the
   instrument itself, and a rising component at lower frequencies due to the
   atmosphere. Noise realizations can also be obtained from the data maps
   themselves by subtracting one map from another.
 
These three components are handled by the python modules ``cmb.py``,
``foregrounds.py``, and ``instrument.py``, respectively.  The latter also
contains a function for obtaining the instrument beam, which is used to smooth
the final simulated map.

The ``scripts`` directory contains scripts that utilize these lower-level
utilities to make simulated skies tailored to specific analyses, e.g. maps for
mock-observing in the power spectrum analyses.

make_3g_sims.py
===============

The script ``make_3g_sims.py`` makes realizations of the quantities mentioned
above (CMB, foregrounds, noise) and produces the combined simulated sky map.
The individual components are saved as alms in .fits files; these are then
added together, smoothed with a beam if one is supplied, and made into a map. 
This map is then masked down to the 3G coverage area, if such a mask is given,
and the temperature and polarization components are saved as CutSkyHealpix maps
in a .fits file, i.e. as expected by ``coordinateutils.fitsio.load_skymap_fits``.

*Running the script*
    There are many parameters that need to be passed to the sim script, all
    of which have defaults set in the argparse section at the top of
    ``make_3g_sims.py``.  View the argument documentation via the help:
    
      python make_3g_sims.py -h
    
    The arguments are separated into groups as follows:
    
    * Configuration parameters : the output directory, the frequencies for
      which to make sims, which components to simulate, whether to apply a beam,
      nside and lmax of healpix maps, path to sky-coverage mask, etc.
    * CMB parameters : which input camb spectra to use
    * Beam parameters : The FWHM of the beam for each frequency band, or path
      to a file containing Bl's
    * Foreground parameters : which components to simulate, which model of
      their power spectra to use, etc.
    * Noise parameters : Map depth, temperature and polarization white-noise
      level and 1/f knee for each frequency band
    
    Many of these parameters may also be set via a config file, as in the
    example ``default_sims.yaml``.  You may use a config file in conjunction
    with command-line options, in which case options set on the command-line
    will override the value set in the config file.  All the settings will
    be recorded in a timestamped yaml file included with the output.
    
    A simple example call (fullsky, no beam) would look like

        python make_3g_sims.py -c default_sims.yaml -o sims_test --cmb --foregrounds

*Outputs*
    The individual simulated components are stored as alms in .fits files, with
    the T, Q, and U components stored in headers 1, 2, and 3, respectively.
    To e.g. read in just the Q alms and plot them as an nside512 map, do

    .. code-block:: python

        qalms = healpy.read_alm('cmb_alms_0000.fits', hdu=2)
        qmap = healpy.alm2map(qalms, 512)
        healpy.mollview(qmap)

    The input to the mock-observing pipeline is in the ``total`` subfolder.
    These .fits contain the temperature and polarization maps as CutSkyHealpix maps.
    To e.g. load these into a dictionary and plot the temperature map, do:

    .. code-block:: python

        m = coordinateutils.fitsio.load_skymap_fits('total_150ghz_map_3g_0000.fits')
        tmap = np.asarray(m['T'])
        healpy.mollview(tmap)
