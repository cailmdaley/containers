-----------
calibration
-----------

The calibration project contains data classes and analysis code for storing the following things:

*Physical properties of the bolometers*
	This includes the relative pointing offsets of the detectors on the focal plane, their polarization angles and efficiences, their bands, and the fabrication name (physical name) of the detector to which the channel corresponds. It does *not* include tuning-dependent parameters like time constants or responsivity. This information is stored in a ``BolometerPropertiesMap``, indexed by logical bolometer ID, in the Calibration frame. The standard name for this object is "BolometerProperties" and it is created from a composite of other calibration information.

*Flux calibration*
	This is stored in a key named SourceFluxCalibration, where "Source" is usually RCW38. The data consist of a mapping from logical detector ID to the measured integral flux of the source within a 4x4 arcminute box in units of the calibrator response. The ``ApplyTCalibration`` module can be used in conjunction with this information to convert timestreams into calibrated surface brightness units. This information, as well as the relative pointing offsets stored in the bolometer properties map, is calculated by the scripts in the ``calibration/scripts/fluxpointcal`` directory.

*Elnod response*
	This stores the response of the detectors to an elevation nod of the telescope as is used for I/Q phasing. See the documentation for ``calibration.elnod.AnalyzeElnod`` for more information.

*Calibrator response*
	This stores the peak-to-peak response of each detector, in Watts, to the chopped thermal calibrator source. This information is used by ``ApplyTCalibration`` to transfer source calibration to another observing field without contamination from loading changes. The signal-to-noise measure provided (usually in CalibratorResponseSN) is also useful as a detector cut.

After calculation, each of these pieces of information is assembled into a single Calibration frame by the script ``build_cal_frame.py``. The resulting files are stored in each observing directory and should be read before any data to ensure that proper calibration is available.

Overview of Calibration Chain
=============================

The standard calibration chain provides temperature, pointing, and polarization calibration from four types of observations. These are processed by the auto-processing into the calibration frames that live with the data such that each calibration frame includes data from the most recent previous calibration observations.

*El nods*
	This provides a responsivity estimate (elnod signal-to-noise) and, on full rate data, the I/Q rotation matrix that maximizes responsivity. When applied in later steps, this rotation provides a small reduction in noise. El-nod processing depends on nothing.

*Calibrator*
	After the el nod, we use the chopped calibrator source as a temperature calibration transfer. This provides a parametrization of the response of every detector on the focal plane in terms of an (approximately) constant reference flux. Calibrator processing depends on the el nod processing, as it needs to operate on the same I/Q combination as later data processing.

*RCW38 Pixel Raster data (aka fast points)*
	RCW38 data from observations in which every detector sees the source (pixel raster) provide the fundamental temperature calibration of each detector in terms of sky fluxes. Results are expressed in terms of the calibrator result from the immediately preceding calibrator observation, allowing the transfer of the calibration from RCW38 to another part of the sky with potentially different loading. These observations thus depend on the calibrator data. Only computes results for detectors with calibrator S/N greater than 20, which removes unresponsive detectors.
	These observations also provide a relative pointing solution for all of the detectors. This pointing solution is *not* used automatically, but instead provides input for manual generation of the BolometerPropertiesMap.

*RCW38 quick observations (aka very fast points)*
	This is a fast observation of RCW38 that does not scan every detector across the source, but does integrate down enough to get high signal-to-noise on RCW38 in a coadd of the full focal plane. This provides an atmospheric opacity correction that allows transfer of the main RCW38 calibration across atmospheric corrections; the stored data are the amplitude of RCW38 in each frequency band expressed in units of what RCW38 should have been with constant sky opacity based on the last long observation. This provides a correction to the temperature calibration and depends on the previous long observation.

*Centaurus A Pixel raster data*
	By fitting the shape of the CenA radio lobes, we can use the polarization structure in the source to fit for bolometer polarization angles. Like relative pointing, this calibration is not time-dependent and is not automatically applied, but is used to generate the BolometerPropertiesMap, which is done manually.

Proper sky calibration and generation of calibration frames depends on all of these succeeding.

Constructing Calibration Frames
===============================

Calibrating an elnod, calibrator stare, pixelraster or very fast point
observation into sky units is not appropriate if the observations are broken up
across fridge cycles or sky regions.  To avoid confusion in these not-so-rare
instances, we exclude the components of the calibration chain from calframes
built for observations that are part of the chain.

Please note that *the following section should be ignored in its entirety if you are looking only at field observations*.

*Elnod Calframe*
    The calframe for an elnod (say) includes only the relevant BolometerProperties
    and the result of that elnod analysis.  Calibrating an elnod to sky units thus
    requires also injecting the calibration analysis frame for an appropriate
    calibrator stare and a very fast point source observation by hand (recall that
    the VFP output frame includes the results of its parent pixelraster observation
    as well, so no need to separately include the pixelraster frame).

*Calibrator Calframe*
    The calframe for a calibrator stare includes:

    * preceding elnod analysis
    * analysis of that calibrator stare
    * bolometer properties

    In this case, sky calibration requires also injecting the appropriate VFP
    analysis frame.

*Pixelraster Calframe*
    The calframe for a calibration source pixelraster includes:

    * preceding elnod analysis
    * preceding calibrator analysis
    * analysis of that pixelraster observation
    * bolometer properties

    Here sky calibration requires also injecting the appropriate VFP analysis frame,
    injected *before* the pixelraster calframe to ensure that flux calibration is
    done off of that pixelraster (and not the one that was used to construct the VFP
    analysis frame).  This is because frames will be read in the order they are
    listed, so the new pixelraster keys will supercede those in the VFP frame.

*Very Fast Point Calframe*
    The calframe for a calibration source very fast point includes:

    * preceding elnod analysis
    * preceding calibrator analysis
    * preceding pixelraster of the same source
    * analysis of that VFP observation
    * bolometer properties

    No additional frames are required to calibrate to sky units.

*All Other Calframes*
    The calframe for any other kind of observation includes:

    * preceding elnod analysis
    * preceding calibrator analysis
    * preceding VFP analysis (including its parent pixelraster results)
    * bolometer properties

    Thus, to calibrate regular (non-calibration) observations into sky units,
    one only needs to inject this calframe into the pipeline.

