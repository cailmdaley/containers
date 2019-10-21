Mapmaker Mapmaker Make Me A Map
-------------------------------

For the math see:
https://pole.uchicago.edu/spt3g/img_auth.php/MapMakingTutorialEarlySpt3g.pdf

Before looking at an example map maker code, it's important to understand a few of the structures in detail that are going to be used.  After reading this document you should check out spt3g_software/std_processing/mapmakers/field_autoproc.py.  For each of the module used in that example there should be detailed documentation.  The process of writing good documentation requires some iteration since it's often hard to know if the information was conveyed clearly.  If any of the documentation isn't clear let me know.

Timestream Flagging
-------------------

The mapmaker will make a map with all the detectors you supply.  The flagging and mapmaking steps are separate.  Flag the data before going on to the map making step

A percentage of our detectors are going to have crap data.  We would like to remove said crap bolometers.  This is accomplished in three steps:

1) Characterize the data (estmate noise, count glitches etc.)
2) Select which detectors to flag based off of the characterized data
3) Remove said crap timestreams.

Steps 1 and 2 are split so that it is easy to modify flagging options based off of the characterized data.  This also facilitates being able to look back at the data and figure out why exactly a bolometer is being flagged.

Flags are stored in a G3MapVectorString that has the form: 
  {bolometer_id: ['BadCalibrator', 'Glitchy', 'ActuallyAKitten']}

The bolometer id maps to a list of reasons why the bolometer is flagged.  

When developing code to flag timestreams we are trying to not force 
flagging decisions on people using intermediate data products.
By that I mean, if we need to compute something in order to flag the timestreams, 
we should compute the thing, store it in the frame, and then use a separate, 
later, Module to actually set the flags. 

As an example, with the glitch finder, we have a module that finds the number of 
glitches above a threshold, and a seperate module for flagging based off of those stored values.

Remember that data products in the frames are supposed to be immutable.  
Do *NOT* do something like frame['Flags'][bolo_id].append("NoisyTimestream").  
BAD BAD BAD, NO. SHAME. BAD

You will need to copy the flags, append the flag you want, 
delete the previous flags in the frame and  then add the new flags.  
Because having you do that every time would be insane please use the add_flag
convenience function (this lives in timestreamflagging/python/flaggingutils.py)

All of the generic timestream flagging code lives in:

timestreamflagging/python/flaggingutils.py

There are some pipe segments with combined flagging routines located in:
std_processing/flagsegments.py

Once we have specified the flags we just call:

timestreamflagging.flaggingutils.RemoveFlaggedTimestreams

and supply it the flags.  


Adding the Map Frame to the Pipeline And Why That's Important
=============================================================
All of the modules in the mapmaker that require having a map (like the map binner or the simulated TOD filler) or require knowing what patch of sky a map covers (the pointing code) get that information from map frames going through the pipeline.  Each map frame in the code has a field 'Id' that stores a unique id for the map.  The modules that need information get passed this Id so they know which map frame to grab information from.  The map frame structure is detailed below.

There is a module in mapmaker.mapmakerutils called MapInjector that will take a list of G3SkyMap objects and the map id.  It packages these into the appropriate map frame format and will add it to the pipeline appropriately.  All of the processing steps (like figuring out which map pixel each detector is pointing at for each TOD sample) that rely on sky map information get that informtion from the Map frame as it goes through the pipeline.

For creating the maps to supply to the map injector, you can either use constructor of the map or call std_processing.CreateSourceMapStub which stores the Ra and Dec of the sources/fields we observe.

Map Frame Structure
===================
Each map frame represents a single patch of sky with a single meaning.  By that I mean, the point source mask will live in a separate map frame from the map we are generating from the map maker.  The sim sky realization will live in a different map frame as well.  

The map frame has a very set structure with specific keys mapping to specific data types.

Required:

"Id" -> G3String:  The string id of the map stored in this frame.  This is used to identify the map to mapmaking modules
"T" -> G3SkyMap:  The intensity map. Or if the data is not physical (like a point source mask)  this is where it is stored.

Optional:
"Q" -> G3SkyMap:  the Q stokes parameter map following IAU conventions (hopefully)
"U" -> G3SkyMap:  the U stokes parameter map following IAU conventions (hopefully)  [IAU not the stupid healpix standard for jerks]

"Wpol" -> G3SkyMapWeights:  The weight map for a polarized sky (if this is included hopefully Q and U are as well)
"Wunpol" -> G3SkyMapWeights:  The weight map for an unpolarized sky 

Anything involving pointing will try to grab the map parameters from the map stored in "T".  The map binner will try to bin the data for every available key in the input map frame.
By that I mean, it "T", "Q", and "U" keys exist in the Map frame, it will bin "T", "Q", and "U", but skip the "Wpol" and "Wunpol" because they weren't present.  If "T" and "Wunpol" are present in the frame it will just bin intensity and unpolarized weight information.

There are two types of sky map objects that can be passed to the mapmaker FlatSkyMap and CutSkyHealpixMap.  The mapmaker should work with either of them.

Detector Pointing
=================
Where each detector is pointing is calculated from the boresight pointing of the telescopes and how far the pointing of each detector is offset from that boresight pointing.  This effectively involves a lot of rotations.  I will write rotation matrices as Rx(theta) which means a rotation about the x-axis by angle theta.  Also, because some of the time euclidean coordinates are useful let's set Ra=Dec=0 to (1,0,0) and use a right hand coordinate system.

If the telescope is pointing at Ra=Dec=0 we can calculate where the detector is pointings as:

Detector Pointing = Rz(x_offset) Ry(y_offset) (1, 0, 0)

Using the definitions of the x and y offsets stored in the bolometer properties (modulo sign)

We then have a rotation that maps (1,0,0) to the pointing on the sky.

We used to estimate this as:
R_trans = Rz(Ra) Ry(Dec)

This has some issues though.  The detector pointing is fixed in local az/el coordinates not RaDec coordinates.  Because Earth's pole has drifted with time we were screwing up this measurement.

So now we estimate R_trans from two sets of points Boresight AzEl->Boresight RaDec and slightly offset Boresight AzEl->slightly offset Boresight RaDec.  

The difference in the two calculations is that the new code also accounts for this rotation about the boresight.

In each frame we store R_trans for each time sample.  This is the rotation that takes (1,0,0) to boresight RaDec.  If you have a map in Galactic coordinates you just need to apply a rotation to the RaDec rotation that converts it to the rotation that takes (1,0,0) to the (l,b) coordinates.

You might also need to account for the polarization angle rotation caused by this rotation around boresight.  The various map making/sim code does not include this estimation of the polarization angle rotation by default.  You need to pass an optional argument to have it handle the polarization angle rotation.  This is important if you are using galactic coordinates.


If you care about implementation details read this, if not, skip it.  The rotations are stored as quaternions.  Quaternions cannot handle improper rotations.  Unfortunately AzEl cooridnates are parity flipped from RaDec coordinates.  When converting AzEl coordinates to the quaternion form we flip the sign of the Elevation.  So Az=0, El = 90 will be at (0,0,-1).  That is so all of our coordinate systems have a consistent handedness.

Example Code Location
----------------------
spt3g_software/std_processing/mapmakers/field_autoproc.py


Working With G3SkyMap objects
-----------------------------
G3SkyMap objects are a base class that has two children: CutSkyHealpixMap and FlatSkyMap.

CutSkyHealpixMaps are maps that use the healpix pixelization, but do not store the entire celestial sphere, they store a cut out.  The maps themselves are 1D arrays that just stores the pixel information.  The CutSkyHealpixMaps have another object stored inside them called HealpixHitPix that handles mapping from the full sky healpix pixelization to the cut sky pixelization.  In general working with the data in the stored buffer form is difficult because the layout of the information does not cleanly map to a sky area.  In general when working with these objects, you either want to use the methods of the CutSkyHealpix objects to get the information you want, or you want to convert the cut sky map into a full sky map and use healpy's routines on them.  healpy is the python library that has the healpix functions in them.  Healpix is the standard pixelization of the sphere used by the CMB community.

FlatSkyMaps are 2d arrays used to store a flat sky projection of a section of the celestial sphere.  When indexing them you can use 2d indexes or 1d indexes where the 1d index is (x + y * n_columns).  Most of the routines just return the 1d index.  The flat sky maps behave almost like numpy arrays, but not quite.  There are complicated implementation reasons for this and it's not easily changeable.  What you can do, however, is make python know you want to treat the flat sky map like a numpy array.  Say flat_sky_map is a FlatSkyMap, then if you use the function: np.asarray(flat_sky_map) the returned object will be a numpy array.  When you use np.asarray it does NOT copy the array, it's pointing to the same section of memory.  Any modifications you make to the array object returned by that call will modify the underlying map.

The FlatSkyMap object stores extra meta information about the patch of sky it covers.  The map projection, the center of the map, the resolution, the dimensions, etc.  The resolution tries to encompass the area of a pixel, but remember that the shape of the pixels is function of the flat sky projection used.

Map Projections
---------------
So, most of the map projections in the software were carried over from SPTsz.  In practice, the map projections that see use are proj0, proj1, proj5 and the healpix pixelization.  Proj0 and Proj1 are both great because they are simple.  Proj5 does a good job of minimally distorting the map for power spectrum estimation.  And healpix is great if you decide you don't want to be lazy and actually use spherical harmonic transforms for n-point estimation.


Proj0/ProjSansonFlamsteed:
https://en.wikipedia.org/wiki/Sinusoidal_projection

Proj1/ProjCAR: 
Projects the Ra/Dec values to the x/y cooridnates directly:
x=(ra-ra0);
y=(-1.*dec+dec0);

Proj2/ProjSIN:
Orthographic https://en.wikipedia.org/wiki/Orthographic_projection_in_cartography

Proj4/ProjStereographic:
https://en.wikipedia.org/wiki/Stereographic_projection

Proj5/ProjLambertAzimuthalEqualArea:
https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection
We technically use the oblique lambert azimuthal equal-area projection, but jeeze, it's already a lot to type.

Proj7/Proj8 are nothing right now

Proj9: is the map projection BICEP uses for no good reason.
