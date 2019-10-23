------------------
timestreamflagging
------------------

The timestream project contains code for flagging timestreams based on the 3 following subsets:

*Physical properties of the bolometers*
	Located in miscflagmodules; includes filtering/grouping by band/wafer, pixel pairs, etc.
  
*Missing Data/ NaNs*
  Located in miscflagmodules; missing keys, bad HK, latched/overbiased/ bolos, NaNs etc.
  
*Noise*
  Located in noiseflagging; variance flagging, High Q noise
  Located in glitchfinding; glitches
  
The project also contains general modules for generating and removing the flags as well as collecting the flag stats; located in flaggingutils.
