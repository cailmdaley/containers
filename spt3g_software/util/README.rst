-----------
Util
-----------

The util folder contains modules with generally useful functions for interacting
with and manipulating G3-specific and generic data. The modules are separated by 
the types of functions they contain, with a brief summary of each as follows:

*extract_data*
	This file contains classes that make interacting with data in G3 files easier.

*files*
	This module contains functions for reading/writing specific kinds of files, such as pickles, camb files, and fridge logs.

*fitting*
	This module contains different fitting and minimization functions that might be useful for maps or timestreams.

*frame_combiner*
	This module contains functionality to combine (add, subtract) frames of the same Id in a G3 File and creates a new frame from that combined data.

*g3curses*
	Wraps some 3G functionality into the curses library.

*genericutils*
	This module contains general functions that extend typical python/numpy operations on basic python data structures.

*healpix_tools*
	A collection of tools for handling healpix maps.

*math*
	This module contains various useful math functions which cannot be found in other packages. These are generic tools which do not use any 3G specific data structures.

*stats*
        This module contains useful statistical functions, like gaussian moments, median absolute deviation, outlier filtering, etc.

*timestream_psd*
        This module contains functionality for plotting the for plotting many detector power spectral densities together.



