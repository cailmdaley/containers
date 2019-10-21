-----------------
GCP Archive Files
-----------------

GCP writes its own streaming file format (which we call "archive" or "arc" files).  Arc files generally contain everything that does not come from the Iceboards, such as telescope pointing.  Scanifaction is the process of stitching arc files and bolometer data together.  However, not all data from the arc files is added to the final G3 files, so it is sometimes necessary to access the arc files directly.  The purpose of this document is to give an overview of some techniques for working with arc files.

An archive file contains a set of sequential data frames (usually 1000).  Each data frame contains all the data (from a variety of sources) for 1 second of data.  Individual data streams are stored in a hierarchy of registers (see the `GCP documentation <https://southpoletelescope.github.io/GCP/registers.html>`_ for a list of available registers).

.. contents:: Contents

Reading Arc Files
-----------------

Arc files can be read frame-by-frame using ``spt3g.gcp.ARFileReader``.  This behaves exactly like ``spt3g.core.G3Reader``, but reads arc files instead.  For example, 

.. code-block:: python

	In [1]: from spt3g import core, gcp

	In [2]: af = gcp.ARCFileReader('/spt/data/arc/20180106_051236.dat')

	In [3]: fr = af(None)[0]
	
	In [4]: print(fr)
	Frame (GcpSlow) [
	"antenna0" (spt3g.core.G3MapFrameObject) => 5 elements
	"array" (spt3g.core.G3MapFrameObject) => 15 elements
	]

	In [5]: az = fr['antenna0']['tracker']['actual'][0]
	
	In [6]: el = fr['antenna0']['tracker']['actual'][1]
	
	In [8]: fr['array']['cryo']['temperature']
	Out[8]: spt3g.core.G3VectorFrameObject([16 elements, 16 elements, 16 elements])
	
	In [9]: fr['array']['cryo']['temperature'][0]
	Out[9]: spt3g.core.G3VectorDouble([0.277949, 0.307301, 1.01168, 1.62593, 3.56612, 3.59573, 3.53362, 23.9017, 26.6354, 24.9882, 0.282698, 0.287977, 0.315117, 3.10634, 2.77282, 36.1949])

Extracting Specific Keys
========================

For a variety of reasons, it is usually more useful to use a higher level data reader.  ``spt3g.util.extractdata.extract_keys`` will read in a set of registers from a set of arcfiles, and return the concatenated results, after applying low-level calibrations (conversions to useful units).  For example,

.. code-block:: python

	In [1]: from spt3g.util import extractdata
	
	In [2]: keys = {'az': ['antenna0', 'tracker', 'actual', 0], 'el': ['antenna0', 'tracker', 'actual', 1], 'azerror': ['antenna0', 'tracker', 'errors', 0]}
	
	In [3]: data = extractdata.extract_keys('/spt/data/arc/20180106_051236.dat', keys)
	
	In [4]: data.keys()
	Out[4]: dict_keys(['az', 'el', 'azerror'])
	
	In [5]: data['az']
	Out[5]: 
	array([ 6.27933628,  6.27933628,  6.27933591, ...,  5.27599024,  # NB: in native G3 units
	        5.27591908,  5.27584942])

This can also be done in a pipeline using ``spt3g.util.extractdata.MultiAccumulator``.  

.. code-block:: python

	In [1]: from spt3g.util import extractdata
	
	In [2]: from spt3g import core, gcp
	
	In [3]: keys = {'az': ['antenna0', 'tracker', 'actual', 0], 'el': ['antenna0', 'tracker', 'actual', 1], 'azerror': ['antenna0', 'tracker', 'errors', 0]}
	
	In [4]: mult_accumulator = extractdata.MultiAccumulator(keys)
	
	In [5]: pipe = core.G3Pipeline()
	
	In [6]: pipe.Add(gcp.ARCFileReader, filename = '/spt/data/arc/20180106_051236.dat')
	
	In [7]: pipe.Add(gcp.ARCExtract)  # For unit conversions.
	
	In [8]: pipe.Add(mult_accumulator)
	
	In [9]: pipe.Run()
	
	In [10]: data = mult_accumulator.extract_values()
	
	In [11]: data['az']
	Out[11]: 
	array([ 6.27933628,  6.27933628,  6.27933591, ...,  5.27599024,  # NB: in native G3 units
                5.27591908,  5.27584942])

Finding Arc Files
-----------------

Arc files are typically 16 minutes (1000 seconds, to be precise), and contain 1000 frames.  They are stored by their start date (YYYYmmdd_hhmmss.dat).  On amundsen and scott, they can be found in ``/spt/data/arc/``.  If you want a simpler interface (or archive data from a specific time range), us ``spt3g.std_processing.ARCTimerangeReader``:

.. code-block:: python

	In [1]: from spt3g import core, std_processing
	
	In [2]: from spt3g.util import extractdata
	
	In [3]: keys = {'az': ['antenna0', 'tracker', 'actual', 0], 'el': ['antenna0', 'tracker', 'actual', 1], 'azerror': ['antenna0', 'tracker', 'errors', 0]}
	
	In [4]: mult_accumulator = extractdata.MultiAccumulator(keys)
	
	In [5]: start = core.G3Time('20180216_001000')
	
	In [6]: stop = core.G3Time('20180216_001200')
	
	In [7]: pipe = core.G3Pipeline()
	
	In [8]: pipe.Add(std_processing.ARCTimerangeReader, start_time = start, stop_time = stop, basedir ='/spt/data/arc')
	
	In [9]: pipe.Add(mult_accumulator)
	
	In [10]: pipe.Run()
	
	In [11]: data = mult_accumulator.extract_values()
	
	In [12]: data['az']
        Out[12]: 
	array([ 4.70857313,  4.70857313,  4.7085735 , ...,  4.70857313,
                4.70857313,  4.70857313])
