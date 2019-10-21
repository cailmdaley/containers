-----
beams
-----
The beams directory contains all the scripts needed to make a beam product. 

The first step in making a beam product is making maps. The only required maps are planet maps, but there are also scripts for making high resolution point source maps (for high ell beam) and for making field maps with corresponding mock observations of Planck (for low ell beam). 

Once the desired maps exist on disk, the top level analysis script, products/beam00.py can be run. This script produces plots of the ell-space beam along with the final text file product with a B_l at each frequency. 
