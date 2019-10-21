-------------------
Running Simulations
-------------------

This document is meant to be an overview of how to start simulating the data.  It assumes
you are vaguely familiar with the map making procedure. 


Overview
========
Running noiseless simulations consists of:

* Generating a simulated map of the sky
* Converting that simulated map to a memory friendly, usable format
* Using that map to create fake TOD.

Useful Resources
================
For most simulations you will want a map that has roughly the same statistical properties of the real sky.  This means that using some distribution of C_ls you will want to convert it to map space using spherical harmonic transforms rather than flat sky FFTs.  The most popular way to do that in the CMB community is to make maps using the Healpix pixelization scheme and the associated libraries.  Here are the libraries/programs we use for that:

* Healpix: An equal area pixelization of the sphere and a library for doing map to alm conversions. 
* healpy:  python bindings for healpix
* lenspix: A program for simulating the CMB and do gravitational lensing of the CMB. 

Generating Simulated Skies
==========================
CR kindly got lenspix compiling and wrote a python routine for running lenspix.  This is an extremely useful resource for generating simulations of the sky.

Not in the future, but for the time being, that repository only compiles on Scott and not
Amundsen.

First clone the repository with:

.. code-block:: shell

        git clone git@github.com:SouthPoleTelescope/Simulations.git  

