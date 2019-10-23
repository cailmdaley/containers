Documentation
-------------

The main documentation for the software is in the docs folder. A pretty, searchable copy can be found at https://scott.grid.uchicago.edu/ndhuang/spt3g_software.  This copy of the documentation is rebuilt nightly.

You can also build your own copy of the html documentation.  Detalis are below.

Dependencies
------------

This depends on Boost and cmake, as well as the usual Python packages. Some additional packages (NetCDF, in particular) will activate optional components of the code if installed. You also need a C++11 compiler. This software is designed to run and work on a variety of operating systems (all Linuxes, Mac OS X, and FreeBSD) and architectures (at least 64-bit x86 and POWER).

Minimum versions:

	- GCC >= 4.7 or clang >= 3.3
	- Boost >= 1.48
	- cmake >= 2.6
	- Python >= 2.7

On Ubuntu/Debian, you can install the non-Python dependencies, including the optional ones, by doing:

.. code-block:: shell

	apt-get install cmake libboost-all-dev libflac-dev libnetcdf-dev libfftw3-dev libgsl0-dev

On RHEL-type systems (SL, CentOS, etc.), do this:

.. code-block:: shell

	yum install cmake netcdf-devel boost-devel flac-devel fftw-devel gsl-devel
	
If your system defaults to Python 2, but you wish to use Python 3, please do the following:

	1. Install Python 3 *from the system package manager*
	2. Make sure the python-3 version of the Boost library is installed (on Ubuntu, this is part of the standard boost-python package referenced above)
	3. When you run cmake below, pass ``-DPYTHON_EXECUTABLE=`which python3```


Setup on sptcloud
-----------------

Note that on any RHEL6 system, you will need a newer compiler than ships with the OS. Please follow whatever directions apply at your site to achieve this. Alternately, if you have OASIS set up on your local system, run this before anything else (also works on several other operating systems, including RHEL7):

.. code-block:: shell

	eval `/cvmfs/spt.opensciencegrid.org/py3-v2/setup.sh`


How to Build
------------

To build:

.. code-block:: shell

	cd spt3g_software
	mkdir build
	cd build
	cmake ..
	make


To build the documentation in the build directory type:

.. code-block:: shell

	./env-shell.sh make docs

This will, hopefully, construct an html version of the documentation.  This builds the documentation in the build/docs folder.  Open build/docs/index.html in your favorite web browser.  You should at least read the quick start portion of the documentation before getting started.

Many of the modules have been parallelized with OpenMP.  OpenMP is a useful tool for speeding up sections of code on machines with multiple processing cores.  By default, OpenMP spawns as many threads as processors available.  When using a computing environment with many processors this is not optimal.  You will want to set the OMP_NUM_THREADS environment variable to be a small integer number.  You will want to add this to your .profile or .bashrc file (if you use cshell, uhh, you can figure it out).

.. code-block:: shell

	export OMP_NUM_THREADS=2

Version Control Hygiene
-----------------------

You can use two mechanisms to access the repository: git and SVN. The following is a brief overview of how to use these in a way that your collaborators will appreciate.

Git
===

To initially check out the repository:

.. code-block:: shell

	git clone https://user@github.com/SouthPoleTelescope/spt3g_software.git

To update your checkout (the --rebase is important, especially if you have local changes):

.. code-block:: shell

	git pull --rebase

To send your changes back:

.. code-block:: shell

	git diff files_to_commit <- Examine this
	git commit files_to_commit
	git push


SVN
===

To initially check out the repository:

.. code-block:: shell

	svn co https://user@github.com/SouthPoleTelescope/spt3g_software/trunk spt3g_software

To update your checkout:

.. code-block:: shell

	svn up

To send your changes back:

.. code-block:: shell

	svn diff files_to_commit <- Examine this
	svn ci files_to_commit

