# - Find HDF5    
# Find the native HDF5 includes and library 
#
# HDF5_INCLUDES - where to find fftw3.h
# HDF5_LIBRARIES - List of libraries when using HDF5.
# HDF5_FOUND - True if HDF5 found.
if (HDF5_INCLUDES)
# Already in cache, be silent 
set (HDF5_FIND_QUIETLY TRUE)
endif (HDF5_INCLUDES)
find_path (HDF5_INCLUDES hdf5.h PATHS $ENV{HDF5_DIR}/include/)
find_library (HDF5_LIBRARIES NAMES hdf5 PATHS $ENV{HDF5_DIR}/lib/)
find_library (HDF5_HL_LIBRARIES NAMES hdf5_hl PATHS $ENV{HDF5_DIR}/lib/)
# handle the QUIETLY and REQUIRED arguments and set HDF5_FOUND to TRUE if
# all listed variables are TRUE 
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (HDF5 DEFAULT_MSG HDF5_LIBRARIES HDF5_INCLUDES)
mark_as_advanced (HDF5_LIBRARIES HDF5_HL_LIBRARIES HDF5_INCLUDES)


