
find_package(GSL)
find_package(FFTW)
find_package(HDF5)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()


if(FFTW_FOUND AND GSL_FOUND)
	add_library(frbutils SHARED
	  src/frbutils.cpp src/python.cpp
	)
	target_link_libraries(frbutils core mapmaker todfilter ${FFTW_LIBRARIES} ${GSL_LIBRARIES})
	include_directories(${GSL_INCLUDES} ${FFTW_INCLUDES})

	add_spt3g_test(testsigest SLOWTEST)
	link_python_dir()
else()
	message(STATUS "One of FFTW or GSL not found. Skipping build of frb code")
endif()
