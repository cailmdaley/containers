
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()


if(FFTW_FOUND AND GSL_FOUND)
	add_library(todfilter SHARED
	    src/filtermodules.cpp src/polyutils.cpp src/python.cpp src/numericutils.cpp
	    src/fftutils.cpp
	)
	target_link_libraries(todfilter core coordinateutils ${FFTW_LIBRARIES} ${GSL_LIBRARIES})
	include_directories(${GSL_INCLUDES} ${FFTW_INCLUDES})

	add_spt3g_test(filter_test_fixture SLOWTEST)
	add_spt3g_test(downsampler_test SLOWTEST)
        add_spt3g_test(timeconstant_filter_test SLOWTEST)

	add_executable(filterperf tests/filtertest.cpp)
	target_link_libraries(filterperf todfilter)
else()
	message(STATUS "One of FFTW or GSL not found. Skipping build of timestream filtering code")
endif()

link_python_dir()
