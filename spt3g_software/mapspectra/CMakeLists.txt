
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()



if(FFTW_FOUND)
	add_library(mapspectra SHARED
	  src/bispec_est.cpp src/inpaint_quick.cpp src/python.cpp 
	  src/apodmask.cpp
	)
	target_link_libraries(mapspectra core coordinateutils ${FFTW_LIBRARIES})
	include_directories(${FFTW_INCLUDES})
else()
	message(STATUS "FFTW not found. Skipping build of spectra code.")
endif()

link_python_dir()

add_spt3g_test(test_powerspectrum SLOWTEST)
add_spt3g_test(test_inpaint SLOWTEST)
