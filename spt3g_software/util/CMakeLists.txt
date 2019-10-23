if(GSL_FOUND AND FFTW_FOUND)

  if (HDF5_FOUND)
    message(STATUS "HDF5 " ${HDF5_INCLUDES} " " ${HDF5_LIBRARIES})
    message(STATUS "HDF5 HL " ${HDF5_HL_LIBRARIES})
    list(APPEND FILTER_EXTRA_SRCS src/hdfutils.cpp)
    list(APPEND FILTER_EXTRA_LIBS ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    include_directories(${HDF5_INCLUDES})
  else()
    message(STATUS "Could not find HDF5 libs, building without HDF5 tools")
  endif()



  add_library(util SHARED
    ${FILTER_EXTRA_SRCS}
    src/python.cpp 
    )

  target_link_libraries(util todfilter ${FILTER_EXTRA_LIBS})
  
  link_python_dir()

  #Tests
  add_spt3g_test(util_test_fixture)
endif()
