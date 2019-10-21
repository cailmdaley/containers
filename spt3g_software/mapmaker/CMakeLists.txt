if(Boost_GEOMETRY_AVAIL)
    set(ABERRATION_SRCS src/aberrationutils.cpp src/aberrationmodules.cpp)
else(Boost_GEOMETRY_AVAIL)
    set(ABERRATION_SRCS )
endif(Boost_GEOMETRY_AVAIL)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()


if(GSL_FOUND)
  add_library(mapmaker SHARED
    src/pointingutils.cpp src/mappingutils.cpp
    src/polarizationutils.cpp src/mappingmodules.cpp
    src/detectormodules.cpp
    src/todmapconverters.cpp
    src/python.cpp
    ${ABERRATION_SRCS}
    )

  include_directories(${GSL_INCLUDES})
  target_link_libraries(mapmaker core coordinateutils calibration ${GSL_LIBRARIES})

else()
  message(STATUS "GSL is missing")
endif()

link_python_dir()
