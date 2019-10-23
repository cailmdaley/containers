find_package(MPI)
if (MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
   add_executable(mpibatch.x src/mpibatch.c)
   target_link_libraries(mpibatch.x ${MPI_LIBRARIES})

   if(MPI_COMPILE_FLAGS)
	set_target_properties(mpibatch.x PROPERTIES
    	COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
   endif()


   if(MPI_LINK_FLAGS)
	set_target_properties(mpibatch.x PROPERTIES
	    LINK_FLAGS "${MPI_LINK_FLAGS}")
   endif()



   add_executable(mpibatch_cl.x src/mpibatch_cl.c)
   target_link_libraries(mpibatch_cl.x ${MPI_LIBRARIES})

   if(MPI_COMPILE_FLAGS)
	set_target_properties(mpibatch_cl.x PROPERTIES
    	COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
   endif()

   if(MPI_LINK_FLAGS)
	set_target_properties(mpibatch_cl.x PROPERTIES
	    LINK_FLAGS "${MPI_LINK_FLAGS}")
   endif()

else()
	message(STATUS "Could not find MPI libraries")
endif()


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
    add_executable(openmp_test.x src/openmp_test_program.c)
endif()

link_python_dir()
