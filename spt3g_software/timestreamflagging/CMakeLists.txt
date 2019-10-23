add_library(timestreamflagging SHARED src/removeflags.cpp src/python.cpp)
target_link_libraries(timestreamflagging core)
link_python_dir()
