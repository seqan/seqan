# Example CMakeLists.txt file that uses the seqan-config.cmake module for
# building a SeqAn-based app.

project (raw_cmake_project CXX)
cmake_minimum_required (VERSION 3.0.0)

# Place binaries into "bin" directory.
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")

# Go to "src" subdirectory.
add_subdirectory (src)
