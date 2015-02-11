# Example CMakeLists.txt file that uses the FindSeqAn.cmake module for
# building a SeqAn-based app.

project (raw_cmake_project)
cmake_minimum_required (VERSION 2.8.2)

# Place binaries into "bin" directory.
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")

# Go to "src" subdirectory.
add_subdirectory (src)
