# Configure SeqAn, enabling features for libbz2 and zlib.
find_package (OpenMP)
find_package (ZLIB)
find_package (BZip2)
find_package (SeqAn REQUIRED)

# Add include directories, defines, and flags for SeqAn (and its dependencies).
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Build the program and link it against the SeqAn dependency libraries.
add_executable (main main.cpp)
target_link_libraries (main ${SEQAN_LIBRARIES})
