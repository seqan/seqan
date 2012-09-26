cmake_minimum_required (VERSION 2.6)
project (core_apps_seqcons)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(seqcons seqcons.cpp)

install(TARGETS seqcons
        DESTINATION bin)
