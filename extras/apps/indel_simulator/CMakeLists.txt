cmake_minimum_required (VERSION 2.6)
project (extras_apps_indel_simulator)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(indel_simulator indel_simulator.cpp indel_simulator.h)
