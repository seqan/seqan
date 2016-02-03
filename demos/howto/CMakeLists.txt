cmake_minimum_required (VERSION 2.6)
project (seqan_demos_howto CXX)

seqan_setup_apps(demos_howto)
seqan_add_all_executables(demos_howto howto_)

seqan_add_all_subdirectories()
