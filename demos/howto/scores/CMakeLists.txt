cmake_minimum_required (VERSION 2.6)
project (seqan_demos_howto_scores CXX)

seqan_setup_apps(demos_howto_scores)
seqan_add_all_executables(demos_howto_scores howto_)

seqan_add_all_subdirectories()
