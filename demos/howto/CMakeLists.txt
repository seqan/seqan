cmake_minimum_required (VERSION 2.6)
project (seqan_core_demos_howto)

seqan_setup_apps(core_demos_howto)
seqan_add_all_executables(core_demos_howto howto_)

seqan_add_all_subdirectories()
