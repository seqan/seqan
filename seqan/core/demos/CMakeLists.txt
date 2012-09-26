cmake_minimum_required (VERSION 2.6)
project (seqan_core_demos)

seqan_setup_apps(core_demos)

seqan_add_all_executables(core_demos)
#seqan_add_all_executables()

seqan_add_all_subdirectories()
