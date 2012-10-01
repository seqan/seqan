cmake_minimum_required (VERSION 2.6)
project (seqan_core)

# We have to register our own include directory here, though.
seqan_setup_includes(include SeqAnCore)

seqan_add_all_subdirectories()  # include_directory(*)
