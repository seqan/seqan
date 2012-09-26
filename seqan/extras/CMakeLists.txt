cmake_minimum_required (VERSION 2.6)
project (seqan_extras)

# We have to register our own include directory here, though.
seqan_setup_includes(include SeqAnExtras)
# SeqAnCore is added include in ../core/CMakeLists.txt.  We only have to
# register that we want it here, too.
seqan_make_seqan_available(SeqAnCore)

seqan_add_all_subdirectories()  # include_directory(*)
