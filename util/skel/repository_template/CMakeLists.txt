cmake_minimum_required (VERSION 2.6)
project (seqan_%(CMAKE_PROJECT_NAME)s)

# We have to register our own include directory here, though.
seqan_setup_includes(include SeqAn%(REPOSITORY_PSEUDO_TARGET_NAME)s)
# SeqAnCore is added include in ../core/CMakeLists.txt.  We only have to
# register that we want it here, too.  Same for SeqAnExtras.
seqan_make_seqan_available(SeqAnCore)
seqan_make_seqan_available(SeqAnExtras)

seqan_add_all_subdirectories()  # include_directory(*)
