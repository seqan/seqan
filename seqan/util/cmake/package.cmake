INCLUDE(InstallRequiredSystemLibraries)

# NOTE that you have to run "make docs" before running cpack.  The reason
# is that we cannot add dependencies to the install target at the moment.
# See: http://public.kitware.com/Bug/view.php?id=8438

# ===========================================================================
# Archive Packages (.tar & .tar.bz2)
# ===========================================================================

SET(CPACK_GENERATOR "ZIP;TBZ2;DEB")
SET(CPACK_PACKAGE_NAME "seqan")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqAn - The C++ library for sequence analysis.")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
SET(CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

seqan_get_version()

SET(CPACK_PACKAGE_VERSION "${SEQAN_VERSION}")
SET(CPACK_PACKAGE_VERSION_MAJOR "${SEQAN_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${SEQAN_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${SEQAN_VERSION_PATCH}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

# Should be the last include.
INCLUDE(CPack)

