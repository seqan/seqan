# If possible, get latest change date from SeqAn SVN.
find_package(Subversion)
message(STATUS "SEQAN_SOURCE_DIR ${SEQAN_SOURCE_DIR}")
if (Subversion_FOUND AND EXISTS ${SEQAN_SOURCE_DIR}/.svn)
  file (TO_CMAKE_PATH "${SEQAN_SOURCE_DIR}" _SEQAN_SOURCE_DIR)
  Subversion_WC_INFO (${_SEQAN_SOURCE_DIR} SEQAN)
  string(REGEX REPLACE "^([0-9]+)-([0-9]+)-([0-9]+) ([0-9]+):([0-9]+).*"
    "\\1\\2\\3\\4\\5" SEQAN_LAST_CHANGE_DATE "${SEQAN_WC_LAST_CHANGED_DATE}")
  set (CF_SEQAN_VERSION ${SEQAN_VERSION_STRING}.${SEQAN_LAST_CHANGE_DATE})
else ()
  set (CF_SEQAN_VERSION "${SEQAN_VERSION_STRING}")
endif ()

# Actually configure the file.
configure_file (${SEQAN_SOURCE_DIR}/util/cmake/ctd/plugin.properties.in ${WORKFLOW_PLUGIN_DIR}/plugin.properties)