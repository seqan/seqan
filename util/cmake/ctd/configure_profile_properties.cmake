# If possible, get latest change date from SeqAn SVN.
message(STATUS "SEQAN_SOURCE_DIR ${SEQAN_SOURCE_DIR}")
message(STATUS "SEQAN_DATE ${SEQAN_DATE}")
if (SEQAN_DATE)
	string(REGEX REPLACE "^([0-9]+)-([0-9]+)-([0-9]+)_([0-9]+):([0-9]+).*"
    "\\1\\2\\3\\4\\5" SEQAN_LAST_CHANGED_DATE "${SEQAN_DATE}")
  set (CF_SEQAN_VERSION ${SEQAN_VERSION_STRING}.${SEQAN_LAST_CHANGED_DATE})
else ()
  set (CF_SEQAN_VERSION "${SEQAN_VERSION_STRING}")
endif ()

if (NOT CTD_PLUGIN_PACKAGE)
  set (CTD_PLUGIN_PACKAGE "de.seqan")
endif ()

if (NOT CTD_PLUGIN_NAME)
  set (CTD_PLUGIN_NAME  "SeqAn")
endif ()

# Actually configure the file.
configure_file (${SEQAN_SOURCE_DIR}/util/cmake/ctd/plugin.properties.in ${WORKFLOW_PLUGIN_DIR}/plugin.properties)
