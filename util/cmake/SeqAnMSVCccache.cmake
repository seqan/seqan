# https://github.com/ccache/ccache/wiki/MS-Visual-Studio#usage-with-cmake

find_program (ccache_exe ccache)

if (ccache_exe)
  message (STATUS "Using ccache: ${ccache_exe}")
  file (COPY_FILE
    ${ccache_exe} ${CMAKE_BINARY_DIR}/cl.exe
    ONLY_IF_DIFFERENT)

  # By default Visual Studio generators will use /Zi which is not compatible
  # with ccache, so tell Visual Studio to use /Z7 instead.
  set (CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")

  set (CMAKE_VS_GLOBALS
    "CLToolExe=cl.exe"
    "CLToolPath=${CMAKE_BINARY_DIR}"
    "UseMultiToolTask=true"
    "DebugInformationFormat=OldStyle"
    CACHE STRING "Visual Studio globals for ccache"
  )
else ()
  message (STATUS "Not using ccache")
endif()
