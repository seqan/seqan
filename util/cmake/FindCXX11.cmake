if(__FIND_CXX11_CMAKE__)
  return()
endif()
set(__FIND_CXX11_CMAKE__ TRUE)

# Visual Studio 2008 (vs9) doesn't seem to support C++11 directly (only as TR1)
if (MSVC AND MSVC_VERSION GREATER 1500)
  set(CXX11_FOUND 1)
  return ()
endif (MSVC AND MSVC_VERSION GREATER 1500)

include(CheckCXXCompilerFlag)
enable_language(CXX)

check_cxx_compiler_flag("-std=c++11" CXX11_FOUND)
if (CXX11_FOUND)
  set (CXX11_CXX_FLAGS "-std=c++11")

  # Tested on Mac OS X 10.8.2 with XCode 4.6 Command Line Tools
  # Clang requires this to find the correct c++11 headers
  if (CMAKE_HOST_APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
     set (CXX11_CXX_FLAGS "${CXX11_CXX_FLAGS} -stdlib=libc++ -Qunused-arguments")
  endif (CMAKE_HOST_APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))

else (CXX11_FOUND)

  check_cxx_compiler_flag("-std=c++0x" CXX11_FOUND)
  if (CXX11_FOUND)
    set (CXX11_CXX_FLAGS "-std=c++0x")
  endif (CXX11_FOUND)

endif (CXX11_FOUND)

