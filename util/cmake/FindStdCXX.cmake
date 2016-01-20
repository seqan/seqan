# Selects the highest available c++ standard
#
# I.e. if c++14 and c++11 is supported, use c++14 standard
#
# === Variables ===
#
# STD_CXX_FLAG     The coresponding compiler flag for the selected c++ standard.
#                  e.g. sets "-std=c++11" for the gcc compiler if only c++11 is
#                  supported.
# CXX11_FOUND      True, if the compiler supports C++11
# CXX14_FOUND      True, if the compiler supports C++14
# CXX11_STL_FOUND  @deprecated

if (STD_CXX_FLAG OR __FIND_STD_CXX_CMAKE__)
    return()
endif()
set(__FIND_STD_CXX_CMAKE__ TRUE)


macro (_seqan_cxx_standard_windows)

    # Visual Studio 2008 (vs9) doesn't seem to support C++11 directly (only as TR1)
    # Visual Studio 2010 (vs10) doesn't support C++11 STL.
    if (MSVC AND MSVC_VERSION GREATER 1600)
        set(CXX11_FOUND 1)
        set(CXX11_STL_FOUND 1)

        return ()
    endif (MSVC AND MSVC_VERSION GREATER 1600)

endmacro(_seqan_cxx_standard_windows)


macro (_seqan_cxx_standard_gcc)

    include(CheckCXXCompilerFlag)
    enable_language(CXX)

    check_cxx_compiler_flag("-std=c++11" CXX11_FOUND)
    check_cxx_compiler_flag("-std=c++14" CXX14_FOUND)

    if (CXX11_FOUND)
        set (STD_CXX_FLAG "-std=c++11")

        # Tested on Mac OS X 10.8.2 with XCode 4.6 Command Line Tools
        # Clang requires this to find the correct c++11 headers
        if (CMAKE_HOST_APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
            set (STD_CXX_FLAG "${STD_CXX_FLAG} -stdlib=libc++ -Qunused-arguments")
        endif (CMAKE_HOST_APPLE AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))

    endif (CXX11_FOUND)

    if (CXX14_FOUND)
        set (STD_CXX_FLAG "-std=c++14")

    endif (CXX14_FOUND)

endmacro(_seqan_cxx_standard_gcc)


_seqan_cxx_standard_windows()
_seqan_cxx_standard_gcc()

# By default, C++11 compiler support implies the C++11 STL.
set(CXX11_STL_FOUND ${CXX11_FOUND})

if (NOT CXX11_FOUND)
    message (FATAL_ERROR "Seqan requires since v2.1.0 C++11.")
    return ()
endif (NOT CXX11_FOUND)
