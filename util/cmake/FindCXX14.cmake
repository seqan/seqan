if(__FIND_CXX14_CMAKE__)
    return()
endif()
set(__FIND_CXX14_CMAKE__ TRUE)

# make sure that c++11 detection happens first, so that the c++11-flag
# doesn't overwrite the c++14 flag
find_package(CXX11)

# Visual Studio is still far away from C++14 (maybe add 2015 or 2015+1 later)
# so no checks for this now

include(CheckCXXCompilerFlag)
enable_language(CXX)

check_cxx_compiler_flag("-std=c++14" CXX14_FOUND)
if (CXX14_FOUND)
    set (CXX11_CXX_FLAGS "-std=c++14")
else (CXX14_FOUND)
    check_cxx_compiler_flag("-std=c++1y" CXX14_FOUND)
    if (CXX14_FOUND)
        set (CXX11_CXX_FLAGS "-std=c++1y")
    endif (CXX14_FOUND)
endif (CXX14_FOUND)

# c++14 implies c++11
if (CXX14_FOUND)
    set (CXX11_FOUND TRUE)
endif (CXX14_FOUND)
