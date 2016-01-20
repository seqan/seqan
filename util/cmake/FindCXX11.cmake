
message(AUTHOR_WARNING "Since Seqan 2.1.0, we require at least C++11 "
        "and select automatically the highest available C++ standard. "
        "Thus, finding C++11 via `find(CXX11)` and including the variable "
        "`CXX11_CXX_FLAGS` into `SEQAN_CXX_FLAGS` or `CMAKE_CXX_FLAGS` are not "
        "necessary anymore. Please remove those instructions.")

include(FindStdCXX)
