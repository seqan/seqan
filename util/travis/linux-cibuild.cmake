# This project name is used for the CDash submission.
SET (CTEST_PROJECT_NAME "SeqAn")

# define build name&co for easier identification on CDash
set(CTEST_BUILD_NAME "travis-$ENV{TRAVIS_BUILD_NUMBER}-$ENV{TRAVIS_REPO_SLUG}-$ENV{TRAVIS_BRANCH}-$ENV{BUILD_NAME}-$ENV{CXX}")
set(CTEST_SITE "travis-ci-build-server")
set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE_DIRECTORY}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/_build")

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

# create cache - mind the newline between variables!
set(INITIAL_CACHE
"CMAKE_BUILD_TYPE=Release
SEQAN_TRAVIS_BUILD:BOOL=ON")
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" ${INITIAL_CACHE})

# customize reporting of errors in CDash
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 1000)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000)

# try to speed up the builds so we don't get killed
if ("$ENV{CXX}" MATCHES ".*(clang\\+\\+-?.*)")
  set(CTEST_BUILD_FLAGS -j6)
else ("$ENV{CXX}" MATCHES ".*(clang\\+\\+-?.*)")
  set(CTEST_BUILD_FLAGS -j4)
endif ("$ENV{CXX}" MATCHES ".*(clang\\+\\+-?.*)")

# we want makefiles
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# Copying the CTestConfig.cmake here is not optimal.  You might have to call
# ctest twice to get an actual build since ctest expects it to be present
# at the first time and will fail.
CONFIGURE_FILE (${CTEST_SOURCE_DIRECTORY}/util/cmake/CTestConfig.cmake
                ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake
                COPYONLY)

# run the classical ctest suite without update
# travis-ci handles this for us
ctest_start     (Continuous)

ctest_configure (BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE _configure_ret)
ctest_build     (BUILD "${CTEST_BINARY_DIRECTORY}" NUMBER_ERRORS _build_errors)
ctest_test      (BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL 8)
ctest_submit()

# indicate errors
if(${_build_errors} GREATER 0 OR NOT ${_configure_ret} EQUAL 0)
  file(WRITE "$ENV{SOURCE_DIRECTORY}/failed" "build_failed")
endif()
