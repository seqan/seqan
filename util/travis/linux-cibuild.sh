#!/bin/bash

# get some infos from git to embed it in the build name
export SOURCE_DIRECTORY=`pwd`
mkdir -p _build

# define the build name
if [ "${TRAVIS_PULL_REQUEST}" != "false" ]; then
  export BUILD_NAME=${TRAVIS_PULL_REQUEST}
elif [[ -n "${TRAVIS_COMMIT_RANGE}" ]]; then
  export BUILD_NAME=${TRAVIS_COMMIT_RANGE}
else
  export BUILD_NAME=${TRAVIS_COMMIT}
fi

# disable OpenMP warnings for clang
if [ "$(echo ${CXX} | cut -c1-5)" = "clang" ]; then
  export CXXFLAGS="${CXXFLAGS} -Qunused-arguments -DSEQAN_IGNORE_MISSING_OPENMP=1"
fi

# compile with c++-17 if g++-7 is used.
if [ "$(echo ${CXX})" = "g++-7" ]; then
  export CXXFLAGS="${CXXFLAGS} -std=c++17"
fi

# Switch version check default to OFF to prevent checks during app tests
CXXFLAGS="${CXXFLAGS} -DSEQAN_VERSION_CHECK_OPT_IN=YES"

ctest -V -S util/travis/linux-cibuild.cmake

# we indicate build failures if ctest experienced any errors
if [ -f ${SOURCE_DIRECTORY}/failed ]; then
  exit -1
fi

FAILED_TEST=$(find _build -name "Test.xml" -type f | xargs grep "<Test Status=\"failed\">" -c)

if [ "${FAILED_TEST}" -gt "0" ]; then
    exit -1
fi

# it seems like everything worked
exit 0
