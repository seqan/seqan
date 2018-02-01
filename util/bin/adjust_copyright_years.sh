#!/bin/bash

# This file replaces in every .h and .cpp file in include, demos and tests
# directories of core and trunk an old copyright line with a new one.
#
# Note: You have to replace the until date in this document.
# Note: Call this file from the seqan root directory.

for i in $( find {include,tests,demos,apps,util,manual} -name "*.h" -o -name "*.cpp" -o -name "LICENSE" -o -name "*.cu" -o -name "*.cmake" -o -name "*.rst" -o -name "COPYRIGHT" ); do
    sed -i.bak -e 's/Copyright (c) 20\(..\)-20..\(.*\)/Copyright (c) 20\1-2018\2/g' $i
done

for i in $( find {include,tests,apps} -name "INFO" ); do
    sed -i.bak -e 's/Copyright: 20\(..\)-20..\(.*\)/Copyright: 20\1-2018\2/g' $i
done
