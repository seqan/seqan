#!/bin/bash

# This file replaces in every .h and .cpp file in include, demos and tests
# directories of core and trunk an old copyright line with a new one.
#
# Note: You have to replace the until date in this document.
# Note: Call this file from the seqan root directory.

for i in $( find {include,tests,demos,apps,util,manual} -name "*.h" -o -name "*.cpp" -o -name "LICENSE" -o -name "*.cu" -o -name "*.cmake" -o -name "*.rst" -o -name "COPYRIGHT" ); do
    sed -i 's/Copyright (c) 20..-20..\(.*\)Knut Reinert, FU Berlin/Copyright (c) 2006-2016\1Knut Reinert, FU Berlin/g' $i
done

for i in $( find {include,tests,apps} -name "INFO" ); do
    sed -i 's/Copyright: 20..-20.*, FU Berlin/Copyright: 2006-2016, FU Berlin/g' $i
done

