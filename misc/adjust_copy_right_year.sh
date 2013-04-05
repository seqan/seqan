#!/bin/bash

#This file replaces in every .h and .cpp file in include and tests directories of core and trunk an old copyright line with a new one. 
#Note: You have to replace the until date in this document
#Note: Call this file from the trunk directory
for i in $( find core/{include,tests, demos} extras/{include,tests, demos} -name "*.h" -o -name "*.cpp" ); do
    sed -i '' 's/Copyright (c) 20..-20.., Knut Reinert, FU Berlin/Copyright (c) 2006-2013, Knut Reinert, FU Berlin/g' $i
done

for i in $( find core/{include,tests} extras/{include,tests} -name "INFO" ); do
    sed -i '' 's/Copyright: 20..-20.*, FU Berlin/Copyright: 2006-2013, FU Berlin/g' $i
done

#find seqan-trunk/core seqan-trunk/extras  -type f -exec sed -i 's/Copyright (c) 2006-201.*, Knut Reinert, FU Berlin/Copyright (c) 2006-2013, Knut Reinert, FU Berlin/g' {} \

