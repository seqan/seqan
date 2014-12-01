#!/bin/bash

# Run doxygen-style documentation system.

rm -f dddoc_cache.bin && ../util/bin/dox.py -b ../core -b ../extras -i ../include/seqan  -i ../include/seqan -i pages --image-dir ../docs2/images/ --image-dir ./images/reduced_aminoacid/ $*
