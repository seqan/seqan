#!/bin/bash

# Run doxygen-style documentation system.

rm -f dddoc_cache.bin && ../util/bin/dox.py -b .. -i ../include/seqan -i pages --image-dir ./images/docs2 --image-dir ./images/reduced_aminoacid/ $*
