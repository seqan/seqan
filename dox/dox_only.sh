#!/bin/bash

# Run doxygen-style documentation system.

rm -f dddoc_cache.bin && ../util/bin/dox.py -b ../core -b ../extras -i ../core/include/seqan  -i ../extras/include/seqan -i pages --image-dir ../docs2/images/ $*
