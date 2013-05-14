#!/bin/bash

# Run doxygen-style documentation system.

rm -f dddoc_cache.bin && ../util/bin/dox.py -b ../ -i ../core/include/seqan
