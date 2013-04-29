#!/bin/bash

# Run doxygen-style documentation system.

rm -f dddoc_cache.bin && ../util/bin/dox.py -ldd ../core/ -ldd ../extras -ldd ../docs2/concepts -ldd ../docs2/pages --debug -b ../
