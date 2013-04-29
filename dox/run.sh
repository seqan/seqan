#!/usr/bin/bash

# Run doxygen-style documentation system.

rm -f dddoc_cache.bin && ../util/bin/dox.py -ldd ../core/ -ldd ../extras -ldd ../docs/concepts -ldd ../docs/pages --debug -b ../
