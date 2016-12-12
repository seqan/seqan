#!/bin/bash
# Run doxygen-style documentation system.
rm -f dddoc_cache.bin && python2 ../util/bin/dox.py -ldd ../ -ldd ../extras -ldd ../docs2/concepts -ldd ../docs2/pages --debug -b ../
