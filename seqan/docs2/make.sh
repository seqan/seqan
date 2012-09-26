#!/bin/sh
rm -f dddoc_cache.bin
rm -rf html/*
../util/bin/dddoc.py -d concepts -d pages ../core/include ../extras/include $@
exit $?
