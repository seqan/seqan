#!/bin/sh
rm -f html/*
./main.py ../core/include ../extras/include -d concepts -d pages $@
exit $?
