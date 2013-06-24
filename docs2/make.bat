@echo off
python ../util/bin/dddoc.py ..\core\include ..\extra\include -d concepts -d pages -I ..\core -I ..\extras %1 %2 %3 %4 %5 %6
