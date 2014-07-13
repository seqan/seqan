****************************************
*                README                *
****************************************

Date: 05.01.2012
Description: Readme file for the usage of the mingw patch.
Version: 1.0

TOC
1. Software information
2. System requirements
3. Installation notes
4. Usage descripton
5. New features
6. Bug fixes
7. Known issues

######
1. SOFTWARE INFORMATION

Software Name: SeqAn MinGW PATH Fix 
Version: 1.1
Binary Files: seqan_mingw_fix.exe
Source Files: seqan_mingw_fix.cpp
Batch Files: seqan_mingw_fix.bat, seqan_mingw_fix_test.bat

Author: Rene Rahn <rene.rahn[at]fu-berlin.de>

######
2. SYSTEM REQUIREMENTS

Minimum supported client: Win XP
Minimum supported server: Windows Server 2003
 
######
3. INSTALLATION NOTES

There is a prebuild binary file, as well as batch scripts, which can be used without further installing.

######
4. USAGE DESCRIPTION

There is an error while trying to genereate MinGW makefiles within cmake, if there exists a sh.exe somewhere
within the PATH environment of Win32 Systems. In order to allow MinGW builds one has to disable all directories
within the PARH environment that include a sh.exe. This often happens if there is a parallel cygwin installation
or the msys environment is within the PATH. Use the seqan_mingw_fix.bat to automatically scan the directories 
and to reset the PATH such that all sh.exe containing directories are disabled. Call the batch script 
in the console before you execute the cmake command. The batch scriptt automatically calls the seqan_mingw_fix.exe
program which reads the PATH and prints the modified PATH to the console. The batch script reads the output 
and sets the new PATH environment, i.e., the PATH is valid for the current session.
Note, that the change of the PATH environment lasts only as long as the current console session.
In order to test the program there is also a script called seqan_mingw_fix_test.bat. This batch script
creates different test scenarios and calls the program in its test mode. All tests schould pass successfully.

######
5. NEW FEATURES
	* added test functionality
	* program accepts different inputs (see help of the program with -h)
	
######
6. BUG FIXES

######
7. KNOWN ISSUES

