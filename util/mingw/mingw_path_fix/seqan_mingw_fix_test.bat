@echo off

REM ==========================================================================
REM                     SeqAn MinGW Fix Test Batch 
REM ==========================================================================
REM Copyright (c) 2006-2016, Knut Reinert, FU Berlin
REM All rights reserved.
REM
REM Redistribution and use in source and binary forms, with or without
REM modification, are permitted provided that the following conditions are met:
Rem
REM 	* Redistributions of source code must retain the above copyright
REM 	  notice, this list of conditions and the following disclaimer.
REM 	* Redistributions in binary form must reproduce the above copyright
REM 	  notice, this list of conditions and the following disclaimer in the
REM 	  documentation and/or other materials provided with the distribution.
REM 	* Neither the name of Knut Reinert or the FU Berlin nor the names of
REM 	  its contributors may be used to endorse or promote products derived
REM 	  from this software without specific prior written permission.
REM
REM THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
REM AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
REM IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
REM ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
REM FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
REM DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
REM SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
REM CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
REM LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
REM OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
REM DAMAGE.
REM
REM ==========================================================================
REM Author: Rene Rahn <rene.rahn@fu-berlin.de>
REM ==========================================================================
REM This script generates different test scenarios and runs the test mode of
REM the corresponding program seqan_mingw_fix.exe.
REM ==========================================================================

REM Test structure generation tool

REM  create test path environment containing one path
set currentDir=%~dp0
echo %currentDir%

mkdir %currentDir%test

REM ### no eplicite test scenario necessary for 1st test case ###

REM ### create test scenario for 2nd test case: one path
mkdir %currentDir%test\test2
echo "" > %currentDir%test\test2\sh.exe
set PATH_TEST_2=%currentDir%test\test2
REM ###

REM ### create test scenario for 3rd test case: two paths, separated by ;
mkdir %currentDir%test\test3
mkdir %currentDir%test\test3\test3_a
mkdir %currentDir%test\test3\test3_b

echo "" > %currentDir%test\test3\test3_b\sh.exe
set PATH_TEST_3=%currentDir%test\test3\test3_a;%currentDir%test\test3\test3_b
REM ###

REM ### create test scenario for 4th test case: two paths, leading ;
mkdir %currentDir%test\test4
mkdir %currentDir%test\test4\test4_a
mkdir %currentDir%test\test4\test4_b

echo "" > %currentDir%test\test4\test4_a\sh.exe
set PATH_TEST_4=;%currentDir%test\test4\test4_a;%currentDir%test\test4\test4_b
REM ### 

REM ### create test scenario for 5th test case: two paths, trailing ;
mkdir %currentDir%test\test5
mkdir %currentDir%test\test5\test5_a
mkdir %currentDir%test\test5\test5_b

echo "" > %currentDir%test\test5\test5_a\sh.exe
set PATH_TEST_5=%currentDir%test\test5\test5_a;%currentDir%test\test5\test5_b;
REM ###

REM ### create test scenario for 6th test case: two paths, leading ;
mkdir "%currentDir%test\test6"
mkdir "%currentDir%test\test6\test6 a"
mkdir "%currentDir%test\test6\test6_b"

echo "" > "%currentDir%test\test6\test6 a\sh.exe"
set PATH_TEST_6=%currentDir%test\test6\test6 a;%currentDir%test\test6\test6_b
REM ###

REM ### create test scenario for 7th test case: two paths, trailing ;
mkdir %currentDir%test\test7
mkdir %currentDir%test\test7\test7_a
mkdir "%currentDir%test\test7\test7 b"

echo "" > %currentDir%test\test7\test7_a\sh.exe
set PATH_TEST_7=%currentDir%test\test7\test7_a;%currentDir%test\test7\test7 b
REM ###

REM ### create test scenario for 8th test case: path without \
mkdir %currentDir%test\test8
mkdir %currentDir%test\test8\test8_a
mkdir %currentDir%test\test8\test8_b

echo "" > %currentDir%test\test8\test8_b\sh.exe
set PATH_TEST_8=;%currentDir%test\test8\test8_a;PathWithoutBackslash;%currentDir%test\test8\test8_b
REM ### 

REM call program in test mode and wait until it is finished.
START /WAIT %currentDir%seqan_mingw_fix.exe -t %currentDir%

REM remove test structures
rmdir %currentDir%test /S /Q
