@echo off

REM ==========================================================================
REM                     SeqAn MinGW Fix Batch
REM ==========================================================================
REM Copyright (c) 2006-2016, Knut Reinert, FU Berlin
REM All rights reserved.
REM
REM Redistribution and use in source and binary forms, with or without
REM modification, are permitted provided that the following conditions are met:
REM
REM		* Redistributions of source code must retain the above copyright
REM		  notice, this list of conditions and the following disclaimer.
REM		* Redistributions in binary form must reproduce the above copyright
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
REM This tool parses the directories of the PATH environment (using the 
REM external program seqan_mingw_fix.exe) and removes all directories that 
REM contains "sh.exe" as one of their childs. The modified PATH is then set 
REM to the current environment.
REM ==========================================================================


echo Disabling sh.exe path...

REM start actual program and redirect stdout to tmp file
REM start cmd /c %program% ^> %tmp_file%

REM ouput old path environment - in case of unexpected corruptions of the PATH variable
echo.
echo old path environment:
echo %path%

REM call seqan_mingw_fix.exe in current working directory: %~dp0 and set path to output of program
for /F "usebackq tokens=*" %%i in (`%~dp0seqan_mingw_fix.exe`) do set PATH=%%i

REM print new PATH
echo.
echo new path environment:
echo %path%

REM done!
echo.
echo done!
