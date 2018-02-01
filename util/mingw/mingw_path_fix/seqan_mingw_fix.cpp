// ==========================================================================
//                         SeqAn MinGW Fix
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//		* Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//		* Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//		* Neither the name of Knut Reinert or the FU Berlin nor the names of
//        its contributors may be used to endorse or promote products derived
//        from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// This is a cpp program used to exclude all directories containing the file 
// "sh.exe" within the PATH environment of win32 systems.
// ==========================================================================

#include <iostream>
#include <string>
#include <vector>

#include <stdlib.h>
#include <windows.h>

using namespace std;

/**
.Function.Internal.extendDirName
..cat:Miscellaneous
..summary:Appends the given file name to the directory while checking for consistent path separators.
..signature:extendDirName(dirName, fileName)
..param.dirName:The directory path.
...type:nolink:string &
..param.fileName:The name of the file to be extended to the path.
...type:nolink:string const &
..return:0 on success, 1 if dirName is invalid.
*/
int extendDirName(string & dirName, string const & fileName)
{
    string const separator = "\\";  // windows style path separator

    if (dirName.length() == 0)
    {
        return 1;
    }
    
    size_t lastOcc = dirName.find_last_of(separator);  // get last occurrence of separator
    // if last char is not a separator, then add it
	if (lastOcc < dirName.length() || lastOcc == dirName.npos)
    {
        dirName.append(separator);
    }
    dirName.append(fileName);  // add file name
    return 0;
}

/**
.Function.Internal.tokenize
..cat:Miscellaneous
..summary:Splits given string into tokens based on the given delimiters.
..signature:tokenize(tokens, str, delimiters)
..param.tokens:A vector storing the tokens.
...type:nolink:vector<string> &
..param.str:The string to be tokenized.
...type:nolink:string const &
..param.delimiters:A list of delimiters.
...type:nolink:string const &
*/
void tokenize(vector<string> & tokens,
            string const & str,
            string const & delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

/**
.Function.Internal.findFile
..cat:Miscellaneous
..summary:Uses WIN32 API to search for the given file with full qualified path.
..signature:findFile(fullFileName)
..param.fullFileName:The absolut path of the file.
...type:nolink:char const *
..return:1 if file not found, otherwise 0.
*/
int findFile(string const & fullFileName)
{
    if (fullFileName.empty())
    {
        return 1;
    }
    WIN32_FIND_DATA findFileData;  // file object
    HANDLE hFind = FindFirstFile(fullFileName.c_str(), &findFileData);  // sore result in HANDLER
    if (hFind == INVALID_HANDLE_VALUE)  // file not found
    {
        return 2;
    }
    return 0;
}

/**
.Function.Internal.updatePathEnv
..cat:Miscellaneous
..summary:Parses the directories of the current PATH environment and removes 
all paths that contain the file "sh.exe".
..signature:updatePathEnv(newPath, currPath)
..param.newPath:A string storing all directories that does not contain sh.exe.
...type:nolink:string &
..param.currPath:The current value of the PATH environment.
...type:nolink:string const &
..return: 0 on success, 1 on failures.
*/
int updatePathEnv(string & newPath, string const & currPath)
{
    if (currPath.empty())
    {	
        return 1;
    }
    
    string const delimiters = ";";  // delimiter for PATH environment (Win32 only)
    string const excludeFile = "sh.exe";  // the name of the file which should be excluded
        
    // this part parses dirs from PATH environment and checks for existing sh.exe

    vector<string> tokens;  // stores each directory contained within the PATH
    tokenize(tokens, currPath, delimiters);
    
    vector<string>::iterator it = tokens.begin();
    
    while (it != tokens.end())
    {
        // append file name to path
		string tmp(*it);
        if (extendDirName(tmp, excludeFile) != 0)
        {
            ++it;
			continue;
        }
        
        if (findFile(tmp) == 2)  // if file is not found, add directory to buffer
        {
            newPath.append(*it);
            newPath.append(delimiters);
        }
        ++it;  // next directory
    }
    return 0;
}

// executes the program in normal mode.
int runNormal()
{
	// run program in normal mode
    string newPath;
	
    string currPath(getenv("PATH"));
    int errorCode = updatePathEnv(newPath, currPath);
    
    if (errorCode != 0)
    {
        cout << currPath << endl;
		return errorCode;
    }
    cout << newPath << endl;  //outputs the modified PATH value
	return errorCode;
}

// Prints error message if testResult is false 
template <typename T>
int printAssert(bool testResult, T expectedValue, T isValue)
{
		if (!testResult)
		{
			cerr << "\nError: Expected value " << expectedValue << " but is " << isValue << "." << endl; 
			return 1;
		}
		return 0;
}

//executes the tests for this program
int runTests(string const & workingDir)
{
	{ // Test 1: empty string
		string test1;
		string newPath;

		int errorCode = updatePathEnv(newPath, test1);
		cerr << "Test 1: empty string  ";
		if (printAssert(errorCode == 1, 1, errorCode) != 0)
			return 1;
		if (printAssert(newPath.empty(), string(""), newPath) != 0)
			return 1;
		cerr << "...Ok!" << endl;
	}
	
	{ // Test 2: one path
		string test2(getenv("PATH_TEST_2"));
		string newPath;
	
		int errorCode = updatePathEnv(newPath, test2);
		cerr << "Test 2: one path  ";
		if (printAssert(errorCode == 0, 0, errorCode) != 0)
			return 1;
		if (printAssert(newPath.empty(), string(""), newPath) != 0)
			return 1;
		cerr << "...Ok!" << endl;
	}
	
	{ // Test 3: two paths, separated by ;
		string test3(getenv("PATH_TEST_3"));
		string newPath;
		
		int errorCode = updatePathEnv(newPath, test3);
		cerr << "Test 3: two paths separated by ;  ";
		if (printAssert(errorCode == 0, 0, errorCode) != 0)
			return 1;
		string testFolder(workingDir);
		testFolder.append("test\\test3\\test3_a;");
		if (printAssert(0 == newPath.compare(testFolder), testFolder, newPath) != 0)
			return 1;
		cerr << "... Ok!" << endl;
	}
	
	{ // Test 4: two paths, leading ;
		string test4(getenv("PATH_TEST_4"));
		string newPath;
		
		int errorCode = updatePathEnv(newPath, test4);
		cerr << "Test 4: two paths, leading ;  ";
		if (printAssert(errorCode == 0, 0, errorCode) != 0)
			return 1;
		string testFolder(workingDir);
		testFolder.append("test\\test4\\test4_b;");
		if (printAssert(0 == newPath.compare(testFolder), testFolder, newPath) != 0)
			return 1;
		cerr << "...Ok!" << endl;
	}
	
	{ // Test 5: two paths, trailing ;
		string test5(getenv("PATH_TEST_5"));
		string newPath;
		
		int errorCode = updatePathEnv(newPath, test5);
		cerr << "Test 5: two paths, trailing ;  ";
		if (printAssert(errorCode == 0, 0, errorCode) != 0)
			return 1;
		string testFolder(workingDir);
		testFolder.append("test\\test5\\test5_b;");
		if (printAssert(0 == newPath.compare(testFolder), testFolder, newPath) != 0)
			return 1;
		cerr << "...Ok!" << endl;
	}

	{ // Test 6: two paths, first one with space
		string test6(getenv("PATH_TEST_6"));
		string newPath;
		
		int errorCode = updatePathEnv(newPath, test6);
		cerr << "Test 6: two paths, first one with space  ";
		if (printAssert(errorCode == 0, 0, errorCode) != 0)
			return 1;
		string testFolder(workingDir);
		testFolder.append("test\\test6\\test6_b;");
		if (printAssert(0 == newPath.compare(testFolder), testFolder, newPath) != 0)
			return 1;
		cerr << "...Ok!" << endl;
	}
	
	{ // Test 7: two paths, second one with space
		string test7(getenv("PATH_TEST_7"));
		string newPath;
		
		int errorCode = updatePathEnv(newPath, test7);
		cerr << "Test 7: two paths, second one with space  ";
		if (printAssert(errorCode == 0, 0, errorCode) != 0)
			return 1;
		string testFolder(workingDir);
		testFolder.append("test\\test7\\test7 b;");
		if (printAssert(0 == newPath.compare(testFolder), testFolder, newPath) != 0)
			return 1;
		cerr << "...Ok!" << endl;
	}
	
	{ // Test 8: two paths, second one with space
		string test8(getenv("PATH_TEST_8"));
		string newPath;
		
		int errorCode = updatePathEnv(newPath, test8);
		cerr << "Test 8: three paths, second one without \\  ";
		if (printAssert(errorCode == 0, 0, errorCode) != 0)
			return 1;
		string testFolder(workingDir);
		testFolder.append("test\\test8\\test8_a;PathWithoutBackslash;");
		if (printAssert(0 == newPath.compare(testFolder), testFolder, newPath) != 0)
			return 1;
		cerr << "...Ok!" << endl;
	}
	return 0;
}

// prints the usage of the program
void printHelp()
{
	cout <<   "\t\t*************************************";
	cout << "\n\t\t*         SeqAn MinGW Fix           *";
	cout << "\n\t\t*************************************";
	cout << "\n\n";
	cout << "SYNOPSIS:";
	cout << "\n\tseqan_mingw_setup.exe";
	cout << "\n\tseqan_mingw_setup.exe -t [path]";
	cout << "\n\tseqan_mingw_setup.exe -h";
	cout << "\n\nDESCRIPTION:";
	cout << "\n\tSeqAn MinGW setup tool is a preinstallation program that removes all";
	cout << "\n\tdirectories from the PATH environment of the current session that";
	cout << "\n\tcontain the program ""sh.exe"" and prints the modified PATH to the";
	cout << "\n\tstandard ouptut. MinGW builds cannot be generated with CMake as long";
	cout << "\n\tas there lives a shell programm within the PATH environment.";
	cout << "\n\tNote, if you want to change the PATH environment immediately you should";
	cout << "\n\tuse the corresponding batch script which calls this program and";
	cout << "\n\tevaluates its output to set the new PATH environment.";
	cout << "\n\nCOMMAND-LINE OPTIONS:";
	cout << "\n\t-t [path]\truns the program in test mode. The path is the current";
	cout << "\n\t\t\tworking directory where the test files are located.\n";
	cout << "\n\t-h\tprints the help of the program.";
}

int main(int argc, char * argv[])
{   
	if (argc == 1)
	{
		return runNormal();
		return 0;
	}
	else if (argc == 3 && argv[1][1] == 't')
	{
		string workingDir(argv[2]);
		int errorCode = runTests(workingDir);
		if (errorCode == 1)
		{
			cerr << "Test Failed!" << endl;
			system("Pause");
			return errorCode;
		}
		cerr << "Test succeded!" << endl;
		system("Pause");
		return 0;
	}
	printHelp();
	system("Pause");
	return 2;
}
