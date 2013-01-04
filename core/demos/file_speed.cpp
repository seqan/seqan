#define SEQAN_PROFILE
//#define SEQAN_DEBUG

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <iostream>

using namespace seqan;
using namespace std;

const unsigned long blockSize = 1 << 12;
const unsigned long repeats = 1 << 19;

CharString block1 = "This a test string";
CharString block2;

template <typename TFile>
void testThroughput(const char *fileName)
{
	TFile myFile;
	typename AsyncRequest<TFile>::Type req1, req2;

	if (!open(myFile, fileName, OPEN_WRONLY | OPEN_CREATE)) {
		cout << "Could not open for writing" << endl;
		return;
	}

	SEQAN_PROTIMESTART(iotime);

	asyncWriteAt(myFile, toCString(block1), blockSize, 0 * blockSize, req1);
	asyncWriteAt(myFile, toCString(block2), blockSize, 1 * blockSize, req2);
	for (unsigned i = 1; i < repeats; ++i)
	{
		waitFor(req1);
		asyncWriteAt(myFile, toCString(block1), blockSize,    2*i  * blockSize, req1);
		waitFor(req2);
		asyncWriteAt(myFile, toCString(block2), blockSize, (2*i+1) * blockSize, req2);
	}
	waitFor(req1);
	waitFor(req2);

	cout << ((repeats*blockSize / (512.0 * 1024.0)) / SEQAN_PROTIMEDIFF(iotime));
	cout << " MB/s" << endl;

	close(myFile);
}

template <typename TFile>
void testExtString(const char *fileName)
{
	String<char, External<ExternalConfig<TFile> > > myString;

	if (!open(myString, fileName, OPEN_WRONLY | OPEN_CREATE)) {
		cout << "Could not open for writing" << endl;
		return;
	}

	SEQAN_PROTIMESTART(iotime);

	for (unsigned i = 0; i < repeats; ++i)
	{
		append(myString, block1);
		append(myString, block2);
	}

	cout << ((repeats*blockSize / (512.0 * 1024.0)) / SEQAN_PROTIMEDIFF(iotime));
	cout << " MB/s" << endl;
}

template <typename TFile>
void testMMapString(const char *fileName)
{
	String<char, MMap<> > myString;

	if (!open(myString, fileName, OPEN_RDWR /*| OPEN_CREATE*/)) {
		cout << "Could not open for writing" << endl;
		return;
	}

	SEQAN_PROTIMESTART(iotime);

	for (unsigned i = 0; i < repeats; ++i)
	{
		append(myString, block1);
		append(myString, block2);
	}

	cout << ((repeats*blockSize / (512.0 * 1024.0)) / SEQAN_PROTIMEDIFF(iotime));
	cout << " MB/s" << endl;
}

int main() 
{
	resize(block1, blockSize);
	resize(block2, blockSize);
	cout << "asyncWrite() using FILE*        ";		testThroughput< FILE* >				("file_speed1.bin");
	cout << "asyncWrite() using sync. File   ";		testThroughput< File< Sync<> > >	("file_speed2.bin");
	cout << "asyncWrite() using async. File  ";		testThroughput< File< Async<> > >	("file_speed3.bin");
	cout << "ExtString using FILE*       ";		testExtString< FILE* >				("file_speed4.bin");
	cout << "ExtString using sync. File  ";		testExtString< File< Sync<> > >		("file_speed5.bin");
	cout << "ExtString using async. File ";		testExtString< File< Async<> > >	("file_speed6.bin");
	cout << "Memory Mapped String        ";		testMMapString< File< Async<> > >	("file_speed7.bin");
	return 0;
}
