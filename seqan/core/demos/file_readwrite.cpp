#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;

char block1[512] = "This a test string";
char block2[512];

int main()
{
///First we create a new / overwrite existing binary file using @Spec.Sync@
	File<Sync<> > myFile1;
	if (!open(myFile1, "file_types.bin", OPEN_WRONLY | OPEN_CREATE)) {
		cout << "Could not open for writing" << endl;
		return 1;
	}
	write(myFile1, block1, sizeof(block1));
	close(myFile1);

///Then we read the binary file using @Spec.Sync@
	File<Sync<> > myFile2;
	if (!open(myFile2, "file_types.bin", OPEN_RDONLY)) {
		cout << "Could not open for reading" << endl;
		return 1;
	}
	read(myFile2, block2, sizeof(block2));
	close(myFile2);

	cout << block1 << endl;
	cout << block2 << endl;

	return 0;
}

