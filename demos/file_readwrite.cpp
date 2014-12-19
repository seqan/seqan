#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

char block1[512] = "This a test string";
char block2[512];

int main()
{
///First we create a new / overwrite existing binary file using @Spec.Sync@
    File<Sync<> > myFile1;
    if (!open(myFile1, "file_types.bin", OPEN_WRONLY | OPEN_CREATE))
    {
        std::cout << "Could not open for writing\n";
        return 1;
    }
    write(myFile1, block1, sizeof(block1));
    close(myFile1);

///Then we read the binary file using @Spec.Sync@
    File<Sync<> > myFile2;
    if (!open(myFile2, "file_types.bin", OPEN_RDONLY))
    {
        std::cout << "Could not open for reading\n";
        return 1;
    }
    read(myFile2, block2, sizeof(block2));
    close(myFile2);

    std::cout << block1 << "\n";
    std::cout << block2 << "\n";

    return 0;
}
