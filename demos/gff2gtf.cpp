///A simple annotation converter. Convert a Gff to Gtf or vice versa.
#include <fstream>
#include <iostream>
#include <string>

#include <seqan/store.h>

using namespace seqan;

int main(int argc, const char * argv[])
{
    if (argc < 3)
    {
        std::cerr << "[OPTION]... <infile.[gff|gtf]> <outfile.[gff|gtf]>" << std::endl;
        return 0;
    }

    FragmentStore<> store;
    GffFileIn inFile;
    if (!open(inFile, argv[1]))
    {
        std::cerr << "Failed to open annotation infile for reading." << std::endl;
        return 1;
    }
    readRecords(store, inFile);

    GffFileOut outFile;
    if (!open(outFile, argv[2]))
    {
        std::cerr << "Failed to open annotation outfile for writing." << std::endl;
        return 1;
    }
    writeRecords(outFile, store);

    return 0;
}
