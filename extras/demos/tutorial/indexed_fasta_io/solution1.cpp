#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: build_fai FILE.fa\n";
        return 1;
    }

    seqan::FaiIndex faiIndex;
    if (build(faiIndex, argv[1]) != 0)
    {
        std::cerr << "ERROR: Could not build FAI index for file " << argv[1] << ".\n";
        return 1;
    }

    seqan::CharString faiFilename = argv[1];
    append(faiFilename, ".fai");

    if (write(faiIndex, toCString(faiFilename)) != 0)
    {
        std::cerr << "ERROR: Could not write the index to file!\n";
        return 1;
    }

    std::cout << "Index file " << faiFilename << " was successfully created.\n";
    return 0;
}
