#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;  // Invalid number of arguments.

    // Open file
    seqan::SeqFileIn inFile(argv[1]);

    // Read file record-wise.
    seqan::CharString id;
    seqan::CharString seq;
    while (!atEnd(inFile))
    {
        readRecord(id, seq, inFile);
        std::cout << id << "\t" << seq << "\n";
    }

    return 0;
}
