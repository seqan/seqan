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
    seqan::Dna5String seq;
    seqan::CharString qual;
    while (!atEnd(inFile))
    {
        readRecord(id, seq, qual, inFile);
        std::cout << id << "\t" << seq << "\t" << qual << "\n";
    }

    return 0;
}
