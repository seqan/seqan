#include <fstream>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main()
{
    seqan::CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/core/demos/input_output/example.fa");

    // Open file
    seqan::SeqFileIn inFile(toCString(path));

    // Read file record-wise.
    seqan::CharString id;
    seqan::Dna5String seq;
    while (!atEnd(inFile))
    {
        readRecord(id, seq, inFile);
        std::cout << id << "\t" << seq << "\n";
    }

    return 0;
}
