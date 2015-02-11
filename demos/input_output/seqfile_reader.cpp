#include <fstream>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/demos/input_output/example.fa");

    // Open file
    SeqFileIn inFile(toCString(path));

    // Read file record-wise.
    CharString id;
    Dna5String seq;
    while (!atEnd(inFile))
    {
        readRecord(id, seq, inFile);
        std::cout << id << "\t" << seq << "\n";
    }

    return 0;
}
