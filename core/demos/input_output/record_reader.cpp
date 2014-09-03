#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

// TODO(singer): rename, delete this file (after the new stream module
// is merged into the develop/master branch!
int main()
{
    seqan::CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/core/demos/input_output/example.fa");

    // Open file
    seqan::SequenceFile<seqan::Input> inFile(toCString(path));

    // Read file record-wise.
    seqan::CharString id;
    seqan::Dna5String seq;
    while (!atEnd(inFile))
    {
        read(inFile, id, seq);

        std::cout << id << "\t" << seq << "\n";
    }

    return 0;
}
