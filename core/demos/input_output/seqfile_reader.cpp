#include <fstream>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/core/demos/input_output/example.fa");

    // Open file
    SeqFileIn inFile(toCString(path));

    // Read file record-wise.
    CharString id;
    Dna5String seq;
    try
    {
        while (!atEnd(inFile))
        {
            readRecord(id, seq, inFile);
            std::cout << id << "\t" << seq << "\n";
        }
    }
    catch (ParseError const & err)
    {
        std::cerr << "Problem reading file: " << err.what() << "\n";
        return 1;
    }

    return 0;
}
