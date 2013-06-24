#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

// USAGE: sequence_read_stream_read FILE
//
// Print the contents of sequence FILE to stdout in tabular format.

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " SEQUENCE.{fa,fq}\n";
        return 1;
    }

    SequenceStream seqStream(argv[1]);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    seqan::CharString id, seq;
    while (!atEnd(seqStream))
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "Problem reading from " << argv[1] << "\n";
            return 1;
        }
        std::cout << id << "\t" << seq << "\n";
    }

    return 0;
}
