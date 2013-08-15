#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/core/demos/seq_io/example.fa");

    // Open file and check for errors.
    SequenceStream seqStream(toCString(path));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open " << path << " for reading.\n";
        return 1;
    }

    // Read from file and print the result to stdout.
    seqan::CharString id, seq;
    while (!atEnd(seqStream))
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "Problem reading from " << path << "\n";
            return 1;
        }
        std::cout << id << "\t" << seq << "\n";
    }

    return 0;
}
