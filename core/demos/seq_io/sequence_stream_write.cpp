#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

// USAGE: sequence_read_stream_write FILE
//
// Print some sequences to the file FILE

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " SEQUENCE.{fa,fq}\n";
        return 1;
    }

    SequenceStream seqStream(argv[1], SequenceStream::WRITE);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for writing.\n";
        return 1;
    }

    if (writeRecord(seqStream, "one", "CGAT") != 0 ||
        writeRecord(seqStream, "two", "ASDF") != 0)
    {
        std::cerr << "ERROR: Problem writing to " << argv[1] << "\n";
        return 1;
    }

    return 0;
}
