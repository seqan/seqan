#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    CharString path = SEQAN_TEMP_FILENAME();
    append(path, ".fa");

    // Open file and check for errors.
    SequenceStream seqStream(toCString(path), SequenceStream::WRITE);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open " << path << " for writing.\n";
        return 1;
    }

    // Write two sequences to the file.
    if (writeRecord(seqStream, "one", "CGAT") != 0 ||
        writeRecord(seqStream, "two", "ASDF") != 0)
    {
        std::cerr << "ERROR: Problem writing to " << path << "\n";
        return 1;
    }

    return 0;
}
