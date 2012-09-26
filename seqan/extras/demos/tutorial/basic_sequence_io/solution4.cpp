#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString qual;

    seqan::SequenceStream seqStream(argv[1]);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    while (!atEnd(seqStream))
    {
        if (readRecord(id, seq, qual, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read from example.fa!\n";
            return 1;
        }

        std::cout << id << '\t' << seq << '\t' << qual << '\n';
    }

    return 0;
}
