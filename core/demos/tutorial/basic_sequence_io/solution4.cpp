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

    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    try
    {
        while (!atEnd(seqFileIn))
        {
            readRecord(id, seq, qual, seqFileIn);
            std::cout << id << '\t' << seq << '\t' << qual << '\n';
        }
    }
    catch (std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
