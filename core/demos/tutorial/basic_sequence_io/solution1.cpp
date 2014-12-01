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

    seqan::SeqFileIn seqFileIn(argv[1]);
    readRecord(id, seq, seqFileIn);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}
