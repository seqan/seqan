// File:   sequence_length.cpp
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    // Check arguments.
    if (argc != 2)
    {
        std::cerr << "Wrong argument count!" << std::endl
                  << "USAGE: sequence_length SEQUENCE.fasta" << std::endl;
        return 1;
    }

    // Open file.
    SeqFileIn file(argv[1]);

    // Read sequence file and print sequence lengths.
    size_t total = 0;
    CharString id;
    CharString seq;
    while (!atEnd(file))
    {
        readRecord(id, seq, file);
        std::cout << id << "\t" << length(seq) << "\n";
        total += length(seq);
    }
    std::cout << "sum\t" << total << std::endl;

    return 0;
}
