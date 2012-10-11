// File:   sequence_length.cpp
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    // Check arguments.
    if (argc != 2) {
        std::cerr << "Wrong argument count!" << std::endl
                  << "USAGE: sequence_length SEQUENCE.fasta" << std::endl;
        return 1;
    }

    // Open file.
    std::fstream inF(argv[1], std::ios::in | std::ios::binary);
    if (!inF.good())
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Create RecordReader.
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(inF);

    // Read sequence file and print sequence lengths.
    size_t total = 0;
    seqan::CharString id;
    seqan::Dna5String seq;
    while (!atEnd(reader))
    {
        if (readRecord(id, seq, reader, seqan::Fasta()) != 0)
        {
            std::cerr << "ERROR: Problem reading file " << argv[1] << ".\n";
            return 1;
        }

        std::cout << id << "\t" << seq << "\n";
        total += length(seq);
    }
    std::cout << "sum\t" << total << std::endl;

    return 0;
}
