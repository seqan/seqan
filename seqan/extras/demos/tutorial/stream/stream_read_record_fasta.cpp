// FRAGMENT(header)
#include <fstream>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/stream.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "ERROR: Invalid argument count." << std::endl
                  << "USAGE: " << argv[0] << " FILE" << std::endl;
        return 1;
    }
    
// FRAGMENT(read-sequences)
    // Open file and create RecordReader.
    std::ifstream fasta(argv[1], std::ios_base::in | std::ios_base::binary);
    if (!fasta.good())
        return 1;
    RecordReader<std::ifstream, SinglePass<> > reader(fasta);
    
    // Define variables for storing the sequences and sequence ids.
    CharString id;
    Dna5String seq;
    // Read FASTA file and output "$id\t$seq".
    while (!atEnd(reader))
    {
        if (readRecord(id, seq, reader, Fasta()) != 0)
        {
            std::cerr << "ERROR reading FASTA." << std::endl;
            return 1;
        }
        std::cout << id << "\t" << seq << "\n";
    }

// FRAGMENT(footer)    
    return 0;
}
