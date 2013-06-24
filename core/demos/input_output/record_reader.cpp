#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;  // Invalid number of arguments.

    // Open file and create RecordReader.
    std::fstream in(argv[1], std::ios::binary | std::ios::in);
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);

    // Read file record-wise.
    seqan::CharString id;
    seqan::Dna5String seq;
    while (!atEnd(reader))
    {
        if (readRecord(id, seq, reader, seqan::Fasta()) != 0)
            return 1;  // Could not read record from file.

        std::cout << id << "\t" << seq << "\n";
    }
    
    return 0;
}
