#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

int main()
{
    seqan::CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/core/demos/input_output/example.fa");

    // Open file and create RecordReader.
    std::fstream in(toCString(path), std::ios::binary | std::ios::in);
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
