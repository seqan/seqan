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

    // Create the AutoSeqStreamFormat object and guess the file format.
    seqan::AutoSeqStreamFormat formatTag;
    if (!guessStreamFormat(reader, formatTag))
        return 1;  // Could not detect file format.

    // Read file record-wise.
    seqan::CharString id;
    seqan::String<seqan::Dna5Q> seq;
    while (!atEnd(reader))
    {
        if (readRecord(id, seq, reader, formatTag) != 0)
            return 1;  // Could not record from file.

        std::cout << id << "\t" << seq << "\n";
    }
    
    return 0;
}
