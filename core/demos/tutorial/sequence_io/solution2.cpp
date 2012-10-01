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

    // Read file in one pass.
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> seqs;
    seqan::StringSet<seqan::CharString> quals;
    if (read2(ids, seqs, quals, reader, seqan::Fastq()) != 0)
        return 1;  // Could not read file.

    for (unsigned i = 0; i < length(ids); ++i)
        std::cout << ids[i] << '\t' << seqs[i] << '\t' << quals[i] << '\n';
    
    return 0;
}
