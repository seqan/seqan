#include <fstream>
#include <iostream>

#include <seqan/file.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;  // Invalid number of arguments.

    // Open memory mapped string.
    seqan::String<char, seqan::MMap<> > mmapString;
    if (!open(mmapString, argv[1], seqan::OPEN_RDONLY))
        return 1;  // Could not open file.

    // Create RecordReader.
    seqan::RecordReader<seqan::String<char, seqan::MMap<> >,
                        seqan::DoublePass<seqan::StringReader> > reader(mmapString);

    // Create the AutoSeqStreamFormat object and guess the file format.
    seqan::AutoSeqStreamFormat formatTag;
    if (!guessStreamFormat(reader, formatTag))
        return 1;  // Could not detect file format.

    // Read file in one pass.
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> seqs;
    seqan::StringSet<seqan::CharString> quals;
    if (read2(ids, seqs, quals, reader, formatTag) != 0)
        return 1;  // Could not read file.

    for (unsigned i = 0; i < length(ids); ++i)
        std::cout << ids[i] << '\t' << seqs[i] << '\t' << quals[i] << '\n';
    
    return 0;
}
