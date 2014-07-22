#include <fstream>
#include <iostream>

#include <seqan/file.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;  // Invalid number of arguments.

    typedef seqan::String<char, seqan::MMap<> > TString;

    // Open memory mapped string.
    TString mmapString;
    if (!open(mmapString, argv[1], seqan::OPEN_RDONLY))
        return 1;  // Could not open file.

    seqan::DirectionIterator<TString, seqan::Input>::Type iter = begin(mmapString);

    // Read file in one pass.
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> seqs;
    seqan::StringSet<seqan::CharString> quals;

    read(ids, seqs, quals, iter, seqan::Fastq());

    for (unsigned i = 0; i < length(ids); ++i)
        std::cout << ids[i] << '\t' << seqs[i] << '\t' << quals[i] << '\n';

    return 0;
}
