#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;  // Invalid number of arguments.

    // Open file
    seqan::SeqFileIn inFile(argv[1]);

    // Read file in one pass.
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> seqs;
    seqan::StringSet<seqan::CharString> quals;
    readRecords(ids, seqs, quals, inFile);

    for (unsigned i = 0; i < length(ids); ++i)
        std::cout << ids[i] << '\t' << seqs[i] << '\t' << quals[i] << '\n';

    return 0;
}
