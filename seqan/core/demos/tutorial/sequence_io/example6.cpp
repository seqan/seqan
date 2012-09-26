#include <seqan/seq_io.h>
#include <iostream>

int main()
{
    seqan::StringSet<seqan::CharString> ids;
    appendValue(ids, "id1");
    appendValue(ids, "id2");
    seqan::StringSet<seqan::Dna5String> seqs;
    appendValue(seqs, "CGATCGATCGAT");
    appendValue(seqs, "AAAAAAAAAAAA");

    for (unsigned i = 0; i < length(ids); ++i)
        if (seqan::writeRecord(std::cout, ids[i], seqs[i], seqan::Fasta()) != 0)
            return 1;  // Error writing.

    return 0;
}
