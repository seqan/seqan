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
    seqan::StringSet<seqan::CharString> quals;
    appendValue(quals, "IIIIIIIIIHII");
    appendValue(quals, "IIIIIIIIIIII");

    seqan::DirectionIterator<std::ostream, seqan::Output>::Type streamIter(std::cout);
    for (unsigned i = 0; i < length(ids); ++i)
        seqan::writeRecord(streamIter, ids[i], seqs[i], quals[i], seqan::Fastq());

    return 0;
}
