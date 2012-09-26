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

    if (seqan::write2(std::cout, ids, seqs, quals, seqan::Fastq()) != 0)
        return 1;  // Error writing.

    return 0;
}
