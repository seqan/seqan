#include <seqan/seq_io.h>
#include <iostream>

int main()
{
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    seqan::StringSet<seqan::CharString> quals;

    seqan::SeqFileIn seqIn(std::cin);
    seqan::readRecords(ids, seqs, quals, seqIn);

    seqan::SeqFileOut seqOut(std::cout, seqan::Fastq());
    seqan::writeRecords(seqOut, ids, seqs, quals);

    return 0;
}
