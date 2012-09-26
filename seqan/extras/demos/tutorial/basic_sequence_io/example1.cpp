#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main()
{
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SequenceStream seqStream("example.fa");
    readRecord(id, seq, seqStream);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}
