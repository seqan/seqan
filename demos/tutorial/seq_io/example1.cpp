#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::SeqFileIn seqFileIn("example.fa");
    readRecord(id, seq, seqFileIn);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}
