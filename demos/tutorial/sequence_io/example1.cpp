#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    CharString id;
    Dna5String seq;

    SeqFileIn seqFileIn("example.fa");
    readRecord(id, seq, seqFileIn);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}
