#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    CharString seqFileName = getAbsolutePath("demos/tutorial/sequence_io/example.fa");
    CharString id;
    Dna5String seq;

    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seq, seqFileIn);
    std::cout << id << '\t' << seq << '\n';

    return 0;
}
