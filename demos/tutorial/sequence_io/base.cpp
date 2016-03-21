#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
//![batch]
    CharString seqFileName = getAbsolutePath("/demos/tutorial/sequence_io/example.fa");

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    SeqFileIn seqFileIn(toCString(seqFileName));

    // Reads up to 10 records.
    readRecords(ids, seqs, seqFileIn, 10);

    // Reads all remaining records.
    readRecords(ids, seqs, seqFileIn);
//![batch]
//![qual]
    CharString seqFqFileName = getAbsolutePath("/demos/tutorial/sequence_io/example.fq");
    CharString id;
    Dna5String seq;
    CharString qual;

    SeqFileIn seqFqFileIn(toCString(seqFqFileName));

    readRecord(id, seq, qual, seqFqFileIn);
//![qual]
    return 0;
}
