//![includes]
#define SEQAN_PROFILE // enable time measurements
#include <seqan/seq_io.h>
#include <iostream>

using namespace seqan;
//![includes]

//![open_file]
int main(int argc, char const * argv[])
{
    SEQAN_PROTIMESTART(loadTime);

    SeqFileIn seqFile;
    if (argc < 2 || !open(seqFile, argv[1]))
        return 1;
//![open_file]

//![read_sequences]
    StringSet<String<Dna5Q> > seqs;
    StringSet<CharString> seqIDs;
    String<Dna5Q> seq;
    CharString qual;
    CharString id;

    size_t seqCount = 0;
    for (; !atEnd(seqFile); ++seqCount)
    {
        readRecord(id, seq, qual, seqFile);     // read record

        // convert ascii to values from 0..62
        // store dna and quality together in Dna5Q
        assignQualities(seq, qual);

        // we use reserve and append, as assign is not supported
        // by StringSet<..., Owner<ConcatDirect<> > >
        appendValue(seqs, seq, Generous());
        appendValue(seqIDs, id, Generous());
    }
//![read_sequences]

//![output]
    std::cout << "Loading " << seqCount << " sequences took " << SEQAN_PROTIMEDIFF(loadTime);
    std::cout << " seconds." << std::endl << std::endl;
    for (unsigned i = 0; i < seqCount && i < 10; ++i)
    {
        std::cout << '>' << seqIDs[i] << std::endl;
        std::cout << seqs[i] << std::endl;
    }

    return 0;
}
//![output]
