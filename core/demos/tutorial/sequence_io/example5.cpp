// FRAGMENT(includes-main)
#include <seqan/file.h>
#include <iostream>

int main (int argc, char const ** argv)
{
// FRAGMENT(open)
    seqan::MultiSeqFile multiSeqFile;
    if (argc < 2 || !open(multiSeqFile.concat, argv[1], seqan::OPEN_RDONLY))
        return 1;

// FRAGMENT(guess)
    seqan::AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
// FRAGMENT(split)
    split(multiSeqFile, format);

// FRAGMENT(load)
    unsigned seqCount = length(multiSeqFile);
    seqan::StringSet<seqan::String<seqan::Dna5Q> > seqs;
    seqan::StringSet<seqan::CharString> seqIDs;

    reserve(seqs, seqCount, seqan::Exact());
    reserve(seqIDs, seqCount, seqan::Exact());

// FRAGMENT(buffers)
    seqan::String<seqan::Dna5Q> seq;
    seqan::CharString qual;
    seqan::CharString id;

// FRAGMENT(output)
    for (unsigned i = 0; i < seqCount; ++i)
    {
        assignSeq(seq, multiSeqFile[i], format);    // read sequence
        assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignSeqId(id, multiSeqFile[i], format);   // read sequence id

        // Convert ascii to values from 0..62.
        // We store DNA and quality together in Dna5Q.
        for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
            assignQualityValue(seq[j], (int)(seqan::ordValue(qual[j]) - 33));

        // We use reserve and append, as assign is not supported by
        // StringSet<..., Owner<ConcatDirect<> > >
        appendValue(seqs, seq);
        appendValue(seqIDs, id);
    }

// FRAGMENT(return)
    return 0;
}
