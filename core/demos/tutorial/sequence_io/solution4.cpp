#include <seqan/file.h>
#include <iostream>

int main (int argc, char const ** argv)
{
    seqan::MultiSeqFile multiSeqFile;
    if (argc < 2 || !open(multiSeqFile.concat, argv[1], seqan::OPEN_RDONLY))
        return 1;

    seqan::AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

    seqan::String<seqan::Dna5> seq;
    seqan::CharString qual;
    seqan::CharString id;

    for (unsigned i = 0; i < length(multiSeqFile); ++i)
    {
        assignSeq(seq, multiSeqFile[i], format);    // read sequence
        assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignSeqId(id, multiSeqFile[i], format);   // read sequence id
        
        std::cout << id << '\t' << seq << '\t' << qual << '\n';
    }

    return 0;
}
