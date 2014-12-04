#include <seqan/bam_io.h>

int main()
{
    seqan::BamFileIn bamFileIn("example.bam");

    seqan::BamHeader header;
    readRecord(header, bamFileIn);

    typedef typename seqan::SmartFileContext<seqan::BamFileIn>::Type TBamContext;

    TBamContext const & bamContext = context(bamFileIn);

    for (unsigned i = 0; i < length(nameStore(bamContext)); ++i)
        std::cout << nameStore(bamContext)[i] << '\t'
                  << sequenceLengths(bamContext)[i] << '\n';

    return 0;
}
