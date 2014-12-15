#include <seqan/bam_io.h>

int main()
{
    seqan::BamFileIn bamFileIn("example.bam");

    seqan::BamHeader header;
    readRecord(header, bamFileIn);

    typedef seqan::SmartFileContext<seqan::BamFileIn, void>::Type TBamContext;

    TBamContext const & bamContext = context(bamFileIn);

    for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
        std::cout << contigNames(bamContext)[i] << '\t'
                  << contigLengths(bamContext)[i] << '\n';

    return 0;
}
