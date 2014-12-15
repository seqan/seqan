#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    seqan::BamFileIn bamFileIn("example.bam");

    seqan::BamHeader header;
    readRecord(header, bamFileIn);

    typedef seqan::SmartFileContext<seqan::BamFileIn, void>::Type TBamContext;

    TBamContext const & bamContext = context(bamFileIn);

    for (unsigned i = 0; i < length(nameStore(bamContext)); ++i)
        std::cout << nameStore(bamContext)[i] << '\t'
                  << sequenceLengths(bamContext)[i] << '\n';

    return 0;
}
