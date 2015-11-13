#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    CharString bamFileName = getAbsolutePath("demos/tutorial/sam_and_bam_io/example.sam");

    BamFileIn bamFileIn(toCString(bamFileName));

    BamHeader header;
    readHeader(header, bamFileIn);

    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    TBamContext const & bamContext = context(bamFileIn);

    for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
        std::cout << contigNames(bamContext)[i] << '\t'
                  << contigLengths(bamContext)[i] << '\n';

    return 0;
}
