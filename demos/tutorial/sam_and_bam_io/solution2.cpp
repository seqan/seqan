#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // Open input file.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, "example.sam"))
    {
        std::cerr << "ERROR: Could not open example.sam!" << std::endl;
        return 1;
    }

    unsigned numUnmappedReads = 0;

    try
    {
        // Read header.
        BamHeader header;
        readHeader(header, bamFileIn);

        // Read records.
        BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            if (hasFlagUnmapped(record))
                numUnmappedReads += 1;
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Number of unmapped reads: " << numUnmappedReads << "\n";

    return 0;
}
