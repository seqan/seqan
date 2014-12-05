#include <seqan/bam_io.h>

int main()
{
    // Open input file.
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, "example.sam"))
    {
        std::cerr << "ERROR: Could not open example.sam!" << std::endl;
        return 1;
    }

    unsigned numUnmappedReads = 0;

    try
    {
        // Read header.
        seqan::BamHeader header;
        readRecord(header, bamFileIn);

        // Read records.
        seqan::BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            if (hasFlagUnmapped(record))
                numUnmappedReads += 1;
        }
    }
    catch (seqan::Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Number of unmapped reads: " << numUnmappedReads << "\n";

    return 0;
}
