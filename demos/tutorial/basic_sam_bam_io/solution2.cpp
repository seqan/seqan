#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, "example.sam"))
    {
        std::cerr << "ERROR: Could not open example.sam!" << std::endl;
        return 1;
    }

    try
    {
        // Copy header.
        seqan::BamHeader header;
        readRecord(header, bamFileIn);

        unsigned numUnmappedReads = 0;
        seqan::BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            if (hasFlagUnmapped(record))
                numUnmappedReads += 1;
        }

        std::cout << "Number of unmapped reads: " << numUnmappedReads << "\n";
    }
    catch (std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
