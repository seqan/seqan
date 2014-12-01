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
    // Open output file, BamFileOut accepts also an ostream and a format tag.
    seqan::BamFileOut bamFileOut(std::cout, seqan::Sam());

    try
    {
        // Copy header.
        seqan::BamHeader header;
        readRecord(header, bamFileIn);
        writeRecord(bamFileOut, header);

        seqan::BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            writeRecord(bamFileOut, record);
        }
    }
    catch (std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
