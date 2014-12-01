#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn("example.sam");
    // Open output file, BamFileOut accepts also an ostream and a format tag.
    seqan::BamFileOut bamFileOut(std::cout, seqan::Sam());
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

    return 0;
}
