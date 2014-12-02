// FRAGMENT(include)
#include <seqan/bam_io.h>

int main()
{
    // FRAGMENT(ctor)
    // Open input SAM file, BamFileIn supports both SAM and BAM files.
    seqan::BamFileIn samFileIn("example.sam");

    // FRAGMENT(open)
    // Open output BAM file by passing the filename to open.
    seqan::BamFileOut bamFileOut;
    open(bamFileOut, "example.bam");

    // FRAGMENT(header)
    // Copy header.
    seqan::BamHeader header;
    readRecord(header, samFileIn);
    writeRecord(bamFileOut, header);

    // FRAGMENT(records)
    // Copy all records.
    seqan::BamAlignmentRecord record;
    while (!atEnd(samFileIn))
    {
        readRecord(record, samFileIn);
        writeRecord(bamFileOut, record);
    }

    return 0;
}
