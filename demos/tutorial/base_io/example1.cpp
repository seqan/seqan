// FRAGMENT(include)
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // FRAGMENT(ctor)
    // Open input BAM file, BamFileIn supports both SAM and BAM files.
    BamFileIn bamFileIn("example.bam");

    // FRAGMENT(open)
    // Open output SAM file by passing the filename to open.
    BamFileOut samFileOut;
    open(samFileOut, "example.sam");

    // FRAGMENT(header)
    // Copy header.
    BamHeader header;
    readRecord(header, bamFileIn);
    writeRecord(samFileOut, header);

    // FRAGMENT(records)
    // Copy all records.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        writeRecord(samFileOut, record);
    }

    return 0;
}
