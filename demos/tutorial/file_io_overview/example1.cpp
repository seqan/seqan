//![include]
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
//![include]
//![ctor]
    CharString bamFileInName = getAbsolutePath("demos/tutorial/file_io_overview/example.bam");
    CharString samFileOutName = getAbsolutePath("demos/tutorial/file_io_overview/example.sam");

    // Open input BAM file, BamFileIn supports both SAM and BAM files.
    BamFileIn bamFileIn(toCString(bamFileInName));

    // Open output SAM file by passing the context of bamFileIn and the filename to open.
    BamFileOut samFileOut(context(bamFileIn), toCString(samFileOutName));
//![ctor]

//![open]
    // Alternative way to open a bam or sam file
    BamFileIn openBamFileIn;
    open(openBamFileIn, toCString(bamFileInName));
//![open]

//![header]
    // Copy header.
    BamHeader header;
    readHeader(header, bamFileIn);
    writeHeader(samFileOut, header);
//![header]

//![records]
    // Copy all records.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        writeRecord(samFileOut, record);
    }

    return 0;
}
//![records]
