#include <iostream>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    std::string pathSam = getAbsolutePath("/demos/bam_io/example.sam");

    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(pathSam)))
    {
        std::cerr << "Can't open the file." << std::endl;
        return 1;
    }

    // Open output stream. The value "-" means reading from stdin or writing to stdout.
    BamFileOut bamFileOut(bamFileIn);
    open(bamFileOut, std::cout, Sam());

    // Copy header. The header is automatically written out before the first record.
    BamHeader header;
    readHeader(header, bamFileIn);
    writeHeader(bamFileOut, header);

    // BamAlignmentRecord stores one record at a time.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        writeRecord(bamFileOut, record);
    }
    return 0;
}
