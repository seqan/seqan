#include <iostream>
#include <seqan/file.h>
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    std::string pathSam = std::string(SEQAN_PATH_TO_ROOT()) + "/core/demos/bam_io/example.sam";
   
    BamFileIn bamFileIn;
    if (!open(bamFileIn, pathSam.c_str()))
    {
        std::cerr << "Can't open the file." << std::endl;
        return 1;
    }

    // Open output stream. The value "-" means reading from stdin or writing to stdout.
    BamFileOut bamFileOut(bamFileIn);
    open(bamFileOut, std::cout, Sam());

    // Copy header. The header is automatically written out before the first record.
    BamHeader header;
    readRecord(header, bamFileIn);
    write(bamFileOut, header);

    // BamAlignmentRecord stores one record at a time.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        write(bamFileOut, record);
    }
    return 0;
}

