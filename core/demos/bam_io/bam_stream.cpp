#include <iostream>
#include <seqan/file.h>
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    std::string pathSam = std::string(SEQAN_PATH_TO_ROOT()) + "/core/demos/bam_io/example.sam";
   
    BamStream bamStreamIn(pathSam.c_str());
    if (!isGood(bamStreamIn))
    {
        std::cerr << "Can't open the file." << std::endl;
        return 1;
    }
    
    // Open output stream. The value "-" means reading from stdin or writing to stdout.
    BamStream bamStreamOut("-", BamStream::WRITE);
    // Copy header. The header is automatically written out before the first record.
    bamStreamOut.header = bamStreamIn.header;

    // BamAlignmentRecord stores one record at a time.
    BamAlignmentRecord record;
    while (!atEnd(bamStreamIn))
    {
        if (readRecord(record, bamStreamIn) != 0)
        {
            std::cerr << "Can't read record!" << std::endl;
            return 1;
        }
        if (writeRecord(bamStreamOut, record) != 0)
        {
            std::cout << "Can't write record!" << std::endl;
            return 1;
        }
    }
    return 0;
}

