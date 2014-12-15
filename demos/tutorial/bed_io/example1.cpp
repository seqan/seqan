#include <seqan/basic.h>
#include <seqan/bed_io.h>

using namespace seqan;

int main()
{
    // Open input bed file.
    BedFileIn bedIn("example.bed");
    // Open output bed file and link to stdout.
    BedFileOut bedOut(std::cout, Bed());

    // Read the file record by record.
    BedRecord<Bed3> record;
    while (!atEnd(bedIn))
    {
        readRecord(record, bedIn);
        writeRecord(bedOut, record);
    }
    
    return 0;
}
