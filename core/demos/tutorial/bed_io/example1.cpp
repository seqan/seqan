#include <seqan/basic.h>
#include <seqan/bed_io.h>

int main()
{
    // Open input bed file.
    seqan::BedFileIn bedIn("example.bed");
    // Open output bed file and link to stdout.
    seqan::BedFileOut bedOut(std::cout, seqan::Bed());

    // Read the file record by record.
    seqan::BedRecord<seqan::Bed3> record;
    while (!atEnd(bedIn))
    {
        readRecord(record, bedIn);
        writeRecord(bedOut, record);
    }
    
    return 0;
}
