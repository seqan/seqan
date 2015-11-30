#include <seqan/bed_io.h>
using namespace seqan;

int main()
{
    // Open input bed file.
    BedFileIn bedIn(toCString(getAbsolutePath("demos/tutorial/bed_io/example.bed")));

    // Attach to standard output.
    BedFileOut bedOut(std::cout, Bed());

    // Copy the file record by record.
    BedRecord<Bed3> record;

    while (!atEnd(bedIn))
    {
        readRecord(record, bedIn);
        writeRecord(bedOut, record);
    }

    return 0;
}
