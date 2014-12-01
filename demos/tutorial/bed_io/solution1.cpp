#include <seqan/basic.h>
#include <seqan/bed_io.h>

int main()
{
    // Open input bed file.
    seqan::BedFileIn bedIn;
    if (!open(bedIn, "example.bed"))
    {
        std::cerr << "ERROR: Could not open example.bed\n";
        return 1;
    }
    // Open output bed file and link to stdout.
    seqan::BedFileOut bedOut(std::cout, seqan::Bed());

    // Read the file record by record.
    seqan::BedRecord<seqan::Bed3> record;

    try
    {
        while (!atEnd(bedIn))
        {
            readRecord(record, bedIn);
            writeRecord(bedOut, record);
        }
    }
    catch (std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
