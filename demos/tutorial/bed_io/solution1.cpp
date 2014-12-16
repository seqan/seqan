#include <seqan/bed_io.h>

int main()
{
    // Open input bed file.
    seqan::BedFileIn bedIn;
    if (!open(bedIn, "example.bed"))
    {
        std::cerr << "ERROR: Could not open example.bed" << std::endl;
        return 1;
    }

    // Attach to standard output.
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
    catch (seqan::Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
