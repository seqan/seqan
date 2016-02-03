#include <seqan/bed_io.h>

using namespace seqan;

int main()
{
    // Open input bed file.
    BedFileIn bedIn;
    if (!open(bedIn, toCString(getAbsolutePath("demos/tutorial/bed_io/example.bed"))))
    {
        std::cerr << "ERROR: Could not open example.bed" << std::endl;
        return 1;
    }
    // Attach to standard output.
    BedFileOut bedOut(std::cout, Bed());

    // Read the file record by record.
    BedRecord<Bed3> record;

    try
    {
        while (!atEnd(bedIn))
        {
            readRecord(record, bedIn);
            writeRecord(bedOut, record);
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
