#include <seqan/gff_io.h>

using namespace seqan;

int main()
{
    // Get path to example file.
    CharString file = getAbsolutePath("demos/tutorial/gff_and_gtf_io/example.gff");

    // Open input file.
    GffFileIn gffIn;
    if (!open(gffIn, toCString(file)))
    {
        std::cerr << "ERROR: Could not open example.gff" << std::endl;
        return 1;
    }

    // Attach to standard output.
    GffFileOut gffOut(std::cout, Gff());

    // Copy the file record by record.
    GffRecord record;

    try
    {
        while (!atEnd(gffIn))
        {
            readRecord(record, gffIn);
            writeRecord(gffOut, record);
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
