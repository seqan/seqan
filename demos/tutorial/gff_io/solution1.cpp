#include <seqan/gff_io.h>

int main()
{
    // Open input file.
    seqan::GffFileIn gffIn;
    if (!open(gffIn, "example.gff"))
    {
        std::cerr << "ERROR: Could not open example.gff" << std::endl;
        return 1;
    }

    // Attach to standard output.
    seqan::GffFileOut gffOut(std::cout, seqan::Gff());

    // Copy the file record by record.
    seqan::GffRecord record;

    try
    {
        while (!atEnd(gffIn))
        {
            readRecord(record, gffIn);
            writeRecord(gffOut, record);
        }
    }
    catch (seqan::Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
