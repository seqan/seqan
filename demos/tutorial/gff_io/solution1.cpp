#include <seqan/basic.h>
#include <seqan/gff_io.h>

using namespace seqan;

int main()
{
    // Open input stream.
    GffFileIn gffIn;
    if (!open(gffIn, "example.gff"))
    {
        std::cerr << "ERROR: Could not open example.gff\n";
        return 1;
    }
    // Open output stream. If target is a ostream we must specify the format.
    GffFileOut gffOut(std::cout, Gff());

    // Read the file record by record.
    GffRecord record;

    try
    {
        while (!atEnd(gffIn))
        {
            readRecord(record, gffIn);
            writeRecord(gffOut, record);
        }
    }
    catch (std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
