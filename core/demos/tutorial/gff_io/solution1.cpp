#include <seqan/basic.h>
#include <seqan/gff_io.h>

int main()
{
    // Open input stream.
    seqan::GffFileIn gffIn;
    if (!open(gffIn, "example.gff"))
    {
        std::cerr << "ERROR: Could not open example.gff\n";
        return 1;
    }
    // Open output stream. If target is a ostream we must specify the format.
    seqan::GffFileOut gffOut(std::cout, seqan::Gff());

    // Read the file record by record.
    seqan::GffRecord record;

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
