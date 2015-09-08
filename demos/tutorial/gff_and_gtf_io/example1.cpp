#include <seqan/gff_io.h>

using namespace seqan;

int main()
{
    // Open input file.
    GffFileIn gffIn("example.gff");

    // Attach to standard output.
    GffFileOut gffOut(std::cout, Gff());

    // Copy the file record by record.
    GffRecord record;
    while (!atEnd(gffIn))
    {
        readRecord(record, gffIn);
        writeRecord(gffOut, record);
    }

    return 0;
}
