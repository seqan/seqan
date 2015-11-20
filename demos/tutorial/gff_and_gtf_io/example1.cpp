#include <seqan/basic.h>
#include <seqan/gff_io.h>

using namespace seqan;

int main()
{
    // Get path to example file.
    CharString file = getAbsolutePath("demos/tutorial/gff_and_gtf_io/example.gff");

    // Open input file.
    GffFileIn gffIn(toCString(file));

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
