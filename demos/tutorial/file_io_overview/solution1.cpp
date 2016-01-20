#include <seqan/bam_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    CharString bamFileInName = getAbsolutePath("demos/tutorial/file_io_overview/example.bam");
    CharString samFileOutName = getAbsolutePath("demos/tutorial/file_io_overview/example.sam");

    // Open input BAM file.
    BamFileIn bamFileIn;
    BamHeader header;
    if (!open(bamFileIn, toCString(bamFileInName)))
    {
        std::cerr << "ERROR: could not open input file " << bamFileInName << ".\n";
        return 1;
    }

    // Open output SAM file.
    BamFileOut samFileOut(context(bamFileIn), toCString(samFileOutName));

    // Copy header.
    try
    {
        readHeader(header, bamFileIn);
        writeHeader(samFileOut, header);
    }
    catch (IOError const & e)
    {
        std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }

    // Copy all records.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(record, bamFileIn);
            writeRecord(samFileOut, record);
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }

    return 0;
}
