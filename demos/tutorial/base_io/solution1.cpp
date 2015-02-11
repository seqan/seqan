#include <seqan/bam_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc < 3)
    {
        std::cerr << "USAGE: " << argv[0] << " INPUT.bam OUTPUT.sam" << "\n";
        return 1;
    }

    // Open input BAM file.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, argv[1]))
    {
        std::cerr << "ERROR: could not open input file " << argv[1] << ".\n";
        return 1;
    }

    // Open output SAM file.
    BamFileOut samFileOut;
    if (!open(samFileOut, argv[2]))
    {
        std::cerr << "ERROR: could not open output file " << argv[2] << ".\n";
        return 1;
    }

    // Copy header.
    BamHeader header;
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
            readHeader(header, bamFileIn);
            writeRecord(samFileOut, record);
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }

    return 0;
}
