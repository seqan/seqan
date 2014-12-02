#include <seqan/bam_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 3)
    {
        std::cerr << "USAGE: " << argv[0] << " INPUT.sam OUTPUT.bam" << "\n";
        return 1;
    }

    // Open input SAM file.
    seqan::BamFileIn samFileIn;
    if (!open(argv[1]))
    {
        std::cerr << "ERROR: could not open file " << argv[1] << ".\n";
        return 1;
    }

    // Open output BAM file.
    seqan::BamFileOut bamFileOut;
    if (!open(argv[2]))
    {
        std::cerr << "ERROR: could not open file " << argv[2] << ".\n";
        return 1;
    }

    // Copy header.
    seqan::BamHeader header;
    try
    {
      readRecord(header, samFileIn);
      writeRecord(bamFileOut, header);
    }
    catch (seqan::IOError const & e)
    {
        std::cerr << "ERROR: could not copy header.\n";
    }

    // Copy all records.
    seqan::BamAlignmentRecord record;
    while (!atEnd(samFileIn))
    {
        try
        {
            readRecord(header, samFileIn);
            writeRecord(bamFileOut, record);
        }
        catch (seqan::IOError const & e)
        {
            std::cerr << "ERROR: could not copy record.\n";
        }
    }

    return 0;
}
