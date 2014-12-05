#include <seqan/bam_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 3)
    {
        std::cerr << "USAGE: " << argv[0] << " INPUT.bam OUTPUT.sam" << "\n";
        return 1;
    }

    // Open input BAM file.
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, argv[1]))
    {
        std::cerr << "ERROR: could not open input file " << argv[1] << ".\n";
        return 1;
    }

    // Open output SAM file.
    seqan::BamFileOut samFileOut;
    if (!open(samFileOut, argv[2]))
    {
        std::cerr << "ERROR: could not open output file " << argv[2] << ".\n";
        return 1;
    }

    // Copy header.
    seqan::BamHeader header;
    try
    {
      readRecord(header, bamFileIn);
      writeRecord(samFileOut, header);
    }
    catch (seqan::ParseError const & e)
    {
        std::cerr << "ERROR: input header is badly formatted. " << e.what() << "\n";
    }
    catch (seqan::IOError const & e)
    {
        std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }

    // Copy all records.
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(header, bamFileIn);
            writeRecord(samFileOut, record);
        }
        catch (seqan::ParseError const & e)
        {
            std::cerr << "ERROR: input record is badly formatted. " << e.what() << "\n";
        }
        catch (seqan::IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }

    return 0;
}
