#include <seqan/basic.h>
#include <seqan/vcf_io.h>

using namespace seqan;

// USAGE: vcf_stream_read IN.vcf
//
// Read VCF file and print the positions (reference, position) to stdout.

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[1] << " IN.vcf\n";
        return 1;
    }

    VcfStream vcfIn(argv[1]);
    if (!isGood(vcfIn))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading!\n";
        return 1;
    }

    VcfRecord record;
    while (!atEnd(vcfIn))
    {
        if (readRecord(record, vcfIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from " << argv[1] << "\n";
            return 1;
        }

        // Note that we print the position 1-based since we use text output
        // whereas it is 0-based in the VcfRecord.
        std::cout << vcfIn.header.sequenceNames[record.rID]
                  << "\t" << (record.beginPos + 1) << "\n";
    }
    
    return 0;
}
