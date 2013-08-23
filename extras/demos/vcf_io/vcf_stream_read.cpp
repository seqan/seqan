#include <seqan/basic.h>
#include <seqan/vcf_io.h>

using namespace seqan;

int main()
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/extras/demos/vcf_io/example.vcf");

    VcfStream vcfIn(toCString(path));
    if (!isGood(vcfIn))
    {
        std::cerr << "ERROR: Could not open " << path << " for reading!\n";
        return 1;
    }

    VcfRecord record;
    while (!atEnd(vcfIn))
    {
        if (readRecord(record, vcfIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from " << path << "\n";
            return 1;
        }

        // Note that we print the position 1-based since we use text output
        // whereas it is 0-based in the VcfRecord.
        std::cout << vcfIn.header.sequenceNames[record.rID]
                  << "\t" << (record.beginPos + 1) << "\n";
    }
    
    return 0;
}
