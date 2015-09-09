#include <seqan/basic.h>
#include <seqan/vcf_io.h>

using namespace seqan;

int main()
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/demos/vcf_io/example.vcf");

    VcfFileIn vcfIn;
    if (!open(vcfIn, toCString(path)))
    {
        std::cerr << "ERROR: Could not open " << path << " for reading!\n";
        return 1;
    }

    VcfHeader header;
    readHeader(header, vcfIn);

    VcfRecord record;
    while (!atEnd(vcfIn))
    {
        readRecord(record, vcfIn);

        // Note that we print the position 1-based since we use text output
        // whereas it is 0-based in the VcfRecord.
        std::cout << contigNames(context(vcfIn))[record.rID]
                  << "\t" << (record.beginPos + 1) << "\n";
    }

    return 0;
}
