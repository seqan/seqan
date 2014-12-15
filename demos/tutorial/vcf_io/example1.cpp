#include <seqan/basic.h>
#include <seqan/vcf_io.h>

using namespace seqan;

int main()
{
    // Open input stream.
    VcfFileIn vcfIn("example.vcf");
    // Open output stream, filename "-" means stdout.
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, Vcf());

    // Copy over header.
    VcfHeader header;
    readRecord(header, vcfIn);
    writeRecord(vcfOut, header);

    // Read the file record by record.
    VcfRecord record;
    while (!atEnd(vcfIn))
    {
        readRecord(record, vcfIn);
        writeRecord(vcfOut, record);
    }
    
    return 0;
}
