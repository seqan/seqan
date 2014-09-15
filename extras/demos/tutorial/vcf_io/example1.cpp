#include <seqan/basic.h>
#include <seqan/vcf_io.h>

int main()
{
    // Open input stream.
    seqan::VcfFileIn vcfIn("example.vcf");
    // Open output stream, filename "-" means stdout.
    seqan::VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, seqan::Vcf());

    // Copy over header.
    seqan::VcfHeader header;
    readRecord(header, vcfIn);
    writeRecord(vcfOut, header);

    // Read the file record by record.
    seqan::VcfRecord record;
    while (!atEnd(vcfIn))
    {
        readRecord(record, vcfIn);
        writeRecord(vcfOut, record);
    }
    
    return 0;
}
