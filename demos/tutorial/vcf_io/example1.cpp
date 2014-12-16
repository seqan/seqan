#include <seqan/vcf_io.h>

int main()
{
    // Open input file.
    seqan::VcfFileIn vcfIn("example.vcf");

    // Attach to standard output.
    seqan::VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, seqan::Vcf());

    // Copy over header.
    seqan::VcfHeader header;
    readRecord(header, vcfIn);
    writeRecord(vcfOut, header);

    // Copy the file record by record.
    seqan::VcfRecord record;
    while (!atEnd(vcfIn))
    {
        readRecord(record, vcfIn);
        writeRecord(vcfOut, record);
    }
    
    return 0;
}
