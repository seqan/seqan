#include <seqan/vcf_io.h>

using namespace seqan;

int main()
{
    // Open input file.
    VcfFileIn vcfIn(toCString(getAbsolutePath("demos/tutorial/vcf_io/example.vcf")));

    // Attach to standard output.
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, Vcf());

    // Copy over header.
    VcfHeader header;
    readHeader(header, vcfIn);
    writeHeader(vcfOut, header);

    // Copy the file record by record.
    VcfRecord record;
    while (!atEnd(vcfIn))
    {
        readRecord(record, vcfIn);
        writeRecord(vcfOut, record);
    }
    
    return 0;
}
