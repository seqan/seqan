#include <seqan/basic.h>
#include <seqan/vcf_io.h>

int main()
{
    // Open input stream.
    seqan::VcfStream vcfIn("example.vcf");
    // Open output stream, filename "-" means stdout.
    seqan::VcfStream vcfOut("-", seqan::VcfStream::WRITE);

    // Copy over header.
    vcfOut.header = vcfIn.header;

    // Read the file record by record.
    seqan::VcfRecord record;
    while (!atEnd(vcfIn))
    {
        readRecord(record, vcfIn);
        writeRecord(vcfOut, record);
    }
    
    return 0;
}
