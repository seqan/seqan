#include <seqan/basic.h>
#include <seqan/vcf_io.h>

int main()
{
    // Open input stream
    seqan::VcfStream vcfIn("example.vcf");
    if (!isGood(vcfIn))
    {
        std::cerr << "ERROR: Could not open example.vcf\n";
        return 1;
    }
    // Open output stream, filename "-" means stdout.
    seqan::VcfStream vcfOut("-", seqan::VcfStream::WRITE);

    // Copy over header.
    vcfOut.header = vcfIn.header;

    // Read the file record by record.
    seqan::VcfRecord record;
    while (!atEnd(vcfIn))
    {
        if (readRecord(record, vcfIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from example.vcf\n";
            return 1;
        }
        if (writeRecord(vcfOut, record) != 0)
        {
            std::cerr << "ERROR: Problem writing to stdout.\n";
            return 1;
        }
    }
    
    return 0;
}
