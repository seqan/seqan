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

    // Get array of counters.
    seqan::String<unsigned> counters;
    resize(counters, length(vcfIn.header.sequenceNames), 0);

    // Read the file record by record.
    seqan::VcfRecord record;
    while (!atEnd(vcfIn))
    {
        if (readRecord(record, vcfIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from example.vcf\n";
            return 1;
        }

        // Register record with counters.
        counters[record.rID] += 1;
    }

    // Print result.
    std::cout << "VARIANTS ON CONTIGS\n";
    for (unsigned i = 0; i < length(vcfIn.header.sequenceNames); ++i)
        std::cout << vcfIn.header.sequenceNames[i] << '\t'
                  << counters[i] << '\n';
    
    return 0;
}
