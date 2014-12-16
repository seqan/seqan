#include <seqan/vcf_io.h>

int main()
{
    try
    {
        // Open input file.
        seqan::VcfFileIn vcfIn("example.vcf");

        // Copy over header.
        seqan::VcfHeader header;
        readRecord(header, vcfIn);

        // Get array of counters.
        seqan::String<unsigned> counters;
        unsigned contigsCount = length(contigNames(context(vcfIn)));
        resize(counters, contigsCount, 0);

        // Read the file record by record.
        seqan::VcfRecord record;
        while (!atEnd(vcfIn))
        {
            readRecord(record, vcfIn);

            // Register record with counters.
            counters[record.rID] += 1;
        }
    }
    catch (seqan::Exceptioon const & e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // Print result.
    std::cout << "VARIANTS ON CONTIGS\n";
    for (unsigned i = 0; i < contigsCount; ++i)
        std::cout << contigNames(context(vcfIn))[i] << '\t'
                  << counters[i] << '\n';

    return 0;
}
