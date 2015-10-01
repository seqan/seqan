#include <seqan/vcf_io.h>

using namespace seqan;

int main()
{
    try
    {
        // Open input file.
        VcfFileIn vcfIn(toCString(getAbsolutePath("demos/tutorial/vcf_io/example.vcf")));

        // Copy over header.
        VcfHeader header;
        readHeader(header, vcfIn);

        // Get array of counters.
        String<unsigned> counters;
        unsigned contigsCount = length(contigNames(context(vcfIn)));
        resize(counters, contigsCount, 0);

        // Read the file record by record.
        VcfRecord record;
        while (!atEnd(vcfIn))
        {
            readRecord(record, vcfIn);

            // Register record with counters.
            counters[record.rID] += 1;
        }

        // Print result.
        std::cout << "VARIANTS ON CONTIGS\n";
        for (unsigned i = 0; i < contigsCount; ++i)
            std::cout << contigNames(context(vcfIn))[i] << '\t'
                      << counters[i] << '\n';
    }
    catch (seqan::Exception const & e)
    {
        std::cerr << "ERROR:" << e.what() << std::endl;
        return 1;
    }

    return 0;
}
