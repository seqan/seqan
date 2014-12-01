#include <seqan/basic.h>
#include <seqan/vcf_io.h>

int main()
{
    try
    {
        // Open input stream.
        seqan::VcfFileIn vcfIn("example.vcf");

        // Copy over header.
        seqan::VcfHeader header;
        readRecord(header, vcfIn);

        // Get array of counters.
        seqan::String<unsigned> counters;
        resize(counters, length(contigNames(context(vcfIn))), 0);

        // Read the file record by record.
        seqan::VcfRecord record;
        while (!atEnd(vcfIn))
        {
            readRecord(record, vcfIn);

            // Register record with counters.
            counters[record.rID] += 1;
        }

        // Print result.
        std::cout << "VARIANTS ON CONTIGS\n";
        for (unsigned i = 0; i < length(contigNames(context(vcfIn))); ++i)
            std::cout << contigNames(context(vcfIn))[i] << '\t'
                      << counters[i] << '\n';
    }
    catch (seqan::IOError &e)
    {
        std::cerr << "=== I/O Error ===\n" << e.what() << std::endl;
        return 1;
    }
    catch (seqan::ParseError &e)
    {
        std::cerr << "=== Parse Error ===\n" << e.what() << std::endl;
        return 1;
    }

    return 0;
}
