#include <seqan/basic.h>
#include <seqan/vcf_io.h>

using namespace seqan;

int main()
{
    try
    {
        // Open input stream.
        VcfFileIn vcfIn("example.vcf");
        // Open output stream
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
    }
    catch (IOError &e)
    {
        std::cerr << "=== I/O Error ===\n" << e.what() << std::endl;
        return 1;
    }
    catch (ParseError &e)
    {
        std::cerr << "=== Parse Error ===\n" << e.what() << std::endl;
        return 1;
    }

    return 0;
}
