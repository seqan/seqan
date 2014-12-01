#include <seqan/basic.h>
#include <seqan/vcf_io.h>

int main()
{
    try
    {
        // Open input stream.
        seqan::VcfFileIn vcfIn("example.vcf");
        // Open output stream
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
