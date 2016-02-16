#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>
#include <seqan/tabix_io.h>

using namespace seqan;

int main(int argc, char const * argv[])
{
    if (argc < 5)
    {
        std::cerr << "USAGE: " << argv[0] << " <variants.vcf.gz> <contig> <begin> <end>\n";
        return 0;
    }

    // Open VCF file
    VcfFileIn vcfFile;
    if (!open(vcfFile, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << " for reading.\n";
        return 1;
    }

    // Read header (to get the contig names)
    seqan::VcfHeader header;
    readHeader(header, vcfFile);
    
    // Open Tabix index
    std::string tbiFileName = (std::string)argv[1] + ".tbi";
    TabixIndex tabixIndex;
    if (!open(tabixIndex, tbiFileName.c_str()))
    {
        std::cerr << "ERROR: Could not read Tabix index file " << tbiFileName << "\n";
        return 1;
    }

    // Search overlapping variants
    bool hasEntries = false;
    int begPos = lexicalCast<int>(argv[3]);
    int endPos = lexicalCast<int>(argv[4]);
    if (!jumpToRegion(vcfFile,
                      hasEntries,
                      argv[2],
                      begPos,
                      endPos,
                      tabixIndex))
    {
        std::cerr << "Contig " << argv[2] << " not found!\n";
        return 1;
    }

    seqan::VcfRecord record;
    while (hasEntries && !atEnd(vcfFile))
    {
        readRecord(record, vcfFile);

        // If we are on the next reference or at the end already then we stop.
        if (record.rID == -1 || contigNames(context(vcfFile))[record.rID] != argv[2] || record.beginPos >= endPos)
            break;
        
        std::cout << record.beginPos << '\t' << record.ref << '\t' << record.alt << std::endl;
    }

    return 0;
}
