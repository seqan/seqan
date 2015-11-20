#include <seqan/gff_io.h>
#include <seqan/misc/name_store_cache.h>

using namespace seqan;

int main()
{
    // Get path to example file.
    CharString file = getAbsolutePath("demos/tutorial/gff_and_gtf_io/example.gff");

    // Open input file.
    GffFileIn gffIn;
    if (!open(gffIn, toCString(file)))
    {
        std::cerr << "ERROR: Could not open example.gff" << std::endl;
        return 1;
    }

    // Array of counters and sequence names.
    String<unsigned> counters;
    StringSet<CharString> seqNames;
    NameStoreCache<StringSet<CharString> > cache(seqNames);

    // Read the file record by record.
    GffRecord record;

    try
    {
        while (!atEnd(gffIn))
        {
            readRecord(record, gffIn);
            unsigned rID = nameToId(cache, record.ref);

            // Resize counters if necessary and increment counter.
            assignValueById(counters, rID, getValueById(counters, rID) + 1);
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    // Print result.
    std::cout << "RECORDS ON CONTIGS\n";
    for (unsigned i = 0; i < length(seqNames); ++i)
        if (counters[i] != 0u)
            std::cout << seqNames[i] << '\t' << counters[i] << '\n';

    return 0;
}
