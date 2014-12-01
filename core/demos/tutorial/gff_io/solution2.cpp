#include <seqan/basic.h>
#include <seqan/gff_io.h>
#include <seqan/misc/misc_name_store_cache.h>

int main()
{
    // Open input gff file.
    seqan::GffFileIn gffIn;
    if (!open(gffIn, "example.gff"))
    {
        std::cerr << "ERROR: Could not open example.gff\n";
        return 1;
    }

    // Array of counters and sequence names.
    seqan::String<unsigned> counters;
    seqan::StringSet<seqan::CharString> seqNames;
    seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > cache(seqNames);

    // Read the file record by record.
    seqan::GffRecord record;

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
    catch (std::runtime_error &e)
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
