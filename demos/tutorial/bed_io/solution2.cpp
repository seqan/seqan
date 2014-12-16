#include <seqan/bed_io.h>
#include <seqan/misc/name_store_cache.h>

int main()
{
    // Open input bed file.
    seqan::BedFileIn bedIn;
    if (!open(bedIn, "example.bed"))
    {
        std::cerr << "ERROR: Could not open example.bed\n";
        return 1;
    }

    // Array of counters and sequence names.
    seqan::String<unsigned> counters;
    seqan::StringSet<seqan::CharString> seqNames;
    seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > cache(seqNames);

    // Read the file record by record.
    seqan::BedRecord<seqan::Bed3> record;

    try
    {
        while (!atEnd(bedIn))
        {
            readRecord(record, bedIn);
            unsigned rID = nameToId(cache, record.ref);

            // Resize counters if necessary and increment counter.
            assignValueById(counters, rID, getValueById(counters, rID) + 1);
        }
    }
    catch (seqan::Exception const & e)
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
