#include <seqan/bed_io.h>
#include <seqan/misc/name_store_cache.h>

using namespace seqan;

int main()
{
    // Open input bed file.
    BedFileIn bedIn;
    if (!open(bedIn, toCString(getAbsolutePath("demos/tutorial/bed_io/example.bed"))))
    {
        std::cerr << "ERROR: Could not open example.bed\n";
        return 1;
    }

    // Array of counters and sequence names.
    String<unsigned> counters;
    StringSet<CharString> seqNames;
    NameStoreCache<StringSet<CharString> > cache(seqNames);

    // Read the file record by record.
    BedRecord<Bed3> record;

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
