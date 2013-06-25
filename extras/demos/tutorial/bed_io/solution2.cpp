#include <seqan/basic.h>
#include <seqan/bed_io.h>

int main()
{
    // Open input stream
    seqan::BedStream bedIn("example.bed");
    if (!isGood(bedIn))
    {
        std::cerr << "ERROR: Could not open example.bed\n";
        return 1;
    }

    // Array of counters and sequence names.
    seqan::String<unsigned> counters;
    seqan::StringSet<seqan::CharString> seqNames;

    // Read the file record by record.
    seqan::BedRecord<seqan::Bed3> record;
    while (!atEnd(bedIn))
    {
        if (readRecord(record, bedIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from example.bed\n";
            return 1;
        }

        // Resize counters and write seqNames if necessary.
        if ((int)length(counters) <= record.rID)
        {
            resize(counters, record.rID + 1, 0);
            resize(seqNames, record.rID + 1);
        }
        if (counters[record.rID] == 0)
            seqNames[record.rID] = record.ref;

        // Register record with counters.
        counters[record.rID] += 1;
    }

    // Print result.
    std::cout << "RECORDS ON CONTIGS\n";
    for (unsigned i = 0; i < length(seqNames); ++i)
        if (counters[i] > 0u)
            std::cout << seqNames[i] << '\t' << counters[i] << '\n';
    
    return 0;
}
