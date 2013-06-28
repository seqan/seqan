#include <seqan/basic.h>
#include <seqan/gff_io.h>

int main()
{
    // Open input stream
    seqan::GffStream gffIn("example.gff");
    if (!isGood(gffIn))
    {
        std::cerr << "ERROR: Could not open example.gff\n";
        return 1;
    }

    // Array of counters and sequence names.
    seqan::String<unsigned> counters;
    seqan::StringSet<seqan::CharString> seqNames;

    // Read the file record by record.
    seqan::GffRecord record;
    while (!atEnd(gffIn))
    {
        if (readRecord(record, gffIn) != 0)
        {
            std::cerr << "ERROR: Problem reading from example.gff\n";
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
