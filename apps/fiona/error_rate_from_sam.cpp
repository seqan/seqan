#include <iostream>
#include <map>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

// USAGE: error_rate_from_sam IN.sam > OUT.txt

// There must only be one alignment per read.
//
// We ignore unaligned reads, and just count them.

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // Checking command line parameters.

    if (argc != 2)
    {
        std::cerr << "USAGE: error_rate IN.sam > OUT.txt\n";
        return 1;
    }

    // Opening file and record reader.

    std::fstream in(argv[1], std::ios::in | std::ios::binary);
    RecordReader<std::fstream, SinglePass<> > reader(in);

    // Allocating BAM header, reference name store, and cache.

    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;

    TNameStore refNameStore;
    TNameStoreCache refNameStoreCache(refNameStore);
    BamIOContext<TNameStore> context(refNameStore, refNameStoreCache);

    // Read header.

    BamHeader header;
    int res = 0;
    res = readRecord(header, context, reader, Sam());
    if (res != 0)
    {
        std::cerr << "Could not read SAM header!\n";
        return 1;
    }

    // Check that the SAM file is sorted by QNAME.
    CharString sortOrder;
    for (unsigned i = 0; i < length(header.records); ++i)
    {
        if (header.records[i].type != BAM_HEADER_FIRST)
            continue;
        unsigned idx = 0;
        if (findTagKey(idx, "SO", header.records[i]))
            sortOrder = header.records[i].tags[idx].i2;
    }
    if (sortOrder != "queryname")
    {
        std::cerr << "SAM file not sorted by 'queryname'!\n";
        return 1;
    }

    // Allocate data structures for counting / histogram building.
    unsigned totalBaseCount = 0;            // Number of bases read.
    unsigned totalErrorCount = 0;           // Number of errors read.
    unsigned totalReadCount = 0;            // Number of reads read.
    unsigned totalErrorneousReadCount = 0;  // Number of reads with errors read, excluding unaligned reads.
    unsigned totalUnalignedReadCount = 0;   // Number of reads without alignments.
    std::map<unsigned, unsigned> histo;     // Histogram error count -> num occurrences.

    // Read records

    BamAlignmentRecord record;
    CharString previousQName;
    bool previousIsFirst = false;

    while (!atEnd(reader))
    {
        res = readRecord(record, context, reader, Sam());
        if (res != 0)
        {
            std::cerr << "Error reading SAM record!\n";
            return 1;
        }

        // Skip secondary records.
        if (hasFlagSecondary(record))
            continue;

        // Skip non-aligning records.
        if (hasFlagUnaligned(record))
            continue;

        // Check that this is not a duplicate non-secondary record.
        if (record.qName == previousQName && hasFlagFirst(record) == previousIsFirst)
        {
            std::cerr << "ERROR: Duplicate non-secondary record for " << record.qName << "\n";
            return 1;
        }

        BamTagsDict bamTags(record.tags);

        // Counter: Total reads, unaligned reads.
        if (record.rId == BamAlignmentRecord::INVALID_REFID)
        {
            totalUnalignedReadCount += 1;
            continue;
        }
        totalReadCount += 1;

        // Get tag with edit distance, must be present for aligned reads.
        unsigned idx = 0;
        if (!findTagKey(idx, bamTags, "NM"))
        {
            std::cerr << "ERROR: Could not find NM tag!\n";
            return 1;
        }
        int editDistance = 0;
        if (!extractTagValue(editDistance, bamTags, idx))
        {
            std::cerr << "ERROR: Could not cast NM tag to int!\n";
            return 1;
        }

        // Count: Reads with errors.
        totalErrorneousReadCount += (editDistance != 0);

        // Count: Bases.
        totalBaseCount += length(record.seq);
        totalErrorCount += editDistance;

        // Register with histogram.
        histo[editDistance] += 1;

        // Update previous QNAME and is-first flag.
        previousQName = record.qName;
        previousIsFirst = hasFlagFirst(record);
    }

    // Print results.

    // The quick stats are what we need for the table in the paper.
    std::cout << "QUICK STATS\n"
              << "by base %\tby read %\n";
    fprintf(stdout, "%5.2f\t\t%5.2f\n\n", 100.0 * totalErrorCount / totalBaseCount, 100.0 * totalErrorneousReadCount / totalReadCount);

    // Print detailed statistics and histogram.
    std::cout << "STATISTICS\n"
              << "total read count      " << totalReadCount << "\t\t(excludes unaligned reads)\n"
              << "unaligned read count  " << totalUnalignedReadCount << "\n"
              << "erroneous read count " << totalErrorneousReadCount << "\n"
              << "per read error rate   " << 100.0 * totalErrorneousReadCount / totalReadCount << "\n"
              << "\n"
              << "total bases           " << totalBaseCount << "\n"
              << "total errors          " << totalErrorCount << "\n"
              << "per base error rate   " << 100.0 * totalErrorCount / totalBaseCount << "\n"
              << "\n"
              << "HISTOGRAM\n"
              << "errors\tcount\t\tpercent\n";
    for (std::map<unsigned, unsigned>::const_iterator it = histo.begin(); it != histo.end(); ++it)
    {
        fprintf(stdout, "%3d\t%12u\t%5.2f\n", it->first, it->second, 100.0 * it->second / totalReadCount);
    }

    return 0;
}
