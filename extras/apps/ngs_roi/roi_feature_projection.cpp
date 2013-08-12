// ==========================================================================
//                   NGS: Regions of Interest Analysis
// ==========================================================================
// Copyright (c) 2012-2013, Bernd Jagla, Institut Pasteur
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Intersection of ROI file with BED and/or GFF files.
// ==========================================================================

// TODO(holtgrew): Actually, union is projection with extension and could be renamed.
// TODO(holtgrew): More stringent checking, i.e. make sure that input is GFF/GTF if grouping/filtering is active.
// TODO(holtgrew): Rename to roi_overlapper?

#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>

#include "project_interval.h"
#include "project_spliced.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class RoiIntersectOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct RoiIntersectOptions
{
    // The combination mode.
    enum CombinationMode
    {
        PROJECTION,
        INTERSECTION,
        UNION,
        DIFF
    };

    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Path to ROI file to read.
    seqan::CharString inputRoiFile;
    // Path to intervals file to read.
    seqan::CharString inputIntervalsFile;
    seqan::CharString inputIntervalsFileExt;  // extension for type

    // Output ROI file.
    seqan::CharString outputRoiFile;

    // The GFF/GTF record type to filter for.  Empty for no filter.
    seqan::CharString gffType;
    // The GFF/GTF key/tag name to use for grouping.  Empty for no grouping.
    seqan::CharString gffGroupBy;

    // The mode to use for the combination.
    CombinationMode mode;

    // Whether or not to run in strand-specific mode.
    bool strandSpecific;

    RoiIntersectOptions() : verbosity(1), mode(PROJECTION), strandSpecific(false)
    {}
};

// --------------------------------------------------------------------------
// Class RoiIntersectApp
// --------------------------------------------------------------------------

// Class RoiIntersectApp

class RoiIntersectApp
{
public:
    typedef seqan::RecordReader<std::ifstream, seqan::SinglePass<> > TRecordReader;

    // The configuration for the application.
    RoiIntersectOptions const & options;

    // Input file streams.
    std::ifstream inIntervals, inRoi;
    std::ofstream outRoi;

    RoiIntersectApp(RoiIntersectOptions const & options) : options(options)
    {}

    // App object's main routine.  Contains mostly console I/O.
    int run();

    // Called by run() for the actual work step.
    int doStreaming();
};

// --------------------------------------------------------------------------
// Class IntersectDriver
// --------------------------------------------------------------------------

// Helper function: Make the given record (must have int rId and int beginPos members) a sentinel (> all others).
template <typename TRecord>
void makeSentinel(TRecord & record)
{
    clear(record.ref);
    record.rID = seqan::maxValue<int>();
    record.beginPos = seqan::maxValue<int>();
}

void makeSentinel(seqan::GffRecord & record)
{
    clear(record.ref);
    record.rID = seqan::maxValue<int>();
    record.beginPos = seqan::maxValue<int>();
}

// Configuration object for directly reading BED records with IntersectDriver.

struct IntersectWithBedConfig
{
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::BedIOContext<TNameStore>     TBedIOContext;

    IntersectWithBedConfig(RoiIntersectOptions const & /*options*/)
    {}

    template <typename TStream, typename TReaderSpec>
    int readRecord(seqan::BedRecord<seqan::Bed6> & record,
                   seqan::RecordReader<TStream, TReaderSpec> & reader,
                   TBedIOContext & bedIOContext)
    {
        return seqan::readRecord(record, reader, bedIOContext, seqan::Bed());
    }
};

// Configuration object for directly reading GFF records with IntersectDriver with implicit conversion to BED.

struct IntersectWithGffConfig
{
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::GffIOContext<TNameStore>     TBedIOContext;

    // running number of read GFF records.
    int gffID;
    // The GFF record type to filter to.
    seqan::CharString gffType;

    // Used for building BED record names.
    std::stringstream ss;

    // Buffer for the current GFF record.
    seqan::GffRecord gffRecord;

    IntersectWithGffConfig(RoiIntersectOptions const & options) : gffID(0), gffType(options.gffType)
    {}

    template <typename TStream, typename TReaderSpec>
    int readRecord(seqan::BedRecord<seqan::Bed6> & bedRecord,
                   seqan::RecordReader<TStream, TReaderSpec> & reader,
                   TBedIOContext & gffIOContext)
    {
        // Read GFF record.  When GFF record type filtering is active then we skip over records with the incorrect type.
        // We make bedRecord a sentinel if there is no GFF record with a valid type left.
        while (!atEnd(reader))
        {
            int res = 0;
            if ((res = seqan::readRecord(gffRecord, reader, gffIOContext, seqan::Gff())) != 0)
                return res;

            if (!empty(gffType) && gffRecord.type != gffType)
                continue;  // Read next.

            // Convert GFF to BED.
            clear(bedRecord);
            bedRecord.ref  = gffRecord.ref;
            bedRecord.rID = gffRecord.rID;
            bedRecord.beginPos = gffRecord.beginPos;
            bedRecord.endPos = gffRecord.endPos;
            bedRecord.score = 0;
            bedRecord.strand = (gffRecord.strand == '-') ? '-' : '+';  // '.' becomes '+'

            // Build BED record name.  We cannot rely on the GFF/GTF record having an ID so we simply construct one.
            ss.str("");
            ss.clear();
            ss << "ggf_" << (gffID++) << "_" << bedRecord.ref << ":" << bedRecord.beginPos << "-" << bedRecord.endPos;
            bedRecord.name = ss.str();

            return 0;
        }

        makeSentinel(bedRecord);
        return 0;
    }
};

// This class template is used for the generic streaming of ROIs against BED records or GFF records being converted into
// BED records.  The names for the I/O contexts and record readers are prefixed for BED even if we actually read from
// GFF and convert to BED afterwards.

template <typename TConfig>
class IntersectDriver
{
public:
    // Names store for the reference names.
    typedef seqan::StringSet<seqan::CharString> TNameStore;

    // TConfig defines some types and we have an instance for reading BED records and converting GFF records to BED
    // records when reading.
    TConfig config;

    // The reference name store and a cache for this store.
    TNameStore refNames;
    seqan::NameStoreCache<TNameStore> refNamesCache;

    // File stream for writing ROI to.
    std::ofstream & outRoi;

    // I/O contexts for translating reference names to ids.
    typename TConfig::TBedIOContext bedIOContext;
    seqan::RoiIOContext<TNameStore> roiIOContext;

    // Record readers for reading the files.
    seqan::RecordReader<std::ifstream, seqan::SinglePass<> > bedReader;
    seqan::RecordReader<std::ifstream, seqan::SinglePass<> > roiReader;

    // BED and ROI records.
    seqan::BedRecord<seqan::Bed6> bedRecord;
    seqan::RoiRecord roiRecord;

    // Options for intersecting.
    RoiIntersectOptions const & options;

    IntersectDriver(std::ofstream & outRoi, std::ifstream & bedStream, std::ifstream & roiStream,
                    RoiIntersectOptions const & options) :
            config(options), refNamesCache(refNames), outRoi(outRoi), bedIOContext(refNames, refNamesCache),
            roiIOContext(refNames, refNamesCache), bedReader(bedStream), roiReader(roiStream),
            options(options)
    {}

    int run()
    {
        // TODO(holtgrew): What happens if there are no records for one contig?

        // Write header.
        outRoi << "##ref\t"
               << "begin_pos\t"
               << "end_pos\t"
               << "region_name\t"
               << "length\t"
               << "strand\t"
               << "max_count\t"
               << "counts\n";

        // Read first records.
        while (!atEnd(bedReader) && value(bedReader) == '#')
            if (skipLine(bedReader) != 0)
            {
                std::cerr << "ERROR: Could not skip header/comment line in BED.\n";
                return 1;
            }
        if (atEnd(bedReader))
        {
            makeSentinel(bedRecord);
        }
        else if (config.readRecord(bedRecord, bedReader, bedIOContext) != 0)
        {
            std::cerr << "ERROR: Problem reading from BED file!\n";
            return 1;
        }
        while (!atEnd(roiReader) && value(roiReader) == '#')
            if (skipLine(roiReader) != 0)
            {
                std::cerr << "ERROR: Could not skip header/comment line in ROI.\n";
                return 1;
            }
        if (atEnd(roiReader))
        {
            makeSentinel(roiRecord);
        }
        else if (readRecord(roiRecord, roiReader, roiIOContext, seqan::Roi()) != 0)
        {
            std::cerr << "ERROR: Problem reading from ROI file!\n";
            return 1;
        }

        // The algorithm objects.
        IntersectBedOptions intersectOptions;
        switch (options.mode)
        {
            case RoiIntersectOptions::PROJECTION:
                intersectOptions.mode = IntersectBedOptions::PROJECTION;
                break;
            case RoiIntersectOptions::INTERSECTION:
                intersectOptions.mode = IntersectBedOptions::INTERSECTION;
                break;
            case RoiIntersectOptions::UNION:
                intersectOptions.mode = IntersectBedOptions::UNION;
                break;
            case RoiIntersectOptions::DIFF:
                intersectOptions.mode = IntersectBedOptions::DIFF;
                break;
            default:
                SEQAN_FAIL("Cannot reach here!");
        }
        intersectOptions.verbosity = options.verbosity;
        // We create two IntersectBed objects, one for each forward and reverse strand.  The forward IntersectBed object is
        // also used in non-strand-specific mode.
        IntersectBed intersectBedF(outRoi, intersectOptions);
        IntersectBed intersectBedR(outRoi, intersectOptions);

        // Stream over all BED and ROI records.
        while (bedRecord.rID != seqan::maxValue<int>() || roiRecord.rID != seqan::maxValue<int>())
        {
            if (options.verbosity >= 3)
            {
                std::cerr << ",--\n"
                          << "| ROI:\t";
                writeRecord(std::cerr, roiRecord, seqan::Roi());
                std::cerr << "| BED:\t";
                writeRecord(std::cerr, bedRecord, seqan::Bed());
                std::cerr << "`--\n";
            }

            // Push smaller, prefering ROI over BED on ties.
            if (bedRecord.rID < roiRecord.rID || (bedRecord.rID == roiRecord.rID && bedRecord.beginPos < roiRecord.beginPos))
            {
                if (!options.strandSpecific || bedRecord.strand == '+')
                    intersectBedF.pushBed(bedRecord);
                else
                    intersectBedR.pushBed(bedRecord);
                clear(bedRecord);

                // Read next records
                if (atEnd(bedReader))
                {
                    makeSentinel(bedRecord);
                    if (options.verbosity >= 2)
                        std::cerr << "BED record is a sentinel now!\n";
                }
                else if (config.readRecord(bedRecord, bedReader, bedIOContext) != 0)
                {
                    std::cerr << "ERROR: Problem reading from BED file!\n";
                    return 1;
                }
            }
            else
            {
                if (!options.strandSpecific || roiRecord.strand == '+')
                    intersectBedF.pushRoi(roiRecord);
                else
                    intersectBedR.pushRoi(roiRecord);
                clear(roiRecord);

                // Read next record.
                if (atEnd(roiReader))
                {
                    makeSentinel(roiRecord);
                    if (options.verbosity >= 2)
                        std::cerr << "ROI record is a sentinel now!\n";
                }
                else if (readRecord(roiRecord, roiReader, roiIOContext, seqan::Roi()) != 0)
                {
                    std::cerr << "ERROR: Problem reading from ROI file!\n";
                    return 1;
                }
            }
        }

        return 0;
    }
};

// --------------------------------------------------------------------------
// Class GroupByDriver
// --------------------------------------------------------------------------

// Helper function.

seqan::Pair<int, int> position(seqan::GffRecord const & record)
{
    return seqan::Pair<int, int>(record.rID, record.beginPos);
}

seqan::Pair<int, int> position(seqan::RoiRecord const & record)
{
    return seqan::Pair<int, int>(record.rID, record.beginPos);
}

// Code for the intersection with grouping of GFF records.

class GroupByDriver
{
public:
    // Names store for the reference names.
    typedef seqan::StringSet<seqan::CharString> TNameStore;

    // The reference name store and a cache for this store.
    TNameStore refNames;
    seqan::NameStoreCache<TNameStore> refNamesCache;

    // File stream for writing ROI to.
    std::ofstream & outRoi;

    // I/O contexts for translating reference names to ids.
    seqan::GffIOContext<TNameStore> gffIOContext;
    seqan::RoiIOContext<TNameStore> roiIOContext;

    // Record readers for reading the files.
    seqan::RecordReader<std::ifstream, seqan::SinglePass<> > gffReader;
    seqan::RecordReader<std::ifstream, seqan::SinglePass<> > roiReader;

    // BED and ROI records.
    seqan::GffRecord gffRecord;
    seqan::RoiRecord roiRecord;

    // Options for intersecting.
    RoiIntersectOptions const & options;

    GroupByDriver(std::ofstream & outRoi, std::ifstream & bedStream, std::ifstream & roiStream,
                  RoiIntersectOptions const & options) :
            refNamesCache(refNames), outRoi(outRoi), gffIOContext(refNames, refNamesCache),
            roiIOContext(refNames, refNamesCache), gffReader(bedStream), roiReader(roiStream),
            options(options)
    {}

    // Actually running the intersection with grouping.
    int run()
    {
        // The current reference name.
        seqan::CharString ref;

        // Write header.
        outRoi << "##ref\t"
               << "begin_pos\t"
               << "end_pos\t"
               << "region_name\t"
               << "length\t"
               << "strand\t"
               << "max_count\t"
               << "counts\n";

        // TODO(holtgrew): What happens if there are no records for one contig?

        // Initialize objects that we will use for the overlapping and output generation.
        ProjectSplicedRoi workerF(outRoi, options.gffGroupBy, options.verbosity);
        ProjectSplicedRoi workerR(outRoi, options.gffGroupBy, options.verbosity);

        // Read first records.
        if (initializeRecords() != 0)
            return 1;
        if (options.verbosity >= 2)
        {
            std::cerr << "FIRST GFF RECORD REF = " << ref << "\n  ";
            writeRecord(std::cerr, gffRecord, seqan::Gff());
        }

        typedef seqan::RecordReader<std::ifstream, seqan::SinglePass<> > TGffReader;
        typedef seqan::Position<TGffReader>::Type TGffReaderPos;
        TGffReaderPos chromBegin = position(gffReader);

        while (gffRecord.rID != seqan::maxValue<int>())
        {
            // First pass: Register all GFF records with worker and reset position to chromosome begin afterwards.
            if (roiRecord.rID != seqan::maxValue<int>())
                ref = std::min(gffRecord.ref, roiRecord.ref);
            else
                ref = gffRecord.ref;  // case where ROI is already at end
            if (options.verbosity >= 2)
                std::cerr << "FIRST PASS REF=" << ref << "\n";
            workerF.beginContig();
            workerR.beginContig();
            chromBegin = position(gffReader);
            seqan::GffRecord firstGffRecord = gffRecord;
            while (gffRecord.ref == ref)
            {
                if (!empty(options.gffType) && gffRecord.type == options.gffType)
                {
                    if (!options.strandSpecific || gffRecord.strand == '+')
                        workerF.updateRanges(gffRecord);
                    else
                        workerR.updateRanges(gffRecord);
                }

                if (atEnd(gffReader))
                {
                    makeSentinel(gffRecord);
                    if (options.verbosity >= 2)
                        std::cerr << "GFF record is a sentinel now.\n";
                }
                else if (readRecord(gffRecord, gffReader, gffIOContext, seqan::Gff()) != 0)
                {
                    std::cerr << "ERROR: Problem reading GFF.\n";
                    return 1;
                }
            }
            if (!atEnd(gffReader))
                SEQAN_ASSERT_GT(gffRecord.ref, ref);
            if (setPosition(gffReader, chromBegin) != 0)
            {
                std::cerr << "ERROR: Could not reset file pointer for second pass!\n";
                return 1;
            }
            gffRecord = firstGffRecord;

            if (options.verbosity >= 2)
                std::cerr << "SECOND PASS REF=" << ref << "\n";
            workerF.beginSecondPass();
            workerR.beginSecondPass();
            while (gffRecord.ref == ref || roiRecord.ref == ref)
            {
                if (options.verbosity >= 3)
                {
                    std::cerr << ",--\n"
                              << "| ROI:\t";
                    // NOTE that the sentinel is given as a negative value because of an overflow.
                    writeRecord(std::cerr, roiRecord, seqan::Roi());
                    std::cerr << "| GFF:\t";
                    writeRecord(std::cerr, gffRecord, seqan::Gff());
                    std::cerr << "`--\n";
                }

                // Push smaller, prefering ROI over GFF on ties.
                if (roiRecord.rID == seqan::maxValue<int>() ||
                    std::make_pair(gffRecord.ref, (int)gffRecord.beginPos) <
                    std::make_pair(roiRecord.ref, roiRecord.beginPos))
                {
                    if (gffRecord.rID == seqan::maxValue<int>())
                        break;  // Break out of outer loop, termination criterion.

                    if (!empty(options.gffType) && gffRecord.type == options.gffType)
                    {
                        if (!options.strandSpecific || gffRecord.strand == '+')
                            workerF.pushGff(gffRecord);
                        else
                            workerR.pushGff(gffRecord);
                    }
                    clear(gffRecord);

                    // Get current position to check for sortedness of GFF.
                    std::pair<seqan::CharString, __uint32> oldPos(gffRecord.ref, gffRecord.beginPos);

                    // Read next records
                    if (atEnd(gffReader))
                    {
                        makeSentinel(gffRecord);
                        if (options.verbosity >= 3)
                            std::cerr << "GFF record is a sentinel now!\n";
                    }
                    else if (readRecord(gffRecord, gffReader, gffIOContext, seqan::Gff()) != 0)
                    {
                        std::cerr << "ERROR: Problem reading from GFF file!\n";
                        return 1;
                    }
                    if (std::make_pair(gffRecord.ref, gffRecord.beginPos) < oldPos)
                    {
                        std::cerr << "ERROR: GFF file is not sorted properly!\n";
                        return 1;
                    }
                }
                else
                {
                    if (!options.strandSpecific || roiRecord.strand == '+')
                        workerF.pushRoi(roiRecord);
                    else
                        workerR.pushRoi(roiRecord);
                    clear(roiRecord);

                    // Get current position to check for sortedness of GFF.
                    std::pair<seqan::CharString, int> oldPos(roiRecord.ref, roiRecord.beginPos);

                    // Read next record.
                    if (atEnd(roiReader))
                    {
                        makeSentinel(roiRecord);
                        if (options.verbosity >= 2)
                            std::cerr << "ROI record is a sentinel now!\n";
                    }
                    else if (readRecord(roiRecord, roiReader, roiIOContext, seqan::Roi()) != 0)
                    {
                        std::cerr << "ERROR: Problem reading from ROI file!\n";
                        return 1;
                    }
                    if (std::make_pair(roiRecord.ref, roiRecord.beginPos) < oldPos)
                    {
                        std::cerr << "ERROR: ROI file is not sorted properly!\n";
                        return 1;
                    }
                }
            }
        }

        // Finish reading ROI file to check for sortedness.
        while (!atEnd(roiReader))
        {
            std::pair<seqan::CharString, int> oldPos(roiRecord.ref, roiRecord.beginPos);

            if (readRecord(roiRecord, roiReader, roiIOContext, seqan::Roi()) != 0)
            {
                std::cerr << "ERROR: Problem reading from ROI file!\n";
                return 1;
            }

            if (std::make_pair(roiRecord.ref, roiRecord.beginPos) < oldPos)
            {
                std::cerr << "ERROR: ROI file is not sorted properly!\n";
                return 1;
            }
        }

        return 0;
    }

    int initializeRecords()
    {
        // TODO(holtgrew): Check for sortedness here as well.

        // Read first records.
        while (!atEnd(gffReader) && value(gffReader) == '#')
            if (skipLine(gffReader) != 0)
            {
                std::cerr << "ERROR: Could not skip header/comment line in GFF.\n";
                return 1;
            }
        if (atEnd(gffReader))
        {
            makeSentinel(gffRecord);
        }
        else
        {
            bool found = false;
            while (!atEnd(gffReader))
            {
                if (readRecord(gffRecord, gffReader, gffIOContext, seqan::Gff()) != 0)
                {
                    std::cerr << "ERROR: Problem reading from GFF file!\n";
                    return 1;
                }
                if (empty(options.gffType) || (options.gffType == gffRecord.type))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
                makeSentinel(gffRecord);
        }
        while (!atEnd(roiReader) && value(roiReader) == '#')
            if (skipLine(roiReader) != 0)
            {
                std::cerr << "ERROR: Could not skip header/comment line in ROI.\n";
                return 1;
            }
        if (atEnd(roiReader))
        {
            makeSentinel(roiRecord);
        }
        else if (readRecord(roiRecord, roiReader, roiIOContext, seqan::Roi()) != 0)
        {
            std::cerr << "ERROR: Problem reading from ROI file!\n";
            return 1;
        }

        return 0;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function combinationModeStr()
// --------------------------------------------------------------------------

char const * combinationModeStr(RoiIntersectOptions::CombinationMode m)
{
    switch (m)
    {
        case RoiIntersectOptions::PROJECTION:
            return "projection";
        case RoiIntersectOptions::INTERSECTION:
            return "intersection";
        case RoiIntersectOptions::UNION:
            return "union";
        case RoiIntersectOptions::DIFF:
            return "difference";
        default:
            return "INVALID";
    }
}

// --------------------------------------------------------------------------
// Function print()
// --------------------------------------------------------------------------

void print(std::ostream & out, RoiIntersectOptions const & options)
{
    // Print the command line arguments back to the user.
    out << "__OPTIONS____________________________________________________________________\n"
        << '\n'
        << "VERBOSITY \t" << options.verbosity << '\n'
        << "\n"
        << "INPUT ROI       \t" << options.inputRoiFile << "\n"
        << "INPUT BED       \t" << options.inputIntervalsFile << "\n"
        << "OUTPUT ROI      \t" << options.outputRoiFile << "\n"
        << "\n"
        << "COMBINATION MODE\t" << combinationModeStr(options.mode) << "\n"
        << "\n"
        << "GFF TYPE        \t" << options.gffType << "\n"
        << "GFF GROUP BY    \t" << options.gffGroupBy << "\n"
        << "\n";
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(RoiIntersectOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("roi_feature_projection");
    setCategory(parser, "NGS ROI Analysis");

    // Set short description, version, and date.
    setShortDescription(parser, "Region Of Interest Projection.");
    setVersion(parser, "0.1");
    setDate(parser, "April 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-ir\\fP \\fIIN.roi\\fP \\fB-if\\fP \\fIIN.{bed,gff,gtf}\\fP \\fB-or\\fP \\fIOUT.roi\\fP");
    addDescription(parser,
                   "Compute the projection of a ROI file to regions from a BED or GFF file.  The result is "
                   "a ROI file where each interval from the BED/GFF/GTF file that overlapped with one input ROI file "
                   "is a region of interest, with the coverage counts projected to the new region of interest.");

    // -----------------------------------------------------------------------
    // General Options
    // -----------------------------------------------------------------------

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // -----------------------------------------------------------------------
    // Input / Output Options
    // -----------------------------------------------------------------------

    addSection(parser, "Input / Output");

    addOption(parser, seqan::ArgParseOption("ir", "in-roi", "ROI file to read.", seqan::ArgParseOption::INPUTFILE, "ROI"));
    setRequired(parser, "in-roi");
    setValidValues(parser, "in-roi", "roi");

    addOption(parser, seqan::ArgParseOption("if", "in-features", "BED, GFF, or GTF file to read.", seqan::ArgParseOption::INPUTFILE, "FILE"));
    setRequired(parser, "in-features");
    setValidValues(parser, "in-features", "bed gff gtf");

    addOption(parser, seqan::ArgParseOption("or", "out-roi", "ROI file to write.", seqan::ArgParseOption::OUTPUTFILE, "ROI"));
    setRequired(parser, "out-roi");
    setValidValues(parser, "out-roi", "roi");

    addOption(parser, seqan::ArgParseOption("g", "genome",
                                            "Path to FASTA file with genome; optional.  When given, this is used for "
                                            "computing the overall region's C+G content.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "genome", "fasta fa");

    // -----------------------------------------------------------------------
    // Combination Options
    // -----------------------------------------------------------------------

    addSection(parser, "Combination Options");

    addOption(parser, seqan::ArgParseOption("m", "mode",
                                            "The mode in which to combine the ROI and BED/GFF file.  See section "
                                            "Combination Modes below for details.", seqan::ArgParseOption::STRING,
                                            "MODE"));
    setDefaultValue(parser, "mode", "projection");
    setValidValues(parser, "mode", "intersection projection union difference");

    addOption(parser, seqan::ArgParseOption("ss", "strand-specific", "Enable strand-specific mode if set."));

    // -----------------------------------------------------------------------
    // GFF Filter Options
    // -----------------------------------------------------------------------

    addSection(parser, "GFF Filters");

    addOption(parser, seqan::ArgParseOption("", "gff-type",
                                            "The GFF/GTF record type (value of third column) to keep.  Keep all if "
                                            "not set or input file type is not GFF/GTF.",
                                            seqan::ArgParseOption::STRING, "STR"));

    addOption(parser, seqan::ArgParseOption("", "gff-group-by",
                                            "The GFF/GTF tag to use for grouping, e.g. \"Parent\", \"transcript_id\". No "
                                            "grouping if empty.  When using the grouping feature, the \\fB--mode\\fP is "
                                            "automatically set to \\fIprojection\\fP.", seqan::ArgParseOption::STRING, "STR"));

    // -----------------------------------------------------------------------
    // Examples
    // -----------------------------------------------------------------------

    addTextSection(parser, "Examples");

    addListItem(parser, "\\fBroi_intersect\\fP \\fB--in-features\\fP \\fIIN.bed\\fP \\fB--in-roi\\fP \\fIIN.roi\\fP \\fB--out-roi\\fP \\fIOUT.roi\\fP",
                "Project the data from IN.roi to the intervals from IN.bed and write out the result to OUT.roi.");

    addListItem(parser, "\\fBroi_intersect\\fP \\fB--mode\\fP \\fIdifference\\fP \\fB--in-features\\fP \\fIIN.bed\\fP \\fB--in-roi\\fP \\fIIN.roi\\fP \\fB--out-roi\\fP \\fIOUT.roi\\fP",
                "Compute symmetric difference of IN.bed and IN.roi and write to OUT.roi.");

    addListItem(parser, "\\fBroi_intersect\\fP \\fB--in-features\\fP \\fIIN.gff\\fP \\fB--gff-type\\fP \\fIexon\\fP \\fB--gff-group-by\\fP \\fIParent\\fP \\fB--in-roi\\fP \\fIIN.roi\\fP \\fB--out-roi\\fP \\fIOUT.roi\\fP",
                "Project the data from IN.ROI to the intervals of type exon in IN.gff.  The resulting projections to the"
                "exons with the same Parent are then joined.  Note that the exons have to be non-overlapping.");

    // -----------------------------------------------------------------------
    // Combination Modes Details
    // -----------------------------------------------------------------------

    addTextSection(parser, "Combination Modes");

    addText(parser,
            "The following combination modes are available.");

    addListItem(parser, "\\fBintersection\\fP",
                "Intersect the intervals  BED/GFF file with the ones from the ROI file.  Only the intersections are "
                "kept and the coverage for each position is taken from the ROI records' coverages.");

    addListItem(parser, "\\fBunion\\fP",
                "Compute union of BED and BED/GFF file intervals, take coverages from ROI records.  This leads "
                "to the joining of ROIs into larger ones.");

    addListItem(parser, "\\fBprojection\\fP",
                "Compute intersection of BED/GFF file intervals with the ROI intervals.  All BED/GFF intervals are "
                "kept and the coverages are taken from the ROI records.");

    addListItem(parser, "\\fBdifference\\fP",
                "Compute symmetric difference of intervals in BED/GFF and ROI file.  BED/GFF intervals without any "
                "overlap in the ROI file are written out as ROI records with max_count=0.  ROI intervals without "
                "overlap in BED/GFF file are written out as they are in the input.");

    // -----------------------------------------------------------------------
    // Parse and Extract Options
    // -----------------------------------------------------------------------

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.inputRoiFile, parser, "in-roi");
    getOptionValue(options.inputIntervalsFile, parser, "in-features");
    options.inputIntervalsFileExt = getOptionFileExtension(parser, "in-features");
    getOptionValue(options.outputRoiFile, parser, "out-roi");

    seqan::CharString tmp;
    getOptionValue(tmp, parser, "mode");
    if (tmp == "projection")
        options.mode = RoiIntersectOptions::PROJECTION;
    else if (tmp == "intersection")
        options.mode = RoiIntersectOptions::INTERSECTION;
    else if (tmp == "union")
        options.mode = RoiIntersectOptions::UNION;
    else if (tmp == "difference")
        options.mode = RoiIntersectOptions::DIFF;
    options.strandSpecific = isSet(parser, "strand-specific");

    getOptionValue(options.gffType, parser, "gff-type");
    getOptionValue(options.gffGroupBy, parser, "gff-group-by");
    if (!empty(options.gffGroupBy))
        options.mode = RoiIntersectOptions::PROJECTION;

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Member Function RoiIntersectApp::doStreaming()
// --------------------------------------------------------------------------

int RoiIntersectApp::doStreaming()
{
    if (empty(options.gffGroupBy))
    {
        // Grouping is not ative, project against BED (or GFF/GTF reduced to BED information).
        if (options.inputIntervalsFileExt == ".bed")
        {
            IntersectDriver<IntersectWithBedConfig> driver(outRoi, inIntervals, inRoi, options);
            return driver.run();
        }
        else
        {
            IntersectDriver<IntersectWithGffConfig> driver(outRoi, inIntervals, inRoi, options);
            return driver.run();
        }
    }
    else
    {
        // Group By is active, project with grouping.
        GroupByDriver driver(outRoi, inIntervals, inRoi, options);
        return driver.run();
    }
}

// --------------------------------------------------------------------------
// Member Function RoiIntersectApp::run()
// --------------------------------------------------------------------------

int RoiIntersectApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "ROI INTERSECT\n"
                  << "=============\n\n";
        print(std::cerr, options);
    }

    if (options.verbosity >= 1)
        std::cerr << "__OPENING FILES______________________________________________________________\n"
                  << "\n";

    if (options.verbosity >= 1)
        std::cerr << "INPUT INTERVALS \t" << options.inputIntervalsFile << " ...";
    inIntervals.open(toCString(options.inputIntervalsFile), std::ios::binary | std::ios::in);
    if (!inIntervals.good())
    {
        std::cerr << "\nERROR: Could not open file!\n";
        return 1;
    }
    if (options.verbosity >= 1)
        std::cerr << " OK\n";

    if (options.verbosity >= 1)
        std::cerr << "INPUT ROI \t" << options.inputRoiFile << " ...";
    inRoi.open(toCString(options.inputRoiFile), std::ios::binary | std::ios::in);
    if (!inRoi.good())
    {
        std::cerr << "\nERROR: Could not open file!\n";
        return 1;
    }
    if (options.verbosity >= 1)
        std::cerr << " OK\n";

    if (options.verbosity >= 1)
        std::cerr << "OUTPUT ROI\t" << options.outputRoiFile << " ...";
    outRoi.open(toCString(options.outputRoiFile), std::ios::binary | std::ios::out);
    if (!inRoi.good())
    {
        std::cerr << "\nERROR: Could not open file!\n";
        return 1;
    }
    if (options.verbosity >= 1)
        std::cerr << " OK\n\n";

    if (options.verbosity >= 1)
        std::cerr << "__PROCESSING FILES___________________________________________________________\n"
                  << "\n"
                  << "STREAMING ...";
    if (doStreaming() != 0)
        return 1;

    if (options.verbosity >= 1)
        std::cerr << " OK\n\n";

    if (options.verbosity >= 1)
        std::cerr << "DONE.\n";

    return 0;

}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    RoiIntersectOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    RoiIntersectApp app(options);
    return app.run();
}
