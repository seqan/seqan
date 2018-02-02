// ==========================================================================
//                   NGS: Regions of Interest Analysis
// ==========================================================================
// Copyright (c) 2012-2018, Bernd Jagla, Institut Pasteur
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
    // The configuration for the application.
    RoiIntersectOptions const & options;

    // Stream for reading GFF and BED files.
    std::ifstream inIntervals;
    // ROI input and output.
    seqan::RoiFileIn  roiFileIn;
    seqan::RoiFileOut roiFileOut;

    RoiIntersectApp(RoiIntersectOptions const & options) : options(options)
    {}

    // App object's main routine.  Contains mostly console I/O.
    void run();

    // Called by run() for the actual work step.
    void doStreaming();
};

// --------------------------------------------------------------------------
// Class IntersectDriver
// --------------------------------------------------------------------------

// Helper function: Make the given record (must have int rId and int beginPos members) a sentinel (> all others).

template <typename TRecord>
void makeSentinel(TRecord & record)
{
    clear(record.ref);
    record.beginPos = std::numeric_limits<int>::max();
}

// Configuration object for directly reading BED records with IntersectDriver.

struct IntersectWithBedConfig
{
    // BedFileIn to use for reading.
    seqan::BedFileIn bedFileIn; 

    IntersectWithBedConfig(std::ifstream & inStream,
                           RoiIntersectOptions const & /*options*/) :
        bedFileIn(inStream)
    {}

    bool atEnd()
    {
        return seqan::atEnd(bedFileIn);
    }

    void skipComments()
    {
        while (!seqan::atEnd(bedFileIn) && *bedFileIn.iter == '#')
            skipLine(bedFileIn.iter);
    }

    void readRecord(seqan::BedRecord<seqan::Bed6> & record)
    {
        if (!seqan::atEnd(bedFileIn))
            seqan::readRecord(record, bedFileIn);
        else
            makeSentinel(record);
    }
};

// Configuration object for directly reading GFF records with IntersectDriver with implicit conversion to BED.

struct IntersectWithGffConfig
{
    // GffFileIn to use for reading.
    seqan::GffFileIn gffFileIn;

    // running number of read GFF records.
    int gffID;
    // The GFF record type to filter to.
    seqan::CharString gffType;

    // Used for building BED record names.
    std::stringstream ss;

    // Buffer for the current GFF record.
    seqan::GffRecord gffRecord;

    IntersectWithGffConfig(std::ifstream & inStream,
                           RoiIntersectOptions const & options) :
        gffFileIn(inStream), gffID(0), gffType(options.gffType)
    {}

    bool atEnd()
    {
        return seqan::atEnd(gffFileIn);
    }

    void skipComments()
    {
        while (!seqan::atEnd(gffFileIn) && *gffFileIn.iter == '#')
            skipLine(gffFileIn.iter);
    }

    void readRecord(seqan::BedRecord<seqan::Bed6> & bedRecord)
    {
        // Read GFF record.  When GFF record type filtering is active then we skip over records with the incorrect type.
        // We make bedRecord a sentinel if there is no GFF record with a valid type left.
        while (!seqan::atEnd(gffFileIn))
        {
            seqan::readRecord(gffRecord, gffFileIn);

            if (!empty(gffType) && gffRecord.type != gffType)
                continue;  // Read next.

            // Convert GFF to BED.
            clear(bedRecord);
            bedRecord.ref  = gffRecord.ref;
            bedRecord.beginPos = gffRecord.beginPos;
            bedRecord.endPos = gffRecord.endPos;
            bedRecord.score = "0";
            bedRecord.strand = (gffRecord.strand == '-') ? '-' : '+';  // '.' becomes '+'

            // Build BED record name.  We cannot rely on the GFF/GTF record having an ID so we simply construct one.
            ss.str("");
            ss.clear();
            ss << "ggf_" << (gffID++) << "_" << bedRecord.ref << ":" << bedRecord.beginPos << "-" << bedRecord.endPos;
            bedRecord.name = ss.str();

            return;
        }

        makeSentinel(bedRecord);
    }
};

// This class template is used for the generic streaming of ROIs against BED records or GFF records being converted into
// BED records.  The names for the I/O contexts and record readers are prefixed for BED even if we actually read from
// GFF and convert to BED afterwards.

template <typename TRecord>
inline bool isSentinel(TRecord const & record)
{
    return empty(record.ref);
}

template <typename TRecordL, typename TRecordR>
inline bool ltRecord(TRecordL const & recordL, TRecordR const & recordR)
{
    if (isSentinel(recordR))
        return true;
    if (isSentinel(recordL))
        return false;
    return ((recordL.ref < recordR.ref) ||
            (recordL.ref == recordR.ref && recordL.beginPos < recordR.beginPos));
}

template <typename TConfig>
class IntersectDriver
{
public:
    // TConfig defines some types and we have an instance for reading BED records and converting GFF records to BED
    // records when reading.
    TConfig config;

    // Files for reading/writing ROI.
    seqan::RoiFileOut & roiFileOut;
    seqan::RoiFileIn & roiFileIn;

    // BED and ROI records.
    seqan::BedRecord<seqan::Bed6> bedRecord;
    seqan::RoiRecord roiRecord;

    // Options for intersecting.
    RoiIntersectOptions const & options;

    IntersectDriver(seqan::RoiFileOut & roiFileOut,
                    std::ifstream & inStream,
                    seqan::RoiFileIn & roiFileIn,
                    RoiIntersectOptions const & options) :
        config(inStream, options), roiFileOut(roiFileOut),
        roiFileIn(roiFileIn), options(options)
    {}

    void run()
    {
        // Read and write header.
        seqan::RoiHeader roiHeader;
        readHeader(roiHeader, roiFileIn);
        clear(roiHeader.extraColumns);  // extra data is removed in projection
        writeHeader(roiFileOut, roiHeader);

        // Skip comments, then read first record through driver.
        config.skipComments();
        config.readRecord(bedRecord);

        // Skip comments, then read first record through ROI file.
        readRoiRecordOrMakeSentinel(roiRecord);

        // Check for both files being empty / returning sentinels.
        if (isSentinel(bedRecord) && isSentinel(roiRecord))
            return;  // both empty, done already

        // Setup the algorithm objects.
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
        IntersectBed intersectBedF(roiFileOut, intersectOptions);
        IntersectBed intersectBedR(roiFileOut, intersectOptions);

        // Stream over all BED and ROI records.
        while (!isSentinel(bedRecord) || !isSentinel(roiRecord))
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
            if (ltRecord(bedRecord, roiRecord))
            {
                if (!options.strandSpecific || bedRecord.strand == '+')
                    intersectBedF.pushBed(bedRecord);
                else
                    intersectBedR.pushBed(bedRecord);
                clear(bedRecord);

                // Read next records
                config.readRecord(bedRecord);
                if (isSentinel(bedRecord) && options.verbosity >= 2)
                    std::cerr << "BED record is a sentinel now!\n";
            }
            else
            {
                if (!options.strandSpecific || roiRecord.strand == '+')
                    intersectBedF.pushRoi(roiRecord);
                else
                    intersectBedR.pushRoi(roiRecord);
                clear(roiRecord);

                // Read next record.
                readRoiRecordOrMakeSentinel(roiRecord);
                if (isSentinel(roiRecord) && options.verbosity >= 2)
                    std::cerr << "ROI record is a sentinel now!\n";
            }
        }
    }

    void readRoiRecordOrMakeSentinel(seqan::RoiRecord & roiRecord)
    {
        if (!atEnd(roiFileIn))
            readRecord(roiRecord, roiFileIn);
        else
            makeSentinel(roiRecord);
    }
};

// --------------------------------------------------------------------------
// Class GroupByDriver
// --------------------------------------------------------------------------

// Code for the intersection with grouping of GFF records.

class GroupByDriver
{
public:
    // File stream for writing ROI to.
    seqan::RoiFileOut & roiFileOut;

    // Record readers for reading the files.
    seqan::GffFileIn gffFileIn;
    seqan::RoiFileIn & roiFileIn;

    // BED and ROI records.
    seqan::GffRecord gffRecord;
    seqan::RoiRecord roiRecord;

    // Options for intersecting.
    RoiIntersectOptions const & options;

    GroupByDriver(seqan::RoiFileOut & roiFileOut, std::ifstream & gffStream,
                  seqan::RoiFileIn & roiFileIn, RoiIntersectOptions const & options) :
            roiFileOut(roiFileOut), gffFileIn(gffStream), roiFileIn(roiFileIn), options(options)
    {}

    // Actually running the intersection with grouping.
    void run()
    {
        // The current reference name.
        seqan::CharString ref;

        // Read Roi header and write out again.
        seqan::RoiHeader roiHeader;
        readHeader(roiHeader, roiFileIn);
        writeHeader(roiFileOut, roiHeader);

        // Initialize objects that we will use for the overlapping and output generation.
        ProjectSplicedRoi workerF(roiFileOut, options.gffGroupBy, options.verbosity);
        ProjectSplicedRoi workerR(roiFileOut, options.gffGroupBy, options.verbosity);

        // Read first records.
        initializeRecords();
        if (options.verbosity >= 2)
        {
            std::cerr << "FIRST GFF RECORD REF = " << ref << "\n  ";
            writeRecord(std::cerr, gffRecord, seqan::Gff());
        }

        typedef seqan::Position<seqan::GffFileIn>::Type TGffFileInPos;
        TGffFileInPos chromBegin = position(gffFileIn);

        while (isSentinel(gffRecord))
        {
            // First pass: Register all GFF records with worker and reset position to chromosome begin afterwards.
            if (isSentinel(roiRecord))
                ref = std::min(gffRecord.ref, roiRecord.ref);
            else
                ref = gffRecord.ref;  // case where ROI is already at end
            if (options.verbosity >= 2)
                std::cerr << "FIRST PASS REF=" << ref << "\n";
            workerF.beginContig();
            workerR.beginContig();
            chromBegin = position(gffFileIn);
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

                readGffRecordOrMakeSentinel(gffRecord);
                if (isSentinel(gffRecord) && options.verbosity >= 2)
                    std::cerr << "GFF record is a sentinel now.\n";
            }
            if (!seqan::atEnd(gffFileIn))
                SEQAN_ASSERT_GT(gffRecord.ref, ref);
            setPosition(gffFileIn, chromBegin);
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
                if (isSentinel(roiRecord) ||
                    std::make_pair(gffRecord.ref, (int)gffRecord.beginPos) <
                    std::make_pair(roiRecord.ref, roiRecord.beginPos))
                {
                    if (isSentinel(gffRecord))
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
                    std::pair<seqan::CharString, uint32_t> oldPos(gffRecord.ref, gffRecord.beginPos);

                    // Read next records
                    readGffRecordOrMakeSentinel(gffRecord);
                    if (isSentinel(gffRecord) && options.verbosity >= 3)
                        std::cerr << "GFF record is a sentinel now!\n";
                    if (std::make_pair(gffRecord.ref, gffRecord.beginPos) < oldPos)
                        throw std::runtime_error("ERROR: GFF file is not sorted properly!");
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
                    readRoiRecordOrMakeSentinel(roiRecord);
                    if (isSentinel(roiRecord) &&options.verbosity >= 2)
                        std::cerr << "ROI record is a sentinel now!\n";
                    if (std::make_pair(roiRecord.ref, roiRecord.beginPos) < oldPos)
                        throw std::runtime_error("ERROR: ROI file is not sorted properly!");
                }
            }
        }

        // Finish reading ROI file to check for sortedness.
        while (!seqan::atEnd(roiFileIn))
        {
            std::pair<seqan::CharString, int> oldPos(roiRecord.ref, roiRecord.beginPos);

            readRecord(roiRecord, roiFileIn);

            if (std::make_pair(roiRecord.ref, roiRecord.beginPos) < oldPos)
                throw std::runtime_error("ERROR: ROI file is not sorted properly!\n");
        }
    }

    void readRoiRecordOrMakeSentinel(seqan::RoiRecord & record)
    {
        if (seqan::atEnd(roiFileIn))
            makeSentinel(record);
        else
            readRecord(record, roiFileIn);
    }

    void readGffRecordOrMakeSentinel(seqan::GffRecord & record)
    {
        if (seqan::atEnd(gffFileIn))
            makeSentinel(record);
        else
            readRecord(record, gffFileIn);
    }

    void initializeRecords()
    {
        // TODO(holtgrew): Check for sortedness here as well.

        // Read first records.
        while (!seqan::atEnd(gffFileIn) && *gffFileIn.iter == '#')
            skipLine(gffFileIn.iter);
        if (atEnd(gffFileIn))
        {
            makeSentinel(gffRecord);
        }
        else
        {
            bool found = false;
            while (!seqan::atEnd(gffFileIn))
            {
                readRecord(gffRecord, gffFileIn);
                if (empty(options.gffType) || (options.gffType == gffRecord.type))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
                makeSentinel(gffRecord);
        }
        while (!seqan::atEnd(roiFileIn) && *roiFileIn.iter == '#')
            skipLine(roiFileIn.iter);
        if (atEnd(roiFileIn))
            makeSentinel(roiRecord);
        else
            readRecord(roiRecord, roiFileIn);
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
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

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

    addOption(parser, seqan::ArgParseOption("ir", "in-roi", "ROI file to read.", seqan::ArgParseOption::INPUT_FILE, "ROI"));
    setRequired(parser, "in-roi");
    setValidValues(parser, "in-roi", seqan::RoiFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("if", "in-features", "BED, GFF, or GTF file to read.", seqan::ArgParseOption::INPUT_FILE, "FILE"));
    setRequired(parser, "in-features");
    std::vector<std::string> extensions = seqan::BedFileIn::getFileExtensions();
    std::vector<std::string> extensionsGff = seqan::GffFileIn::getFileExtensions();
    extensions.insert(extensions.end(), extensionsGff.begin(), extensionsGff.end());
    setValidValues(parser, "in-features", extensions);

    addOption(parser, seqan::ArgParseOption("or", "out-roi", "ROI file to write.", seqan::ArgParseOption::OUTPUT_FILE, "ROI"));
    setRequired(parser, "out-roi");
    setValidValues(parser, "out-roi", "roi");

    addOption(parser, seqan::ArgParseOption("g", "genome",
                                            "Path to FASTA file with genome; optional.  When given, this is used for "
                                            "computing the overall region's C+G content.",
                                            seqan::ArgParseOption::INPUT_FILE, "FASTA"));
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

void RoiIntersectApp::doStreaming()
{
    if (empty(options.gffGroupBy))
    {
        // Grouping is not ative, project against BED (or GFF/GTF reduced to BED information).
        if (options.inputIntervalsFileExt == "bed" || options.inputIntervalsFileExt == ".bed")
        {
            IntersectDriver<IntersectWithBedConfig> driver(roiFileOut, inIntervals, roiFileIn, options);
            driver.run();
        }
        else
        {
            IntersectDriver<IntersectWithGffConfig> driver(roiFileOut, inIntervals, roiFileIn, options);
            driver.run();
        }
    }
    else
    {
        // Group By is active, project with grouping.
        GroupByDriver driver(roiFileOut, inIntervals, roiFileIn, options);
        driver.run();
    }
}

// --------------------------------------------------------------------------
// Member Function RoiIntersectApp::run()
// --------------------------------------------------------------------------

void RoiIntersectApp::run()
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
        throw std::runtime_error("ERROR: Could not open file!");
    if (options.verbosity >= 1)
        std::cerr << " OK\n";

    if (options.verbosity >= 1)
        std::cerr << "INPUT ROI \t" << options.inputRoiFile << " ...";
    if (!open(roiFileIn, toCString(options.inputRoiFile)))
        throw std::runtime_error("ERROR: Could not open file!");
    if (options.verbosity >= 1)
        std::cerr << " OK\n";

    if (options.verbosity >= 1)
        std::cerr << "OUTPUT ROI\t" << options.outputRoiFile << " ...";
    if (!open(roiFileOut, toCString(options.outputRoiFile)))
        throw std::runtime_error("\nERROR: Could not open file!");
    if (options.verbosity >= 1)
        std::cerr << " OK\n\n";

    if (options.verbosity >= 1)
        std::cerr << "__PROCESSING FILES___________________________________________________________\n"
                  << "\n"
                  << "STREAMING ...";
    doStreaming();

    if (options.verbosity >= 1)
        std::cerr << " OK\n\n";

    if (options.verbosity >= 1)
        std::cerr << "DONE.\n";
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
    try
    {
        app.run();
    }
    catch (std::runtime_error const & err)
    {
        return 1;
    }
    return 0;
}
