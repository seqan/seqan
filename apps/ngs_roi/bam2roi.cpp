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
// Author: Bernd Jagla <bernd.jagla@pasteur.fr>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================


#include <fstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

#include "roi_builder.h"

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

// Options for the bam2roi Application.

struct Options
{
    // Verbosity of output: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // -----------------------------------------------------------------------
    // Input / Output Options
    // -----------------------------------------------------------------------

    // Paths to input read files.
    seqan::CharString inputFileName;

    // Paths to output read files of accepted reads.
    seqan::CharString outputFileName;

    // -----------------------------------------------------------------------
    // ROI Creation Options
    // -----------------------------------------------------------------------

	// Whether or not the experiment is strand-specific.
	bool strandSpecific;

    // Whether or not to use paired information.
    bool usePairing;

    // Whether or not to link over skipped bases.
    bool linkOverSkipped;

    Options() : verbosity(0), strandSpecific(false), usePairing(false), linkOverSkipped(false)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function yesNo()                                                 [Options]
// --------------------------------------------------------------------------

char const * yesNo(bool b)
{
    return b ? "YES" : "NO";
}

// --------------------------------------------------------------------------
// Function print()                                                 [Options]
// --------------------------------------------------------------------------

void print(std::ostream & out, Options const & options)
{
    out << "__OPTIONS_____________________________________________________________________\n"
        << "\n"
        << "INPUT FILE       \t" << options.inputFileName << "\n"
        << "OUTPUT FILE      \t" << options.outputFileName << "\n"
        << "\n"
        << "STRAND SPECIFIC  \t" << yesNo(options.strandSpecific) << "\n"
        << "USE PAIRED INFO  \t" << yesNo(options.usePairing) << "\n"
        << "LINK OVER SKIPPED\t" << yesNo(options.linkOverSkipped) << "\n"
        << "\n";
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

// Parse command line with ArgumentParser class.

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bam2roi");
    setCategory(parser, "NGS ROI Analysis");

    // Set short description, version, and date.
    setShortDescription(parser, "Create ROI from BAM file.");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "\\fB-if\\fP \\fIIN.bam\\fP \\fB-of\\fP \\fIOUT.roi\\fP");

    addDescription(parser,
                   "Calculated consecutive regions of coverage from alignment file \\fIIN.bam\\fP "
                   "and write regions of interst to file \\fIOUT.roi\\fP. "
				   "Counting is performed over the entire region (including intron and N-regions) "
				   "based on the CIGAR string of the alignment record.");
	// it is in the scope of this application to make any assumtions over splice junctions
	// since there are many boundary issues concerning splice junctions of many reads we
	// consider this an alignment artefact. Splice junctions can be analyzed separately by
	// using special GFF annotations in conjunction with these roi files.

    // -----------------------------------------------------------------------
    // General Options
    // -----------------------------------------------------------------------

    addSection(parser, "General Options");

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose mode."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose mode."));

    // -----------------------------------------------------------------------
    // Input / Output Options
    // -----------------------------------------------------------------------

    addSection(parser, "Input / Output Parameters");

    addOption(parser, seqan::ArgParseOption("if", "input-file", "SAM/BAM formatted file.  Must be sorted by coordinate.",
                                            seqan::ArgParseOption::INPUT_FILE));
    setValidValues(parser, "input-file", seqan::BamFileIn::getFileExtensions());
    setRequired(parser, "input-file");

    addOption(parser, seqan::ArgParseOption("of", "output-file", "Output file with regions of interest.",
                                            seqan::ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-file", seqan::RoiFileIn::getFileExtensions());
    setRequired(parser, "output-file");

    // -----------------------------------------------------------------------
    // Optional Parameters
    // -----------------------------------------------------------------------

    addSection(parser, "ROI Construction Options");

	addOption(parser, seqan::ArgParseOption("ss", "strand-specific",
                                            "Calculate strand-specific ROIs (see section Strand Specificness below."));

	addOption(parser, seqan::ArgParseOption("ip", "ignore-pairing",
                                            "Ignore paired information.  Also see Section ROI Creation Details."));

	addOption(parser, seqan::ArgParseOption("ls", "link-over-skipped", "Link over skipped bases in the read alignment."));

    // -----------------------------------------------------------------------
    // Documentation of ROI Creation Details.
    // -----------------------------------------------------------------------

    addTextSection(parser, "ROI Creation Details");

    addText(parser,
            "By default, ROIs are created from single-end data by looking at the alignment only.  Matches, mismatches, "
            "and deletions from the CIGAR string are counted as consecutive bases with aligned read bases.  Skipped bases "
            "(CIGAR operation \"N\") are counted as uncovered and do not link two regions.");

    addText(parser,
            "This behaviour can be changed by specifying \\fB--link-over-skipped\\fP where these skipped bases are not "
            "counted as being covered but they link two otherwise not connected ROIs into one with possibly 0 counts in "
            "between.");

    addText(parser,
            "For paired data, we by default link two ROIs that have no interconnecting bases if there is a pair of reads "
            "where one mate maps to either ROI.  You can switch off this behaviour using \\fP--ignore-pairing\\fP.  The "
            "forward and reverse strand of the paired read will still be used when using strand specific information.");

    // -----------------------------------------------------------------------
    // Documentation of ROI file format.
    // -----------------------------------------------------------------------

    addTextSection(parser, "ROI File Format");

    addText(parser,
            "The file begins with a sequence of comments, followed by one column header, followed by the records.");

    addText(parser, "\\fIComments\\fP start with a single hash (#) and are ignored in analysis.");

    addText(parser,
            "The \\fIcolumn header\\fP starts with two hashes (##) and is followed by the tab separated column names. "
            "The column names must not contain spaces.");

    addText(parser,
            "The columns \\fI1-7\\fP are fixed.  The last column is always \"counts\".  The columns are as follows");

    addListItem(parser, "1. ref", "The reference name.");
    addListItem(parser, "2. begin_pos", "The begin position.");
    addListItem(parser, "3. end_pos", "The end position.");
    addListItem(parser, "4. region_name", "The name of the region.");
    addListItem(parser, "5. length", "The length of the region.");
    addListItem(parser, "6. strand", "The strand, one of \"+\" and \"-\". \"+\" in case of being not strand-specific.");
    addListItem(parser, "7. max_count", "The largest value of the \\fIcounts\\fP columns.");
    addListItem(parser, "additional columns", "Additional annotation data for the ROI.");
    addListItem(parser, "N. counts", "Comma-separated list of \\fIlength\\fP unsignd integers.");

    // TODO(holtgrew): Write me!

    // -----------------------------------------------------------------------
    // Parsing and Value Extraction
    // -----------------------------------------------------------------------

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.

    options.verbosity = 1;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.inputFileName, parser, "input-file");
    getOptionValue(options.outputFileName, parser, "output-file");

    options.strandSpecific = isSet(parser, "strand-specific");
    options.usePairing = !isSet(parser, "ignore-pairing");
    options.linkOverSkipped = isSet(parser, "link-over-skipped");

	return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
	// Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0 (e.g. help
    // was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Print program name and options.

    if (options.verbosity >= 1)
    {
        std::cerr << "BAM TO ROI\n"
                  << "==========\n"
                  << "\n";
        print(std::cerr, options);
    }

    // -----------------------------------------------------------------------
    // Open input and output file
    // -----------------------------------------------------------------------

    if (options.verbosity >= 1)
        std::cerr << "__OPENING FILES_______________________________________________________________\n"
                  << "\n";

    if (options.verbosity >= 1)
         std::cerr << "Opening " << options.inputFileName << " ...";
	seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(options.inputFileName)))
	{
		std::cerr << "ERROR: Could not open " << options.inputFileName << "\n";
        return 1;
	}
    if (options.verbosity >= 1)
        std::cerr << "OK\n";

    if (options.verbosity >= 1)
         std::cerr << "Opening " << options.outputFileName << " ...";
    seqan::RoiFileOut roiFileOut;
    if (!open(roiFileOut, toCString(options.outputFileName)))
	{
		std::cerr << "ERROR: Could not open " << options.outputFileName << "\n";
        return 1;
	}
    if (options.verbosity >= 1)
        std::cerr << "OK\n";

    // -----------------------------------------------------------------------
    // Build ROIs
    // -----------------------------------------------------------------------

    if (options.verbosity >= 1)
        std::cerr << "__BUILD ROI___________________________________________________________________\n"
                  << "\n"
                  << "Working ...";

    // We have two RoiBuilder objects, one for ROIs on the forward strand and one for ROIs on the reverse strand.  In
    // the case of not being strand specific, we only use the forward builder.

    RoiBuilderOptions roiBuilderOptions(options.verbosity, options.strandSpecific,
                                        options.usePairing, options.linkOverSkipped);
    RoiBuilder roiBuilderF(roiFileOut, roiBuilderOptions);
    roiBuilderF.writeHeader();  // only once
    RoiBuilder roiBuilderR(roiFileOut, roiBuilderOptions);
    // Set the reference sequence names.
    seqan::BamHeader header;
    readHeader(header, bamFileIn);
    for (unsigned i = 0; i < length(contigNames(context(bamFileIn))); ++i)
    {
        appendValue(roiBuilderF.refNames, contigNames(context(bamFileIn))[i]);
        appendValue(roiBuilderR.refNames, contigNames(context(bamFileIn))[i]);
    }

    // TODO(holtgrew): This is only suited for the Illumina mate pair protocol at the moment (--> <--).
    int oldRId = 0;
    int oldPos = 0;
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);

        // Break if record is unmapped.
        if (hasFlagUnmapped(record))
            break;  // Found first unmapped record!

        // Check sorting.
        if (record.rID < oldRId || (record.rID == oldRId && record.beginPos < oldPos))
        {
            std::cerr << "\nERROR: The BAM file is not sorted properly!\n";
            return 1;
        }
        oldRId = record.rID;
        oldPos = record.beginPos;

        // Process record, different cases, depending on strandedness and pairedness.
        if (hasFlagMultiple(record))
        {
            if (hasFlagFirst(record) && !hasFlagRC(record))
                roiBuilderF.pushRecord(record);
            else if (hasFlagFirst(record) && hasFlagRC(record))
                roiBuilderR.pushRecord(record);
            else if (hasFlagLast(record) && !hasFlagRC(record))
                roiBuilderF.pushRecord(record);
            else  // (hasFlagLast(record) && hasFlagRC(record))
                roiBuilderR.pushRecord(record);
        }
        else
        {
            if (options.strandSpecific)
            {
                if (hasFlagRC(record))
                    roiBuilderR.pushRecord(record);
                else
                    roiBuilderF.pushRecord(record);
            }
            else
            {
                roiBuilderF.pushRecord(record);
            }
        }
    }

    if (options.verbosity >= 1)
        std::cerr << " OK\n"
                  << "\n"
                  << "Done.\n";

    return 0;
}
