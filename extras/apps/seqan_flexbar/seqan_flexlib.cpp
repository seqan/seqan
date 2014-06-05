// ==========================================================================
//                                SeqAn-Flexbar
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// Author: Benjamin Strauch <b.strauch@fu-berlin.de>
// ==========================================================================
// This file provides the argument parsing functionality is of
// seqan-flexbar which is based in the implementation of the original
// flexbar program in [1].  In addition, the file contains the main body
// of the program, selecting the sub-routines.
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================

#include "seqan_flexbar.h"

#define SEQAN_PROFILE
#ifdef SEQAN_ENABLE_DEBUG
#define DEBUG_MSG(str) do { std::cerr << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#include <iostream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#ifdef _WIN32
#include <direct.h>
#endif

// Headers for creating subdirectories.
#include <errno.h>
// For setting the number of threads used by OpenMP.
#ifdef _OPENMP
#include <omp.h>
#endif

#include "read_trimming.h"
#include "adapter_trimming.h"
#include "demultiplex.h"
#include "general_processing.h"

typedef seqan::String<seqan::Dna5Q> Dna5QString;

// Global variables are evil, this is for adaption and should be removed
// after refactorization.
FlexiProgram flexiProgram;

// --------------------------------------------------------------------------
// Class ArgumentParserBuilder
// --------------------------------------------------------------------------

class ArgumentParserBuilder
{
public:
    virtual seqan::ArgumentParser build() = 0;

protected:
    void addGeneralOptions(seqan::ArgumentParser & parser);
    void addFilteringOptions(seqan::ArgumentParser & parser);
    void addDemultiplexingOptions(seqan::ArgumentParser & parser);
    void addAdapterTrimmingOptions(seqan::ArgumentParser & parser);
    void addReadTrimmingOptions(seqan::ArgumentParser & parser);
};

void ArgumentParserBuilder::addGeneralOptions(seqan::ArgumentParser & parser)
{
    // GENERAL OPTIONS -----------------------------------------
    addSection(parser, "General Options");
    
    seqan::ArgParseOption recordOpt = seqan::ArgParseOption(
        "r", "records", "Number of records to be read in one run.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(recordOpt, 10000);
    setMinValue(recordOpt, "10");
    addOption(parser, recordOpt);

    seqan::ArgParseOption noQualOpt = seqan::ArgParseOption(
        "nq", "noQualities", "Force .fasta format for output files.");
        addOption(parser, noQualOpt);

    seqan::ArgParseOption threadOpt = seqan::ArgParseOption(
        "tnum", "threads", "Number of threads used.",
        seqan::ArgParseOption::INTEGER, "THREADS");
    setDefaultValue(threadOpt, 1);
    setMinValue(threadOpt, "1");
    addOption(parser, threadOpt);

    if(flexiProgram == ADAPTER_REMOVAL || flexiProgram == FILTERING || flexiProgram == QUALITY_CONTROL)
    {
        seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
            "o", "output", "Name of the output file.",
            seqan::ArgParseOption::OUTPUTFILE, "OUTPUT");
        setValidValues(outputOpt, ".fasta .fa .fasta.gz .fa.gz .fastq .fq .fastq.gz .fq.gz");
        addOption(parser, outputOpt);
    }
    else
    {
        seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
            "o", "output", "Folder for output (must already exist).",
            seqan::ArgParseOption::OUTPUTPREFIX, "OUTPUT");
        setValidValues(outputOpt, ".fasta .fa .fasta.gz .fa.gz .fastq .fq .fastq.gz .fq.gz");
        addOption(parser, outputOpt);
    }

    if (flexiProgram == ALL_STEPS)
    {
        seqan::ArgParseOption adTagOpt = seqan::ArgParseOption(
            "t", "tag", "Tags IDs of sequences which had adapters removed and/or were quality-trimmed.");
        addOption(parser, adTagOpt);
    }

    seqan::ArgParseOption adInfoOpt = seqan::ArgParseOption(
        "ni", "noInfo", "Don't print paramter overwiev to console.");
    addOption(parser, adInfoOpt);
}

void ArgumentParserBuilder::addFilteringOptions(seqan::ArgumentParser & parser)
{
    // Preprocessing and Filtering
    addSection(parser, "Preprocessing and Filtering");

    seqan::ArgParseOption leftTrimOpt = seqan::ArgParseOption(
        "tl", "trimLeft", "Number of Bases to be trimmed from the 5'end(s) before further processing.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(leftTrimOpt, 0);
    setMinValue(leftTrimOpt, "0");
    addOption(parser, leftTrimOpt);

    seqan::ArgParseOption rigthTrimOpt = seqan::ArgParseOption(
        "tr", "trimRight", "Number of Bases to be trimmed from the 3'end(s) before further processing.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(rigthTrimOpt, 0);
    setMinValue(rigthTrimOpt, "0");
    addOption(parser, rigthTrimOpt);

    seqan::ArgParseOption minLenOpt = seqan::ArgParseOption(
        "ml", "minLength", "Required minimal length of reads after all PREprocessing steps.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(minLenOpt, 0);
    setMinValue(minLenOpt, "0");
    addOption(parser, minLenOpt);

    seqan::ArgParseOption uncalledOpt = seqan::ArgParseOption(
        "u", "uncalled", "Number of allowed uncalled bases per sequence.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(uncalledOpt, 0);
    setMinValue(uncalledOpt, "0");
    addOption(parser, uncalledOpt);

    seqan::ArgParseOption substituteOption = seqan::ArgParseOption(
        "s", "substitute", "Substitue Dna-character for N's.",
        seqan::ArgParseOption::STRING, "BASE");
    setDefaultValue(substituteOption, "A");
    setValidValues(substituteOption, "A C G T");
    addOption(parser, substituteOption);

     seqan::ArgParseOption finMinLenOpt = seqan::ArgParseOption(
        "fm", "finalMinLength", "Deletes read (and mate)"
        " if on of them is shorter than the given value after the complete worflow.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setMinValue(finMinLenOpt, "1");
    addOption(parser, finMinLenOpt);

    seqan::ArgParseOption finLenOpt = seqan::ArgParseOption(
        "fl", "finalLength", "Trims reads to desired length after the complete workflow.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(finLenOpt, 0);
    setMinValue(finLenOpt, "1");
    addOption(parser, finLenOpt);
}

void ArgumentParserBuilder::addDemultiplexingOptions(seqan::ArgumentParser & parser)
{
    // Barcode Demultiplexing
    addSection(parser, "Demultiplexing Options");
    
    seqan::ArgParseOption barcodeFileOpt = seqan::ArgParseOption(
        "b", "barcodes", "FastA file containing the used barcodes and their IDs. Necessary for demutiplexing.",
        seqan::ArgParseArgument::INPUTFILE, "BARCODE_FILE");
    setValidValues(barcodeFileOpt, "fasta fasta.gz fastq fastq.gz");
    addOption(parser, barcodeFileOpt);

    seqan::ArgParseOption multiplexFileOpt = seqan::ArgParseOption(
        "x", "multiplex", "FastA/FastQ file containing the barcode for each read.",
        seqan::ArgParseArgument::INPUTFILE, "MULTIPLEX_FILE");
    setValidValues(multiplexFileOpt, "fasta fasta.gz fastq fastq.gz");
    addOption(parser, multiplexFileOpt);
    
    addOption(parser, seqan::ArgParseOption(
        "app", "approximate", "Select approximate barcode demultiplexing, allowing one mismatch."));

    addOption(parser, seqan::ArgParseOption(
        "hc", "hardClip", "Select hardClip option, clipping the first length(barcode) bases in any case."));

    addOption(parser, seqan::ArgParseOption(
        "ex", "exclude", "Exclude unidentified reads from further processing."));
}

void ArgumentParserBuilder::addAdapterTrimmingOptions(seqan::ArgumentParser & parser)
{
    // ADAPTER TRIMMING
    addSection(parser, "Adapter removal options");
    seqan::ArgParseOption adapterFileOpt = seqan::ArgParseOption(
        "a", "adapters", "FastA file containing the two adapter sequences. "
        "The adapters according to the layout: 5'-adapter1-read-adapter2-3'.",
        seqan::ArgParseArgument::INPUTFILE, "ADAPTER_FILE");
    setValidValues(adapterFileOpt, "fasta fasta.gz fastq fastq.gz");
    addOption(parser, adapterFileOpt);

    seqan::ArgParseOption noAdapterOpt = seqan::ArgParseOption(
        "na", "no-adapter", "Trim adapters from paired-end reads without using reference adapters.");
    addOption(parser, noAdapterOpt);

    seqan::ArgParseOption pairedModeOpt = seqan::ArgParseOption(
        "np", "no-paired", "Trim paired-end input with single-end trimming method.");
    addOption(parser, pairedModeOpt);

    seqan::ArgParseOption rateOpt = seqan::ArgParseOption(
        "e", "errors", "Allowed errors in adapter detection.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(rateOpt, 0);
    addOption(parser, rateOpt);

    seqan::ArgParseOption overlapOpt = seqan::ArgParseOption(
        "ol", "overlap", "Minimum length of overlap for a significant adapter match.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(overlapOpt, 0);
    addOption(parser, overlapOpt);

    if (flexiProgram != ALL_STEPS)
    {
        seqan::ArgParseOption adTagOpt = seqan::ArgParseOption(
            "t", "tag", "Tags IDs of sequences which had adapters removed.");
        addOption(parser, adTagOpt);
    }

}

void ArgumentParserBuilder::addReadTrimmingOptions(seqan::ArgumentParser & parser)
{
    // READ TRIMMING
    addSection(parser, "Quality trimming options");
    seqan::ArgParseOption qualOpt = seqan::ArgParseOption(
        "q", "quality", "Quality threshold for read trimming.",
        seqan::ArgParseArgument::INTEGER, "PHRED");
    setMinValue(qualOpt, "0");
    setMaxValue(qualOpt, "40");
    addOption(parser, qualOpt);

    seqan::ArgParseOption lenOpt = seqan::ArgParseOption(
        "l", "length", 
        "Minimum read length after trimming. " 
        "Shorter reads will substitutet by a single N or removed if the paired read also is too short.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(lenOpt, 1);
    setMinValue(lenOpt, "1");
    addOption(parser, lenOpt);

    seqan::ArgParseOption trimOpt = seqan::ArgParseOption(
        "m", "method", "Method for trimming reads.",
        seqan::ArgParseArgument::STRING, "METHOD");
    setDefaultValue(trimOpt, "WIN");
    setValidValues(trimOpt, "WIN BWA TAIL");
    addOption(parser, trimOpt);

    if (flexiProgram != ALL_STEPS)
    {
        seqan::ArgParseOption adTagOpt = seqan::ArgParseOption(
            "t", "tag", "Tags IDs of sequences which were quality-trimmed.");
        addOption(parser, adTagOpt);
    }

    addTextSection(parser, "EXAMPLES");
    std::string appName;
    assign(appName, getAppName(parser));
    addText(parser, appName + " INPUT -b BARCODES -app -o FOLDER\\Prefix");
    addText(parser, appName + " INPUT -a ADAPTERS -q 20 -o FOLDER");
    addText(parser, appName + " INPUT1 INPUT2 -b BARCODES -na -ex -q 20 -tnum 2 -l 50 -o FOLDER");
}

// --------------------------------------------------------------------------
// Class FilteringParserBuilder
// --------------------------------------------------------------------------

class FilteringParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build();

private:
    void addHeader(seqan::ArgumentParser & parser);
};

void FilteringParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Filtering Toolkit of seqan_flexbar.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
       "This program is a sub-routine of SeqAn-Flexbar (a reimplementation of"
       " the original flexbar[1]) and can be used to filter reads and apply "
       "sequence independent trimming options");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
            "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
            "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    std::string version = "1.0.1";
#ifdef SEQAN_REVISION
    version += std::string(" [") + std::string(SEQAN_REVISION) + "]";
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#endif
    seqan::setVersion(parser, version);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUTFILE, "READS", true);
    setValidValues(fileArg, ".fasta .fa .fasta.gz .fa.gz .fastq .fq .fastq.gz .fq.gz");
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

seqan::ArgumentParser FilteringParserBuilder::build()
{
    seqan::ArgumentParser parser("sflexFilter");

    addHeader(parser);
    addGeneralOptions(parser);
    addFilteringOptions(parser);

    return parser;
}

// --------------------------------------------------------------------------
// Class AdapterRemovalParserBuilder
// --------------------------------------------------------------------------

class AdapterRemovalParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build();

private:
    void addHeader(seqan::ArgumentParser & parser);
};

void AdapterRemovalParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Adapter Removal Toolkit of seqan_flexbar.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of SeqAn-Flexbar (a reimplementation of"
        " the original flexbar[1]) and removes adapter sequences from reads.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
            "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
            "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    std::string version = "1.0.1";
#ifdef SEQAN_REVISION
    version += std::string(" [") + std::string(SEQAN_REVISION) + "]";
#endif
#ifdef SEQAN_DATE
    seqan::setDate(parser, SEQAN_DATE);
#endif
    seqan::setVersion(parser, version);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUTFILE, "READS", true);
    setValidValues(fileArg, ".fasta .fa .fasta.gz .fa.gz .fastq .fq .fastq.gz .fq.gz");
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

seqan::ArgumentParser AdapterRemovalParserBuilder::build()
{
    seqan::ArgumentParser parser("sflexAR");

    addHeader(parser);
    addGeneralOptions(parser);
    addAdapterTrimmingOptions(parser);

    return parser;
}



// --------------------------------------------------------------------------
// Class DemultiplexingParserBuilder
// --------------------------------------------------------------------------

class DemultiplexingParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build();

private:
    void addHeader(seqan::ArgumentParser & parser);
};

void DemultiplexingParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Demultiplexing Toolkit of seqan_flexbar.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of SeqAn-Flexbar (a reimplementation of"
        " the original flexbar[1]) and can be used for demultiplexing of reads.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
            "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
            "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");


    std::string version = "1.0.1";
#ifdef SEQAN_REVISION
    version += std::string(" [") + std::string(SEQAN_REVISION) + "]";
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#endif
    seqan::setVersion(parser, version);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUTFILE, "READS", true);
    setValidValues(fileArg, ".fasta .fa .fasta.gz .fa.gz .fastq .fq .fastq.gz .fq.gz");
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

seqan::ArgumentParser DemultiplexingParserBuilder::build()
{
    seqan::ArgumentParser parser("sflexDMulti");

    addHeader(parser);
    addGeneralOptions(parser);
    addDemultiplexingOptions(parser);

    return parser;
}


// --------------------------------------------------------------------------
// Class QualityControlParserBuilder
// --------------------------------------------------------------------------

class QualityControlParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build();

private:
    void addHeader(seqan::ArgumentParser & parser);
};

void QualityControlParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Quality Control Toolkit of seqan_flexbar.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of SeqAn-Flexbar (a reimplementation of"
        " the original flexbar [1]) and can be used for quality controlling of reads.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
            "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
            "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    std::string version = "1.0.1";
#ifdef SEQAN_REVISION
    version += std::string(" [") + std::string(SEQAN_REVISION) + "]";
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#endif
    seqan::setVersion(parser, version);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUTFILE, "READS", true);
    setValidValues(fileArg, ".fasta .fa .fasta.gz .fa.gz .fastq .fq .fastq.gz .fq.gz");
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

seqan::ArgumentParser QualityControlParserBuilder::build()
{
    seqan::ArgumentParser parser("sflexQC");

    addHeader(parser);
    addGeneralOptions(parser);
    addReadTrimmingOptions(parser);

    return parser;
}

// --------------------------------------------------------------------------
// Class AllStepsParserBuilder
// --------------------------------------------------------------------------

class AllStepsParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build();

private:
    void addHeader(seqan::ArgumentParser & parser);
};

void AllStepsParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn NGS-Data Processing Toolkit");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "SeqAn-Flexbar is a toolkit for the processing of sequenced NGS reads "
        "and based on the original Flexbar implementation of Dodt [1]. It is "
        "possible to demultiplex the reads and order them according to "
        "different kind of barcodes, to remove adapter contamination from "
        "reads, to trim low quality bases, filter N's or trim whole reads. The "
        "different tools are controlled through command line parameters and can "
        "operate on both single- and paired-end read data.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
            "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
            "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    addDescription(parser, "(c) Copyright 2008-2013 by Sebastian Roskosch.");

    std::string version = "1.0.1";
#ifdef SEQAN_REVISION
    version += std::string(" [") + std::string(SEQAN_REVISION) + "]";
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#endif
    seqan::setVersion(parser, version);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUTFILE, "READS", true);
    setValidValues(fileArg, ".fasta .fa .fasta.gz .fa.gz .fastq .fq .fastq.gz .fq.gz");
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

seqan::ArgumentParser AllStepsParserBuilder::build()
{
    seqan::ArgumentParser parser("seqan_flexbar");

    addHeader(parser);
    addGeneralOptions(parser);
    addFilteringOptions(parser);
    addDemultiplexingOptions(parser);
    addAdapterTrimmingOptions(parser);
    addReadTrimmingOptions(parser);

    return parser;
}

// --------------------------------------------------------------------------
// Function initParser()
// --------------------------------------------------------------------------

//Defining the the argument parser
seqan::ArgumentParser initParser()
{
    std::auto_ptr<ArgumentParserBuilder> argParseBuilder;

    switch (flexiProgram)
    {
        case ADAPTER_REMOVAL:
            argParseBuilder.reset(new AdapterRemovalParserBuilder());
            break;
        case DEMULTIPLEXING:
            argParseBuilder.reset(new DemultiplexingParserBuilder());
            break;
        case FILTERING:
            argParseBuilder.reset(new FilteringParserBuilder());
            break;
        case QUALITY_CONTROL:
            argParseBuilder.reset(new QualityControlParserBuilder());
            break;
        case ALL_STEPS:
            argParseBuilder.reset(new AllStepsParserBuilder());
            break;
        default:
            SEQAN_FAIL("Invalid program type.");
    }

    return argParseBuilder->build();
}

// --------------------------------------------------------------------------
// Class ProcessingParams
// --------------------------------------------------------------------------

/**
.Class.ProcessingParams:
..summary:Struct holding all general pre- and postprocessing parameters.
.Memvar.ProcessingParams#substitute
..class:Class.ProcessingParams
..summary:Substitute Base for N's
.Memvar.ProcessingParams#uncalled
..class:Class.ProcessingParams
..summary:Number of allowed N's per sequence.
.Memvar.ProcessingParams#trimLeft
..class:Class.ProcessingParams
..summary:Number of bases to be cut from 5'-end before further processing.
.Memvar.ProcessingParams#trimRight
..class:Class.ProcessingParams
..summary:Number of bases to be cut from 3'-end before further processing.
.Memvar.ProcessingParams#minLen
..class:Class.ProcessingParams
..summary:Minimum length sequence length AFTER trimming at 5'- and 3'-end.
.Memvar.ProcessingParams#finalMinLength
..class:Class.ProcessingParams
..summary:Minimum sequence length AFTER COMPLETE workflow. If a sequence (or it's mate)
 is shorter, it/they will be removed.
.Memvar.ProcessingParams#finalLength
..class:Class.ProcessingParams
..summary:Sequence length all sequences will be trimmed to (cutting from 3'-end) to AFTER COMPLETE workflow.
.Memvar.ProcessingParams#runPre
..class:Class.ProcessingParams
..summary:bool indicating that preprocessing actions shall be performed (true).
.Memvar.ProcessingParams#runPost
..class:Class.ProcessingParams
..summary:bool indicating that postprocessing actions shall be performed (true).

.Memfunc.ProcessingParams#ProcessingParams:
..class:Class.ProcessingParams
..summary:Default Constructor
..signature:ProcessingParams(void)
..remarks: finalLength is initialised with 0, because 0 indicates that no action shall be performed.
*/
struct ProcessingParams
{
    seqan::Dna substitute;
    unsigned uncalled;
    unsigned trimLeft;
    unsigned trimRight;
    unsigned minLen;
    unsigned finalMinLength;
    unsigned finalLength;
    bool runPre;
    bool runPost;

    ProcessingParams() :
        substitute('A'), 
        uncalled(0),
        trimLeft(0),
        trimRight(0),
        minLen(0),
        finalMinLength(0),
        finalLength(0),
        runPre(false),
        runPost(false) {};
};

/**
.Class.TrimmingMode:
..summary:Enum for storing the trimming mode selected by the user
*/
enum TrimmingMode
{
    E_WINDOW,
    E_BWA,
    E_TAIL
};

/**
.Class.MatchMode:
..summary:Enum for storing the match options used for adapter trimming.
*/
enum MatchMode
{
    E_AUTO,
    E_USER
};

/**
.Class.DemultiplexingParams:
..summary:Struct holding all demultiplexing parameters.
.Memvar.DemultiplexingParams#barcodeFile
..class:Class.DemultiplexingParams
..summary:String holding the path to the barcode-file.
.Memvar.DemultiplexParams#barcodes
..class:Class.DemultiplexingParams
..summary:StringSet holding the barcodes.
.Memvar.DemultiplexParams#multiplexFile
..class:Class.DemultiplexingParams
..summary:String hilding the path to the multiplex-file.
.Memvar.DemultiplexParams#multiplex
..class:Class.DemultiplexingParams
..summary:StringSet holding the multiplex barcodes.
.Memvar.DemultiplexParams#approximate
..class:Class.DemultiplexingParams
..summary:bool indicating if approximate search shall be used (true).
.Memvar.DemultiplexParams#hardClip
..class:Class.DemultiplexingParams
..summary:bool indicating if hardClip shall be used (true).
.Memvar.DemultiplexParams#run
..class:Class.DemultiplexingParams
..summary:bool indicating if demultiplexing operations shall be performed (true).
.Memvar.DemultiplexParams#runx
..class:Class.DemultiplexingParams
..summary:bool indicating if multiplex barcodes are used (true).
.Memvar.DemultiplexParams#stats
..class:Class.DemultiplexingParams
..summary:DemultiplexStats-object holding interesting numbers about the demultiplexing process.
.Memvar.DemultiplexParams#exclude
..class:Class.DemultiplexingParams
..summary:bool indicating if unidentified sequences shall be excluded from further processing (true).
*/
struct DemultiplexingParams
{
    seqan::String<char> barcodeFile;
    seqan::StringSet<seqan::String<seqan::Dna5> > barcodes;
    seqan::StringSet<seqan::String<char> > barcodeIds;
    seqan::String<char> multiplexFile;
    seqan::StringSet<seqan::String<seqan::Dna5Q> > multiplex;
    bool approximate;
    bool hardClip;
    bool run;
    bool runx;
    bool exclude;
    DemultiplexStats stats;

    DemultiplexingParams() :
        approximate(false),
        hardClip(false),
        run(false),
        runx(false)
    {};
};
/**
.Class.AdapterTrimmingParams:
..summary:Struct holding all adapter trimming parameters.
.Memvar.AdapterTrimmingParams#paired
..class:Class.AdapterTrimmingParams
..summary:bool indicating if paired-end data is used (true).
.Memvar.AdapterTrimmingParams#noAdapter
..class:Class.AdapterTrimmingParams
..summary:bool indicating if no adapter file is provided (true).
.Memvar.AdapterTrimmingParams#run
..class:Class.AdapterTrimmingParams
..summary:bool indicating if trimming operations shall be performed.
.Memvar.AdapterTrimmingParams#adapter1
..class:Class.AdapterTrimmingParams
..summary:String holding the first adapter sequence.
.Memvar.AdapterTrimmingParams#adapter2
..class:Class.AdapterTrimmingParams
..summary:String holding the second adapter sequence.
.Memvar.AdapterTrimmingParams#mode
..class:Class.AdapterTrimmingParams
..summary:The Mode-object determining which trimming mode shall be used.
.Memvar.AdapterTrimmingParams#mmode
..class:Class.AdapterTrimmingParams
..summary:Enum indicating which trimming mode is used. Necessary for casting mode.
.Memvar.AdapterTrimmingParams#stats
..class:Class.AdapterTrimmingParams
..summary:AdapterTrimmingStats object holding interesting numbers about the trimmed process.
..type:Class.AdapterTrimmingStats
 */
struct AdapterTrimmingParams
{
    bool paired;
    bool noAdapter;
    bool run;
    Dna5QString adapter1;
    Dna5QString adapter2;
    Mode mode;
    MatchMode mmode;
    AdapterTrimmingStats stats;
    AdapterTrimmingParams() : paired(false), noAdapter(false), run(false), mmode(E_AUTO) {};
};

/**
.Class.QualityTrimmingParams:
..summary:Struct holding all quality trimming  parameters.
.Memvar.QualityTrimmingParams#trim_mode
..class:Class.QualityTrimmingParams
..summary:The TrimmingMode-object, determining which algorithm shall be used.
..type:Class.TrimmingMode
.Memvar.QualityTrimmingParams#cutoff
..class:Class.QualityTrimmingParams
..summary:int holding the cutoff score.
.Memvar.QualityTrimmingParams#min_length
..class:Class.QualityTrimmingParams
..summary:int holding the minimal length of a sequence after trimming.
.Memvar.QualityTrimmingParams#run
..class:Class.QualityTrimmingParams
..summary:bool indicating if quality trimming operations shall be performed.
.Memvar.QualityTrimmingParams#stats
..class:Class.QualityTrimmingParams
..summary:QualityTrimmingsStats object holding interesting numbers about the quality trimming process.
..type:Class.QualityTrimmingStats
 */
struct QualityTrimmingParams
{
    TrimmingMode trim_mode;
    int cutoff;
    int min_length;
    bool run;
    QualityTrimmingStats stats;
       QualityTrimmingParams() : trim_mode(E_WINDOW), cutoff(-1), min_length(1), run(false) {};
};

/**
.Class.ProgramParams:
..summary:Struct holding generall program paramters.
.Memvar.ProgramParams#fileCount
..class:Class.ProgramParams
..summary:int holding the number of read files.
.Memvar.ProgramParams#readCount
..class:Class.ProgramParams
..summary:int holding the total number of reads.
.Memvar.ProgramParams#processTime
..class:Class.ProgramParams
..summary:double for storing the process time.
.Memvar.ProgramParams#ioTime
..class:Class.ProgramParams
..summary:double for storing the I/O time.
.Memvar.ProgramParams#fileStream1
..class:Class.ProgramParams
..summary:SequenceSteam object for the forward-reads file.
..type:Class.SequenceStream
.Memvar.ProgramParams#fileStream2
..class:Class.ProgramParams
..summary:SequenceSteam object for the backward-reads file.
..type:Class.SequenceStream
*/
struct ProgramParams
{
    int fileCount;
    int readCount;
    double processTime, ioTime;
    seqan::SequenceStream fileStream1, fileStream2;

    ProgramParams() : fileCount(0), readCount(0), processTime(0), ioTime(0) {};
};

/**
.Class.OutputStreams:
..summary:Class that dynamically manages output streams that write out sets of sequences.
..see:Class.SequenceStream
.Memvar.OutputStreams#pairedFileStreams
..class:Class.OutputStreams
..summary:map of one int and a pair of SequenceStream Pointers. Used for paired-end data.
.Memvar.OutputStreams#fileStreams
..class:Class.OutputStreams
..summary:map of one int and a SequenceStream Pointer. Used for single-end data.
.Memvar.OutputStreams#basePath
..class:Class.OutputStreams
..summary:CharString for storing the path.
.Memvar.OutputStreams#extension
..class:Class.OutputStreams
..summary:CharString for storing the file extension.


.Memfunc.OutputStreams#OutputStreams:
..class:Class.OutputStreams
..summary:Constructor
..signature:OutputStreams(base, format, compress)
..param.base:CharString holding of the directory path.
...type:nolink:CharString
..param.format:File format of the output files.
...type:nolink:SeqIOFileFormat_::Type
..param.compress:bool, indicating if the output shall be in gzip format (true)
...type:nolink:bool

.Memfunc.OutputStreams#exists:
..class:Class.OutputStreams
..summary:Checks if a key already exists in the map
..signature:exists(key, map)
..param.key:The key.
...type:nolink:int
..param.map:The map.
...type:nolink:std::map
..returns:A bool, true if the key exists, false otherwise.
...type:nolink:bool

.Memfunc.OutputStreams#addStream:
..class:Class.OutputStreams
..summary:Adds a new output stream to the collection of streams. Used for single-end data.
..signature:addStream(filename, id)
..param.filename:The name of the file to be created.
...type:nolink:CharString
..param.id:The associated id used to identify the stream.
...type:nolink:int
..returns:void

.Memfunc.OutputStreams#addStreams:
..class:Class.OutputStreams
..summary:Adds two new output streams to the collection of streams. Used for paired-end data.
..signature:addStreams(filename1, filename2, id)
..param.filename1:The name of the forward-read file to be created.
...type:nolink:CharString
..param.filename1:The name of the backward-read file to be created.
...type:nolink:CharString
..param.id:The associated id used to identify the streams.
...type:nolink:int
..returns:void

.Memfunc.OutputStreams#updateStreams:
..class:Class.OutputStreams
..summary:This method takes a String of integers and checks if these integers are
 already associated with a stream. If not, a new stream is added and the opened
 file is named according to the list of names. One or two files are created.
..signature:updateStreams(map, names, pair)
..param.map:The list of IDs for which the existence of a file shall be checked.
...type:Class.StringSet
..param.names:The list of names to be used when creating new files.
...type:Class.StringSet
..param.pair:Indicates whether one or two (a pair of files) shall be created.
...type:nolink:bool
..returns:void

.Memfunc.OutputStreams#writeSeqs:
..class:Class.OutputStreams
..summary:Writes the sets of ids and sequences to their corresponding files.
..signature:writeSeqs(ids, seqs, map, names)
..signature:writeSeqs(ids, seqs, ids2, seqs2, map, names)
..param.ids:A list of sets of sequence IDs for the forward-reads (as returned by readRecord etc.).
...type:Class.StringSet
..param.seqs:A StringSet of sets of forward-reads.
...type:Class.StringSet
..param.ids2:A list of sets of sequence IDs for the bacward-reads (as returned by readRecord etc.).
...type:Class.StringSet
..param.seqs2:A StringSet of sets of backward-reads.
...type:Class.StringSet
..param.map:A map of the sets of sequences to their corresponding output streams.
...type:nolink:std::map
..param.names:Names to be used when creating new streams.
...type:Class.StringSet
..returns:void    
*/
class OutputStreams
{
    typedef seqan::SequenceStream * PSeqStream;
    typedef seqan::Pair<PSeqStream, PSeqStream> TStreamPair;
    std::map<int, TStreamPair> pairedFileStreams;
    std::map<int, PSeqStream> fileStreams;
    const seqan::CharString basePath;
    seqan::CharString extension;
public:
//Constructor prepares the file extension which will be used for all streams created 
//by this object and saves a base directory path.
    OutputStreams(seqan::CharString base, bool /*noQuality*/) : basePath(base)
    {
        /*
        seqan::CharString fileExt("");
        if (noQuality)
        {
            seqan::append(fileExt, seqan::CharString(".fasta"));
        }
        else
        {
            seqan::append(fileExt, seqan::CharString(".fastq"));
        }
        if (compress)
        {
            seqan::append(fileExt, seqan::CharString(".gz"));
        }
        extension = fileExt;
        */
    }    
     //Checks whether a key exists in a map. 
    template <typename TKey, typename TMap>
    bool exists(TKey& key, TMap& map)
    {
        return map.find(key) == map.end();
    }
    //Adds a new output streams to the collection of streams.
    void addStream(seqan::CharString fileName, int id)
    {
        //Prepend basePath and append file extension to the filename.
        seqan::CharString path = basePath;
        seqan::append(path, fileName);
        seqan::append(path, this->extension);
        char* file = seqan::toCString(path);
        PSeqStream stream = new SequenceStream(file, seqan::SequenceStream::WRITE);
        fileStreams[id] = stream;
    }
     //Adds two new output streams to the collection of streams.
    void addStreams(seqan::CharString fileName1, seqan::CharString fileName2, int id)
    {
        //Prepend basePath and append file extension to the filename.
        seqan::CharString path1 = basePath, path2 = basePath;
        seqan::append(path1, fileName1);
        seqan::append(path1, this->extension);
        seqan::append(path2, fileName2);
        seqan::append(path2, this->extension);
        char* file1 = seqan::toCString(path1);
        char* file2 = seqan::toCString(path2);
        PSeqStream stream1 = new SequenceStream(file1, seqan::SequenceStream::WRITE);
        PSeqStream stream2 = new SequenceStream(file2, seqan::SequenceStream::WRITE);
        pairedFileStreams[id] = TStreamPair(stream1, stream2);
    }
    //This method takes a String of integers and checks if these integers are
    //already associated with a stream. If not, a new stream is added and the opened
    //file is named according to the list of names. One or two files are created.
    template <typename TMap, typename TNames>
    void updateStreams(TMap& map, TNames& names, bool pair)
    {
        for (unsigned i=0; i < length(map); ++i)
        {
            unsigned streamIndex = map[i];
            bool missing = pair ? exists(streamIndex, pairedFileStreams) : exists(streamIndex, fileStreams);
            // If no stream for this id exists, create one.
            if (missing)
            {
                // If the index is 0 (unidentified) create special stream.
                // Otherwise use index to get appropriate name for output file.
                seqan::CharString file;
                if (streamIndex > 0)
                {
                    file = names[streamIndex-1];
                }
                else
                {
                    file = seqan::CharString("unidentified");
                }
                // Add file extension to stream and create it.
                if (pair)
                {
                    // Create a new subfolder at basePath/[barcodeID, unidentified].
                    seqan::CharString folderPath(basePath);
                    seqan::append(folderPath, file);
                    int result = 0;
                    #ifdef __linux__
                        result = mkdir(seqan::toCString(folderPath), 0777);
                    #elif _WIN32
                        result = _mkdir(seqan::toCString(folderPath));
                    #endif
                        if (result==-1)
                    {
                        if (errno == EEXIST)
                        {
                            std::cerr << "Warning: Directory " << folderPath << " already exists.\n";
                        }
                    }
                    // Turn file target from [barcodeID,unidentified]
                    // to subfolder [barcodeID, unidentified]/[barcodeID, unidentified]
                    seqan::CharString filePath(file);
                    seqan::append(file, "/");
                    seqan::append(file, filePath);
                    // To differentiate between the paired reads, add index to the file name.
                    seqan::CharString file2 = file;
                    seqan::append(file, seqan::CharString("_1"));
                    seqan::append(file2, seqan::CharString("_2"));
                    this->addStreams(file, file2, streamIndex);
                }
                else
                {
                    this->addStream(file, streamIndex);
                }
            }
        }
    }
    //Writes the sets of ids and sequences to their corresponding files. Used for single-end data.
    template <typename TIds, typename TSeqs, typename TMap, typename TNames>
    void writeSeqs(TIds& ids, TSeqs& seqs, TMap& map, TNames& names)
    {
        updateStreams(map, names, false);
        for (unsigned i = 0; i < length(seqs); ++i)
        {
            unsigned streamIndex = map[i];
            seqan::writeAll(*fileStreams[streamIndex], ids[i], seqs[i]);
        }
    }
    //Overload of writeSeqs for paired-end data.
    template <typename TIds, typename TSeqs, typename TMap, typename TNames>
    void writeSeqs(TIds& ids1, TSeqs& seqs1, TIds& ids2, TSeqs& seqs2, TMap& map, TNames& names)
    {
        updateStreams(map, names, true);
        for (unsigned i = 0; i < length(seqs1); ++i)
        {
            unsigned streamIndex = map[i];
            TStreamPair tmp = pairedFileStreams[streamIndex];
            seqan::writeAll(*tmp.i1, ids1[i], seqs1[i]);
            seqan::writeAll(*tmp.i2, ids2[i], seqs2[i]);
        }
    }
    //Destructor of the object holding the output streams. Needed
     //to destroy the streams properly after they have been created with new.
    ~OutputStreams()
    {
        // Delete all created file streams.
        for (unsigned i=0; i < length(fileStreams); ++i)
        {
            delete fileStreams[i];
        }
        for (unsigned i=0; i < length(pairedFileStreams); ++i)
        {
            TStreamPair tmp = pairedFileStreams[i];
            delete tmp.i1;
            delete tmp.i2;
        }
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Function lastOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TToken>
typename Iterator<TString const, Standard>::Type
lastOf(TString const & string, TToken const & token)
{
    typedef typename Iterator<TString const, Standard>::Type TIterator;

    TIterator it = end(string, Standard()) - length(token);

    for (TIterator itBegin = begin(string, Standard());
         it != itBegin && !isEqual(infix(string, it, it + length(token)), token);
         goPrevious(it)) ;

    return it;
}

// ----------------------------------------------------------------------------
// Helper Function trimExtension()
// ----------------------------------------------------------------------------

template <typename TString>
Segment<TString, PrefixSegment>
trimExtension(TString & string)
{
    return prefix(string, lastOf(string, '.'));
}

// ----------------------------------------------------------------------------
// Helper Function getExtension()
// ----------------------------------------------------------------------------

template <typename TString>
Segment<TString, SuffixSegment>
getExtension(TString & string)
{
    return suffix(string, lastOf(string, '.') + 1);
}


/**
.Function.loadSeqs:
..summary:Function for loading the sequence files.
..signature:loadSeqs(seqStream, ids, seqs, records)
..param.ids:StringSet of String<Char> the IDs shall be stored in.
...type:Class.StringSet
..param.seqStream:SequenceStream object of the file.
...type:Class.SequenceStream
..param.seqs:StringSet of Dna5Q-Strings the sequences shall be stored in.
...type:Class.StringSet
..param.records:Unsigned int holding the number of records to be read in one block.
...type:nolink:unsigned
..returns:An integer: 1 on errors, 0 otherwise.
...type:nolink:int
..see:Class..SequenceStream
*/
int loadSeqs(seqan::SequenceStream& seqStream, seqan::StringSet<seqan::String<char> >& ids,
    seqan::StringSet<seqan::String<seqan::Dna5Q> >& seqs, unsigned records)
{
    resize(ids, records, Exact());
    resize(seqs, records, Exact());
    if (!seqan::isGood(seqStream))
    {
        std::cerr << "Error while opening the sequence-file.\n";
        return 1;
    }

    if (seqan::readBatch(ids, seqs, seqStream, records) != 0)
    {
        std::cerr << "Error while reading the sequences.\n";
        return 1;
    }
    return 0;
}
/**
.Function.loadBarcodes:
..summary:Function for loading the barcode file.
..signature:loadBarcodes(path, params)
..param.path:Path Pointer to the barcode file.
...type:nolink:char const *
..param.params:The DemultiplexingParams object the barcodes shall be written to.
...type:Class.DemultiplexingParams
..remarks:If approximate search has been selected, the loaded barcodes will be modified,
 resulting in a a bigger set, holding all barcode variations with 1 missmatch.
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int loadBarcodes(char const * path, DemultiplexingParams& params)
{
    seqan::SequenceStream bcStream(path, seqan::SequenceStream::READ);
    if (!seqan::isGood (bcStream))
    {
        std::cerr << "Error while opening file'" <<  params.barcodeFile << "'.\n";
        return 1;
    }

    if (seqan::readAll(params.barcodeIds, params.barcodes, bcStream) != 0)
    {
        std::cerr << "Error while reading barcodes from '" << params.barcodeFile << "'.\n";
        return 1;
    }
    if (params.approximate)                            //modifies the barcodes for approximate matching
    {
        buildAllVariations(params.barcodes);
    }
    resize(params.stats.groups, length(params.barcodes)+1);
    for (unsigned i = 0; i < length(params.stats.groups); ++i)
    {
        params.stats.groups[i] = 0;             //Sets the right size for the stats String and fills it with 0.
    }
    return 0;
}

/**
.Function.loadMultiplex:
..summary:Function for loading the multiplex file.
..signature:loadMultiplex(multiplexStream, params, records)
..param.multiplexStream:SequenceStream object of the file.
...type:Class.SequenceStream
..param.params:The DemultiplexingParams-object the multiplex barcodes shall be written to.
...type:Class.DemultiplexingParams
..param.records:Unsigned int holding the number of records to be read in one block.
...type:nolink:unsigned
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int loadMultiplex (seqan::SequenceStream& multiplexStream, DemultiplexingParams& params, unsigned records)
{
    resize(params.multiplex, records, Exact());
    seqan::StringSet<seqan::String<char> > ids;
    resize(ids, records);
    if (!seqan::isGood(multiplexStream))
    {
        std::cerr << "Error while opening file '" <<  params.multiplexFile << "'.\n";
        return 1;
    }
    if (seqan::readBatch(ids, params.multiplex, multiplexStream, records) != 0)
    {
        std::cerr << "Error while reading barcodes from '" << params.multiplexFile << "'.\n";
        return 1;
    }
    return 0;
}

/**
.Function.openStream:
..summary:Function for opening a new sequenceStream.
..signature:openStream(file, stream)
..param.file:CharString of the file's path.
...type:nolink:CharString
..param.stream:sequenceStream to be opened.
...type:Class.SequenceStream
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int openStream(seqan::CharString const & file, seqan::SequenceStream & stream)
{
    open(stream, seqan::toCString(file));
    if (!isGood(stream))
    {
        std::cerr << "Error while opening input file '" << file << "'.\n";
        return 1;
    }
    return 0;
}

/**
.Function.loadDemultiplexingParams:
..summary:Function for loading all demultiplexing parameters.
..signature:loadDemultiplexingParams(parser, params)
..param.parser:The SeqAn ArgumentParser object.
...type:Class.ArgumentParser
..param.params:The DemultiplexingParams object the paramters shall be written to.
...type:Class.DemultiplexingParams
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int loadDemultiplexingParams(seqan::ArgumentParser const& parser, DemultiplexingParams& params)
{
    // APPROXIMATE/EXACT MATCHING---------------------    
    params.approximate = seqan::isSet(parser, "app");
    // HARD CLIP MODE --------------------------------
    params.hardClip = seqan::isSet(parser, "hc") && !(isSet(parser, "ex"));
    // EXCLUDE UNIDENTIFIED --------------------------
    params.exclude = isSet(parser, "ex") && isSet(parser, "b");
    // RUN -------------------------------------------
    params.run = isSet(parser, "b");
    if (isSet(parser, "x"))
        params.runx = true;
    else
        params.runx = false;
    // BARCODES --------------------------------------
    if (isSet(parser, "b"))
    {
        getOptionValue(params.barcodeFile, parser, "b");
        if (loadBarcodes(toCString(params.barcodeFile), params) != 0)
            return 1;
    }
    else
    {
        resize(params.stats.groups , 1);
        params.stats.groups[0] = 0;
    }
    return 0;
}

/**
.Function.loadAdapterTrimmingParams:
..summary:Function for loading all adapter trimming parameters.
..signature:loadAdapterTrimmingParams(parser, params)
..param.parser:The SeqAn ArgumentParser object.
...type:Class.ArgumentParser
..param.params:The AdapterTrimmingParams object the paramters shall be written to.
...type:Class.AdapterTrimmingParams
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int loadAdapterTrimmingParams(seqan::ArgumentParser const& parser, AdapterTrimmingParams & params)
{
    // PAIRED-END ------------------------------
    int fileCount =  getArgumentValueCount(parser, 0);
    // Only consider paired-end mode if two files are given and user wants paired mode.
    params.paired = fileCount == 2 && !isSet(parser, "np");
    params.noAdapter = isSet(parser, "na");
    // Set run flag, depending on essential parameters.
    params.run = isSet(parser,"a") || (params.noAdapter && fileCount == 2);
    // ADAPTER SEQUENCES ----------------------------
    seqan::CharString adapterFile, id;
    // If adapter sequences are given, we read them in any case.
    if (isSet(parser, "a"))
    {
        getOptionValue(adapterFile, parser, "a");
        seqan::SequenceStream adapterStream(toCString(adapterFile));
        if (!seqan::isGood(adapterStream))
        {
            std::cerr << "Error while opening file'" << adapterFile << "'.\n";
            return 1;
        }
        int err = seqan::readRecord(id, params.adapter1, adapterStream);
            err += seqan::readRecord(id, params.adapter2, adapterStream);
        if (err != 0)
        {
            std::cerr << "Error while reading adapters from '" << adapterFile << "'.\n";
            return 1;
        }
    }
    // If they are not given, but we would need them (single-end trimming), output error.
    else if (isSet(parser, "np") || (params.noAdapter && fileCount == 1))
    {
        std::cerr << "Unpaired adapter removal requires adapter sequences.\n";
        return 1;
    }
    // TRIMMING MODE ----------------------------
    // User must fully specify mode, if he wants to. (But not specifying both is ok too.)
    if (seqan::isSet(parser, "e") != seqan::isSet(parser, "overlap"))
    {
        std::cerr << "User must define both error rate (-e) and minimum overlap (-o).\n";
        return 1;
    }
    // Read both options (if one is given, the above check guarantees that the other is too.)
    if (seqan::isSet(parser, "e")){
        // If user tried to specify alignment parameters, but didn't activate adapter trimming, warn and exit.
        if (!params.run)
        {
            std::cerr << "Adapter removal alignment parameters require adapters or --no-adapter flag.\n";
            return 1;
        }
        int o;
        int e;
        getOptionValue(o, parser, "overlap");
        getOptionValue(e, parser, "e");
        params.mode = User(o, e);
        params.mmode = E_USER;
    }
    // Otherwise use the automatic configuration.
    else
    {
        params.mode = Auto();
        params.mmode = E_AUTO;
    }
    return 0;
}

/**
.Function.loadQualityTrimmingParams:
..summary:Function for loading all quality trimming parameters.
..signature:loadQualityTrimmingParams(parser, params)
..param.parser:The SeqAn ArgumentParser object.
...type:Class.ArgumentParser
..param.params:The QualityTrimmingParams object the paramters shall be written to.
...type:Class.QualityTrimmingParams
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int loadQualityTrimmingParams(seqan::ArgumentParser const & parser, QualityTrimmingParams & params)
{
    // TRIMMING METHOD ----------------------------
    std::string method;
    getOptionValue(method, parser, "m");
    if (method == "WIN")
    {
        params.trim_mode = E_WINDOW;
    }
    else if (method == "BWA")
    {
        params.trim_mode = E_BWA;
    }
    else
    {
        params.trim_mode = E_TAIL;
    }
    // QUALITY CUTOFF ----------------------------
    if (isSet(parser, "q"))
    {
        getOptionValue(params.cutoff, parser, "q");
    }
    // MINIMUM SEQUENCE LENGTH -------------------
    getOptionValue(params.min_length, parser, "l");
    // Set run flag, depending on essential parameters (which are in a valid state at this point).
    params.run = isSet(parser, "q");
    return 0;
}

/**
.Function.loadProgramParams:
..summary:Function for loading general parameters.
..signature:loadProgramParams(parser, params)
..param.parser:The SeqAn ArgumentParser object.
...type:Class.ArgumentParser
..param.params:The ProgramParams object the paramters shall be written to.
...type:Class.ProgramParams
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int loadProgramParams(seqan::ArgumentParser const & parser, ProgramParams & params)
{
    params.fileCount = getArgumentValueCount(parser, 0);
    // Load files.
    seqan::CharString fileName1, fileName2;
    getArgumentValue(fileName1, parser, 0, 0);
    if (openStream(fileName1, params.fileStream1) != 0)
    {
        return 1;
    }
    if (params.fileCount == 2)
    {
        getArgumentValue(fileName2, parser, 0, 1);
        if (openStream(fileName2, params.fileStream2) != 0)
        {
            return 1;
        }
        if (params.fileStream1._fileFormat != params.fileStream2._fileFormat)
        {
            std::cerr << "Input files must have the same file format.\n";
            return 1;
        }
    }
    return 0;
}

/**
.Function.checkParams:
..summary:Function controlling the input paramters.
..signature:checkParams(programParams, processingParams, demultiplexingParams, adapterTrimmingParams,
 qualityTrimmingParams)
..param.programParams:The ProgramParams object.
...type:Class.ProgramParams
..param.programParams:The ProcessingParams object.
...type:Class.ProcessingParams
..param.demultiplexingParams:The DemultiplexingParams object.
...type:Class.DemultiplexingParams
..param.adapterTrimmingParams:The AdapterTrimmingParams object.
...type:Class.AdapterTrimmingParams
..param.qualityTrimmingParams:The QualityTrimmingParams object.
...type:Class.QualityTrimmingParams
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
int checkParams(ProgramParams const & programParams, ProcessingParams const & processingParams,
    DemultiplexingParams const & demultiplexingParams, AdapterTrimmingParams const & adapterTrimmingParams,
    QualityTrimmingParams & qualityTrimmingParams)
{
    // Were there options that activated at least one processing stage?

    if (!(adapterTrimmingParams.run || qualityTrimmingParams.run || demultiplexingParams.run || processingParams.runPre
        || processingParams.runPost))
    {
        std::cerr << "\nNo processing stage was specified.\n";
        return 1;
    }
    // If quality trimming was selected, check if file format includes qualities.
    if (qualityTrimmingParams.run)
    {
        if (programParams.fileStream1._fileFormat != seqan::SeqIOFileFormat_::FILE_FORMAT_FASTQ
          || ((programParams.fileCount == 2)
          && programParams.fileStream2._fileFormat != seqan::SeqIOFileFormat_::FILE_FORMAT_FASTQ))
        {
            std::cerr << "\nQuality trimming requires quality information, please specify fastq files." << std::endl;
            return 1;
        }
    }
    // If we don't demultiplex (and therefore take file names from the barcode IDs), set file names.
    if (!demultiplexingParams.run)
    {
        //std::cout << basename(toCString(fileName1)) << "\n\n";
    }
    return 0;
}

// PROGRAM STAGES ---------------------

//Preprocessing Stage
/**
.Function.preprocessingStage
..summary:Function for calling all preprocessing operations.
..signature:preprocessingStage(seqs, ids, seqsRev, idsRev, demultiplexingParams,
 processingParams, parser, generalStats)
..param.seqs:The forward reads.
...type:Class.StringSet
..param.ids:The ids associated with the forward reads.
...type:Class.StringSet
..param.seqsRev:The backward reads.
...type:Class.StringSet
..param.idsRev:The ids associated with the backward reads.
...type:Class.StringSet
..param.demultiplexingParams:The DemultiplexingParams object.
...type:Class.DemultiplexingParams
..param.processingParams:The processingParams object.
...type:Class.ProcessingParams
..param.parser:The parser object
...type:Class.ArgumentParser
..param.generalStats:The GeneralStats object holding general statistic data.
...type:Class.GeneralStats
..remarks:This function calls the desired instances of preTrim and processN.
..returns:void
..see:Function.processN
..see:Function.preTrim
*/
template<typename TSeqs, typename TIds>
void preprocessingStage(TSeqs& seqs, TIds& ids, DemultiplexingParams& demultiplexingParams,
    ProcessingParams& processingParams, seqan::ArgumentParser const & parser, GeneralStats& generalStats)
{
    if (processingParams.runPre)
    {
        //Trimming and filtering
        if (processingParams.trimLeft + processingParams.trimRight + processingParams.minLen != 0)
        {
            if (demultiplexingParams.multiplexFile != "")
            {
                preTrim(seqs, ids, demultiplexingParams.multiplex, processingParams.trimLeft,
                    processingParams.trimRight, processingParams.minLen, generalStats);
            }
            else
            {
                preTrim(seqs, ids, processingParams.trimLeft, processingParams.trimRight,
                    processingParams.minLen, generalStats);
            }
        }
        //Detecting uncalled Bases
        if (seqan::isSet(parser, "u"))
        {
            if (seqan::isSet(parser, "s"))
            {
                if (demultiplexingParams.multiplexFile != "")
                {
                    processN(seqs, ids, demultiplexingParams.multiplex, processingParams.uncalled,
                        processingParams.substitute, generalStats);
                }
                else
                {
                    processN(seqs, ids, processingParams.uncalled, processingParams.substitute, generalStats);
                }
            }
            else
            {
                if (demultiplexingParams.multiplexFile != "")
                {
                    processN(seqs, ids, demultiplexingParams.multiplex, processingParams.uncalled, generalStats);
                }
                else
                {
                    processN(seqs, ids, processingParams.uncalled, generalStats);
                }
            }
        }
    }
}

//Overload for paired-end data
template<typename TSeqs, typename TIds>
void preprocessingStage(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev,
    DemultiplexingParams& demultiplexingParams, ProcessingParams& processingParams,
    seqan::ArgumentParser const & parser, GeneralStats& generalStats)
{
    if (processingParams.runPre)
    {
        //Trimming and filtering
        if (processingParams.trimLeft + processingParams.trimRight + processingParams.minLen != 0)
        {
            if (demultiplexingParams.multiplexFile != "")
            {
                preTrim(seqs, ids, seqsRev, idsRev, demultiplexingParams.multiplex, processingParams.trimLeft,
                    processingParams.trimRight, processingParams.minLen, generalStats);
            }
            else
            {
                preTrim(seqs, ids, seqsRev, idsRev, processingParams.trimLeft, processingParams.trimRight,
                    processingParams.minLen, generalStats);
            }
        }
        //Detecting uncalled Bases
        if (seqan::isSet(parser, "u"))
        {
            if (seqan::isSet(parser, "s"))
            {
                if (demultiplexingParams.multiplexFile != "")
                {
                    processN(seqs, ids, seqsRev, idsRev, demultiplexingParams.multiplex, processingParams.uncalled,
                        processingParams.substitute, generalStats);
                }
                else
                {
                    processN(seqs, ids, seqsRev, idsRev, processingParams.uncalled, processingParams.substitute,
                        generalStats);
                }
            }
            else
            {
                if (demultiplexingParams.multiplexFile != "")
                {
                    processN(seqs, ids, seqsRev, idsRev, demultiplexingParams.multiplex, processingParams.uncalled,
                        generalStats);
                }
                else
                {
                    processN(seqs, ids, seqsRev, idsRev, processingParams.uncalled, generalStats);
                }
            }
        }
    }
}
// DEMULTIPLEXING
/**
.Function.demultiplexingStage:
..summary:Function for checking the parameters and calling the desired demultiplexing functions.
..signature:demultiplexingStage(params, seqs, ids, esaFinder, map)
..signature:demultiplexingStage(params, seqs, seqsRev, ids, idsRev, esaFinder,  map)
..param.params:The DemultiplexingParams object.
...type:Class.DemultiplexingParams
..param.seqs:StringSet of forward-reads the operations shall be performed on.
...type:Class.StringSet
..param.seqsRev:StringSet of backward-reads the operations shall be performed on.
...type:Class.StringSet
..param.ids:StringSet of IDs of the forward-reads.
...type:Class.StringSet
..param.idsRev:StringSet of IDs of the backward-reads.
...type:Class.StringSet
..param.map:Map the information about the barcode groups shall be stored in.
...type:nolink:std::map
..returns:An interger: 1 on errors, 0 otherwise.
...type:nolink:int
*/
//Version for single-end data
template <typename TSeqsVec, typename TIdsVec, typename TFinder,typename TMap>
int demultiplexingStage(DemultiplexingParams& params, TSeqsVec& seqs, TIdsVec& ids, TFinder& esaFinder, 
    TMap& map, GeneralStats& generalStats)
{
    if (!params.run)
    {
        return 0;
    }
    seqan::StringSet<seqan::String<int> > groups;
    if (params.runx &&! params.approximate)
    {
        //DEBUG_MSG(std::cout << "Demultiplexing exact multiplex single-end reads.\n");
        doAll(groups, params.multiplex, params.barcodes, esaFinder, params.stats, params.exclude);
    }
    else if (!params.approximate)
    {
        if (!check(seqs[0], ids[0], params.barcodes, generalStats))        // On Errors with barcodes return 1;
        {
            return 1;
        }
        //DEBUG_MSG("Demultiplexing exact inline single-end reads.\n");
        doAll(groups, seqs[0], params.barcodes, esaFinder, params.hardClip, params.stats, params.exclude);
    }
    else if (params.runx && params.approximate)
    {
        //DEBUG_MSG("Demultiplexing approximate multiplex single-end reads.\n");
        doAll(groups, params.multiplex, params.barcodes, esaFinder, params.stats, params.approximate, params.exclude);
    }
    else
    {
        if (!check(seqs[0], ids[0], params.barcodes, generalStats))            // On Errors with barcodes return 1;
        {
            return 1;
        }
        //DEBUG_MSG("Demultiplexing approximate inline single-end reads.\n");
        doAll(groups, seqs[0], params.barcodes, esaFinder, params.hardClip, params.stats,
            params.approximate, params.exclude);
    }
    // Saves information on how groups correspond to barcodes.
    clear(map);
    for (unsigned i = 0; i < length(groups); ++i)
    {
        if (length(groups[i]) != 0)
        {
            appendValue(map,i);
        }
    }
    // Sorting the results into the sequence- and ID Strings
    seqan::String<seqan::StringSet<seqan::String<seqan::Dna5Q> > > sortedSeqs; 
    seqan::String<seqan::StringSet<seqan::String<char> > > sortedIds;
    buildSets(seqs[0], ids[0], groups, sortedSeqs, sortedIds);
    resize(seqs, length(sortedSeqs));
    resize(ids, length(sortedIds));
    seqs = sortedSeqs;
    ids = sortedIds;
    return 0;
}
//Version for paired-end data
template <typename TSeqsVec, typename TIdsVec, typename TFinder, typename TMap>
int demultiplexingStage(DemultiplexingParams& params, TSeqsVec& seqs, TSeqsVec& seqsRev, TIdsVec& ids,
    TIdsVec& idsRev, TFinder& esaFinder, TMap& map, GeneralStats& generalStats)
{
    if (!params.run)
    {
        return 0;
    }
    seqan::StringSet<seqan::String<int> > groups;
    if (params.runx && !params.approximate)
    {
        //DEBUG_MSG("Demultiplexing exact multiplex paired-end reads.\n");
        doAll(groups, params.multiplex, params.barcodes, esaFinder, params.stats, params.exclude);
    }
    else if (!params.approximate)
    {
        if (!check(seqs[0], ids[0], seqsRev[0], idsRev[0] , params.barcodes, generalStats))
        {
            return 1;            // On Errors with barcodes return 1;            
        }
        //DEBUG_MSG("Demultiplexing exact inline paired-end reads.\n");
        doAll(groups, seqs[0], params.barcodes, esaFinder, params.hardClip, params.stats, params.exclude);
    }
    else if (params.runx && params.approximate)
    {
        //DEBUG_MSG("Demultiplexing approximate multiplex paired-end reads.\n");
        doAll(groups, params.multiplex, params.barcodes, esaFinder, params.stats, params.approximate, params.exclude);
    }
    else
    {
        if (!check(seqs[0], ids[0], seqsRev[0], idsRev[0], params.barcodes, generalStats))
        {
            return 1;           // On Erros with barcodes return 1;
        }
        //DEBUG_MSG("Demultiplexing approximate inline paired-end reads.\n");
        doAll(groups, seqs[0], params.barcodes, esaFinder, params.hardClip, params.stats,
            params.approximate, params.exclude);
    }
    // Saves information on how groups correspond to barcodes.
    clear(map);
    for (unsigned i = 0; i < length(groups); ++i)
    {
        if (length(groups[i]) != 0) 
        {
           appendValue(map,i);
        }
    }
    // Sorting the results into the sequence- and ID Strings
    seqan::String<seqan::StringSet<seqan::String<seqan::Dna5Q> > > sortedSeqs; 
    seqan::String<seqan::StringSet<seqan::String<seqan::Dna5Q> > > sortedSeqsRev; 
    seqan::String<seqan::StringSet<seqan::String<char> > > sortedIds;
    seqan::String<seqan::StringSet<seqan::String<char> > > sortedIdsRev;
    buildSets(seqs[0], seqsRev[0] ,ids[0], idsRev[0], groups, sortedSeqs, sortedSeqsRev, sortedIds, sortedIdsRev);
    resize(seqs, length(sortedSeqs));
    resize(seqsRev, length(sortedSeqs));
    resize(ids, length(sortedIds));
    resize(idsRev, length(sortedIds));
    seqs = sortedSeqs;
    seqsRev = sortedSeqsRev;
    ids = sortedIds;
    idsRev = sortedIdsRev;
    return 0;
}

// ADAPTER TRIMMING
/**
.Function.adapterTrimmingStage:
..summary:Function for checking the parameters and calling the desired adapter trimming functions.
..signature:adapterTrimmingStage(params, seqSet1, idSet1, tagOpt)
..signature:adapterTrimmingStage(params, seqSet1, idSet1, seqSet2, idSet2, tagOpt)
..param.params:The AdapterTrimmingParams object.
...type:Class.AdapterTrimmingParams
..param.seqSet1:The set of forward reads.
...type:Class.StringSet
..param.idSet1:The set of IDs associated with the forward reads.
...type:Class.StringSet
..param.seqSet2:The set of backward reads.
...type:Class.StringSet
..param.idSet2:The set of IDs associated with the backward reads.
...type:Class.StringSet
..param.tagOpt:Bool incidating that reads with removed adapter shall be tagged (i.e their IDs).
...type:nolink:bool
..returns:void
 */
//Version for single-end data
template <typename TSeqs, typename TIds>
void adapterTrimmingStage(AdapterTrimmingParams& params, TSeqs& seqSet, TIds& idSet, bool tagOpt)
{
    if (!params.run)
    {
        return;
    }
    //DEBUG_MSG("Trimming single-end adapters.\n");
    switch(params.mmode)
    {
        case E_USER:
        {
            for (unsigned i = 0; i < length(seqSet); ++i)
            {
                stripAdapterBatch(seqSet[i], idSet[i], params.adapter2, (User&) params.mode, params.stats, false, tagOpt);
            }
            break;
        }
        case E_AUTO:
        {
            for (unsigned i = 0; i < length(seqSet); ++i)
            {
                stripAdapterBatch(seqSet[i], idSet[i], params.adapter2, (Auto&) params.mode, params.stats, false, tagOpt);
            }
            break;
        }
    }
}

//Overload for paired-end data
template <typename TSeqs, typename TIds>
void adapterTrimmingStage(AdapterTrimmingParams& params, TSeqs& seqSet1, TIds& idSet1, TSeqs& seqSet2, TIds& idSet2, bool tagOpt)
{
    if (!params.run)
    {
        return;
    }
    typename seqan::Iterator<TSeqs, Rooted>::Type it1, it2;
    for (unsigned i = 0; i < length(seqSet1); ++i)
    {
        if (!params.paired)
        {
            //DEBUG_MSG("Trimming paired-end adapters in single-end mode.\n");
            switch(params.mmode)
            {
                case E_USER:
                {
                    stripAdapterBatch(seqSet1[i], idSet1[i], params.adapter2, (User&)params.mode, params.stats, tagOpt);
                    stripReverseAdapterBatch(seqSet2[i], idSet2[i], params.adapter1, (User&)params.mode, params.stats, tagOpt);
                    break;
                }
                case E_AUTO:
                {
                    stripAdapterBatch(seqSet1[i], idSet1[i], params.adapter2, (Auto&)params.mode, params.stats, tagOpt);
                    stripReverseAdapterBatch(seqSet2[i], idSet2[i], params.adapter1, (Auto&)params.mode, params.stats, tagOpt);
                    break;
                }
            }
        }
        else
        {
            //DEBUG_MSG("Trimming paired-end adapters.\n");
            stripPairBatch(seqSet1[i], idSet1[i], seqSet2[i], idSet2[i], params.stats, tagOpt);
        }
    }
}

// QUALITY TRIMMING
/**
.Function.qualityTrimmingStage:
..summary:Function for checking the parameters and calling the desired quality trimming functions.
..signature:qualityTrimmingStage(params, idSet, seqSet)
..signature:qualityTrimmingStage(params, idSet, seqSet, idSet2, seqSet2)
..param.params:The QualityTrimmingParams object.
...type:Class.QualityTrimmingParams
..param.idSet:The StringSet of IDs of forward-reads to be trimmed.
...type:Class.StringSet
..param.seqSet:The StringSet of forward-reads to be trimmed.
...type:Class.StringSet
..param.idSet2:The StringSet of IDs of backward-reads to be trimmed.
...type:Class.StringSet
..param.seqSet2:The StringSet of backward-reads to be trimmed.
...type:Class.StringSet
..returns:void
 */
//Version for single-ende data
template <typename TIds, typename TSeqs>
void qualityTrimmingStage(QualityTrimmingParams& params, TIds& idSet, TSeqs& seqSet, bool tagOpt)
{
    if (params.run)
    {
        //DEBUG_MSG("Trimming qualities.\n");
        // Templates don't support runtime polymorphism, so code paths for all possibilities.
        switch (params.trim_mode)
        {
            case E_WINDOW:
            {
                for (unsigned i=0; i < length(idSet); ++i)
                {
                    trimBatch(seqSet[i], idSet[i], params.cutoff, Mean(5), tagOpt);
                }
                break;
            }
            case E_BWA:
                {
                for (unsigned i=0; i < length(idSet); ++i)
                {
                    trimBatch(seqSet[i], idSet[i],params.cutoff, BWA(), tagOpt);
                }
                break;
                }
            case E_TAIL:
            {
                for (unsigned i=0; i < length(idSet); ++i)
                {
                    trimBatch(seqSet[i], idSet[i], params.cutoff, Tail(), tagOpt);
                }
                break;
            }
        }
    }
    for (unsigned i=0; i < length(idSet); ++i)
    {
        dropReads(idSet[i], seqSet[i], params.min_length, params.stats);
    }
}

//Overload for paired-end data
template <typename TIds, typename TSeqs>
void qualityTrimmingStage(QualityTrimmingParams& params, TIds& idSet1, TSeqs& seqSet1, TIds& idSet2, TSeqs& seqSet2, bool tagOpt)
{
    if (params.run)
    {
        //DEBUG_MSG("Trimming (pair) qualities.\n");
        // Templates don't support runtime polymorphism, so code paths for all possibilities.
        switch (params.trim_mode)
        {
            case E_WINDOW:
            {
                for (unsigned i=0; i < length(idSet1); ++i)
                {
                    trimPairBatch(seqSet1[i], idSet1[i], seqSet2[i], idSet2[i], params.cutoff, Mean(5), tagOpt);
                }
                break;
            }
            case E_BWA:
            {
                for (unsigned i=0; i < length(idSet1); ++i)
                {
                    trimPairBatch(seqSet1[i], idSet1[i], seqSet2[i], idSet2[i], params.cutoff, BWA(), tagOpt);
                }
                break;
            }
            case E_TAIL:
            {
                for (unsigned i=0; i < length(idSet1); ++i)
                {
                    trimPairBatch(seqSet1[i], idSet1[i], seqSet2[i], idSet2[i], params.cutoff, Tail(), tagOpt);
                }
                break;
            }
        }
    }
    for (unsigned i = 0; i < length(idSet1); ++i)
    {
        dropReads(idSet1[i], seqSet1[i], idSet2[i], seqSet2[i], params.min_length, params.stats);
    }
}
//Postprocessing
template<typename TSeqSet, typename TIdSet>
void postprocessingStage(TSeqSet& seqSet, TIdSet& idSet, ProcessingParams& params, GeneralStats& stats)
{
    if (params.runPost)
    {
        if ((params.finalMinLength != 0) && (params.finalLength == 0))
        {
            for (unsigned i = 0; i < length(seqSet); ++i)
            {
                preTrim(seqSet[i], idSet[i], 0, 0, params.finalMinLength, stats);
            }
        }
        else if (params.finalLength != 0)
        {
            for (unsigned i = 0; i < length(seqSet); ++i)
            {
                trimTo(seqSet[i], idSet[i], params.finalLength, stats); 
            }
        }
    }
}
//Overload for paired-end data
template<typename TSeqSet, typename TIdSet>
void postprocessingStage(TSeqSet& seqSet, TIdSet& idSet, TSeqSet& seqSet2, TIdSet& idSet2, ProcessingParams& params, GeneralStats& stats)
{
    if (params.runPost)
    {
        if ((params.finalMinLength != 0) && (params.finalLength == 0))
        {
            for (unsigned i = 0; i < length(seqSet); ++i)
            {
                preTrim(seqSet[i], idSet[i], seqSet2[i], idSet2[i], 0, 0, params.finalMinLength, stats);
            }
        }
        else if (params.finalLength != 0)
        {
            for (unsigned i = 0; i < length(seqSet); ++i)
            {
                trimTo(seqSet[i], idSet[i], seqSet2[i], idSet2[i], params.finalLength, stats);
            }
        }
    }    
}

// END PROGRAM STAGES ---------------------
/**
.Function.printStatistics:
..summary:Function for calculating and printing the statistics.
..signature:printStatistics(programParams, demultiplexingParams, adapterTrimmingParams, qualityTrimmingParams)
..param.programParams:The ProgramParams object.
...type:Class.ProgramParams
..param.demultiplexingParams:The DemultiplexingParams object.
...type:Class.DemultiplexingParams
..param.adapterTrimmingParams:The AdapterTrimmingParams object.
...type:Class.AdapterTrimmingParams
..param.qualityTrimmingParams:The QualityTrimmingParams object.
...type:Class.QualityTrimmingParams
..returns:void
 */
void printStatistics(ProgramParams& programParams, GeneralStats& generalStats, DemultiplexingParams& demultiplexParams,
                AdapterTrimmingParams& adapterParams, QualityTrimmingParams& qualityParams, bool timing)
{
    bool paired = programParams.fileCount == 2;
    bool adapter = adapterParams.run;
    int read_factor = (1+paired);
    std::cout << std::endl;
    std::cout << "\r\rRead statistics\n";
    std::cout << "===============\n";
    std::cout << "Reads processed:\t" << read_factor * programParams.readCount;
    if (paired) 
    {
        std::cout << " (2 * " << programParams.readCount << ")";
    }
    std::cout << std::endl;
    double dropped = qualityParams.stats.dropped_1 + qualityParams.stats.dropped_2 + generalStats.removedSeqs
        + generalStats.removedSeqsShort + (demultiplexParams.exclude * demultiplexParams.stats.groups[0]);
    double dropped_d = (generalStats.removedSeqs + generalStats.removedSeqsShort) / read_factor;
    std::cout << "  Reads dropped:\t" << dropped << "\t(" << std::setprecision(3) 
        << dropped / (double(programParams.readCount * read_factor)) * 100 << "%)\n";
    if (dropped != 0.0)
    {
        if (generalStats.removedSeqs != 0)
        {
            std::cout << "  Due to N content:\t" << generalStats.removedSeqs << "\t(" << std::setprecision(3)
                << double(generalStats.removedSeqs) / dropped * 100 << "%)\n";
        }
        if (generalStats.removedSeqsShort != 0)
        {
            std::cout << "  Due to shortness:\t" << dropped - generalStats.removedSeqs << "\t(" 
                << std::setprecision(3) << double(generalStats.removedSeqsShort) / dropped * 100 << "%)\n";
        }
    }
    if (generalStats.uncalledBases != 0)
        {
            std::cout << "  Remaining uncalled (or masked) bases: " << generalStats.uncalledBases << "\n";
        }
    //Statistics for Demultiplexing
    if (demultiplexParams.run)    
    {
        std::cout << "\nBarcode Demultiplexing statistics\n";
        std::cout << "=================================\n";
        if (demultiplexParams.approximate)                //Calculates the original barcode ID from the modified barcodes
        {
            int val = length(demultiplexParams.barcodes[0])*5; //Number of barcode variations for one barcode
            seqan::String<unsigned> newgroups;
            resize(newgroups, (length(demultiplexParams.stats.groups)-1)/val+1);
            newgroups[0] = demultiplexParams.stats.groups[0];
            for (unsigned i = 1; i < length(newgroups); ++i)
            {
                unsigned sum = 0;
                for (unsigned j = i*val-val+1; j <= i*val; ++j)
                {
                    sum += demultiplexParams.stats.groups[j];
                }
                newgroups[i] = sum;
            }
            clear(demultiplexParams.stats.groups);
            demultiplexParams.stats.groups = newgroups;
            resize(newgroups,0);
        }
        unsigned usedBarcodes = 0;
        for (unsigned i = 0; i < length(demultiplexParams.stats.groups); ++i)
        {
            if (demultiplexParams.stats.groups[i] != 0)
            {
                ++usedBarcodes;
            }
        }
        if (usedBarcodes != 0) //in case 0 barcodes were used
        {
            --usedBarcodes;        //correction for unidentified group    
        }
        unsigned barcodesTotal = length(demultiplexParams.barcodes);
        if (demultiplexParams.approximate)
        {
            barcodesTotal = barcodesTotal/(length(demultiplexParams.barcodes[0])*5);
        }
        std::cout << "Barcodes used: " << usedBarcodes << "/" << barcodesTotal << "\t\t(" 
            << std::setprecision(3) << (double)usedBarcodes/(double)barcodesTotal*100 << "%)\n";
        std::cout << "Reads per barcode:\n";
        std::cout << "Unidentified:\t" << read_factor * demultiplexParams.stats.groups[0];
        if (programParams.readCount - dropped_d != 0)
        {
            std::cout  << "\t\t(" << std::setprecision(3) << (double)demultiplexParams.stats.groups[0] /
                ((double)programParams.readCount - dropped_d) * 100 << "%)";
        }  
        std::cout << "\n";
        for (unsigned i = 1; i < length(demultiplexParams.stats.groups); ++i)
        {
            std::cout << demultiplexParams.barcodeIds[i-1]<<":\t" << read_factor * demultiplexParams.stats.groups[i];
            if (programParams.readCount - dropped_d != 0)
            {
                std::cout  << "\t\t(" << std::setprecision(3) << (double)demultiplexParams.stats.groups[i] /
                    ((double)programParams.readCount - dropped_d) * 100 << "%)";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }
    std::cout << "File statistics\n";
    std::cout << "===============\n";
    // How many reads are left.
    int survived1 = programParams.readCount - qualityParams.stats.dropped_1
        - ((generalStats.removedSeqs +  generalStats.removedSeqsShort
        + (demultiplexParams.exclude * demultiplexParams.stats.groups[0])) / read_factor);
    int survived2 = programParams.readCount - qualityParams.stats.dropped_2
        - ((generalStats.removedSeqs + generalStats.removedSeqsShort
        + (demultiplexParams.exclude * demultiplexParams.stats.groups[0])) / read_factor);
    // In percentage points.
    double surv_proc_1 = (double)survived1 / (double)programParams.readCount * 100;
    double surv_proc_2 = (double)survived2 / (double)programParams.readCount * 100;
    std::cout << "File 1:\n";
    std::cout << "-------\n";
    std::cout << "  Surviving: " << survived1 << "/" << programParams.readCount
              << " (" << std::setprecision(3) << surv_proc_1 << "%)\n";
    if (adapter)
    {
            std::cout << "   Adapters: " << adapterParams.stats.a1count << "\n";
    }
    std::cout << std::endl;
    if (paired)
    {
        std::cout << "File 2:\n";
        std::cout << "-------\n";
        std::cout << "  Surviving: " << survived2 << "/"
                  << programParams.readCount << " (" << surv_proc_2 << "%)\n";
        if (adapter)
        {
            std::cout << "   Adapters: " << adapterParams.stats.a2count << "\n";
        }
        std::cout << std::endl;
    }
    if (adapter && (adapterParams.stats.a2count + adapterParams.stats.a2count != 0))
    {
        int mean = adapterParams.stats.overlapSum/(read_factor*programParams.readCount);
        std::cout << "Adapter sizes:\n";
        std::cout << "Min: " << adapterParams.stats.minOverlap << ",  Mean: " << mean
                << ", Max: " << adapterParams.stats.maxOverlap << "\n\n";
    }
    // Print processing and IO time. IO is (approx.) the whole loop without the processing part.
    if (timing)
    {
        std::cout << "Time statistics:\n";
        std::cout << "==================\n";
        std::cout << "Processing time: " << std::setw(5) << programParams.processTime << " seconds.\n";
        std::cout << "       I/O time: " << std::setw(5) << programParams.ioTime << " seconds.\n";
        std::cout << std::endl;
    }
}

// END FUNCTION DEFINITIONS ---------------------------------------------

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
// Program entry point.

int flexbarMain(int argc, char const ** argv)
{
    SEQAN_PROTIMESTART(loopTime);
    seqan::ArgumentParser parser = initParser();

    // Additional checks
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if input was successfully parsed.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check if one or two input files (single or paired-end) were given.
    int fileCount = getArgumentValueCount(parser, 0);
    if (!(fileCount == 1 || fileCount == 2)){
        printShortHelp(parser);
        return 1;
    }

    //--------------------------------------------------
    // Parse general parameters.
    //--------------------------------------------------
    GeneralStats generalStats;

    unsigned records;
    getOptionValue(records, parser, "r");

    seqan::CharString output;
    getOptionValue(output, parser, "output");

    int threads;
    getOptionValue(threads, parser, "tnum");
    omp_set_num_threads(threads);

    //--------------------------------------------------
    // Parse pre- and postprocessing parameters.
    //--------------------------------------------------
    ProcessingParams processingParams;

    if(flexiProgram == FILTERING || flexiProgram == ALL_STEPS)
    {
        seqan::CharString substituteString;
        getOptionValue(substituteString, parser, "s");
        processingParams.substitute = substituteString[0];
        getOptionValue(processingParams.uncalled, parser, "u");
        if (seqan::isSet(parser, "s") && (!seqan::isSet(parser, "u") || processingParams.uncalled == 0))
        {
            std::cout << "\nWarning: Substitute for uncalled bases is set, but will be ineffective unless number of allowed"
                << " uncalled bases is set to a value higher than 0 by using the paramter -u [VALUE]\n";
            return 1;
        }
        getOptionValue(processingParams.trimLeft, parser, "tl");
        getOptionValue(processingParams.trimRight, parser, "tr");
        getOptionValue(processingParams.minLen, parser, "ml");
        if (seqan::isSet(parser, "fl"))
        {
            getOptionValue(processingParams.finalLength , parser, "fl");
        }
        if (seqan::isSet(parser, "fm"))
        {
            getOptionValue(processingParams.finalMinLength , parser, "fm");
        }
        processingParams.runPre = ((processingParams.minLen + processingParams.trimLeft + processingParams.trimRight != 0)
            || isSet(parser, "u"));
        processingParams.runPost = (processingParams.finalLength + processingParams.finalMinLength != 0);
        if(flexiProgram == FILTERING)
        {
            processingParams.runPre = true;
            processingParams.runPost = true;
        }

    }
    //--------------------------------------------------
    // Parse demultiplexing parameters.
    //--------------------------------------------------

    DemultiplexingParams demultiplexingParams;
    seqan::SequenceStream multiplexStream;    //Initialising the SequenceStream for the multiplex file

    if(flexiProgram == DEMULTIPLEXING || flexiProgram == ALL_STEPS)
    {
        getOptionValue(demultiplexingParams.multiplexFile, parser, "x");
        if (isSet(parser, "x"))
            open(multiplexStream, seqan::toCString(demultiplexingParams.multiplexFile));

        if (loadDemultiplexingParams(parser, demultiplexingParams) != 0)
            return 1;

        if(flexiProgram == DEMULTIPLEXING)
            demultiplexingParams.run = true;
    }
    else
    {
        resize(demultiplexingParams.stats.groups , 1, 0u);
    }

    //--------------------------------------------------
    // Process Barcodes
    //--------------------------------------------------

    seqan::Index<seqan::StringSet<seqan::String<seqan::Dna5> >,
        seqan::IndexEsa<> > indexSet(demultiplexingParams.barcodes);
    seqan::Finder<seqan::Index<seqan::StringSet<seqan::String<seqan::Dna5> >,
        seqan::IndexEsa<> > > esaFinder(indexSet);
    if(flexiProgram == DEMULTIPLEXING || flexiProgram == ALL_STEPS)
    {
        if (isSet(parser, "b"))
        {
            indexRequire(indexSet, FibreSA());
        }
        else
        {
            if(flexiProgram == DEMULTIPLEXING)
            {
               std::cerr << "No Barcodefile was provided." << std::endl;
               return 1;
            }
        }

    }

    //--------------------------------------------------
    // Parse quality trimming parameters.
    //--------------------------------------------------

    QualityTrimmingParams qualityTrimmingParams;
    if(flexiProgram == QUALITY_CONTROL || flexiProgram == ALL_STEPS)
    {
        if ( loadQualityTrimmingParams(parser, qualityTrimmingParams) != 0 )
            return 1;

        if(flexiProgram == QUALITY_CONTROL)
            qualityTrimmingParams.run = true;
    }

    //--------------------------------------------------
    // Parse adapter trimming parameters.
    //--------------------------------------------------

    AdapterTrimmingParams adapterTrimmingParams;
    bool tagOpt = false;

    if(flexiProgram == ADAPTER_REMOVAL || flexiProgram == QUALITY_CONTROL|| flexiProgram == ALL_STEPS)
    {
        if(flexiProgram == ADAPTER_REMOVAL || flexiProgram == ALL_STEPS)
        {
            if (loadAdapterTrimmingParams(parser, adapterTrimmingParams) != 0)
            {
                return 1;
            }
        }
        tagOpt = seqan::isSet(parser, "t");
        if (tagOpt && !(qualityTrimmingParams.run || adapterTrimmingParams.run))
        {
            std::cout << "\nWarning: Tag option was selected without selecting adapter removal or quality-trimming stage.\n";
                return 1;
        }
        if(flexiProgram == ADAPTER_REMOVAL)
            adapterTrimmingParams.run = true;
    }

    //--------------------------------------------------
    // General program parameters and additional checks.
    //--------------------------------------------------

    ProgramParams programParams;
    if (loadProgramParams(parser, programParams) != 0)
        return 1;

    if (checkParams(programParams, processingParams, demultiplexingParams, adapterTrimmingParams, qualityTrimmingParams) != 0)
        return 1;

    //--------------------------------------------------
    // Processing
    //--------------------------------------------------

    // Prepare output stream object and initial mapping from StringSets to files.
    bool noQuality = (programParams.fileStream1._fileFormat == seqan::SeqIOFileFormat_::FILE_FORMAT_FASTA);
    if (isSet(parser, "nq"))
    {
        noQuality = true;
    }
    OutputStreams outputStreams(output, noQuality);
    seqan::String<unsigned> map;
    resize(map,1);
    map[0] = 0;
    // Output additional Information on selected stages:
    if (!isSet(parser, "ni"))
    {
        seqan::CharString filename1, filename2, out;
        getOptionValue(out, parser, "output");
        std::cout << "Overview:\n";
        std::cout << "=========\n";
        getArgumentValue(filename1, parser, 0, 0);
        std::cout << "Forward-read file: " << filename1 << "\n";
        if (programParams.fileCount == 2)
        {
            getArgumentValue(filename2, parser, 0, 1);
            std::cout << "Backward-read file: " << filename2 << "\n";
        }
        else
        {
            std::cout << "Backward-read file: NONE\n";
        }
        if (isSet(parser, "output"))
        {
            std::cout << "Output-path: " << out <<"\n";
        }
        else
        {
            std::cout << "Output-path: Working Directory\n";
        }
        std:: cout << "\n"; 
        std::cout << "General Options:\n";
        std::cout << "\tReads per block: " << records * ((programParams.fileCount == 2) + 1)<< "\n";
        /*
        if (isSet(parser, "c"))
        {
            std::cout << "\tCompress output: YES" << "\n";
        }
        else
        {
            std::cout << "\tCompress output: NO" << "\n";
        }
        */
        if (isSet(parser, "nq") && (programParams.fileStream1._fileFormat != seqan::SeqIOFileFormat_::FILE_FORMAT_FASTA))
        {
            std::cout << "\tForce no-quality output: YES\n";
        }
        else if (programParams.fileStream1._fileFormat != seqan::SeqIOFileFormat_::FILE_FORMAT_FASTA)
        {
            std::cout << "\tForce no-quality output: NO\n";
        }
        std::cout << "\tNumber of threads: " << threads << "\n";
        if(flexiProgram == ADAPTER_REMOVAL || flexiProgram == QUALITY_CONTROL|| flexiProgram == ALL_STEPS)
        {
            if (isSet(parser, "t"))
            {
                std::cout << "\tTag quality-trimmed or adapter-removed reads: YES\n";
            }
            else if (adapterTrimmingParams.run||qualityTrimmingParams.run)
            {
                std::cout << "\tTag quality-trimmed or adapter-removed reads: NO\n";
            }
            std::cout << "\n";
        }
    
        if(flexiProgram == FILTERING || flexiProgram == ALL_STEPS)
        {
            std::cout << "Pre-, Postprocessing and Filtering:\n";
            std::cout << "\tPre-trim 5'-end length: " << processingParams.trimLeft << "\n";
            std::cout << "\tPre-trim 3'-end length: " << processingParams.trimRight << "\n";
            std::cout << "\tExclude reads shorter than: " << processingParams.minLen << "\n";
            if (isSet(parser, "u"))
            {
                std::cout << "\tAllowed uncalled bases per sequence: " << processingParams.uncalled << "\n";
            }
            if (isSet(parser, "s"))
            {
                std::cout << "\tSubstitute for uncalled bases: " << processingParams.substitute << "\n";
            }
            if (isSet(parser, "fm") && (processingParams.finalLength == 0))
            {
                std::cout << "\tMinimum sequence length after COMPLETE workflow: " << processingParams.finalMinLength << "\n";
            }
            else if (processingParams.finalLength != 0)
            {
                    std::cout << "\tTrim sequences after COMPLETE workflow to length: " <<processingParams.finalLength<< "\n";
            }
            std::cout << "\n";   
        }
        if (demultiplexingParams.run)
        {
            std::cout << "Barcode Demultiplexing:\n";
            std::cout << "\tBarcode file: " << demultiplexingParams.barcodeFile << "\n";
            if (demultiplexingParams.runx)
            {
                std::cout << "\tMultiplex barcode file: " << demultiplexingParams.multiplexFile << "\n";
            }
            else
            {
                std::cout << "\tMultiplex barcodes:  NO" << demultiplexingParams.multiplexFile << "\n";
            }
            if (demultiplexingParams.approximate)
            {
                std::cout << "\tApproximate matching: YES\n";
            }
            else
            {
                std::cout << "\tApproximate matching: NO\n";
            }
            if (demultiplexingParams.hardClip && !demultiplexingParams.runx)
            {
                std::cout << "\tHardClip mode: YES\n";
            }
            else
            {
                std::cout << "\tHardClip mode: NO\n";
            }
            if (demultiplexingParams.exclude)
            {
                std::cout << "\tExclude unidentified sequences: YES\n";
            }
            else
            {
                std::cout << "\tExclude unidentified sequences: NO\n";
            }
            std::cout << "\n";
        }
        if (adapterTrimmingParams.run)
        {
            std::cout << "Adapter Removal:\n";
            if (isSet(parser, "a"))
            {
                seqan::CharString adapter;
                getOptionValue(adapter, parser, "a");
                std::cout << "\tAdapter file: " << adapter << "\n";
            }
            else
            {
                std::cout << "\tAdapter file: NONE\n";
            }
            if (adapterTrimmingParams.noAdapter)
            {
                std::cout << "\tDon't use Adapter file: YES\n";
            }
            else
            {
                std::cout << "\tDon't use Adapter file: NO\n";
            }
            if (isSet(parser, "np"))
            {
                std::cout << "\tForce single-end method: YES\n";
            }
            else
            {
                std::cout << "\tForce single-end method: NO\n";
            }
            if (isSet(parser, "e"))
            {
                unsigned e, o;
                getOptionValue(e, parser, "e");
                getOptionValue(o, parser, "overlap");
                std::cout << "\tAllowed mismatches: " << e << "\n";
                std::cout << "\tMinimum overlap " << o << "\n";
            }
            std::cout << "\n";
        }
        if (qualityTrimmingParams.run)
        {
            std::cout << "Quality Trimming:\n";
            std::cout << "\tMinumum PHRED-Quality: " << qualityTrimmingParams.cutoff << "\n";
            std::string method;
            getOptionValue(method, parser, "m");
            std::cout << "\tMethod: " << method << "\n";
            if (isSet(parser, "l"))
            {
                std::cout << "\tMinimum length after trimming: " << qualityTrimmingParams.min_length << "\n";
            }
            std::cout << "\n";
        }
    }
    // Start processing. Different functions are needed for one or two input files.
    std::cout << "\nProcessing reads...\n" << std::endl;


    if (fileCount == 1)
    {
        if (!demultiplexingParams.run)
            outputStreams.addStream("", 0);

        seqan::String<seqan::StringSet<seqan::CharString> > idSet;
        seqan::String<seqan::StringSet<Dna5QString> > seqSet;

        while (!atEnd(programParams.fileStream1))
        {
            clear(idSet);
            clear(seqSet);

            appendValue(idSet, seqan::StringSet<seqan::CharString>());
            appendValue(seqSet, seqan::StringSet<Dna5QString>());

            if (loadSeqs(programParams.fileStream1, idSet[0], seqSet[0], records) == 0)
            {
                programParams.readCount += length(idSet[0]);
                SEQAN_PROTIMESTART(processTime);            // START of processing time.

                if (demultiplexingParams.runx && (loadMultiplex(multiplexStream, demultiplexingParams , records) !=0))
                {
                    return 1;
                }


                // Preprocessing and Filtering
                preprocessingStage(seqSet[0], idSet[0], demultiplexingParams,
                    processingParams, parser, generalStats);

                // Demultiplexing
                if (demultiplexingStage(demultiplexingParams, seqSet, idSet,
                                    esaFinder, map, generalStats) != 0)
                {
                    return 1;
                }

                // Adapter trimming
                adapterTrimmingStage(adapterTrimmingParams, seqSet, idSet, tagOpt);

                // Quality trimming
                qualityTrimmingStage(qualityTrimmingParams, idSet, seqSet, tagOpt);

                // Postprocessing
                postprocessingStage(seqSet, idSet, processingParams, generalStats);
                programParams.processTime += SEQAN_PROTIMEDIFF(processTime);    // END of processing time.

                // Append to output file.
                outputStreams.writeSeqs(idSet, seqSet, map, demultiplexingParams.barcodeIds);
                // Information
                std::cout << "\r" << programParams.readCount;
            }
            else
            {
                return 1;
            }
        }
     }
     else
     {
        if (!demultiplexingParams.run)
            outputStreams.addStreams("", "", 0);

        seqan::String<seqan::StringSet<seqan::CharString> > idSet1, idSet2;
        seqan::String<seqan::StringSet<Dna5QString> > seqSet1, seqSet2;

        while (!(atEnd(programParams.fileStream1) || atEnd(programParams.fileStream2)))
        {
            clear(idSet1); 
            clear(idSet2);
            clear(seqSet1); 
            clear(seqSet2);

            appendValue(idSet1, seqan::StringSet<seqan::CharString>());
            appendValue(seqSet1, seqan::StringSet<Dna5QString>());
            appendValue(idSet2, seqan::StringSet<seqan::CharString>());
            appendValue(seqSet2, seqan::StringSet<Dna5QString>());

            if (loadSeqs(programParams.fileStream1, idSet1[0], seqSet1[0], records) == 0 &&
                loadSeqs(programParams.fileStream2, idSet2[0], seqSet2[0], records) == 0)
            {
                programParams.readCount += length(idSet1[0]);
                SEQAN_PROTIMESTART(processTime); // START of processing time.

                if (demultiplexingParams.runx && (loadMultiplex(multiplexStream, demultiplexingParams , records) !=0))
                {
                    return 1;
                }
                // Generall Processing
                preprocessingStage(seqSet1[0], idSet1[0], seqSet2[0], idSet2[0], demultiplexingParams,
                    processingParams, parser, generalStats);

                // Demultiplexing.
                if (demultiplexingStage(demultiplexingParams, seqSet1, seqSet2, idSet1, idSet2,
                                    esaFinder, map, generalStats) != 0)
                    return 1;

                // Adapter trimming.
                adapterTrimmingStage(adapterTrimmingParams, seqSet1, idSet1, seqSet2, idSet2, tagOpt);

                // Quality trimming.
                qualityTrimmingStage(qualityTrimmingParams, idSet1, seqSet1, idSet2, seqSet2, tagOpt);

                // Postprocessing
                postprocessingStage(seqSet1, idSet1, seqSet2, idSet2, processingParams, generalStats);
                programParams.processTime += SEQAN_PROTIMEDIFF(processTime); // END of processing time.

                // Append to output file.
                outputStreams.writeSeqs(idSet1, seqSet1, idSet2, seqSet2, map, demultiplexingParams.barcodeIds);
                // Information
                std::cout << "\r" << 2*programParams.readCount;
            }
            else return 1;
        }
    }
    double loop = SEQAN_PROTIMEDIFF(loopTime);
    programParams.ioTime = loop - programParams.processTime;

    printStatistics(programParams, generalStats, demultiplexingParams, adapterTrimmingParams, qualityTrimmingParams, !isSet(parser, "ni"));

    return 0;
}
