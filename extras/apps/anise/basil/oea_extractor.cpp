// ==========================================================================
//                              OEA EXTRACTOR
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Given a coordinate sorted BAM file, extract OEA records in the format
// required by NovelSeq.
// ==========================================================================

#include <stdexcept>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>

#include "progress_indicator.h"
#include "library_info.h"
#include "bam_filter_pipeline.h"

#include "shared/orphan_extractor.h"

typedef std::runtime_error OeaExtractorAppException;

struct OeaExtractorOptions
{
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity { 1 };
    // Path to input SAM/BAM file.
    seqan::CharString inputFile;
    // Path to output OEA SAM/BAM file.
    seqan::CharString outputOeaFile;
    // Path to output orphan SAM file.
    seqan::CharString outputOrphansFile;

    // Number of records to use for library estimation.
    unsigned autoLibraryNumRecords { 100*1000 };

    // Information about the BAM library.
    BamLibraryInfo libraryInfo;

    OeaExtractorOptions() = default;
};

class OeaExtractorApp
{
public:
    // Run extractor app.
    int run(int argc, char const ** argv);

private:
    // Parse command line.
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv);
    // Run the filtration pipeline and write out resulting mapping file.
    void runPipeline();

    // OEA Extraction Options.
    OeaExtractorOptions options;

    // Input and output BAM stream.
    seqan::BamStream outStream;
};

seqan::ArgumentParser::ParseResult OeaExtractorApp::parseCommandLine(int argc, char const ** argv)
{
    // Declare ArgumentParser variable.
    seqan::ArgumentParser parser("oea_extractor");

    // Set short description, version, and date.
    setShortDescription(parser, "OEA Extractor");
#ifdef SEQAN_REVISION
        setVersion(parser, "0.2.0-beta.1 [" + std::string(SEQAN_REVISION) + "]");
#else
        setVersion(parser, "0.2.0-beta.1");
#endif
#ifdef SEQAN_DATE
        setDate(parser, SEQAN_DATE);
#endif
    // Set usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-im\\fP \\fIIN.{sam,bam}\\fP \\fB-om\\fP \\fIOUT.{sam,bam}\\fP \\fB-of\\fP \\fIOUT.fa\\fP");
    addDescription(parser, "Extract OEA alignments and orphans from coordinate sorted BAM file.");

    // Define OeaExtractorOptions -- General
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Only print on errors."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    // Define OeaExtractorOptions -- Input / Output Related
    addSection(parser, "Input / Output Options");
    addOption(parser, seqan::ArgParseOption("im", "input-mapping", "SAM/BAM file to use as the input.", seqan::ArgParseOption::INPUTFILE, "IN"));
    setRequired(parser, "input-mapping");
    setValidValues(parser, "input-mapping", "sam bam");
    addOption(parser, seqan::ArgParseOption("om", "output-oeas", "SAM/BAM file to write to", seqan::ArgParseOption::OUTPUTFILE, "OUT"));
    setRequired(parser, "output-oeas");
    setValidValues(parser, "output-oeas", "sam bam");
    addOption(parser, seqan::ArgParseOption("of", "output-orphans", "FASTA file to write to", seqan::ArgParseOption::OUTPUTFILE, "OUT"));
    setRequired(parser, "output-orphans");
    setValidValues(parser, "output-orphans", "fq fastq");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // General Options
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    // Input / Output Options.
    getOptionValue(options.inputFile, parser, "input-mapping");
    getOptionValue(options.outputOeaFile, parser, "output-oeas");
    getOptionValue(options.outputOrphansFile, parser, "output-orphans");

    return res;
}

void OeaExtractorApp::runPipeline()
{
    std::unique_ptr<ProgressBar> progress;

    BamFilterOptions bamFilterOptions;
    bamFilterOptions.maxCoverage = 200;
    bamFilterOptions.minAlignmentQuality = 0;
    bamFilterOptions.clippingMinLength = 15;
    bamFilterOptions.disableRealignment = true;
    bamFilterOptions.realignmentNumThreads = 1;
    bamFilterOptions.realignmentChunkSize = 1024;
    bamFilterOptions.maxAlignmentLength = 2000;
    bamFilterOptions.libraryInfo = options.libraryInfo;
    bamFilterOptions.fragmentSizeFactor = 8;
    // bamFilterOptions.refFileName = options.inputFile;

    BamReaderOptions bamReaderOptions;
    bamReaderOptions.bamFileName = options.inputFile;

    typedef ThreadSafeQueue<std::vector<seqan::BamAlignmentRecord *>> TQueue;

    // Setup BamReader with progress callback.
    BamReader bamReader(bamReaderOptions);
    outStream.header = bamReader.bamStream().header;
    // Setup progress indicator.
    unsigned long long const MIB = 1024 * 1024;
    progress.reset(new ProgressBar(std::cerr, 0, fileSize(bamReader.bamStream()) / MIB,
                                   (options.verbosity == 1)));
    progress->setLabel(toCString(bamReaderOptions.bamFileName));
    progress->updateDisplay();
    auto callback = [&bamReader,&progress,MIB] { progress->advanceTo(positionInFile(bamReader.bamStream()) / MIB); };
    bamReader.setProgressCallback(callback);
    int BUFFER_SIZE = 2;
    BamFilter bamFilter(bamReader.bamIOContext(), bamFilterOptions);
    Pipeline<std::vector<seqan::BamAlignmentRecord *>> pipeline = bamFilter.makePipeline(BUFFER_SIZE);
    pipeline.start();

    // Pipeline of threads.
    std::vector<std::thread> threads;

    // Read records from BAM File and push them to the queue back.
    threads.push_back(std::thread(
            [&]() {
                std::vector<seqan::BamAlignmentRecord *> buffer;
                while (!bamReader.atEnd())
                {
                    bamReader.read(buffer);
                    pipeline.frontQueue().push(std::move(buffer));  // one copy is enough
                }
                // Close the next pipeline steps.
                pipeline.frontQueue().close();
            }));

    // Write one record per aligned pair, can appear in any order.
    auto writePairRecord = [&](seqan::BamAlignmentRecord record) {
        static seqan::BamAlignmentRecord buffer;
        static bool first = true;
        if (first)
        {
            first = false;
            buffer = record;
            return;
        }

        seqan::BamAlignmentRecord * anchor = &record;
        seqan::BamAlignmentRecord * shadow = &buffer;
        if (hasFlagUnmapped(*anchor))
            std::swap(anchor, shadow);

        SEQAN_ASSERT_EQ(anchor->qName, shadow->qName);
        SEQAN_ASSERT(hasFlagUnmapped(*anchor) || hasFlagUnmapped(*shadow));
        SEQAN_ASSERT_NEQ(hasFlagUnmapped(*anchor), hasFlagUnmapped(*shadow));

        seqan::BamTagsDict tagsDict(anchor->tags);
        setTagValue(tagsDict, "NS", toCString(shadow->seq));
        setTagValue(tagsDict, "NQ", toCString(shadow->qual));

        if (writeRecord(outStream, *anchor) != 0)
            throw OeaExtractorAppException("Problem writing SAM record.");
        first = true;
    };

    // Consume records and write them to outStream.
    threads.push_back(std::thread(
            [&]() {
                TQueue::QueueResult res;
                std::vector<seqan::BamAlignmentRecord *> buffer;
                while ((res = pipeline.backQueue().pop(buffer)) == TQueue::OK)
                {
                    // Write to outStream.
                    for (auto ptr : buffer)
                    {
                        writePairRecord(*ptr);
                        delete ptr;
                    }
                }
            }));


    // Wait for all threads and the pipeline to finish.
    for (auto & thread : threads)
        thread.join();
    pipeline.join();
    
    // Finalize progress display.
    progress->finish();
}

int OeaExtractorApp::run(int argc, char const ** argv)
{
    // Parse Options.
    if (options.verbosity >= 1)
        std::cerr << "Parsing command line...\n";
    seqan::ArgumentParser::ParseResult res = parseCommandLine(argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return (res == seqan::ArgumentParser::PARSE_ERROR);

    // Estimate library size.
    if (options.verbosity >= 1)
        std::cerr << "Estimating library size...\n";
    BamLibraryEstimator estimator(options.autoLibraryNumRecords);
    if (estimator.run(options.libraryInfo, toCString(options.inputFile)) != 0)
        throw OeaExtractorAppException("Could not estimate library info!");
    if (options.verbosity >= 1)
        std::cerr << "Insert Size\n"
                  << "\tmedian\t" << options.libraryInfo.median << "\n"
                  << "\tstd dev\t" << options.libraryInfo.stdDev << "\n"
                  << "\tmax normal\t" << options.libraryInfo.maxNormalISize << "\n"
                  << "\torientation\t" << getOrientationStr(options.libraryInfo.defaultOrient) << "\n";

    if (options.verbosity >= 1)
        std::cerr << "Extracting orphans...\n";
    OrphanExtractor extractor(toCString(options.inputFile), toCString(options.outputOrphansFile),
                              options.verbosity);
    extractor.run();

    if (options.verbosity >= 1)
        std::cerr << "Extracting OEA alignments...\n";
    open(outStream, toCString(options.outputOeaFile), seqan::BamStream::WRITE);
    if (!isGood(outStream))
        throw OeaExtractorAppException("Could not open output OEA file for writing.");
    runPipeline();

    if (options.verbosity >= 1)
        std::cerr << "All done. Have a nice day!\n";

    return 0;
}

int main(int argc, char const ** argv)
{
    OeaExtractorApp app;
    return app.run(argc, argv);
}
