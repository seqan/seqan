// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Compute transcripts from a genome and a GFF/GTF file.  Optionally, you
// can apply a VCF file to the genome before splicing.
//
// Transcripts must not span structural variants.
// ==========================================================================

// Note: We treat all given variants as phased.

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>
#include <seqan/gff_io.h>

#include "vcf_materialization.h"
#include "mason_options.h"
#include "mason_types.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MyGffRecord
// --------------------------------------------------------------------------

// Subclass of GffRecord that has an rID member.

class MyGffRecord : public seqan::GffRecord
{
public:
    int rID;

    static const int INVALID_IDX;

    MyGffRecord() : seqan::GffRecord(), rID(std::numeric_limits<int>::max())
    {}
};

const int MyGffRecord::INVALID_IDX = std::numeric_limits<int>::max();

// --------------------------------------------------------------------------
// Class SplicingInstruction
// --------------------------------------------------------------------------

// Represents one exon.

struct SplicingInstruction
{
    // ID of the transcript that this instruction belongs to.
    int transcriptID;
    // Begin and end position of the exon.
    int beginPos, endPos;
    // The strand of the instruction '-' or '+'.
    char strand;

    SplicingInstruction() : transcriptID(-1), beginPos(-1), endPos(-1), strand('.')
    {}

    SplicingInstruction(int transcriptID, int beginPos, int endPos, char strand) :
            transcriptID(transcriptID), beginPos(beginPos), endPos(endPos), strand(strand)
    {}

    bool operator<(SplicingInstruction const & other) const
    {
        return std::make_pair(transcriptID, std::make_pair(beginPos, std::make_pair(endPos, strand))) <
                std::make_pair(other.transcriptID,
                               std::make_pair(other.beginPos, std::make_pair(other.endPos, other.strand)));
    }
};

bool differentTranscript(SplicingInstruction const & lhs, SplicingInstruction const & rhs)
{
    return lhs.transcriptID != rhs.transcriptID;
}

// --------------------------------------------------------------------------
// Class MasonSplicingApp
// --------------------------------------------------------------------------

class MasonSplicingApp
{
public:
    // The configuration to use.
    MasonSplicingOptions const & options;

    // The random number generation.
    TRng rng;

    // Materialization of VCF.
    VcfMaterializer vcfMat;

    // Input GFF/GTF stream.
    seqan::GffFileIn gffFileIn;

    // Output sequence stream.
    seqan::SeqFileOut seqFileOut;

    MasonSplicingApp(MasonSplicingOptions const & _options) :
            options(_options), rng(options.seed),
            vcfMat(rng, toCString(options.matOptions.fastaFileName), toCString(options.matOptions.vcfFileName))
    {}

    int run()
    {
        // Intialization
        std::cerr << "__INITIALIZATION_____________________________________________________________\n"
                  << "\n";

        std::cerr << "Opening files...";
        try
        {
            vcfMat.init();

            if (!open(seqFileOut, toCString(options.outputFileName)))
                throw MasonIOException("Could not open output file.");

            if (!open(gffFileIn, toCString(options.inputGffFile)))
                throw MasonIOException("Could not open GFF/GTF file.");
        }
        catch (MasonIOException & e)
        {
            std::cerr << "\nERROR: " << e.what() << "\n";
            return 1;
        }
        std::cerr << " OK\n";

        // Perform genome simulation.
        std::cerr << "\n__COMPUTING TRANSCRIPTS______________________________________________________\n"
                  << "\n";

        // Read first GFF record.
        MyGffRecord record;
        _readFirstRecord(record);
        if (record.rID == std::numeric_limits<int>::max())
            return 0;  // at end, could not read any, done

        // Transcript names.
        typedef seqan::StringSet<seqan::CharString> TNameStore;
        typedef seqan::NameStoreCache<TNameStore> TNameStoreCache;
        TNameStore transcriptNames;
        TNameStoreCache transcriptNamesCache(transcriptNames);

        // The splicing instructions for the current contig.
        std::vector<SplicingInstruction> splicingInstructions;

        // Materialized sequence.
        seqan::Dna5String seq;
        // Tanscript ids, used as a buffer below.
        seqan::String<unsigned> transcriptIDs;

        // Read GFF/GTF file contig by contig (must be sorted by reference name).  For each contig, we all recors,
        // create simulation instructions and then build the transcripts for each haplotype.
        while (record.rID != std::numeric_limits<int>::max())  // sentinel, at end
        {
            seqan::CharString refName = record.ref;
            std::cerr << "Splicing for " << refName << " ...";

            // Read GFF records for this contig.
            MyGffRecord firstGffRecord = record;
            while (record.rID == firstGffRecord.rID)
            {
                if (empty(options.gffType) || (record.type == options.gffType))
                {
                    // Make transcript names known to the record.
                    _appendTranscriptNames(transcriptIDs, transcriptNames, transcriptNamesCache, record);
                    // Add the splicing instructions for this record to the list for this contig.
                    for (unsigned i = 0; i < length(transcriptIDs); ++i)
                        splicingInstructions.push_back(SplicingInstruction(transcriptIDs[i], record.beginPos,
                                                                           record.endPos, record.strand));
                }

                if (atEnd(gffFileIn))
                {
                    record.rID = std::numeric_limits<int>::max();
                    break;
                }

                readRecord(record, gffFileIn);
                // Translate ref to idx from VCF.
                unsigned idx = 0;
                if (!getIdByName(idx, vcfMat.faiIndex, record.ref))
                    throw MasonIOException("Reference name from GFF/GTF not in VCF!");
                record.rID = idx;
            }

            // ---------------------------------------------------------------
            // Process the splicing instructions.
            // ---------------------------------------------------------------

            // First, sort them.
            std::sort(splicingInstructions.begin(), splicingInstructions.end());

            // Materialize all haplotypes of this contig
            int rID = 0, hID = 0;  // reference and haplotype id
            // Get index of the gff record's reference in the VCF file.
            unsigned idx = 0;
            if (!getIdByName(idx, vcfMat.faiIndex, refName))
            {
                std::stringstream ss;
                ss << "Reference from GFF file " << refName << " unknown in FASTA/FAI file.";
                throw MasonIOException(ss.str());
            }
            rID = idx;

            vcfMat.currRID = rID - 1;
            std::vector<SmallVarInfo> varInfos;  // small variants for counting in read alignments
            std::vector<std::pair<int, int> > breakpoints;  // unused/ignored
            while (vcfMat.materializeNext(seq, varInfos, breakpoints, rID, hID))
            {
                std::cerr << " (allele " << (hID + 1) << ")";
                if (rID != (int)idx)
                    break;  // no more haplotypes for this reference
                _performSplicing(splicingInstructions, seq, transcriptNames, hID, vcfMat);
            }

            std::cerr << " DONE.\n";

            // ---------------------------------------------------------------
            // Handle contig switching.
            // ---------------------------------------------------------------

            // Check that the input GFF file is clustered (weaker than sorted) by reference name.
            if (record.rID < firstGffRecord.rID)
                throw MasonIOException("GFF file not sorted or clustered by reference.");
            // Reset transcript names and cache.
            clear(transcriptNames);
            refresh(transcriptNamesCache);
            // Flush splicing instructions.
            splicingInstructions.clear();
        }

        std::cerr << "\nDone splicing FASTA.\n";

        return 0;
    }

    // Perform splicing of transcripts.
    void _performSplicing(std::vector<SplicingInstruction> const & instructions,
                          seqan::Dna5String const & seq,
                          seqan::StringSet<seqan::CharString> const & tNames,
                          int hID,  // -1 in case of no variants
                          VcfMaterializer const & vcfMat)
    {
        typedef std::vector<SplicingInstruction>::const_iterator TIter;
        TIter it = instructions.begin();
        TIter itEnd = std::adjacent_find(it, instructions.end(), differentTranscript);
        if (itEnd != instructions.end())
            ++itEnd;

        seqan::Dna5String transcript, buffer;

        do
        {
            clear(transcript);

            bool onBreakpoint = false;
            int tID = it->transcriptID;
            for (; it != itEnd; ++it)
            {
                // Convert from original coordinate system to coordinate system with SVs.
                std::pair<int, int> smallVarInt = vcfMat.posMap.originalToSmallVarInterval(
                        it->beginPos, it->endPos);
                GenomicInterval gi = vcfMat.posMap.getGenomicIntervalSmallVarPos(smallVarInt.first);
                SEQAN_ASSERT_GT(smallVarInt.second, 0);
                GenomicInterval giR = vcfMat.posMap.getGenomicIntervalSmallVarPos(smallVarInt.second - 1);
                bool overlapsWithBreakpoint = (gi != giR);
                std::pair<int, int> largeVarInt = vcfMat.posMap.smallVarToLargeVarInterval(
                        smallVarInt.first, smallVarInt.second);

                // Transcripts with exons overlapping breakpoints are not written out.
                if (overlapsWithBreakpoint)
                {
                    onBreakpoint = true;
                    break;
                }

                // Append buffer to transcript in original state or reverse-complemented.
                buffer = infix(seq, largeVarInt.first, largeVarInt.second);
                if (it->strand != gi.strand)
                    reverseComplement(buffer);
                append(transcript, buffer);
            }

            if (onBreakpoint)
            {
                std::cerr << "\nWARNING: Exon lies on breakpoint!\n";
                while (it != instructions.end() && it->transcriptID == tID)
                    ++it;
            }
            else
            {
                std::stringstream ss;
                ss << tNames[tID];
                if (!empty(options.matOptions.vcfFileName))
                    ss << options.haplotypeNameSep << (hID + 1);
                writeRecord(seqFileOut, ss.str(), transcript);
            }

            // Search next range.
            itEnd = std::adjacent_find(it, instructions.end(), differentTranscript);
            if (itEnd != instructions.end())
                ++itEnd;
        }
        while (it != instructions.end());
    }

    // Append the transcript names for the given record.
    void _appendTranscriptNames(seqan::String<unsigned> & tIDs,  // transcript ids to write out
                                seqan::StringSet<seqan::CharString> & contigNames,
                                seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > & cache,
                                MyGffRecord const & record)
    {
        clear(tIDs);

        seqan::CharString groupNames;
        for (unsigned i = 0; i < length(record.tagNames); ++i)
            if (record.tagNames[i] == options.gffGroupBy)
                groupNames = record.tagValues[i];
        if (empty(groupNames))
            return;  // Record has no group names.

        // Write out the ids of the transcripts that the record belongs to as indices in contigNames.
        unsigned idx = 0;
        seqan::StringSet<seqan::CharString> ss;
        strSplit(ss, groupNames, seqan::EqualsChar<','>());
        for (unsigned i = 0; i < length(ss); ++i)
        {
            if (empty(ss[i]))
                continue;
            if (!getIdByName(idx, cache, ss[i]))
            {
                appendValue(tIDs, length(contigNames));
                appendName(cache, ss[i]);
            }
            else
            {
                appendValue(tIDs, idx);
            }
        }
    }

    void _readFirstRecord(MyGffRecord & record)
    {
        record.rID = record.INVALID_IDX;  // uninitialized

        bool found = false;
        while (!found && !atEnd(gffFileIn))
        {
            readRecord(record, gffFileIn);

            // Translate ref to idx from VCF.
            unsigned idx = 0;
            if (!getIdByName(idx, vcfMat.faiIndex, record.ref))
                throw MasonIOException("Reference name from GFF/GTF not in VCF!");
            record.rID = idx;

            if (empty(options.gffType) || (options.gffType == record.type))
            {
                found = true;
                break;
            }
        }
        if (!found)
            record.rID = std::numeric_limits<int>::max();
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonSplicingOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_splicing");
    // Set short description, version, and date.
    setShortDescription(parser, "Generating Transcripts");
    setDateAndVersion(parser);
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[OPTIONS] \\fB-ir\\fP \\fIIN.fa\\fP \\fB-ig\\fP \\fIIN.gff\\fP [\\fB-iv\\fP \\fIIN.vcf\\fP] \\fB-o\\fP \\fIOUT.fa\\fP");
    addDescription(parser,
                   "Create transcripts from \\fIIN.fa\\fP using the annotations from \\fIIN.gff\\fP.  The resulting "
                   "transcripts are written to \\fIOUT.fa\\fP.");
    addDescription(parser,
                   "You can pass an optional VCF file \\fIIN.vcf\\fP and the transcripts will be created from the "
                   "haplotypes stored in the VCF file.");

    // Add option and text sections.
    options.addOptions(parser);
    options.addTextSections(parser);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.getOptionValues(parser);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    MasonSplicingOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "MASON SPLICING\n"
              << "==============\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cerr << "__OPTIONS____________________________________________________________________\n"
                  << "\n";
        options.print(std::cerr);
    }

    MasonSplicingApp app(options);
    return app.run();
}
