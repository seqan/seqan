// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#include "finishing_step.h"

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include "shared/progress_indicator.h"
#include "shared/streaming_exception.h"

#include "anise/anise_options.h"
#include "anise/app_state.h"
#include "anise/file_name_tokens.h"
#include "anise/library_info.h"
#include "anise/site_data.h"
#include "anise/site_state.h"
#include "anise/temporary_file_manager.h"
#include "anise/time_log.h"
#include "anise/bam_tag_names.h"

namespace  // anonymous namespace
{

// ---------------------------------------------------------------------------
// Function trimAfterSpace()
// ---------------------------------------------------------------------------

// Trim after the first whitespace.
void trimAfterSpace(seqan::CharString & s)
{
    unsigned i = 0;
    for (; i < length(s); ++i)
        if (isspace(s[i]))
            break;
    resize(s, i);
}

// ----------------------------------------------------------------------------
// Function eraseIf()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPred>
void eraseIf(std::vector<TValue> & v, TPred pred)
{
    v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
}

// --------------------------------------------------------------------------
// HappyMateCounts
// --------------------------------------------------------------------------

struct HappyMateCounts
{
    // Total number of mates.
    int numTotal;
    // Number of mates with link to other.
    int numSingle;
    // Number of OEA mates.
    int numOea;
    // Number of happy mates.
    int numHappy;
    // Number of unhappy mates.
    int numUnhappy;

    HappyMateCounts() : numTotal(0), numSingle(0), numOea(0), numHappy(0), numUnhappy(0)
    {}

    void addSingle()
    {
        numTotal += 1;
        numSingle += 1;
    }

    void addOea()
    {
        numTotal += 1;
        numOea += 1;
    }

    void addHappy()
    {
        numTotal += 1;
        numHappy += 1;
    }

    void addUnhappy()
    {
        numTotal += 1;
        numUnhappy += 1;
    }

    void print(std::ostream & out) const
    {
        out << "TOTAL=" << numTotal << " SINGLE=" << numSingle << " OEA=" << numOea
            << " HAPPY=" << numHappy << " UNHAPPY=" << numUnhappy;
    }
};

// --------------------------------------------------------------------------
// HappyMateCounter
// --------------------------------------------------------------------------

// Counts happyness read-wise, not link wise.

class HappyMateCounter
{
public:

    HappyMateCounter(int numRefs,
                     BamLibraryInfo const & libraryInfo,
                     AniseOptions const & options) :
            numRefs(numRefs), counts(numRefs), options(options), libraryInfo(libraryInfo)
    {}

    void operator()(seqan::BamAlignmentRecord const & record)
    {
        if (hasFlagUnmapped(record))
            return;  // orphan or shadow
        if (hasFlagNextUnmapped(record))
            counts[record.rID].addOea();
        else if (record.rID != record.rNextId)
            counts[record.rID].addSingle();
        else if (abs(abs(record.tLen) - libraryInfo.median) <=
                 libraryInfo.stdDev * options.fragmentSizeFactor)
            counts[record.rID].addHappy();
        else
            counts[record.rID].addUnhappy();
    }

    int numRefs;
    std::vector<HappyMateCounts> counts;
    AniseOptions const & options;
    BamLibraryInfo const & libraryInfo;
};

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class FinishingStepImpl
// --------------------------------------------------------------------------

class FinishingStepImpl
{
public:

    FinishingStepImpl(TemporaryFileManager & tmpMgr, AniseOptions const & options) :
            tmpMgr(tmpMgr), options(options)
    {}

    unsigned limitRID(unsigned num) const
    {
        unsigned result = num;
        if (options.onlyWriteOutBest)
            result = std::min(result, 1u);
        return result;
    }

    void run()
    {
        // Load application state.
        appState.load(tmpMgr);

        // Open output stream.
        seqan::SequenceStream outSeqs;
        open(outSeqs, toCString(options.outputFasta), seqan::SequenceStream::WRITE);
        if (!isGood(outSeqs))
            throw StreamingException() << "Could not open output FASTA file " << toCString(options.outputFasta);

        // Open output SAM stream and build the header.
        seqan::BamStream outMapping;
        if (!empty(options.outputMapping))
        {
            open(outMapping, toCString(options.outputMapping), seqan::BamStream::WRITE);
            if (!isGood(outMapping))
                throw StreamingException() << "Could not open output mapping file " << toCString(options.outputMapping);
            outMapping.header = buildBamHeader();
        }

        // Merge the output.
        ProgressBar pb(std::cerr, 0, appState.numSites, (options.verbosity == AniseOptions::NORMAL));
        pb.setLabel("  building result");
        pb.updateDisplay();

        // Merge each site's scaffold seqs and alignments to output.
        for (int siteID = 0, contigOffset = 0; siteID < appState.numSites; ++siteID)
        {
            if (options.debugSiteID != -1 && options.debugSiteID != siteID)
                continue;  // Skip, not selected.

            // Load site state to get last active step no.
            AssemblySiteState siteState;
            siteState.load(tmpMgr, siteID);
            if (siteState.stepNo == 0)
                continue;  // Skip if no usable assembly in first step.

            // Load site's data.
            SiteData siteData;
            siteData.load(tmpMgr, siteID, siteState.stepNo);
            if (empty(siteData.scaffold.seqs))
                continue;  // Skip if no usable scaffold sequence.

            // Count happy/unhappy mates.
            HappyMateCounter hmc(length(siteData.scaffold.seqs), appState.libraryInfos[0], options);
            std::for_each(siteData.readSet.bamRecords.begin(), siteData.readSet.bamRecords.end(),
                          [&](seqan::BamAlignmentRecord const & record) { hmc(record); });

            // Soft-mask bases in scaffold from original alignment.
            seqan::StringSet<seqan::CharString> maskedSeqs = softMaskScaffold(siteData.scaffold.seqs, siteData.readSet.bamRecords);
            seqan::StringSet<seqan::CharString> refNames;
            unsigned maxRID = limitRID(length(siteData.scaffold.refNames));
            char const * label[] = { "no", "yes" };
            for (unsigned rID = 0; rID < maxRID; ++rID)
            {
                // Append site info to scaffold seq names.
                std::stringstream ss;
                ss << siteData.scaffold.refNames[rID] << " REF=" << siteState.refName << " POS="
                   << siteState.pos + 1 << " STEPS=" << siteState.stepNo << " ";
                //hmc.counts[rID].print(ss);  // TODO(holtgrew): Add back happiness!

                // std::cerr << siteData.state.comment << "\n";

                ss << " ANCHORED_LEFT=" << label[!!siteData.scaffold.scaffoldInfos[rID].anchoredLeft]
                   << " ANCHORED_RIGHT=" << label[!!siteData.scaffold.scaffoldInfos[rID].anchoredRight]
                   << " SPANNING=" << label[!!siteData.scaffold.scaffoldInfos[rID].spanning];

                // Tag for why we stopped.
                if (siteData.state.comment.find("too many") != std::string::npos)
                    ss << " STOPPED=too_many_reads";
                else if (siteData.state.comment.find("too few new") != std::string::npos)
                    ss << " STOPPED=no_more_reads";
                else
                    ss << " STOPPED=round_limit";

                appendValue(refNames, ss.str());
            }
            resize(maskedSeqs, length(refNames));
            if (writeAll(outSeqs, refNames, maskedSeqs) != 0)
                throw std::runtime_error("Problem writing to sequences out file.");

            // Remove orphan BAM pairs.
            eraseIf(siteData.readSet.bamRecords, [](seqan::BamAlignmentRecord const & arg) {
                    return (hasFlagUnmapped(arg) && hasFlagNextUnmapped(arg)); });
            // Remove BAM records on wrong assembly.
            eraseIf(siteData.readSet.bamRecords, [maxRID](seqan::BamAlignmentRecord const & arg) {
                    return (arg.rID >= (int)maxRID || arg.rNextId >= (int)maxRID); });

            // Sort BAM records by coordinate.
            auto ltCoordinate = [](seqan::BamAlignmentRecord const & lhs, seqan::BamAlignmentRecord const & rhs)
            {
                auto coord = [](seqan::BamAlignmentRecord const & arg) {
                    return std::pair<unsigned, int>(arg.rID, arg.beginPos);
                };
                return coord(lhs) < coord(rhs);
            };
            std::stable_sort(siteData.readSet.bamRecords.begin(), siteData.readSet.bamRecords.end(), ltCoordinate);

            // Write BAM records.
            if (!empty(options.outputMapping))
                writeBamRecords(outMapping, siteData.readSet.bamRecords, contigOffset);

            // Update contig offset.
            contigOffset += limitRID(length(siteData.scaffold.seqs));

            pb.advanceTo(siteID);
        }
        pb.finish();
    }

private:

    // Convert CIGAR string, e.g. "1M3I20M" into a string of CigarElement<> object.
    seqan::String<seqan::CigarElement<> > parseCigar(seqan::CharString str) const
    {
        seqan::String<seqan::CigarElement<> > result;
        seqan::CigarElement<> element;

        auto it = begin(str, seqan::Standard());
        auto itEnd = end(str, seqan::Standard());

        seqan::CharString buffer;
        while (it != itEnd && *it != '\0')
        {
            clear(buffer);
            for (; isdigit(*it); ++it)
                appendValue(buffer, *it);
            if (it == itEnd || *it == '\0')
                break;  // invalid, but ignore
            element.count = seqan::lexicalCast<__uint32>(buffer);
            element.operation = *it++;
            appendValue(result, element);
        }

        return result;
    }

    // Get number of soft-clipped characters at beginning/end of CIGAR string.
    std::pair<int, int> clippingSize(seqan::String<seqan::CigarElement<> > const & str) const
    {
        auto clippedCount = [](seqan::CigarElement<> const & el) {
            return el.operation == 'S' ? el.count : 0;
        };

        if (empty(str))
            return std::make_pair(0, 0);
        else if (length(str) == 1u)
            return std::make_pair(clippedCount(front(str)), 0);
        else
            return std::make_pair(clippedCount(front(str)), clippedCount(back(str)));
    }

    // Soft-mask bases in scaffold from original alignments.
    seqan::StringSet<seqan::CharString> softMaskScaffold(seqan::StringSet<seqan::Dna5String> const & contigs,
                                                         std::vector<seqan::BamAlignmentRecord> /*const*/ & records)
    {
        // Build (yet) unmasked sequences.
        seqan::StringSet<seqan::CharString> seqs;
        for (auto const & contig : contigs)
            appendValue(seqs, contig);

        // Soft-masked bases that were already in the input mapping.
        for (auto /*const*/ & record : records)
        {
            seqan::BamTagsDict tagsDict(record.tags);
            unsigned ignored = 0, cigarIdx;
            // Ignore if currently not aligned.
            if (hasFlagUnmapped(record))
                continue;
            // Ignore if no original alignment RID/pos or original cigar.
            if (!findTagKey(ignored, tagsDict, "oR") || !findTagKey(cigarIdx, tagsDict, "oP") ||
                !findTagKey(cigarIdx, tagsDict, "oC"))
                continue;
            seqan::CharString cigarStr;
            if ((getTagType(tagsDict, cigarIdx) != 'Z' && !extractTagValue(cigarStr, tagsDict, cigarIdx)) || empty(cigarStr))
                continue;

            auto clipped = clippingSize(parseCigar(suffix(cigarStr, 1)));
            if (clipped.first || clipped.second)
                continue;  // Skip if had clipping, conservative masking.

            int beginPos = record.beginPos, endPos = record.beginPos + getAlignmentLengthInRef(record);
            // std::cerr << "CLIPPING FOR " << record.qName << "(" << record.seq << ") : " << beginPos << " - " << endPos << "\n";
            auto itBegin = iter(seqs[record.rID], beginPos, seqan::Standard());
            auto itEnd = iter(seqs[record.rID], endPos, seqan::Standard());
            std::transform(itBegin, itEnd, itBegin, tolower);
        }

        return seqs;
    }

    // Write out BAM records after updating reference id.
    void writeBamRecords(seqan::BamStream & out,
                         std::vector<seqan::BamAlignmentRecord> & records,
                         int contigOffset)
    {
        // Update BAM records before writing out.
        std::for_each(records.begin(), records.end(),
                      [contigOffset](seqan::BamAlignmentRecord & record) {
                          record.rID += contigOffset;
                          record.rNextId += contigOffset;
                      });
        // Restore original SEQ and QUAL of record.
        std::for_each(records.begin(), records.end(),
                      [](seqan::BamAlignmentRecord & record) {
                            seqan::BamTagsDict tagsDict(record.tags);
                            unsigned idx = 0;
                            if (!findTagKey(idx, tagsDict, KEY_ORIG_SEQ))
                                return;  // no oS tag
                            extractTagValue(record.seq, tagsDict, idx);
                            resize(record.cigar, 1);
                            front(record.cigar).operation = 'M';
                            front(record.cigar).count = length(record.seq);
                            resize(record.qual, length(record.seq), 'I');
                            idx = 0;
                            if (!findTagKey(idx, tagsDict, KEY_ORIG_QUAL))
                                return;  // no oQ tag
                            extractTagValue(record.qual, tagsDict, idx);
                      });

        for (auto const & record : records)
            if (writeRecord(out, record) != 0)
                throw std::runtime_error("Problem writing to output mapping file.");
    }

    // Load the scaffolds for each site and write them to the output file.
    seqan::BamHeader buildBamHeader()
    {
        ProgressBar pb(std::cerr, 0, appState.numSites, (options.verbosity == AniseOptions::NORMAL));
        pb.setLabel("  building BAM header");
        pb.updateDisplay();

        // Initialize BAM header.
        seqan::BamHeader bamHeader;
        seqan::BamHeaderRecord hd;
        hd.type = seqan::BAM_HEADER_FIRST;
        setTagValue("VN", "1.4", hd);
        setTagValue("SO", "coordinate", hd);
        appendValue(bamHeader.records, hd);
        seqan::BamHeaderRecord pg;
        pg.type = seqan::BAM_HEADER_PROGRAM;
        setTagValue("ID", "anise", pg);
        setTagValue("PN", "anise", pg);
        setTagValue("CL", options.commandLine.c_str(), pg);
        appendValue(bamHeader.records, pg);

        seqan::StringSet<seqan::CharString> ids, seqs;
        for (int siteID = 0; siteID < appState.numSites; ++siteID)
        {
            // Load site state to get last active step no.
            AssemblySiteState siteState;
            siteState.load(tmpMgr, siteID);
            if (siteState.stepNo == 0)
                continue;  // Skip if no usable assembly in first step.

            // Open scaffold seqs file.
            std::string path = tmpMgr.fileName(SCAFFOLD_SEQS_TOKEN, SCAFFOLD_SEQS_EXT, siteID, siteState.stepNo);
            seqan::SequenceStream inSeqs;
            open(inSeqs, path.c_str());
            if (atEnd(inSeqs))
                continue;  // Ignore if no usable assembly in first place.
            if (!isGood(inSeqs))
                throw StreamingException() << "Could not open temporary scaffold seqs file " << path.c_str();
            if (readAll(ids, seqs, inSeqs) != 0)
                throw StreamingException() << "Could not read from temporary scaffold seqs file " << path.c_str();

            // Fill sequenceInfos and @SQ headers.
            for (unsigned i = 0; i < limitRID(length(ids)); ++i)
            {
                trimAfterSpace(ids[i]);

                resize(bamHeader.sequenceInfos, length(bamHeader.sequenceInfos) + 1);
                back(bamHeader.sequenceInfos).i1 = ids[i];
                back(bamHeader.sequenceInfos).i2 = length(seqs[i]);

                seqan::BamHeaderRecord sq;
                sq.type = seqan::BAM_HEADER_REFERENCE;
                setTagValue("SN", ids[i], sq);
                std::stringstream ss;
                ss << length(seqs[i]);
                setTagValue("LN", ss.str().c_str(), sq);
                appendValue(bamHeader.records, sq);
            }

            // Advance progress bar.
            pb.advanceBy(1);
            pb.updateDisplay();
        }
        pb.finish();

        return bamHeader;
    }

    // Application state object to use.
    AppState appState;
    // Manager object to use for opening temporary files.
    TemporaryFileManager & tmpMgr;
    // Application configuration.
    AniseOptions const & options;
};

// --------------------------------------------------------------------------
// Class FinishingStep
// --------------------------------------------------------------------------

FinishingStep::FinishingStep(TemporaryFileManager & manager, AniseOptions const & options) :
        impl(new FinishingStepImpl(manager, options))
{}

FinishingStep::~FinishingStep()
{}

void FinishingStep::run()
{
    impl->run();
}
