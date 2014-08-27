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

#include "preparation_step.h"

#include <atomic>
#include <mutex>
#include <iostream>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>

#include "shared/progress_indicator.h"

#include "anise/anise_options.h"
#include "anise/app_state.h"
#include "anise/bam_tag_names.h"
#include "anise/assembler.h"
#include "anise/file_name_tokens.h"
#include "anise/file_name_tokens.h"
#include "anise/parallel_utils.h"
#include "anise/realignment.h"
#include "anise/site_data.h"
#include "anise/site_state.h"
#include "anise/store_utils.h"
#include "anise/temporary_file_manager.h"
#include "anise/time_log.h"

namespace  // anonymous namespace
{

// ----------------------------------------------------------------------------
// Function startsWithClipping().
// ----------------------------------------------------------------------------

// Returns whether the CIGAR string starts with a clipping character.

template <typename TOperation, typename TCount, typename TSpec>
inline bool startsWithClipping(seqan::String<seqan::CigarElement<TOperation, TCount>, TSpec> const & cigarString,
                               unsigned minCount = 0)
{
    if (empty(cigarString))
        return false;
    return ((front(cigarString).operation == 'S') || (front(cigarString).operation == 'H')) &&
            (front(cigarString).count >= minCount);
}

// ----------------------------------------------------------------------------
// Function startsWithClipping().
// ----------------------------------------------------------------------------

// Returns whether the CIGAR string ends with a clipping character.

template <typename TOperation, typename TCount, typename TSpec>
inline bool endsWithClipping(seqan::String<seqan::CigarElement<TOperation, TCount>, TSpec> const & cigarString,
                             unsigned minCount = 0)
{
    if (empty(cigarString))
        return false;
    return ((back(cigarString).operation == 'S') || (back(cigarString).operation == 'H')) &&
            (back(cigarString).count >= minCount);
}

// --------------------------------------------------------------------------
// Function hasClipping()
// --------------------------------------------------------------------------

inline bool hasClipping(seqan::BamAlignmentRecord const & record, unsigned minCount = 0)
{
    return startsWithClipping(record.cigar, minCount) || endsWithClipping(record.cigar, minCount);
}

// ----------------------------------------------------------------------------
// Function clippingPosition()
// ----------------------------------------------------------------------------

// Returns the position of the clipping in the reference.

inline int clippingPosition(seqan::BamAlignmentRecord const & record, unsigned minCount = 0)
{
    if (startsWithClipping(record.cigar, minCount))
        return record.beginPos;
    else if (endsWithClipping(record.cigar, minCount))
        return record.beginPos + getAlignmentLengthInRef(record);
    else
        SEQAN_FAIL("Record's CIGAR string has no clipping!");
    return -1;
}

// --------------------------------------------------------------------------
// Class OrphanExtractor
// --------------------------------------------------------------------------

// Extract the orphans from the input mapping file.  The orphans are written to a FASTQ file.  Additionally, we create a
// binary file with bool flags that determine whether an orphan read has been used already and consequently is not going
// to be mapped again in a further step.

class OrphanExtractor
{
public:

    OrphanExtractor(AniseOptions const & options) : options(options)
    {}

    // Returns the number of orphans.
    size_t run(TemporaryFileManager & tmpMgr)
    {
        size_t numOrphans = createFastq(tmpMgr);
        return numOrphans;
    }

private:

    // Create the orphans.fq file.  Returns seqan::maxValue<size_t>() if exists already.
    size_t createFastq(TemporaryFileManager & tmpMgr);

    // Returns true if the orphans exist.
    bool orphansFastqExists(std::string const & path);

    // The BamStream to use for reading the reads and the BamIndex to use for jumping.
    seqan::BamStream bamStream;
    seqan::BamIndex<seqan::Bai> bamIndex;
    // The SequenceStream to use for writing orphans.
    seqan::SequenceStream orphansStream;
    // The options to use.
    AniseOptions const & options;
};

bool OrphanExtractor::orphansFastqExists(std::string const & path)
{
    open(orphansStream, path.c_str(), seqan::SequenceStream::READ);
    return isGood(orphansStream);
}

size_t OrphanExtractor::createFastq(TemporaryFileManager & tmpMgr)
{
    size_t numOrphans = 0;

    std::string path = tmpMgr.fileName(ORPHANS_TOKEN, ORPHANS_EXT);
    if (orphansFastqExists(path))
        return seqan::maxValue<size_t>();  // OK, exists

    // Open orphans file for writing.
    open(orphansStream, path.c_str(), seqan::SequenceStream::WRITE);
    if (!isGood(orphansStream))
        throw std::runtime_error("Could not open orphans FASTQ file for writing!");

    // Open BAM and BAI file.
    open(bamStream, toCString(options.inputMapping));
    if (!isGood(bamStream))
        throw std::runtime_error("Could not open input mapping for orphans extraction!");
    seqan::CharString baiFilename = options.inputMapping;
    append(baiFilename, ".bai");
    seqan::BamIndex<seqan::Bai> bamIndex;
    if (read(bamIndex, toCString(baiFilename)) != 0)
        throw std::runtime_error("Could not open BAI file for orphans extraction!");

    // Jump to orphans.
    if (!jumpToOrphans(bamStream, bamIndex))
        throw std::runtime_error("Could not jump to orphans in BAM file.");

    // Computation for progress display and progress bar for this.
    unsigned const MIB = 1024 * 1024;
    __int64 pos = positionInFile(bamStream) / MIB;
    __int64 size = fileSize(bamStream) / MIB;
    unsigned const MAX_BATCH_SIZE = 1000;
    unsigned batchSize = 0;
    ProgressBar pb(std::cerr, pos, size, (options.verbosity == AniseOptions::NORMAL));
    pb.setLabel("  orphan extraction [MiB of BAM]");
    pb.updateDisplay();

    // Actually read the orphans.
    seqan::BamAlignmentRecord recordL, recordR;
    while (!atEnd(bamStream))
    {
        // Read pair of orphans.
        if (readRecord(recordL, bamStream) != 0)
            throw std::runtime_error("Could not read first mate record from BAM file.");
        if (!hasFlagUnmapped(recordL) || !hasFlagNextUnmapped(recordL))
            throw std::runtime_error("First mate record is not an orphans!");
        if (readRecord(recordR, bamStream) != 0)
            throw std::runtime_error("Could not read second mate record from BAM file.");
        if (!hasFlagUnmapped(recordR) || !hasFlagNextUnmapped(recordR))
            throw std::runtime_error("Second mate record is not an orphans!");

        // Update progress bar.
        if (++batchSize > MAX_BATCH_SIZE)
        {
            batchSize = 0;
            pos = positionInFile(bamStream);
            pb.advanceTo(pos / MIB);
        }

        // Perform sanity checks.
        if (recordL.qName != recordR.qName)
        {
            std::stringstream ss;
            ss << "\nERROR: Invalid pairing for reads " << recordL.qName << ", " << recordR.qName << ". "
               << "Are the records sorted by QNAME?";
            throw std::runtime_error(ss.str());
        }

        // Enforce ordering.
        if (hasFlagFirst(recordR))
        {
            using std::swap;
            swap(recordL, recordR);
        }

        // Write out left and right mate of orphan pair.
        if (writeRecord(orphansStream, recordL.qName, recordL.seq, recordL.qual) != 0)
            throw std::runtime_error("Could not write first mate to orphans FASTQ file.");
        if (writeRecord(orphansStream, recordR.qName, recordR.seq, recordR.qual) != 0)
            throw std::runtime_error("Could not write second mate to orphans FASTQ file.");

        numOrphans += 2;
    }

    // Finish progress bar display.
    pb.finish();

    // Manually close orphans stream to protect against crashes later that might have kep the stream unflushed.
    close(orphansStream);
    
    return numOrphans;
}

// --------------------------------------------------------------------------
// Class AlignmentExtractor
// --------------------------------------------------------------------------

class AlignmentExtractor
{
public:

    AlignmentExtractor(AniseOptions const & options);

    // Read the alignments in the given range in siteState with respect to scaffold and store it in readSet.
    void run(ReadSet & readSet, AssemblySiteState const & siteState, std::pair<int, int> range);

private:

    // Jump to region around the site (lengt of the maximal assumed fragment left of the position of the site in siteState).
    void jumpToRegion(AssemblySiteState const & siteState);
    // Perform the first scanning pass.  This will collect all true OEA reads and collect the name of the reads that
    // have a mate with clipping.
    void firstPass(ReadSet & readSet, std::set<seqan::CharString> & clippedNames, AssemblySiteState const & siteState,
                   std::pair<int, int> range);
    // Perform the second scanning pass.  This will collect the records that have one clipped mate and memoized in the
    // first pass.
    void secondPass(ReadSet & readSet, std::set<seqan::CharString> const & clippedNames, AssemblySiteState const & siteState,
                    std::pair<int, int> range);

    // Write tags with original rID, position, and cigar string.
    void updateRecordTags(seqan::BamAlignmentRecord & record) const;

    // The BamStream to use for reading the records.
    seqan::BamStream bamStream;
    // The class to use for reading the BAI files.
    seqan::BamIndex<seqan::Bai> bamIndex;
    // The configuration, used for the BAM path.
    AniseOptions const & options;
};

AlignmentExtractor::AlignmentExtractor(AniseOptions const & options) : options(options)
{
    // Open BamStream.
    open(bamStream, toCString(options.inputMapping));
    if (!isGood(bamStream))
        throw std::runtime_error("Could not open BAM mapping file for reading.");

    // Read BAI file.
    std::string baiPath = toCString(options.inputMapping);
    baiPath.append(".bai");
    if (read(bamIndex, baiPath.c_str()) != 0)
        throw std::runtime_error("Could not read BAI file.");
}

void AlignmentExtractor::run(ReadSet & readSet, AssemblySiteState const & siteState, std::pair<int, int> range)
{
    // Jump to the region beginning for the first pass.
    jumpToRegion(siteState);
    // Perform first pass.
    std::set<seqan::CharString> clippedNames;
    firstPass(readSet, clippedNames, siteState, range);
    // Jump to the region beginning for the second pass.
    jumpToRegion(siteState);
    // Perform the second pass.
    secondPass(readSet, clippedNames, siteState, range);
}

void AlignmentExtractor::jumpToRegion(AssemblySiteState const & siteState)
{
    int pos = std::max(0, siteState.pos - options.maxFragmentSize());
    int endPos = bamStream.header.sequenceInfos[siteState.rID].i2;;
    if (siteState.pos + 2 * options.maxFragmentSize() < endPos)
        endPos = siteState.pos + 2 * options.maxFragmentSize();
    bool hasAlignments = false;
    if (!seqan::jumpToRegion(bamStream, hasAlignments, siteState.rID, pos, endPos, bamIndex))
        throw std::runtime_error("Could not jump in alignment file.");
    if (!hasAlignments)
    {
        std::stringstream ss;
        ss << "The BAM file does not have any alignments for insert site (siteID=" << siteState.siteID << ")";
        throw std::runtime_error(ss.str().c_str());
    }
}

void AlignmentExtractor::firstPass(ReadSet & readSet, std::set<seqan::CharString> & clippedNames,
                                   AssemblySiteState const & siteState, std::pair<int, int> range)
{
    (void)range;  // TODO(holtgrew): Will be needed again when we update beginpos
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStream))
    {
        if (readRecord(record, bamStream) != 0)
            throw std::runtime_error("Problem reading from BAM file.");
        updateRecordTags(record);
        if (record.rID == seqan::BamAlignmentRecord::INVALID_REFID || record.rID != siteState.rID ||
            record.beginPos > siteState.pos + options.maxFragmentSize())
            break;  // Done
        // If we reach here, the record cannot be an orphan.
        if (record.beginPos < siteState.pos - options.maxFragmentSize())
            continue;  // Skip if too far away from pos.

        // If both reads are mapped, then record the name if it is clipped.  Regardless of the clipping, do not record
        // the alignment info here if it has been mapped.
        if (hasFlagUnmapped(record) == hasFlagNextUnmapped(record))  // both are mapped
        {
            // TODO(holtgrew): Only consider clipping at the border of the fragment.
            if (hasClipping(record))
            {
                int MIN_CLIPPING = 15;
                // Consider all with sufficient clipping or all that are sufficiently closely clipped.
                if (options.assemblySiteFringeRadius == -1 && hasClipping(record, MIN_CLIPPING))
                    clippedNames.insert(record.qName);
                else if (abs(clippingPosition(record) - siteState.pos) <= options.assemblySiteFringeRadius)
                    clippedNames.insert(record.qName);
            }

            continue;  // Skip if both are mapped.
        }

        // If we reach here then the record is for an OEA pair.

        // Ignore OEAs that point into the wrong direction.
        if ((!hasFlagUnmapped(record) && record.beginPos > siteState.pos && !hasFlagRC(record)) ||
            (hasFlagUnmapped(record) && record.beginPos > siteState.pos && !hasFlagNextRC(record)))
            continue;  // Is right of site and points right or corresponding shadow.
        if ((!hasFlagUnmapped(record) && record.beginPos < siteState.pos && hasFlagRC(record)) ||
            (hasFlagUnmapped(record) && record.beginPos < siteState.pos && hasFlagNextRC(record)))
            continue;  // Is left of site and points left or corresponding shadow.

        // TODO(holtgrew): We assume good separation here and above but this may not be the case.
        // Update rID and position of the BamAlignmentRecord.
        enum { ALIGNED, SHADOW } type = hasFlagUnmapped(record) ? SHADOW : ALIGNED;
        int rID = (record.beginPos >= siteState.pos);
        // int pos = record.beginPos - (rID ? range.second : range.first);
        record.rID = rID;
        // record.beginPos = pos;
        // Update tLen to reflect OEA mapping.
        record.tLen = 0;

        record.rNextId = (record.pNext >= siteState.pos);
        // record.pNext -= (record.rNextId ? range.second : range.first);

        // Store the read sequence, this is also used for reverse-complementing below.
        appendValue(readSet.seqs, record.seq);

        // In case of a shadow records, fix the sequence orientation such that it is correct for concordant pairs.  Fix
        // flag of anchor record as well.
        if (type == SHADOW && hasFlagRC(record) == hasFlagNextRC(record))
        {
            record.flag ^= seqan::BAM_FLAG_RC;
            reverseComplement(back(readSet.seqs));
            record.seq = back(readSet.seqs);
            reverse(record.qual);
        }
        else if (type == ALIGNED && hasFlagRC(record) == hasFlagNextRC(record))
        {
            record.flag ^= seqan::BAM_FLAG_NEXT_RC;
        }

        // Store the BAM alignment record itself.
        readSet.bamRecords.emplace_back(record);
    }

    // Cleanup readSet.bamRecords, in case of flag problems (observed with BWA 0.7.5a-r405).
    unsigned offset = 0;
    for (unsigned i = 0; i + 2 < readSet.bamRecords.size();)
    {
        if (readSet.bamRecords[i].qName == readSet.bamRecords[i + 1].qName)
        {
            readSet.seqs[offset] = std::move(readSet.seqs[i]);
            readSet.seqs[offset + 1] = std::move(readSet.seqs[i + 1]);
            readSet.bamRecords[offset] = std::move(readSet.bamRecords[i]);
            readSet.bamRecords[offset + 1] = std::move(readSet.bamRecords[i + 1]);
            i += 2;
            offset += 2;
        }
        else
        {
            i += 1;
        }
    }
    resize(readSet.seqs, offset);
    readSet.bamRecords.resize(offset);
}

void AlignmentExtractor::secondPass(ReadSet & readSet, std::set<seqan::CharString> const & clippedNames,
                                    AssemblySiteState const & siteState, std::pair<int, int> range)
{
    (void)range;  // TODO(holtgrew): Required when we correct begin pos again.
    // We will collect the pairs with clipping in the std::map clippingPairs.  The first entry will be the record that
    // is not clipped (will be stored as anchored) and the second one is the record that is clipped (will be stored as
    // orphan).  We ignore incomplete pairs (could happen if the assembly window radius is too small.
    std::map<seqan::CharString, seqan::Pair<seqan::BamAlignmentRecord> > clippingPairs;

    // TODO(holtgrew): Require minimal clipping size?
    // TODO(holtgrew): Allow only one clipped mate?
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStream))
    {
        if (readRecord(record, bamStream) != 0)
            throw std::runtime_error("Problem reading from BAM file.");
        updateRecordTags(record);
        if (record.rID == seqan::BamAlignmentRecord::INVALID_REFID || record.rID != siteState.rID ||
            record.beginPos > siteState.pos + options.maxFragmentSize())
            break;  // Done
        // If we reach here, the record cannot be an orphan.
        if (hasFlagUnmapped(record) != hasFlagNextUnmapped(record))
            continue;  // skip OEA pairs
        if (record.beginPos < siteState.pos - options.maxFragmentSize())
            continue;  // Skip if too far away from pos.
        if (clippedNames.count(record.qName) == 0u)
            continue;  // Has no clipped mate.

        // Store record.
        seqan::BamAlignmentRecord * info = hasClipping(record) ?
                &clippingPairs[record.qName].i2 : &clippingPairs[record.qName].i1;
        *info = record;

        seqan::BamAlignmentRecord & infoAligned = clippingPairs[record.qName].i1;
        seqan::BamAlignmentRecord & infoShadow = clippingPairs[record.qName].i2;
        if (infoAligned.rID == seqan::BamAlignmentRecord::INVALID_REFID ||
            infoShadow.rID == seqan::BamAlignmentRecord::INVALID_REFID)
            continue;  // We have not seen both yet.

        // If we have now seen alignments for both mates of the clipping pair then we can compute the coordinates.
        // Instead of checking whether we are left or right of the called insert site we assume that the clipped record
        // is at the insert site or spans over it.  We use this to orient ourselves.
        if (infoAligned.beginPos < infoShadow.beginPos)
        {
            infoAligned.rID = 0;
            // TODO(holtgrew): Fix the begin position stuff here.
            // infoAligned.beginPos -= range.first;
        }
        else
        {
            infoAligned.rID = 1;
            // infoAligned.beginPos -= range.second;
        }
        infoShadow.rID = infoShadow.rNextId = infoAligned.rNextId = infoAligned.rID;
        infoShadow.beginPos = infoAligned.beginPos;
        infoShadow.mapQ = 0;

        infoShadow.beginPos = infoShadow.pNext = infoAligned.pNext = infoAligned.beginPos;

        // Shadow cannot have CIGAR string.
        clear(infoShadow.cigar);

        // Update tLen to reflect OEA mapping.
        infoAligned.tLen = 0;
        infoShadow.tLen = 0;

        // Update flags.
        bool alignedFirst = hasFlagFirst(infoAligned);
        bool alignedReverse = hasFlagRC(infoAligned);
        bool shadowFirst = hasFlagFirst(infoShadow);
        bool shadowReverse = hasFlagRC(infoShadow);

        infoAligned.flag = seqan::BAM_FLAG_MULTIPLE | seqan::BAM_FLAG_NEXT_UNMAPPED;
        infoShadow.flag = seqan::BAM_FLAG_MULTIPLE | seqan::BAM_FLAG_UNMAPPED;
        if (alignedFirst)
        {
            infoAligned.flag |= seqan::BAM_FLAG_FIRST;
            infoShadow.flag |= seqan::BAM_FLAG_LAST;
        }
        if (alignedReverse)
        {
            infoAligned.flag |= seqan::BAM_FLAG_RC;
            infoShadow.flag |= seqan::BAM_FLAG_NEXT_RC;
        }
        if (shadowFirst)
        {
            infoAligned.flag |= seqan::BAM_FLAG_LAST;
            infoShadow.flag |= seqan::BAM_FLAG_FIRST;
        }
        if (shadowReverse)
        {
            infoAligned.flag |= seqan::BAM_FLAG_NEXT_RC;
            infoShadow.flag |= seqan::BAM_FLAG_RC;
        }

        // Append to alignmentInfos if both are complete.
        if (infoAligned.rID != seqan::BamAlignmentRecord::INVALID_REFID &&
            infoShadow.rID != seqan::BamAlignmentRecord::INVALID_REFID)
        {
            appendValue(readSet.seqs, infoAligned.seq);
            appendValue(readSet.seqs, infoShadow.seq);
            readSet.bamRecords.emplace_back(infoAligned);
            readSet.bamRecords.emplace_back(infoShadow);
        }
    }
}

void AlignmentExtractor::updateRecordTags(seqan::BamAlignmentRecord & record) const
{
    if (hasFlagUnmapped(record))
        return;

    std::stringstream ss;
    std::for_each(begin(record.cigar, seqan::Standard()), end(record.cigar, seqan::Standard()),
                  [&](seqan::CigarElement<> el) { ss << el.count << el.operation; });

    static char const * const label[] = {"fwd", "rev"};

    seqan::BamTagsDict tagsDict(record.tags);
    setTagValue(tagsDict, KEY_BIRTH_STEP, 0);  // mapped in step 0
    setTagValue(tagsDict, KEY_ORIG_REF, toCString(bamStream.header.sequenceInfos[record.rID].i1));  // original ref
    setTagValue(tagsDict, KEY_ORIG_STRAND, label[!!hasFlagRC(record)]);
    setTagValue(tagsDict, KEY_ORIG_POS, record.beginPos + 1);  // original pos
    setTagValue(tagsDict, KEY_ORIG_CIGAR, ss.str().c_str());  // original CIGAR string
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class PreparationStepImpl
// --------------------------------------------------------------------------

class PreparationStepImpl
{
public:

    PreparationStepImpl(TemporaryFileManager & manager, AniseOptions const & options, AppState const & appState) :
            manager(manager), options(options), appState(appState), orphanExtractor(options)
    {}

    void run();

private:

    // Create the initial site data information.
    SiteData createInitialSiteData(AlignmentExtractor & alignmentExtractor,
                                   seqan::VcfRecord const & record,
                                   int siteID) const;

    // Create the initial assembly site state.
    AssemblySiteState createInitialSiteState(seqan::VcfRecord const & record, int siteID) const;
    // Prepare extraction around the sites and return the number of sites.
    unsigned prepareExtraction();
    // Extraction loop for ref seqs and read sets.
    void extractSites(unsigned numSites);
    // Extract initial reference sequences.
    //
    // This extracts a window each left and right of the position of siteState.  In the middle, some bases are left out
    // in the inner fringe.
    //
    // Returns the computed begin and end position.
    std::pair<int, int> extractInitialRefSeqs(Scaffold & scaffold,
                                              AssemblySiteState const & siteState) const;
    // Extract initial read set.  Reads overlapping with the inner fringe are marked as shadows.
    //
    // siteState can be marked inactive if no reads could be extracted.
    void extractInitialReadSet(ReadSet & readSet, AssemblySiteState & siteState,
                               AlignmentExtractor & alignmentExtractor,
                               std::pair<int, int> range) const;

    // The manager for temporary files and options.
    TemporaryFileManager & manager;
    AniseOptions const & options;
    // Application state.
    AppState const & appState;
    // We use an OrphanExtractor object to extract the orphan reads.
    OrphanExtractor orphanExtractor;
    // Used for reading the reference sequence.
    seqan::FaiIndex faiIndex;
    // The VcfStream to use for reading VCF records.
    seqan::VcfStream vcfStream;
};

void PreparationStepImpl::run()
{
    // Load application state.
    AppState appState;
    appState.load(manager);

    // Extract the orphans.
    if (options.debugSiteID == -1)
    {
        appState.numOrphans = orphanExtractor.run(manager);
        appState.save(manager);
    }

    // Run extraction for all sites.
    unsigned numSites = prepareExtraction();
    extractSites(numSites);

    // Advance state in app state.
    appState.load(manager);
    appState.numSites = numSites;
    appState.superStep = AppState::SuperStep::ASSEMBLING;
    appState.stepNo = 1;
    appState.save(manager);
}

unsigned PreparationStepImpl::prepareExtraction()
{
    // Count VCF records.
    open(vcfStream, toCString(options.inputVcf));
    if (!isGood(vcfStream))
        throw std::runtime_error("Could not open input VCF file for site extraction.");
    unsigned numSites = 0;
    seqan::VcfRecord record;
    for (; !atEnd(vcfStream); ++numSites)
        if (readRecord(record, vcfStream) != 0)
            throw std::runtime_error("Problem reading from VCF file.");

    // Open VCF file.
    open(vcfStream, toCString(options.inputVcf));
    if (!isGood(vcfStream))
        throw std::runtime_error("Could not open input VCF file for site extraction.");

    // Build index for FASTA file and write out to disk.
    if (read(faiIndex, toCString(options.inputReference)) != 0)
    {
        if (build(faiIndex, toCString(options.inputReference)) != 0)
            throw std::runtime_error("Could read FASTA file.");
        if (write(faiIndex) != 0)
            throw std::runtime_error("Could not write out FAI file.");
    }

    return numSites;
}

void PreparationStepImpl::extractSites(unsigned numSites)
{
    // Progress bar setup for site data extraction.
    ProgressBar pb(std::cerr, 0, numSites, (options.verbosity == AniseOptions::NORMAL));
    pb.setLabel("  extracted sites");
    pb.updateDisplay();

    // Mutex to use for access to vcfStream and index of current site.
    std::mutex mutex;
    std::atomic<int> nextSiteID(0);

    // This loop is executed by each thread.
    auto extractorLoop = [&]() {
        // The alignment extractor to use in this thread.
        AlignmentExtractor alignmentExtractor(options);
        // This VcfRecord is used by this thread.
        seqan::VcfRecord record;

        while (!atEnd(vcfStream))
        {
            int siteID = -1;
            // Read the record and get site ID.
            {
                std::lock_guard<std::mutex> lock(mutex);
                if (atEnd(vcfStream))
                    return;
                if (readRecord(record, vcfStream) != 0)
                    throw std::runtime_error("Problem reading from VCF file.");
                siteID = nextSiteID.fetch_add(1);
            }

            if (options.debugSiteID != -1 && options.debugSiteID != siteID)
            {
                pb.advanceBy(1);
                continue;  // Skip, not selected.
            }

            // Perform the initial assembly step.  The input is the initial site data that is directly generated from
            // the input BAM and FASTA (+.fai) files.  The resulting assembly is stored as step 0's result.
            SiteData extractedSiteData = createInitialSiteData(alignmentExtractor, record, siteID);
            if (!extractedSiteData.state.active)
            {
                extractedSiteData.save(manager);
                pb.advanceBy(1);
                continue;
            }
            SiteData siteData;
            withTimeLog("ASSEMBLY", appState.stepNo, siteID, [&]() {
                    performAssemblyStep(siteData, extractedSiteData, options, appState, manager);
                });
            siteData.save(manager);

            // Advance in progess bar.
            pb.advanceBy(1);
        }
    };

    // Run extraction loop with multiple threads.
    forkJoin(options.numThreads, extractorLoop);

    // Finish progress bar display.
    pb.finish();
}

SiteData PreparationStepImpl::createInitialSiteData(
        AlignmentExtractor & alignmentExtractor,
        seqan::VcfRecord const & record,
        int siteID) const
{
    SiteData result;

    std::stringstream ss;
    ss << "site_" << siteID;

    // Create the site states from the VCF file.
    result.state = createInitialSiteState(record, siteID);
    // Extract the reference sequence around the tentative insertion site.
    result.scaffold.siteName = ss.str();
    auto range = extractInitialRefSeqs(result.scaffold, result.state);
    // Extract the read sets around the tentative insertion sites.
    extractInitialReadSet(result.readSet, result.state, alignmentExtractor, range);

    return result;
}

AssemblySiteState PreparationStepImpl::createInitialSiteState(seqan::VcfRecord const & record, int siteID) const
{
    AssemblySiteState result(siteID);
    result.rID = record.rID;
    result.refName = toCString(vcfStream.header.sequenceNames[record.rID]);
    result.pos = record.beginPos;
    result.stepNo = 0;
    result.active = true;

    return result;
}

std::pair<int, int> PreparationStepImpl::extractInitialRefSeqs(Scaffold & scaffold,
                                                               AssemblySiteState const & siteState) const
{
    // Radius around insert site to cut out of contigs.  Any record overlapping with the fringe will be put into the MSA
    // assembly as shadow sequence.
    int const INNER_FRINGE_RADIUS = 20;

    // Cut out contigs left and right of insert site and get the names.
    //
    // First, get begin and end position of contigs to cut out.
    int contigLength = sequenceLength(faiIndex, siteState.rID);
    int leftBeginPos = std::max(siteState.pos - options.assemblySiteWindowRadius, 0);
    int leftEndPos = std::max(leftBeginPos, siteState.pos - INNER_FRINGE_RADIUS);
    std::stringstream ssLeftContigName;
    ssLeftContigName << sequenceName(faiIndex, siteState.rID) << ":" << leftBeginPos+1 << "-" << leftEndPos;
    appendValue(scaffold.refNames, ssLeftContigName.str());
    int rightBeginPos = std::min(siteState.pos + INNER_FRINGE_RADIUS, contigLength);
    int rightEndPos = std::max(rightBeginPos, siteState.pos + options.assemblySiteWindowRadius);
    std::stringstream ssRightContigName;
    ssRightContigName << sequenceName(faiIndex, siteState.rID) << ":" << rightBeginPos+1 << "-" << rightEndPos;
    appendValue(scaffold.refNames, ssRightContigName.str());

    // Get left and right contig sequence.
    resize(scaffold.seqs, 2);
    int res = 0;
    if ((res = readRegion(scaffold.seqs[0], faiIndex, siteState.rID, leftBeginPos, leftEndPos)) != 0)
        throw std::runtime_error("Problem reading left scaffold seq.");
    if ((res = readRegion(scaffold.seqs[1], faiIndex, siteState.rID, rightBeginPos, rightEndPos)) != 0)
        throw std::runtime_error("Problem reading right scaffold seq.");

    return std::make_pair(leftBeginPos, rightEndPos);
}

void PreparationStepImpl::extractInitialReadSet(ReadSet & readSet,
                                                AssemblySiteState & siteState,
                                                AlignmentExtractor & alignmentExtractor,
                                                std::pair<int, int> range) const
{
    // Extract the read alignments using an AlignmentExtractor.
    alignmentExtractor.run(readSet, siteState, range);
    readSet.computeStats();

    // If there are no reads in the inital read set then disable the site from the very beginning.  Also disable in case
    // of too many reads in initial read set.
    if (readSet.bamRecords.empty() || readSet.bamRecords.size() > (unsigned)options.stopInitialReadCount)
    {
        TimeLog::instance().log("STOPPING SITE", 0, siteState.siteID);
        siteState.active = false;
        std::stringstream ss;
        ss << "deactivating: too many alignments (" << readSet.bamRecords.size() << " > "
           << (unsigned)options.stopInitialReadCount << ")";
        siteState.comment = ss.str();
    }
}

// --------------------------------------------------------------------------
// Class PreparationStep
// --------------------------------------------------------------------------

PreparationStep::PreparationStep(TemporaryFileManager & manager, AniseOptions const & options, AppState const & appState) :
        impl(new PreparationStepImpl(manager, options, appState))
{}

PreparationStep::~PreparationStep()
{}

void PreparationStep::run()
{
    impl->run();
}

