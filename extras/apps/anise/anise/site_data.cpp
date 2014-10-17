// ==========================================================================
//                                  ANISE
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

#include "site_data.h"

#include <fstream>
#include <set>
#include <sstream>

#include <seqan/bam_io.h>
#include <seqan/misc/misc_name_store_cache.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include "asm/overlapper.h"

#include "shared/streaming_exception.h"

#include "anise/file_name_tokens.h"
#include "anise/bam_tag_names.h"

namespace {  // anonymous namespace

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

// Helper class, helps wrapping name store and caching stuff.

class ScaffoldLoader
{
public:
    ScaffoldLoader(TemporaryFileManager & tmpMgr, int siteID, int stepNo) :
            tmpMgr(tmpMgr), siteID(siteID), stepNo(stepNo), refNamesCache(refNames)
    {}

    void run(Scaffold & scaffold)
    {
        // Set scaffold.siteName.
        std::stringstream ns;
        ns << "site_" << siteID;
        scaffold.siteName = ns.str();

        // Load sequences and build refNames and refNamesCache.
        loadSequences(scaffold);
    }

private:
    // Load sequences and names and build the name store cache.
    void loadSequences(Scaffold & scaffold);

    // TemporaryFileManager used for loading the scaffold.
    TemporaryFileManager & tmpMgr;
    // Current site id.
    int siteID;
    // Current step number.
    int stepNo;

    // Name store for reference names and the cache.
    seqan::StringSet<seqan::CharString> refNames;
    seqan::NameStoreCache<seqan::StringSet<seqan::CharString>> refNamesCache;
};

void ScaffoldLoader::loadSequences(Scaffold & scaffold)
{
    std::string path = tmpMgr.fileName(SCAFFOLD_SEQS_TOKEN, SCAFFOLD_SEQS_EXT, siteID, stepNo);
    seqan::SeqFileIn ss;
    if (!open(ss, path.c_str()))
        throw AniseIOException() << "Could not open scaffold seqs file " << path << " for reading.";
    if (atEnd(ss))
        return;  // empty scaffold
    try
    {
        readRecords(scaffold.refNames, scaffold.seqs, ss);
    }
    catch (seqan::ParseError const & e)
    {
        throw AniseIOException() << "Problem reading from " << path;
    }

    // Parse out spans_insert info.
    resize(scaffold.scaffoldInfos, length(scaffold.refNames));
    for (unsigned i = 0; i < length(scaffold.refNames); ++i)
    {
        std::string s = toCString(scaffold.refNames[i]);
        scaffold.scaffoldInfos[i].anchoredLeft = (s.find("anchored_left=yes") != std::string::npos);
        scaffold.scaffoldInfos[i].anchoredRight = (s.find("anchored_right=yes") != std::string::npos);
        scaffold.scaffoldInfos[i].spanning = (s.find("spanning=yes") != std::string::npos);
        trimAfterSpace(scaffold.refNames[i]);
    }

    // Build the name store cache.
    refNames = scaffold.refNames;
    refresh(refNamesCache);
}

// ----------------------------------------------------------------------------
// Class ScaffoldSaver
// ----------------------------------------------------------------------------

class ScaffoldSaver
{
public:
    ScaffoldSaver(TemporaryFileManager & tmpMgr, int siteID, int stepNo) :
            tmpMgr(tmpMgr), siteID(siteID), stepNo(stepNo)
    {}

    void run(Scaffold const & scaffold) const
    {
        saveSequences(scaffold);
    }

private:
    void saveSequences(Scaffold const & scaffold) const;

    // TemporaryFileManager used for loading the scaffold.
    TemporaryFileManager & tmpMgr;
    // Current site id.
    int siteID;
    // Current step number.
    int stepNo;
};

void ScaffoldSaver::saveSequences(Scaffold const & scaffold) const
{
    // Create masked sequences.
    seqan::StringSet<seqan::CharString> seqs;
    for (unsigned i = 0; i < length(scaffold.seqs); ++i)
    {
        appendValue(seqs, scaffold.seqs[i]);
        if (i < length(scaffold.masks))
            for (unsigned j = 0, jEnd = std::min(length(scaffold.masks[i]), length(scaffold.seqs[i])); j < jEnd; ++j)
                if (!scaffold.masks[i][j])
                    seqs[i][j] = tolower(seqs[i][j]);
    }

    // Create tags for anchoring and spanning.
    seqan::StringSet<seqan::CharString> refNames = scaffold.refNames;
    char const * label[] = { "no", "yes" };
    for (unsigned i = 0; (i < length(refNames)) && (i < length(scaffold.scaffoldInfos)); ++i)
    {
        append(refNames[i], " anchored_left=");
        append(refNames[i], label[!!scaffold.scaffoldInfos[i].anchoredLeft]);
        append(refNames[i], " anchored_right=");
        append(refNames[i], label[!!scaffold.scaffoldInfos[i].anchoredRight]);
        append(refNames[i], " spanning=");
        append(refNames[i], label[!!scaffold.scaffoldInfos[i].spanning]);
    }

    std::string path = tmpMgr.fileName(SCAFFOLD_SEQS_TOKEN, SCAFFOLD_SEQS_EXT, siteID, stepNo);
    seqan::SeqFileOut ss;
    if (!open(ss, path.c_str()))
        throw AniseIOException() << "Could not open scaffold seqs file " << path << " for writing.";
    try
    {
        writeRecords(ss, refNames, seqs);
    }
    catch (seqan::IOError const & e)
    {
        throw AniseIOException() << "Problem writing to " << path;
    }
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class Scaffold
// ----------------------------------------------------------------------------

void Scaffold::load(TemporaryFileManager & tmpMgr, int siteID, int stepNo)
{
    ScaffoldLoader loader(tmpMgr, siteID, stepNo);
    loader.run(*this);
}

void Scaffold::save(TemporaryFileManager & tmpMgr, int siteID, int stepNo) const
{
    ScaffoldSaver saver(tmpMgr, siteID, stepNo);
    saver.run(*this);
}

// ----------------------------------------------------------------------------
// Class ReadSet
// ----------------------------------------------------------------------------

void ReadSet::load(TemporaryFileManager & manager, int siteID, int stepNo, Scaffold const & scaffold)
{
    std::string path = manager.fileName(READS_TOKEN, READS_EXT, siteID, stepNo);
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, path.c_str()))
        throw AniseIOException() << "Could not open read set SAM file " << path << " for reading.";

    // TODO(holtgrew): Check that the reference name are consistent with scaffold.
    (void)scaffold;

    // Memoize read names for left and right so we do not load reads twice.  Read mapping can map the same read multiple
    // times to a site and we ignore all but the first match.
    std::set<seqan::CharString> seenL, seenR;

    seqan::BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(record, bamFileIn);
        }
        catch (seqan::ParseError const & e)
        {
            throw AniseIOException() << "Problem reading record from " << path;
        }
        if (hasFlagSecondary(record))
            continue;  // Skip secondary records.
        if ((hasFlagFirst(record) && seenL.count(record.qName)) ||
            (hasFlagLast(record) && seenR.count(record.qName)))
            continue;  // Skip, already seen before.
        // In any case, mark read as seen now.
        if (hasFlagFirst(record))
        {
            seenL.insert(record.qName);
        }
        else if (hasFlagLast(record))
        {
            seenR.insert(record.qName);
        }
        else
        {
            seenL.insert(record.qName);
            seenR.insert(record.qName);
        }

        // Ignore if it has "aD" tag set to 1; was disabled in disableReadAndMate().  Store in removedReads, though, so
        // we can write out a "blocker" record with "aD" tag set again.
        seqan::BamTagsDict tagsDict(record.tags);
        unsigned idx = 0;
        int value = 0;
        if (findTagKey(idx, tagsDict, KEY_ANISE_REMOVED) && extractTagValue(value, tagsDict, idx) && value)
        {
            removedReads.insert(record.qName);
            continue;
        }

        // Also ignore if orphan map, we require at least one anchor.
        if (hasFlagUnmapped(record) && hasFlagNextUnmapped(record))
            continue;

        appendValue(seqs, record.seq);
        bamRecords.emplace_back(record);
    }

    SEQAN_ASSERT_EQ(bamRecords.size() % 2u, 0u);
    SEQAN_ASSERT_EQ(length(seqs) % 2u, 0u);

    computeStats();
}

void ReadSet::save(TemporaryFileManager & manager, int siteID, int stepNo, Scaffold const & scaffold) const
{
    SEQAN_ASSERT_EQ(bamRecords.size() % 2u, 0u);
    SEQAN_ASSERT_EQ(length(seqs) % 2u, 0u);

    std::string path = manager.fileName(READS_TOKEN, READS_EXT, siteID, stepNo);
    seqan::BamFileOut bamFileOut;
    if (!open(bamFileOut, path.c_str()))
        throw AniseIOException() << "Could not open read set SAM file " << path << " for writing.";

    auto trimmedNames = scaffold.refNames;
    for (unsigned i = 0; i < length(trimmedNames); ++i)
        trimAfterSpace(trimmedNames[i]);

    // Build header.
    seqan::BamHeader bamHeader;

    // Fill header records.
    seqan::BamHeaderRecord hd;
    hd.type = seqan::BAM_HEADER_FIRST;
    setTagValue("VN", "1.5", hd);
    appendValue(bamHeader, hd);

    for (unsigned i = 0; i < length(scaffold.seqs); ++i)
    {
        seqan::BamHeaderRecord sq;
        sq.type = seqan::BAM_HEADER_REFERENCE;
        setTagValue("SN", trimmedNames[i], sq);
        std::stringstream ss;
        ss << length(scaffold.seqs[i]);
        setTagValue("LN", ss.str().c_str(), sq);
        appendValue(bamHeader, sq);
    }

    seqan::BamHeaderRecord co;
    co.type = seqan::BAM_HEADER_COMMENT;
    setTagValue("  ", "Temporary read set file file created by ANISE.", co);
    appendValue(bamHeader, co);

    writeRecord(bamFileOut, bamHeader);

    // Do not write out records of mates with one or more reads that are unaligned and have been mapped too many
    // steps before.  We first collect their names.
    int const MAX_STEPS = 2;
    std::set<seqan::CharString> toRemove;
    for (auto const & record : bamRecords)
    {
        if (!hasFlagUnmapped(record) && !hasFlagNextUnmapped(record))
            continue;
        seqan::BamTagsDict tagsDict(const_cast<seqan::CharString &>(record.tags));
        unsigned idx = 0;
        int mappedInStep = 0;
        if (findTagKey(idx, tagsDict, KEY_BIRTH_STEP) && !extractTagValue(mappedInStep, tagsDict, idx) &&
            (stepNo - mappedInStep > MAX_STEPS))
            toRemove.insert(record.qName);
    }

    // Write out alignment records.
    for (auto const & record : bamRecords)
        if (!toRemove.count(record.qName))
        {
            try
            {
                writeRecord(bamFileOut, record);
            }
            catch (seqan::IOError const & e)
            {
                throw AniseIOException() << "Problem writing record to " << path;
            }
        }

    // Write out records for records for removed reads.  By keeping these in, subsequent mappings will be ignored.
    seqan::BamAlignmentRecord removedRecord;
    seqan::BamTagsDict tagsDict(removedRecord.tags);
    setTagValue(tagsDict, KEY_ANISE_REMOVED, 1);
    removedRecord.flag = seqan::BAM_FLAG_UNMAPPED;
    for (auto const & name : removedReads)
    {
        removedRecord.qName = name;
        try
        {
            writeRecord(bamFileOut, removedRecord);
        }
        catch (seqan::IOError const & e)
        {
            throw AniseIOException() << "Problem writing record to " << path;
        }
    }
}

void ReadSet::disableReadAndMate(unsigned readID)
{
    // Will be inserted into "removedReads" when read and written out from this the next time.
    for (unsigned i = 2 * (readID / 2); i <= 2 * (readID / 2) + 1; ++i)
    {
        seqan::BamTagsDict tagsDict(bamRecords[i].tags);
        setTagValue(tagsDict, KEY_ANISE_REMOVED, 1);  // "aD" == "disabled by ANISE"
    }
}

void ReadSet::computeStats()
{
    numOrphans = 0;
    long long sum = 0;
    for (auto const & record : bamRecords)
    {
        numOrphans += (hasFlagFirst(record) && hasFlagUnmapped(record) && hasFlagNextUnmapped(record));
        sum += length(record.seq);
    }
    avgReadLength = (int)(1.0 * sum / bamRecords.size() + 0.5);
}

// ----------------------------------------------------------------------------
// Class SiteData
// ----------------------------------------------------------------------------

void SiteData::load(TemporaryFileManager & tmpMgr, int siteID, int stepNo)
{
    state.load(tmpMgr, siteID);
    if (stepNo != -1)
        state.stepNo = stepNo;

    scaffold.load(tmpMgr, state.siteID, state.stepNo);
    readSet.load(tmpMgr, state.siteID, state.stepNo, scaffold);
}

void SiteData::save(TemporaryFileManager & tmpMgr)
{
    state.save(tmpMgr);
    scaffold.save(tmpMgr, state.siteID, state.stepNo);
    readSet.save(tmpMgr, state.siteID, state.stepNo, scaffold);
}
