// ==========================================================================
//                                   ANISE
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

#include "read_mapping_substep.h"

#include <atomic>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include "shared/progress_indicator.h"
#include "shared/streaming_exception.h"

#include "anise/anise_options.h"
#include "anise/app_state.h"
#include "anise/bam_tag_names.h"
#include "anise/file_name_tokens.h"
#include "anise/parallel_utils.h"
#include "anise/read_mapping_impl.h"
#include "anise/site_data.h"
#include "anise/site_state.h"
#include "anise/temporary_file_manager.h"

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
// Class ReadMappingData
// ----------------------------------------------------------------------------

// Data structure for input/output of read mapping.
//
// Is aware of disabled reads.

class ReadMappingData
{
public:
    // The read sequences.
    seqan::StringSet<seqan::Dna5String> reads;
    // The qualities.
    seqan::StringSet<seqan::CharString> quals;
    // The reads for read mapping.
    seqan::StringSet<seqan::CharString> readNames;
    // Mapping from index in reads/readNames to index in input FASTQ file and activeReads/nextActiveReads.
    seqan::String<unsigned> globalReadIds;
    // The result from read mapping.
    std::vector<seqan::BamAlignmentRecord> bamRecords;

    void mapBatch(seqan::StringSet<seqan::Dna5String> & contigs,
                  AniseOptions const & options);
};

void ReadMappingData::mapBatch(seqan::StringSet<seqan::Dna5String> & contigs,
                               AniseOptions const & options)
{
    clear(bamRecords);
    FilteringAllMapper mapper(FilteringAllMapper::Options(
            /*errorRate=*/0.01 * options.readMappingErrorRate,
            /*maxMatches=*/100,
            /*debug=*/(options.verbosity >= 3)));
    mapper.run(bamRecords, contigs, reads, quals, readNames);
}

// ----------------------------------------------------------------------------
// Class ReadLoader
// ----------------------------------------------------------------------------

// Responsible for loading reads in a thread-safe fashion.

class ReadLoader
{
public:
    // Reference to the SequenceStream to load sequences from.
    seqan::SequenceStream & readsIn;

    // The number of reads read from readsIn and the previous number of read reads.
    int numReadReads;
    int prevNumReadReads;

    // The mutex to use for locking the loader.
    std::mutex mutex;

    // Whether or not we are done here.
    std::atomic<bool> _done;

    ReadLoader(seqan::SequenceStream & readsIn) :
            readsIn(readsIn), numReadReads(0), prevNumReadReads(0), _done(false)
    {}

    // Load next batch into result.
    std::unique_ptr<ReadMappingData> loadNextBatch(int batchSize);

    bool done() const
    {
        return _done;
    }
};

std::unique_ptr<ReadMappingData> ReadLoader::loadNextBatch(int batchSize)
{
    std::unique_ptr<ReadMappingData> result(new ReadMappingData);

    std::lock_guard<std::mutex> lock(mutex);
    if (done())
        return result;  // Already done, short-circuit.

    seqan::CharString id, qual;
    seqan::Dna5String seq;
    prevNumReadReads = numReadReads;
    for (int i = 0, j = 0; j < batchSize * 2 && !atEnd(readsIn); ++i, ++numReadReads)
    {
        if (readRecord(id, seq, qual, readsIn) != 0)
            throw AniseIOException() << "Problem reading reads.";
        appendValue(result->globalReadIds, prevNumReadReads + i);
        appendValue(result->reads, seq);
        appendValue(result->quals, qual);
        appendValue(result->readNames, id);
        j++;
    }

    if (atEnd(readsIn))
        _done = true;

    return result;
}

// --------------------------------------------------------------------------
// Class SiteContextData
// --------------------------------------------------------------------------

// Data with BamIOContext for one site.

class SiteContextData
{
    // The contig names.
    seqan::StringSet<seqan::CharString> nameStore;
    // Name store cache for BamIOContext.
    seqan::NameStoreCache<decltype(nameStore)> nameStoreCache;

public:

    // The BamIOContext that is used when writing out the record.
    seqan::BamIOContext<decltype(nameStore)> context;

    SiteContextData() : nameStoreCache(nameStore), context(nameStore, nameStoreCache)
    {}

    // Register contig with the given name.
    void addContig(seqan::CharString const & name)
    {
        appendName(nameStore, name, nameStoreCache);
    }
};

// --------------------------------------------------------------------------
// Class ReadMappingContigsData
// --------------------------------------------------------------------------

// The contigs together with mapping between contig and site.

class ReadMappingContigsData
{
public:

    ReadMappingContigsData(AppState const & appState) : appState(appState)
    {}

    // Load the contigs from temporary files.  Builds the SiteContextData record at the same time.
    void load(std::vector<std::unique_ptr<SiteContextData>> & contextData,
              TemporaryFileManager & tmpMgr,
              AniseOptions const & options);

    // Map global contig ID to site number.
    int siteForGlobalContig(int contigID) const
    {
        if (contigID == -1)
            return -1;
        return contigToSite[contigID];
    }

    // Map global contig ID to site number.
    int globalRefIDToLocal(int contigID) const
    {
        if (contigID == -1)
            return -1;
        return globalToLocal[contigID];
    }

    // Return the contig sequences.
    seqan::StringSet<seqan::Dna5String> const & contigs() const
    {
        return _contigs;
    }

private:

    // Mapping from global contig id (index in contigs) to site-local contig id.
    std::vector<int> globalToLocal;
    // Mapping from global contig id (index in contigs) to site (e.g. index in contextData).
    std::vector<int> contigToSite;
    // The contig sequences.
    seqan::StringSet<seqan::Dna5String> _contigs;
    // The application state.
    AppState const & appState;
};

void ReadMappingContigsData::load(std::vector<std::unique_ptr<SiteContextData>> & contextData,
                                  TemporaryFileManager & tmpMgr,
                                  AniseOptions const & options)
{
    seqan::CharString id, seq;

    for (int siteID = 0; siteID < appState.numSites; ++siteID)
    {
        (void)options;
        // if (options.debugSiteID != -1 && options.debugSiteID != siteID)
        //     continue;  // Skip, not selected for debugging.

        // Load site state.
        AssemblySiteState siteState;
        siteState.load(tmpMgr, siteID);
        if (!siteState.active)
            continue;  // Skip if site not active any more.

        // Create new ReadMappingContigsData.
        contextData[siteID].reset(new SiteContextData);

        // Load contigs.
        std::string path = tmpMgr.fileName(SCAFFOLD_SEQS_TOKEN, SCAFFOLD_SEQS_EXT, siteID, appState.stepNo);
        seqan::SequenceStream contigsIn(path.c_str());
        if (!isGood(contigsIn))
            throw AniseIOException() << "Problem opening scaffold seqs file " << path;

        for (int contigID = 0; !atEnd(contigsIn); ++contigID)
        {
            if (readRecord(id, seq, contigsIn) != 0)
                throw AniseIOException() << "Problem reading from scaffold seqs file " << path;
            trimAfterSpace(id);

            // std::cerr << "SEQ\t" << seq << "\n";

            // Convert lower case sequence characters to Ns (soft to hard masking).
            std::transform(begin(seq, seqan::Standard()), end(seq, seqan::Standard()),
                           begin(seq, seqan::Standard()),
                           [](char c) { return islower(c) ? 'N' : c; });

            // std::cerr << "SEQ\t" << seq << "\n";

            contextData[siteID]->addContig(id);
            appendValue(_contigs, seq);
            globalToLocal.push_back(contigID);
            contigToSite.push_back(siteID);
        }
    }
}

// ----------------------------------------------------------------------------
// Class ResultDistributor
// ----------------------------------------------------------------------------

// Is responsible for thread-safe distribution of read mapping results to the reads TSV data files.

class ResultDistributor
{
    // The context for writing out.
    std::vector<std::unique_ptr<SiteContextData>> const & contexts;
    // Information about the contigs, includes mapping from global contig id to site id etc.
    ReadMappingContigsData const & contigsData;
    // Global options.
    AniseOptions const & options;
    // The application state.
    AppState const & appState;

    // The mutex to use for locking the distributor.
    std::mutex mutex;

public:

    // How many reads mapped on the given site.
    std::vector<int> newMappedOnSite;
    // The number of freshly mapped and distributed records.
    int numNewMapped;

    ResultDistributor(std::vector<std::unique_ptr<SiteContextData>> const & contexts,
                      ReadMappingContigsData const & contigsData,
                      AniseOptions const & options,
                      AppState const & appState) :
            contexts(contexts), contigsData(contigsData), options(options), appState(appState),
            newMappedOnSite(appState.numSites, 0), numNewMapped(0)
    {}

    // Distribute mapping results and mark as inactive if necessary.
    void run(ReadMappingData & results, TemporaryFileManager & tmpMgr);
};

void ResultDistributor::run(ReadMappingData & results,
                            TemporaryFileManager & tmpMgr)
{
    // The records from results are sorted by readID.

    std::lock_guard<std::mutex> lock(mutex);

    if (options.verbosity >= AniseOptions::VERBOSE)
        std::cerr << "Distributing results for step " << appState.stepNo << "\n";

    // Distribute records to buffers for each site, records will remain sorted by read ID.  We will translate the contig
    // ID from global to local on the fly.
    std::vector<std::vector<seqan::BamAlignmentRecord>> siteRecords(appState.numSites);
    for (auto & record : results.bamRecords)
    {
        int siteID = contigsData.siteForGlobalContig(record.rID);
        record.rID = contigsData.globalRefIDToLocal(record.rID);
        siteRecords[siteID].resize(siteRecords[siteID].size() + 1);
        using std::swap;
        swap(siteRecords[siteID].back(), record);
    }

    // Process records site-wise.
    for (unsigned siteID = 0; siteID < siteRecords.size(); ++siteID)
    {
        auto const & records = siteRecords[siteID];
        if (records.empty())
            continue;  // skip if empty

        std::vector<seqan::BamAlignmentRecord> out;  // records for this site

        auto it = records.begin(), itNext = std::next(records.begin());
        while (it != records.end())
        {
            bool itIsFirst = !(it->_qID % 2u);  // whether it points to first
            for (; itNext != records.end() && it->_qID == itNext->_qID; ++itNext)
                continue;  // fast forward itNext to match of next record

            if (itNext != records.end() && itIsFirst && it->_qID + 1 == itNext->_qID &&
                hasFlagRC(*it) != hasFlagRC(*itNext))
            {
                // Allocate result records and get shortcuts.
                out.resize(out.size() + 2);
                auto & record1 = out[out.size() - 2];
                auto & record2 = out[out.size() - 1];

                // Start out with *it and *itNext;
                record1 = *it;
                record2 = *itNext;

                // Fix pointers to next positions.
                record1.rNextId = record2.rID;
                record1.pNext = record2.beginPos;
                record2.rNextId = record1.rID;
                record2.pNext = record1.beginPos;

                // Update flags.
                record1.flag |= seqan::BAM_FLAG_MULTIPLE | seqan::BAM_FLAG_FIRST;
                record2.flag |= seqan::BAM_FLAG_MULTIPLE | seqan::BAM_FLAG_LAST;
                if (hasFlagRC(*it))
                    record2.flag |= seqan::BAM_FLAG_NEXT_RC;
                else
                    record1.flag |= seqan::BAM_FLAG_NEXT_RC;

                // Set TLEN
                if (record1.rID == record2.rID)
                {
                    auto itL = &record1, itR = &record2;
                    if (itL->beginPos > itR->beginPos)
                        std::swap(itL, itR);
                    itL->tLen = itR->beginPos + getAlignmentLengthInRef(*itR) - itL->beginPos;
                    itR->tLen = -itL->tLen;
                }
            }
            else
            {
                // We have found only one match for it.
                unsigned readID = it->_qID;
                unsigned otherID = (readID % 2u) ? (readID - 1) : (readID + 1);

                // Allocate result records and get shortcuts.
                out.resize(out.size() + 2);
                auto & record1 = out[out.size() - 2];
                auto & record2 = out[out.size() - 1];

                // Fill first (aligned) record.
                record1.seq = results.reads[readID];
                record1.qual = results.quals[readID];
                record1._qID = readID;
                record1.qName = results.readNames[readID];
                SEQAN_ASSERT_EQ(record1.qName, it->qName);
                record1.flag = seqan::BAM_FLAG_MULTIPLE | seqan::BAM_FLAG_NEXT_UNMAPPED;
                if (hasFlagRC(*it))
                {
                    record1.flag |= seqan::BAM_FLAG_RC;
                    reverseComplement(record1.seq);
                    reverse(record1.qual);
                }
                else
                {
                    record1.flag |= seqan::BAM_FLAG_NEXT_RC;
                }
                if (hasFlagFirst(*it))
                    record1.flag |= seqan::BAM_FLAG_FIRST;
                else
                    record1.flag |= seqan::BAM_FLAG_LAST;
                record1.rID = record1.rNextId = it->rID;
                record1.beginPos = record1.pNext = it->beginPos;
                record1.cigar = it->cigar;
                record1.tLen = 0;
                record1.tags = it->tags;

                // Fill second (unaligned record).
                record2.seq = results.reads[otherID];
                record2.qual = results.quals[otherID];
                record2._qID = otherID;
                record2.qName = results.readNames[otherID];
                SEQAN_ASSERT_EQ(record2.qName, it->qName);
                record2.flag = seqan::BAM_FLAG_MULTIPLE | seqan::BAM_FLAG_UNMAPPED;
                if (!hasFlagRC(*it))
                {
                    record2.flag |= seqan::BAM_FLAG_RC;
                    reverseComplement(record2.seq);
                    reverse(record2.qual);
                }
                else
                {
                    record2.flag |= seqan::BAM_FLAG_NEXT_RC;
                }
                if (hasFlagFirst(*it))
                    record2.flag |= seqan::BAM_FLAG_LAST;
                else
                    record2.flag |= seqan::BAM_FLAG_FIRST;
                record2.rID = record2.rNextId = it->rID;
                record2.beginPos = record2.pNext = it->beginPos;
                record2.tLen = 0;

                // Swap records if first not flagged first.
                using std::swap;
                if (!hasFlagFirst(record1))
                    swap(record1, record2);
            }

            // Fast forward it to next match for next id if was a pair match (regardless of whether we interpreted it as
            // a pair match or not.
            if (itNext != records.end() && it->_qID == itNext->_qID + 1)
                for (it = itNext; it != records.end() && it->_qID == itNext->_qID; ++it)
                    continue;

            // Set it to itNext, will cause itNext to search for next record in next round.
            it = itNext;
        }

        // Set "mS" tag to the step that the reads were mapped in.
        int stepNo = appState.stepNo;
        std::for_each(out.begin(), out.end(),
                      [stepNo](seqan::BamAlignmentRecord & record) {
                          seqan::BamTagsDict tagsDict(record.tags);
                          setTagValue(tagsDict, KEY_BIRTH_STEP, stepNo);
                      });

        // Write out records for the current site.
        if (!out.empty())
        {
            // Open SAM file.
            std::fstream f;
            tmpMgr.open(f, std::ios::binary | std::ios::out | std::ios::app, READS_TOKEN, READS_EXT,
                        siteID, appState.stepNo);
            // Write to SAM file.
            for (auto const & record : out)
                if (write2(f, record, contexts[siteID]->context, seqan::Sam()) != 0)
                    throw std::runtime_error("Problem writing to reads SAM file.");
        }
    }
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class ReadMappingSubstepImpl
// --------------------------------------------------------------------------

class ReadMappingSubstepImpl
{
public:

    ReadMappingSubstepImpl(TemporaryFileManager & manager, int stepNo, AniseOptions const & options) :
            tmpMgr(manager), options(options)
    {
        SEQAN_CHECK(stepNo > 0, "Cannot map in zero-th step.");
    }

    // Run read mapping state.
    void run()
    {
        // This function performs all the loading and initialization of the required data structures and hands the
        // actual mapping off to mapReads().

        // Load application state.
        appState.load(tmpMgr);

        // The site contexts.  Constructed here since the size mut be known in advance since it contains
        // std::unique_ptr<> objects and we cannot do std::vector<>.resize() here.
        std::vector<std::unique_ptr<SiteContextData>> contextData(appState.numSites);
        // Load the read mapping contigs and make BamIOContext for the reference names.
        ReadMappingContigsData contigsData(appState);
        contigsData.load(contextData, tmpMgr, options);
        // Create result distributor.
        ResultDistributor distributor(contextData, contigsData, options, appState);

        if (options.verbosity >= AniseOptions::NORMAL)
            std::cerr << "Mapping " << appState.numOrphans << " against " << lengthSum(contigsData.contigs())
                      << " reference bases\n";

        // Perform the actual mapping.
        mapReads(distributor, contigsData);
    }

private:

    // Disable reads that do not fulfill quality requirements (too many Ns, too short).
    void disableBadReads(unsigned siteID, unsigned stepNo)
    {
        SiteData siteData;
        siteData.load(tmpMgr, siteID, stepNo);

        for (unsigned readID = 0; readID < siteData.readSet.bamRecords.size(); ++readID)
        {
            auto & record = siteData.readSet.bamRecords[readID];
            unsigned const MAX_NS = 3u;
            unsigned const MIN_LEN = 60u;
            unsigned numNs = 0;
            for (auto c : record.seq)
                numNs += (c == 'N');
            if (length(record.seq) < MIN_LEN || numNs > MAX_NS)
                siteData.readSet.disableReadAndMate(readID);
        }

        siteData.save(tmpMgr);
    }

    // Determine number of reads mapped in step stepNo on site with siteID.
    int newReadsOnSite(int siteID, int stepNo)
    {
        int result = 0;

        SiteData siteData;
        siteData.load(tmpMgr, siteID, stepNo);

        std::set<seqan::CharString> seen;

        for (unsigned recordID = 0; recordID < siteData.readSet.bamRecords.size(); recordID += 2)
        {
            auto & record = siteData.readSet.bamRecords[recordID];
            if (hasFlagSecondary(record) || seen.count(record.qName))
                continue;  // skip, seen already
            seen.insert(record.qName);

            seqan::BamTagsDict tagsDict(record.tags);
            unsigned idx = 0;
            int mappedStep = 0;
            if (!findTagKey(idx, tagsDict, KEY_BIRTH_STEP) || !extractTagValue(mappedStep, tagsDict, idx))
                continue;

            result += (mappedStep >= stepNo);
        }

        return result;
    }

    void mapReads(ResultDistributor & distributor, ReadMappingContigsData const & contigsData)
    {
        // Open orphans input file.
        std::string inPath = tmpMgr.fileName("orphans", ".fq");
        seqan::SequenceStream inReads(inPath.c_str());
        if (!isGood(inReads))
            throw AniseIOException() << "Could not open orphans file for reading " << inPath;
        ReadLoader loader(inReads);

        // Create progress bar.
        // TODO(holtgrew): Using real active reads.
        ProgressBar pb(std::cerr, 0, appState.numOrphans, (options.verbosity == AniseOptions::NORMAL));
        pb.setLabel("  mapping reads");
        pb.updateDisplay();

        // This lambda implements the batch-wise read mapping loop
        auto readMappingLoop = [&](){
            // Copy contigs since each thread needs to reverse-complement it.
            auto contigsCopy = contigsData.contigs();

            // TODO(holtgrew): Parallelize and lock-ize.
            while (!atEnd(inReads))
            {
                auto batch = loader.loadNextBatch(options.readMappingBatchSize);
                batch->mapBatch(contigsCopy, options);
                distributor.run(*batch, tmpMgr);

                // Update progress bar.
                pb.advanceBy(length(batch->reads));
            }
        };

        // Execute read mapping loop in fork-join manner.
        forkJoin(options.numThreads, readMappingLoop);
        pb.finish();

        // Check reads.
        {
            // Progress bar.
            ProgressBar pb(std::cerr, 0, appState.numSites, (options.verbosity == AniseOptions::NORMAL));
            pb.setLabel("  check alignment counts");
            pb.updateDisplay();

            // Job generator to atomically generate jobs for each site.
            JobGenerator gen(appState.numSites);

            // This lambda contains the code for on thread, grabs site ids from gen.
            std::mutex mutex;
            int const MIN_NEW_MAPPED = 5;
            int maxNewMapped = 0, sumNewMapped = 0;
            auto checkLoop = [&]() {
                int siteID = -1;
                while ((siteID = gen()) != -1)
                {
                    AssemblySiteState siteState;
                    siteState.load(tmpMgr, siteID);
                    if (siteState.active)  // ignore inactive sites
                    {
                        disableBadReads(siteState.siteID, siteState.stepNo);
                        int count = newReadsOnSite(siteState.siteID, siteState.stepNo);
                        {
                            std::lock_guard<std::mutex> lock(mutex);
                            maxNewMapped = std::max(maxNewMapped, count);
                            sumNewMapped += count;
                        }
                        if (siteState.active && count < MIN_NEW_MAPPED)
                        {
                            std::stringstream ss;
                            ss << "deactivating: too few new new alignments: (" << count << " < "
                               << MIN_NEW_MAPPED << ")";
                            siteState.comment = ss.str();
                            siteState.active = false;
                            siteState.save(tmpMgr);
                        }
                    }
                    pb.advanceBy(1);
                }
            };

            forkJoin(options.numThreads, checkLoop);

            pb.finish();

            // Get largest number of newly mapped reads.  If there is no site with a certain minimum of new reads then
            // switch app state to FINISHING.
            std::cerr << "  => found " << sumNewMapped << " new pair mappings (up to "
                      << maxNewMapped << " per site)\n";
            if (maxNewMapped < MIN_NEW_MAPPED)
            {
                appState.superStep = AppState::SuperStep::FINISHING;
                appState.save(tmpMgr);
            }
        }
    }

    // The current global app state.
    AppState appState;
    // Temporary file manager to use.
    TemporaryFileManager & tmpMgr;
    // The application configuration.
    AniseOptions const & options;
};

// --------------------------------------------------------------------------
// Class ReadMappingSubstep
// --------------------------------------------------------------------------

ReadMappingSubstep::ReadMappingSubstep(TemporaryFileManager & manager, int stepNo,AniseOptions const & options) :
        impl(new ReadMappingSubstepImpl(manager, stepNo, options))
{}

ReadMappingSubstep::~ReadMappingSubstep()
{}

void ReadMappingSubstep::run()
{
    impl->run();
}
