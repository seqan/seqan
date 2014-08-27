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
// Data types related to a site (scaffold and read set).
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SITE_DATA_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SITE_DATA_H_

#include <utility>
#include <set>
#include <vector>
#include <string>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include "anise/site_state.h"
#include "anise/temporary_file_manager.h"

#include "asm/overlapper.h"

// ============================================================================
// Forwards
// ============================================================================

class ConflictStore;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Scaffod
// ----------------------------------------------------------------------------

struct Scaffold
{
    struct ScaffoldInfo
    {
        // Whether the sequence was anchored left (step 0 forward read on the left of the assembly).
        bool anchoredLeft  { false };
        // Whether the sequence was anchored left (step 0 reverse read on the right of the assembly).
        bool anchoredRight { false };
        // Whether or not an s-t path was created with links and overlaps only.
        bool spanning { false };

        ScaffoldInfo() = default;
        ScaffoldInfo(bool anchoredLeft, bool anchoredRight, bool spanning) :
                anchoredLeft(anchoredLeft), anchoredRight(anchoredRight), spanning(spanning)
        {}
    };

    // The name of the site, generated from siteID on loading.
    std::string siteName;

    // Reference sequence names.  Used for identifying sequences in files and output to the user but unused in the
    // program itself.
    seqan::StringSet<seqan::CharString> refNames;
    // The actual scaffold sequences.
    seqan::StringSet<seqan::Dna5String> seqs;
    // Optionally, a masking for the sequences, the characters in seqs will be converted to lower case if the
    // corresponding mask value is false.
    seqan::StringSet<seqan::String<bool>> masks;
    // Optionally, a flag for each sequence whether or not the assembly met from both ends for this sequence.
    seqan::String<ScaffoldInfo> scaffoldInfos;

    // Load given the temporary file manager, site id and step no.
    void load(TemporaryFileManager & manager, int siteID, int stepNo);
    // Save given the temporary file manager, site id and step no.
    void save(TemporaryFileManager & manager, int siteID, int stepNo) const;
};

// ----------------------------------------------------------------------------
// Class ReadSet
// ----------------------------------------------------------------------------

// Stores alignments of the reads against a scaffold.  Read sets are serialized as SAM files and we keep the
// BamAlignmentRecord objects in the ReadSet.  Besides, we store the read sequences separately from the records in a
// Dna5String StringSet for easier access in alignments etc.
//
// We keep the invariant that we only store read pairs and store them consecutively.
//
// We append records to site read sets in the read mapping step.  These are appended as orphan record pairs

struct ReadSet
{
    ReadSet() : avgReadLength(-1), numOrphans(-1)
    {}

    // The extract read sequences.
    seqan::StringSet<seqan::Dna5String> seqs;
    // The BAM alignment records, also include the pairwise alignment to the scaffold reference sequences.
    std::vector<seqan::BamAlignmentRecord> bamRecords;
    // A list of read (pair) names that are ignored.  We will write out a dummy record so they remain removed from this
    // site.
    std::set<seqan::CharString> removedReads;

    // Average read length, computed on loading.
    int avgReadLength;
    // Number of orphan records, computed upon loading.
    int numOrphans;

    // Load given the temporary file manager, site id and step no.
    //
    // scaffold is passed for the reference names in the SAM file.
    void load(TemporaryFileManager & manager, int siteID, int stepNo, Scaffold const & scaffold);
    // Save given the temporary file manager, site id and step no.
    //
    // scaffold is passed for the reference names in the SAM file.
    void save(TemporaryFileManager & manager, int siteID, int stepNo, Scaffold const & scaffold) const;

    // Mark a read and its mate as invalid using the BAM tags.  This will only be interpreted after loading the read set
    // the next time since these reads will be skipped.
    void disableReadAndMate(unsigned readID);
    // TODO(holtgrew): Implement skipping of ignored reads.

    // Compute avgReadLength and numOrphans.
    void computeStats();
};

// ----------------------------------------------------------------------------
// Class SiteData
// ----------------------------------------------------------------------------

// Aggregation of scaffold and read set.

struct SiteData
{
    AssemblySiteState state;
    Scaffold scaffold;
    ReadSet readSet;

    SiteData(int siteID = -1) : state(siteID)
    {}

    // Load site data.  If stepNo is not given then load for the step as indicated in site state.
    void load(TemporaryFileManager & manager, int siteID, int stepNo = -1);
    void save(TemporaryFileManager & manager);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SITE_DATA_H_
