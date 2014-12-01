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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_SEPARATION_SITES_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_SEPARATION_SITES_H_

#include <seqan/sequence.h>

#include "rep_sep/local_variation_store.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

struct ReadSeparatorOptions;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeparationSite
// ----------------------------------------------------------------------------

class SeparationSite
{
public:
    // The read ids that are part of this site.
    seqan::String<unsigned> readIDs;
    // The columns that are part of this site.
    seqan::String<unsigned> columnIDs;

    void print()
    {
        std::cerr << "reads  ";
        for (unsigned i = 0; i < length(readIDs); ++i)
            std::cerr << " " << readIDs[i];
        std::cerr << "\ncolumns";
        for (unsigned i = 0; i < length(columnIDs); ++i)
            std::cerr << " " << columnIDs[i];
        std::cerr << "\n";
    }

    unsigned numColumns() const
    {
        return length(columnIDs);
    }

    void clear()
    {
        seqan::clear(readIDs);
        seqan::clear(columnIDs);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

inline unsigned length(SeparationSite const & site)
{
    return length(site.columnIDs);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

inline bool empty(SeparationSite const & site)
{
    return !length(site.columnIDs);
}

// ----------------------------------------------------------------------------
// Function enumerateSeparationSites()
// ----------------------------------------------------------------------------

void enumerateSeparationSites(seqan::String<SeparationSite> & sites,
                              LocalVariationStore const & varStore,
                              ReadSeparatorOptions const & options);

void enumerateSeparationSites(seqan::String<SeparationSite> & sites,
                              LocalVariationStore const & varStore);

// ----------------------------------------------------------------------------
// Function checkSeparationSite()
// ----------------------------------------------------------------------------

void checkSeparationSite(SeparationSite const & site,
                         LocalVariationStore const & varStore,
                         int siteID = 0);

// ----------------------------------------------------------------------------
// Function checkSeparationSites()
// ----------------------------------------------------------------------------

void checkSeparationSites(seqan::String<SeparationSite> const & sites,
                          LocalVariationStore const & varStore);

// ----------------------------------------------------------------------------
// Function extractSiteStore()
// ----------------------------------------------------------------------------

void extractSiteStore(LocalVariationStore & siteStore,
                      SeparationSite & site,
                      LocalVariationStore const & varStore);

// ----------------------------------------------------------------------------
// Function extractSiteStore()
// ----------------------------------------------------------------------------

void buildSiteStrings(seqan::StringSet<seqan::String<LocalVariationStore::TConsensusAlphabet> > & result,
                      LocalVariationStore const & store,
                      SeparationSite const & site,
                      bool debug = false);

// ----------------------------------------------------------------------------
// Function compressRedundantReads()
// ----------------------------------------------------------------------------

// Replace the redundant reads (i.e. duplicated sequence over the deviating columns) by one read and build
// compressionMap.
//
// The store must already be filtered down to a site that we will also reconstruct here.
//
// At the same time a compressed separation site is created.
//
// TODO(holtgrew): Re-number reads?

void compressRedundantReads(LocalVariationStore & smallStore,
                            SeparationSite & smallSite,
                            LocalVariationStore const & varStore);

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_SEPARATION_SITES_H_
