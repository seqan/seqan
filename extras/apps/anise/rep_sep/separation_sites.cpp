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

#include "separation_sites.h"

#include "rep_sep/rep_sep_options.h"

namespace rep_sep {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function includes()
// ----------------------------------------------------------------------------

// A wrapper for std::includes() that accepts two strings (rhs in lhs).

template <typename TValue, typename TSpec>
inline bool includes(seqan::String<TValue, TSpec> const & lhs,
                     seqan::String<TValue, TSpec> const & rhs)
{
    return std::includes(begin(lhs, seqan::Standard()), end(lhs, seqan::Standard()),
                         begin(rhs, seqan::Standard()), end(rhs, seqan::Standard()));
}

// ----------------------------------------------------------------------------
// Function setIntersection()
// ----------------------------------------------------------------------------

// A wrapper for std::set_intersection() that accepts containers.

template <typename TValue, typename TSpec>
inline void setIntersection(seqan::String<TValue, TSpec> & result,
                            seqan::String<TValue, TSpec> const & lhs,
                            seqan::String<TValue, TSpec> const & rhs)
{
    resize(result, std::min(length(lhs), length(rhs)), 0);
    resize(result, std::set_intersection(begin(lhs, seqan::Standard()), end(lhs, seqan::Standard()),
                                         begin(rhs, seqan::Standard()), end(rhs, seqan::Standard()),
                                         begin(result, seqan::Standard())) - begin(result, seqan::Standard()));
}

// ----------------------------------------------------------------------------
// Function setIntersectionSize()
// ----------------------------------------------------------------------------

// A wrapper for std::set_intersection() that computes the length of the intersection.

template <typename TValue, typename TSpec>
inline size_t setIntersectionSize(seqan::String<TValue, TSpec> const & lhs,
                                  seqan::String<TValue, TSpec> const & rhs)
{
    seqan::String<TValue, TSpec> tmp;
    resize(tmp, std::min(length(lhs), length(rhs)), 0);
    return std::set_intersection(begin(lhs, seqan::Standard()), end(lhs, seqan::Standard()),
                                 begin(rhs, seqan::Standard()), end(rhs, seqan::Standard()),
                                 begin(tmp, seqan::Standard())) - begin(tmp, seqan::Standard());
}

// ---------------------------------------------------------------------------
// Class SeparationSiteEnumerator
// ---------------------------------------------------------------------------

class SeparationSiteEnumerator
{
public:
    LocalVariationStore const * varStore;
    ReadSeparatorOptions const & options;

    SeparationSiteEnumerator(ReadSeparatorOptions const & options) : varStore(), options(options)
    {}

    SeparationSiteEnumerator(LocalVariationStore const & varStore,
                             ReadSeparatorOptions const & options) :
            varStore(&varStore), options(options)
    {}

    void run(seqan::String<SeparationSite> & sites)
    {
        using namespace seqan;
        for (unsigned i = 0; i < length(*varStore); ++i)
        {
            SeparationSite newSite;  // the new site to construct

            if (empty(sites))
            {
                bootstrapSite(newSite, i);
                reverse(newSite.columnIDs);
                checkSeparationSite(newSite, *varStore, -1);
                if (options.verbosity >= 3)
                {
                    std::cerr << "BOOTSTRAPPED\n";
                    newSite.print();
                }
            }
            else
            {
                if (options.verbosity >= 3)
                {
                    std::cerr << "EXTENDING\n";
                    newSite.print();
                }
                if (extendSites(sites, newSite, i))
                    continue;
                if (options.verbosity >= 3)
                {
                    std::cerr << "EXTENDED TO\n";
                    newSite.print();
                }
            }

            // Reverse, column ids are in reverse order by algorithm.
            reverse(newSite.columnIDs);

            // Check that the new site is not a subset of the previous subset i.e. that columns ore reads are added.
            //
            // This is a simple heuristic.
            // TODO(holtgrew): Extend sPre before expansion in extendSites?
            if (!empty(sites) && includes(back(sites).columnIDs, newSite.columnIDs) &&
                includes(back(sites).readIDs, newSite.readIDs))
                continue;

            if (options.verbosity >= 3)
            {
                std::cerr << "APPENDING NEW SITE\n";
                newSite.print();
            }

            // TODO(holtgrew): Why is the following necessary? Bug in algorithm? Ask Leon!
            if (length(newSite.readIDs) <= 1u)
                continue;  // Skip site with less than two reads, makes no sense to cluster!

            // Append new site to set of selected site.
            checkSeparationSite(newSite, *varStore, -1);
            appendValue(sites, newSite);
        }
    }

    // Construct initial site.
    void bootstrapSite(SeparationSite & site, unsigned i)
    {
        using namespace seqan;
        appendValue(site.columnIDs, i);
        site.readIDs = varStore->coveringReads[i];

        for (unsigned j = i + 1; j < length(*varStore); ++j)
        {
            seqan::Tag<seqan::Standard_> tag;
            seqan::StringSet<seqan::String<unsigned> > const & covReads = varStore->coveringReads;

            // Whether or not varStore->coveringReads[j] is a subset of varStore->coveringReads[i].
            bool isSubset = includes(covReads[i], covReads[j]);
            if (options.verbosity >= 3)
            {
                std::cerr << "covReads[i] ==";
                for (unsigned k = 0; k < length(covReads[i]); ++k)
                    std::cerr << " " << covReads[i][k];
                std::cerr << "\n";
                std::cerr << "covReads[j] ==";
                for (unsigned k = 0; k < length(covReads[j]); ++k)
                    std::cerr << " " << covReads[j][k];
                std::cerr << "\n";
            }
            if (isSubset)
            {
                appendValue(site.columnIDs, j);
                unsigned n = std::set_intersection(
                        begin(site.readIDs, tag), end(site.readIDs, tag),
                        begin(covReads[j], tag), end(covReads[j], tag),
                        begin(site.readIDs, tag)) - begin(site.readIDs, tag);
                resize(site.readIDs, n);
            }
            else
            {
                if (setIntersectionSize(site.readIDs, covReads[j]) == 0u)
                    return;  // No further overlap possible.
            }
        }

        if (options.verbosity >= 3)
        {
            std::cerr << "BOOTSTRAPPED\n";
            site.print();
        }
    }

    // Extend the current site.
    //
    // Returns true if we want to call "continue" a level above.
    bool extendSites(seqan::String<SeparationSite> & sites, SeparationSite & site, unsigned i)
    {
        using namespace seqan;
        // Compute read overlap of last added site and next column.  If this is larger than tauMin then the window added
        // by this i would is completely contained in the next one.
        if (i + 1 < length(*varStore))  // Computation makes no sense otherwise.
        {
            int sharedCount = setIntersectionSize(varStore->coveringReads[i + 1], back(sites).readIDs);
            if (sharedCount >= options.tauMin)
                return true;
        }

        if (options.verbosity >= 3)
        {
            std::cerr << "EXTENDING SITE\n";
        }

        // Add columns c_j to new site, starting with column c_i as long as the reads from c_j are a subset of the reads
        // of c_i or the intersection of the reads of c_i and c_j contains more than rMin elements.
        appendValue(site.columnIDs, i);
        for (int j = i - 1; j >= 0; --j)
        {
            bool isSubset = includes(varStore->coveringReads[i], varStore->coveringReads[j]);
            int sharedCount = setIntersectionSize(varStore->coveringReads[i], varStore->coveringReads[j]);

            if (isSubset || sharedCount >= options.rMin)
                appendValue(site.columnIDs, j);
            else
                break;
        }

        // Possibly expand s^{pre}.
        unsigned c1Pre = front(back(sites).columnIDs);
        int sharedCount = setIntersectionSize(varStore->coveringReads[c1Pre], varStore->coveringReads[i]);
        if (sharedCount < options.tauMin)
            expand(sites);

        // Add reads to new site.
        if (!empty(site) && !empty(site.columnIDs))
        {
            site.readIDs = varStore->coveringReads[front(site.columnIDs)];
            for (unsigned i = 1; i < length(site.columnIDs); ++i)
                setIntersection(site.readIDs, site.readIDs, varStore->coveringReads[site.columnIDs[i]]);
        }
        if (empty(site.readIDs))
            return true;

        return false;
    }

    // Successively create sites by removing the leftmost column of sPre = back(sites).  Each time the number of
    // participating reads is increased, add this as a new site.
    void expand(seqan::String<SeparationSite> & sites)
    {
        using namespace seqan;
        SeparationSite sPre = back(sites);

        // Number of reads participating in the previously added site (initialize with number of reads in sPre).
        unsigned numReads = length(sPre.readIDs);

        String<unsigned> tmp;  // buffer with ids of read intersection
        for (unsigned i = 1; i < length(sPre.columnIDs); ++i)
        {
            clear(tmp);
            // Intersect reads i..(length(sPre.columnIDs) - 1).
            if (i + 1 == length(sPre.columnIDs))
                tmp = varStore->coveringReads[sPre.columnIDs[i]];
            else
                setIntersection(tmp, varStore->coveringReads[sPre.columnIDs[i]],
                                varStore->coveringReads[sPre.columnIDs[i + 1]]);
            for (unsigned j = i + 2; j < length(sPre.columnIDs); ++j)
                setIntersection(tmp, tmp, varStore->coveringReads[sPre.columnIDs[j]]);

            if (length(tmp) > numReads)
            {
                SeparationSite site;
                site.readIDs = tmp;
                site.columnIDs = suffix(sPre.columnIDs, i);
                numReads = length(tmp);
                if (options.verbosity >= 3)
                {
                    std::cerr << "EXPANDED TO\n";
                    site.print();
                }
                appendValue(sites, site);
            }
        }
    }
};

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function checkSeparationSite()
// ----------------------------------------------------------------------------

void checkSeparationSite(SeparationSite const & site,
                         LocalVariationStore const & varStore,
                         int siteID)
{
    for (unsigned j = 0; j < length(site.columnIDs); ++j)
    {
        unsigned c = site.columnIDs[j];
#if SEQAN_READ_SEP_DEBUG
        if (!includes(varStore.coveringReads[c], site.readIDs))
        {
            std::cerr << "covering reads\t";
            for (unsigned i = 0; i < length(varStore.coveringReads[c]); ++i)
                std::cerr << "\t" << varStore.coveringReads[c][i];
            std::cerr << "\n";
            std::cerr << "read ids\t";
            for (unsigned i = 0; i < length(site.readIDs); ++i)
                std::cerr << "\t" << site.readIDs[i];
            std::cerr << "\n";
        }
#endif  // #if SEQAN_READ_SEP_DEBUG
        SEQAN_CHECK(includes(varStore.coveringReads[c], site.readIDs),
                    "site id = %u, column id = %u", siteID, c);
    }
}

void checkSeparationSites(seqan::String<SeparationSite> const & sites,
                          LocalVariationStore const & varStore)
{
#if SEQAN_READ_SEP_DEBUG
    std::cerr << "length(sites) == " << length(sites) << "\n";
#endif  // #if SEQAN_READ_SEP_DEBUG
    for (unsigned i = 0; i < length(sites); ++i)
        checkSeparationSite(sites[i], varStore, i);
}

// ----------------------------------------------------------------------------
// Function enumerateSeparationSites()
// ----------------------------------------------------------------------------

// Enumerate separation sites using Kuchenbecker's method.

void enumerateSeparationSites(seqan::String<SeparationSite> & sites,
                              LocalVariationStore const & varStore,
                              ReadSeparatorOptions const & options)
{
    SeparationSiteEnumerator sse(varStore, options);
    sse.run(sites);
#if SEQAN_ENABLE_DEBUG
    checkSeparationSites(sites, varStore);
#endif  // #if SEQAN_ENABLE_DEBUG
}

void enumerateSeparationSites(seqan::String<SeparationSite> & sites,
                              LocalVariationStore const & varStore)
{
    ReadSeparatorOptions options;
    enumerateSeparationSites(sites, varStore, options);
}

}  // namespace rep_sep
