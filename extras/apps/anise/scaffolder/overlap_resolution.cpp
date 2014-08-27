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

#include "overlap_resolution.h"

#include <iterator>
#include <sstream>
#include <string>

#include <seqan/sequence.h>
#include <seqan/file.h>

#include "asm/overlapper.h"

#include "scaffolder/scaffolding_result.h"

namespace scaffolder {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class simultaneous_iterator
// ----------------------------------------------------------------------------

template <typename FwdIt> class simultaneous_iterator {
public:
    simultaneous_iterator(FwdIt left_first, FwdIt right_first)
            : m_left_first(left_first), m_right_first(right_first) {}

    bool operator!=(const simultaneous_iterator& other) const {
        return std::make_pair(m_left_first, m_right_first) != std::make_pair(other.m_left_first, other.m_right_first);
    }

    simultaneous_iterator& operator++() {
        ++m_left_first;
        ++m_right_first;
        return *this;
    }

    typedef typename std::iterator_traits<FwdIt>::reference Ref;
    typedef std::pair<Ref, Ref> Pair;

    Pair operator*() const {
        return Pair(*m_left_first, *m_right_first); // NOT std::make_pair()!
    }

private:
    FwdIt m_left_first;
    FwdIt m_right_first;
};

// ----------------------------------------------------------------------------
// Class simultaneous_range
// ----------------------------------------------------------------------------

template <typename FwdIt> class simultaneous_range {
public:
    simultaneous_range(FwdIt left_first, FwdIt right_first, FwdIt left_last, FwdIt right_last)
            : m_left_first(left_first), m_right_first(right_first), m_left_last(left_last), m_right_last(right_last)
    {}

    simultaneous_iterator<FwdIt> begin() const {
        return simultaneous_iterator<FwdIt>(m_left_first, m_right_first);
    }

    simultaneous_iterator<FwdIt> end() const {
        return simultaneous_iterator<FwdIt>(m_left_last, m_right_last);
    }

private:
    FwdIt m_left_first;
    FwdIt m_right_first;
    FwdIt m_left_last;
    FwdIt m_right_last;
};

// ----------------------------------------------------------------------------
// Function make_simultaneous_range()
// ----------------------------------------------------------------------------

template <typename C> auto make_simultaneous_range(C& c, C& d) -> simultaneous_range<decltype(c.begin())> {
    return simultaneous_range<decltype(c.begin())>(c.begin(), d.begin(), c.end(), d.end());
}

// --------------------------------------------------------------------------
// Class Coordinate
// --------------------------------------------------------------------------

// Small helper type.

struct Coordinate
{
    int contigID { 0 };
    int pos { 0 };

    Coordinate() = default;
    Coordinate(int contigID, int pos) : contigID(contigID), pos(pos) {}
};

// --------------------------------------------------------------------------
// Class OverlapResolver
// --------------------------------------------------------------------------

// Once we have a good scaffold, we should look at potential overlaps between the contigs.  The contigs with overlaps
// can then be merged.
//
// We then linearize along the scaffold, fill the gaps with Ns.

template <typename TContigSeq>
class OverlapResolver
{
public:
    // Configuration of the OverlapResolver.
    struct Options
    {
        int k { 3 };
        int verbosity { 0 };
        int minOverlapMergeLength { 10 };
        std::string prefix { "prefix" };

        Options() = default;
    };

private:

    // Input / Output.

    // Resulting scaffold sequences.
    seqan::StringSet<TContigSeq> & scaffolds;
    // ScaffoldingResult to materialize.
    ScaffoldingResult scaffoldingResult;  // copied since updated
    // Input contigs that are to be scaffolded.
    seqan::StringSet<TContigSeq> const & contigs;
    // Configuration of the overlap resolver.
    Options options;

public:

    OverlapResolver(seqan::StringSet<TContigSeq> & scaffolds,
                    ScaffoldingResult const & scaffoldingResult,
                    seqan::StringSet<TContigSeq> const & contigs,
                    Options options) :
            scaffolds(scaffolds), scaffoldingResult(scaffoldingResult), contigs(contigs), options(options)
    {}

    // Return true if there were overlaps to be resolved.
    bool run();

private:

    // Compute overlaps.
    std::unique_ptr<assembler::OverlapStore> computeOverlaps() const;

    // Compute overlap candidates from inferred positions.
    std::vector<assembler::OverlapCandidate> computeOverlapCandidates() const;

    // Use overlaps to make the positions in scaffoldingResult more precise or separate presumably overlapping contigs
    // that do not actually overlap in their sequence by one position.
    void updateScaffoldingResult(assembler::OverlapStore const & ovlStore);

    // Create scaffold sequences from scaffoldingResult and contigs.
    void linearizeAndFillGaps();
};

template <typename TContigSeq>
bool OverlapResolver<TContigSeq>::run()
{
    // Generate contig overlap candidates and verify.
    auto ovlStore = computeOverlaps();

    if (options.verbosity >= 3)
    {
        std::cerr << "Before update\n";
        scaffoldingResult.print(std::cerr);
    }

    // Merge overlapping contigs guided by the verified overlaps.
    updateScaffoldingResult(*ovlStore);

    if (options.verbosity >= 3)
    {
        std::cerr << "After update\n";
        scaffoldingResult.print(std::cerr);
    }
    
    // Project the contigs on the line for each scaffold and fill gaps with Ns.
    linearizeAndFillGaps();

    return (length(scaffolds) != length(contigs));
}

template <typename TContigSeq>
std::unique_ptr<assembler::OverlapStore> OverlapResolver<TContigSeq>::computeOverlaps() const
{
    assembler::OverlapperOptions ovlOptions;
    ovlOptions.logging = (options.verbosity >= 3);
    ovlOptions.overlapMinLength = 10;
    ovlOptions.overlapErrorRate = 0.10;  // higher error rate such that we can join even different haplotypes
    assembler::Overlapper overlapper(ovlOptions);  // default config is fine?

    auto candidates = computeOverlapCandidates();

    std::unique_ptr<assembler::OverlapStore> ovlStore(new assembler::OverlapStore);
    for (auto const & cand : candidates)
        overlapper.computeOverlap(*ovlStore, contigs[cand.seq0], contigs[cand.seq1], cand);
    ovlStore->refresh();
    if (options.verbosity >= 3)
        ovlStore->print(std::cerr);
    return ovlStore;
}

template <typename TContigSeq>
std::vector<assembler::OverlapCandidate> OverlapResolver<TContigSeq>::computeOverlapCandidates() const
{
    std::vector<assembler::OverlapCandidate> result;

    int const MIN_BAND = 20;

    // Iterate over all scaffolds and adjacent pairs in the scaffolds.
    for (auto const & scaffold : scaffoldingResult.scaffolds)
        for (auto itL = scaffold.begin(); std::next(itL) != scaffold.end(); ++itL)
            for (auto itR = std::next(itL); itR != scaffold.end(); ++itR)
            {
                if (itL->pos + options.k * itL->posSD + (int)itL->length < itR->pos - options.k * itR->posSD)
                    break;  // no tentative overlap

                int overlap0 = itL->pos + itL->length - itR->pos;  // overlap in itL
                int diag = itR->pos - itL->pos;
                int band = options.k * std::max(2, std::min(itL->posSD, (int)itL->length) + std::min(itR->posSD, (int)itR->length));
                if (band < MIN_BAND)
                    band = MIN_BAND;
                if (overlap0 >= options.minOverlapMergeLength)
                {
                    int lDiag = diag - band;
                    int uDiag = diag + band;
                    result.push_back(assembler::OverlapCandidate(itL->id, itR->id, lDiag, uDiag));
                    if (options.verbosity >= 3)
                        std::cerr << "ADDING OVERLAP CANDIDATE " << result.back() << "\n";
                }
            }

    return result;
}

template <typename TContigSeq>
void OverlapResolver<TContigSeq>::updateScaffoldingResult(assembler::OverlapStore const & ovlStore)
{
    ScaffoldingResult tmp;

    for (auto const & scaffold : scaffoldingResult.scaffolds)
    {
        tmp.scaffolds.push_back(ScaffoldingResult::TScaffold());
        auto & out = tmp.scaffolds.back();
        for (auto posContig : scaffold)
        {
            if (out.empty())  // always add first
            {
                out.push_back(PositionedContig(0, 0, posContig.id, posContig.length));
                continue;
            }
            auto & last = out.back();
            if (!last.tentativeOverlap(posContig))  // no tentative overlap, keep current pos
            {
                out.push_back(posContig);
                continue;
            }
            // If we reach here, there is a tentative overlap by the positions.
            auto it = ovlStore.find(last.id, posContig.id);
            if (it == ovlStore.end())
                it = ovlStore.find(posContig.id, last.id);
            if (it == ovlStore.end())
            {
                // Does not match an overlap in the store; update if indicating overlap or touching
                if (posContig.pos <= last.pos + (int)last.length)
                    posContig.pos = last.pos + last.length + 3;
            }
            else
            {
                // An overlap was indicated and we can make the position exact.
                int ovlLen = it->len0 - it->begin1;  // overlap length
                posContig.pos = last.pos + last.length - ovlLen;
                posContig.posSD = 0;  // exact now
            }
            out.push_back(posContig);
        }
    }

    using std::swap;
    tmp.shiftScaffolds();
    swap(tmp, scaffoldingResult);
}

template <typename TContigSeq>
void OverlapResolver<TContigSeq>::linearizeAndFillGaps()
{
    // Below, we have to obtain a mapping from reference id in siteData to the id of the target scaffold and an
    // offset therein.
    std::vector<Coordinate> mapping;
    mapping.resize(length(contigs));

    // Build new refNames and refs for siteData as well as the mapping.
    seqan::StringSet<seqan::CharString> refNames;
    seqan::StringSet<TContigSeq> refs;
    int scaffoldID = 0;
    for (auto const & scaffold : scaffoldingResult.scaffolds)
    {
        // Fill the sequence that we know of and create mapping.
        TContigSeq seq;
        SEQAN_ASSERT_NOT(scaffold.empty());
        unsigned newSize = 0;
        for (auto const & posContig : scaffold)
            newSize = std::max((int)newSize, (int)(posContig.pos + posContig.length));
        resize(seq, newSize, 'N');
        for (auto posContig : scaffold)
        {
            SEQAN_ASSERT_EQ(length(contigs[posContig.id]), (unsigned)posContig.length);
            std::copy(begin(contigs[posContig.id], seqan::Standard()),
                      end(contigs[posContig.id], seqan::Standard()),
                      iter(seq, posContig.pos, seqan::Standard()));
            mapping[posContig.id] = Coordinate(scaffoldID, posContig.pos);
        }
        appendValue(refs, seq);
        if (options.verbosity >= 3)
            std::cerr << "appendValue(refs, {len=" << length(seq) << "}" << seq << ")\n";

        // Build the scaffold ref name.
        std::stringstream ss;
        ss << options.prefix << "_scaffold_" << scaffoldID++;
        appendValue(refNames, ss.str());
    }

    if (options.verbosity >= 3)
    {
        std::cerr << "# of scaffolds: " << scaffoldID << "\n"
                  << "\nMAPPING\n";
        unsigned idx = 0;
        for (auto coord : mapping)
            std::cerr << idx++ << "\t" << coord.contigID << "\t" << coord.pos << "\n";
    }

    // Write out result.
    using std::swap;
    swap(scaffolds, refs);
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function resolveOverlaps()
// ----------------------------------------------------------------------------

void resolveOverlaps(seqan::StringSet<seqan::Dna5String> & scaffolds,
                     ScaffoldingResult const & scaffoldingResult,
                     seqan::StringSet<seqan::Dna5String> const & contigs)
{
    typedef OverlapResolver<seqan::Dna5String> TOverlapResolver;
    TOverlapResolver::Options options;
    // options.verbosity = 1000;
    TOverlapResolver resolver(scaffolds, scaffoldingResult, contigs, options);
    resolver.run();
}

void resolveOverlaps(seqan::StringSet<seqan::String<seqan::Dna5Q>> & scaffolds,
                     ScaffoldingResult const & scaffoldingResult,
                     seqan::StringSet<seqan::String<seqan::Dna5Q>> const & contigs)
{
    typedef OverlapResolver<seqan::String<seqan::Dna5Q>> TOverlapResolver;

    TOverlapResolver::Options options;
    // options.verbosity = 1000;
    TOverlapResolver resolver(scaffolds, scaffoldingResult, contigs, options);
    resolver.run();
}

}  // namespace scaffolder
