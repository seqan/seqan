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

#include "overlapper.h"

#include <numeric>

namespace {  // anonymous namespace
}  // anonymous namespace

namespace assembler {

// ----------------------------------------------------------------------------
// Class OverlapStoreImpl
// ----------------------------------------------------------------------------

class OverlapStoreImpl
{
public:
    typedef OverlapStore::TIterator TIterator;

    OverlapStoreImpl() = default;

    // Add the given alignment to the store.
    void insert(Overlap const & ovl, TFragments const & alignment);

    // Refresh the store such that the query functions work.
    void refresh();

    // Filter out alignments if they do not fit the predicate.
    void eraseIf(std::function<bool(Overlap)> pred);

    // Retrieve iterator to first / last element.
    TIterator begin() const { return overlaps_.begin(); }
    TIterator end() const { return overlaps_.end(); }

    // Return overlaps.
    std::vector<Overlap> const & overlaps() const { return overlaps_; }

    // Find alignment between the ordered (!) sequence pair.
    TIterator find(unsigned seq0, unsigned seq1) const;

    // Retrieve alignment fragments for the given iterator, must be valid.
    TFragments const & alignment(TIterator it) const;

    // Print store to stream.
    void print(std::ostream & out) const;

    // Returns store size.
    size_t size() const { return overlaps_.size(); }

private:

    // Check whether reading is allowed.
    void checkRead() const
    {
        if (needRefresh)
            throw std::runtime_error("OverlapStore requires refresh()\n");
        SEQAN_CHECK(overlaps_.size() == length(alignments), "Must have same size.");
    }

    // Whether or not we need a refresh.
    bool needRefresh { false };
    
    // The overlaps and the alignments.
    std::vector<Overlap> overlaps_;
    seqan::StringSet<TFragments> alignments;
};

void OverlapStoreImpl::insert(Overlap const & ovl, TFragments const & ali)
{
    overlaps_.push_back(ovl.normalize());
    appendValue(alignments, ali);

    // Mark store as requiring refresh.
    needRefresh = true;
}

void OverlapStoreImpl::refresh()
{
    if (overlaps_.empty())
        return;

    // Generate set of indices in overlaps/alignments sorted by alignment key.
    std::vector<unsigned> idxs(overlaps_.size());
    std::iota(idxs.begin(), idxs.end(), 0);
    std::stable_sort(idxs.begin(), idxs.end(), [&](unsigned lhs, unsigned rhs) {
            return (overlaps_[lhs].key() < overlaps_[rhs].key());
        });

    // Build output.
    unsigned numNew = 0;
    std::vector<Overlap> tmpOverlaps(overlaps_.size());
    seqan::StringSet<TFragments> tmpAlignments;
    resize(tmpAlignments, overlaps_.size());

    // Always write out first.
    tmpOverlaps[0] = overlaps_[idxs[0]];
    tmpAlignments[0] = alignments[idxs[0]];
    numNew = 1;

    for (auto it = std::next(idxs.begin()); it != idxs.end(); ++it)
        if (tmpOverlaps[numNew - 1].key() != overlaps_[*it].key())
        {
            tmpOverlaps[numNew] = overlaps_[*it];
            tmpAlignments[numNew] = alignments[*it];
            numNew += 1;
        }

    // Write out, resize in case of duplicates.
    tmpOverlaps.resize(numNew);
    resize(tmpAlignments, numNew);
    using std::swap;
    swap(tmpOverlaps, overlaps_);
    swap(tmpAlignments, alignments);

    needRefresh = false;  // mark as fresh
}

void OverlapStoreImpl::eraseIf(std::function<bool(Overlap)> pred)
{
    auto itO = overlaps_.begin();
    auto itA = seqan::begin(alignments, seqan::Standard());

    auto it = overlaps_.begin();
    auto it2 = seqan::begin(alignments, seqan::Standard());
    while (it != overlaps_.end())
    {
        if (!pred(*it))
        {
            *itO = *it;
            *itA = *it2;
            ++itO;
            ++itA;
        }
        ++it;
        ++it2;
    }

    unsigned numNew = it - overlaps_.begin();
    overlaps_.resize(numNew);
    resize(alignments, numNew);
}

OverlapStoreImpl::TIterator OverlapStoreImpl::find(unsigned seq0, unsigned seq1) const
{
    checkRead();

    auto cmp = [](Overlap const & ovl, std::pair<unsigned, unsigned> const & val) {
        return (ovl.key() < val);
    };

    auto it = std::lower_bound(overlaps_.begin(), overlaps_.end(), std::make_pair(seq0, seq1), cmp);
    if (it != overlaps_.end() && it->key() == std::make_pair(seq0, seq1))
        return it;
    else
        return overlaps_.end();
}

TFragments const & OverlapStoreImpl::alignment(TIterator it) const
{
    return alignments[it - overlaps_.begin()];
}

void OverlapStoreImpl::print(std::ostream & out) const
{
    char const * labels[2] = { "FRESH", "NEED_REFRESH" };
    
    out << "OverlapStore\n"
        << "state\t" << labels[needRefresh] << "\n"
        << "  Overlaps\n";
    for (auto const & ovl : overlaps_)
        out << "  " << ovl << "\n";
}

// ----------------------------------------------------------------------------
// Class OverlapStore
// ----------------------------------------------------------------------------

OverlapStore::OverlapStore() : impl(new OverlapStoreImpl())
{}

OverlapStore::~OverlapStore()
{}

void OverlapStore::insert(Overlap const & ovl, TFragments const & alignments)
{
    impl->insert(ovl, alignments);
}

void OverlapStore::eraseIf(std::function<bool(Overlap)> pred)
{
    impl->eraseIf(pred);
}

void OverlapStore::refresh()
{
    impl->refresh();
}

std::vector<Overlap> const & OverlapStore::overlaps() const
{
    return impl->overlaps();
}

OverlapStore::TIterator OverlapStore::begin() const
{
    return impl->begin();
}

OverlapStore::TIterator OverlapStore::end() const
{
    return impl->end();
}

OverlapStore::TIterator OverlapStore::find(unsigned seq0, unsigned seq1) const
{
    return impl->find(seq0, seq1);
}

TFragments const & OverlapStore::alignment(TIterator it) const
{
    return impl->alignment(it);
}

void OverlapStore::print(std::ostream & out) const
{
    impl->print(out);
}

size_t OverlapStore::size() const
{
    return impl->size();
}

}  // namespace assembler
