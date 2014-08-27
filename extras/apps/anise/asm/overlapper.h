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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASSEMBLER_OVERLAPPER_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ASSEMBLER_OVERLAPPER_H_

#include <functional>
#include <iosfwd>
#include <memory>

// TODO(holtgrew): Forwards for seqan::stuff
#include <seqan/sequence.h>
#include <seqan/graph_types/graph_impl_fragment.h>

#include "asm/options.h"

namespace assembler {

// ============================================================================
// Forwards
// ============================================================================

class OverlapperImpl;
class OverlapStoreImpl;
class OverlapCandidateGeneratorImpl;

// ============================================================================
// Tags, Classes, Enums, Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Class OverlapCandidate
// ----------------------------------------------------------------------------

class OverlapCandidate
{
public:
    unsigned seq0 { 0 };
    unsigned seq1 { 0 };
    int lDiag { 0 };
    int uDiag { 0 };

    // Constructors.
    OverlapCandidate() = default;
    OverlapCandidate(unsigned seq0, unsigned seq1, int lDiag, int uDiag) :
            seq0(seq0), seq1(seq1), lDiag(lDiag), uDiag(uDiag)
    {}

    // Compares two OverlapCanidate objects lexicographically.
    bool operator<(OverlapCandidate const & other) const
    { return std::make_pair(seq0, seq1) < std::make_pair(other.seq0, other.seq1); }
};

inline std::ostream & operator<<(std::ostream & out, OverlapCandidate const & cand)
{
    return out << "OverlapCandidate(seq0=" << cand.seq0 << ", seq1=" << cand.seq1
               << ", lDiag=" << cand.lDiag << ", uDiag=" << cand.uDiag << ")";
}

// ----------------------------------------------------------------------------
// Class OverlapCandidateGeneratorOptions
// ----------------------------------------------------------------------------

struct OverlapCandidateGeneratorOptions
{
    // k-mers with more than this number of occurences are ignored
    unsigned kMerMaxOcc { 200 };
    // k-mer length for overlap candidate generation.
    unsigned k { 20 };  // 5% error rate
    // Wheter or not to enable logging.
    bool logging { false };
    
    OverlapCandidateGeneratorOptions() = default;
};

// ----------------------------------------------------------------------------
// Class OverlapCandidateGenerator
// ----------------------------------------------------------------------------

class OverlapCandidateGenerator
{
public:

    OverlapCandidateGenerator(OverlapCandidateGeneratorOptions options = OverlapCandidateGeneratorOptions());
    ~OverlapCandidateGenerator();  // for pimpl

    // Perform candidate generation and pass to functor
    void run(std::function<void(OverlapCandidate)> func,
             seqan::StringSet<seqan::Dna5String> const & seqs) const;

private:
    std::unique_ptr<OverlapCandidateGeneratorImpl> impl;
};

// ----------------------------------------------------------------------------
// Typedef TFragments
// ----------------------------------------------------------------------------

typedef seqan::String<seqan::Fragment<>> TFragments;

// ----------------------------------------------------------------------------
// Class Overlap
// ----------------------------------------------------------------------------

// Stores an overlap, usually such that seq0 is left of seq1 or if both start at the same position then seq0 is the id
// of the longer sequence (container).  In case of stacking, seq0 < seq1.

class Overlap
{
public:
    // Value for invalid entry below.
    static const unsigned INVALID = (unsigned)-1;
    
    // Identifiers of the sequence.
    unsigned seq0 { INVALID };
    unsigned seq1 { INVALID };
    // Lengths of the sequences.
    unsigned len0 { INVALID };
    unsigned len1 { INVALID };
    // Begin positions of the sequences.
    unsigned begin0 { INVALID };
    unsigned begin1 { INVALID };
    // Edit distance errors of the alignment.
    unsigned errors { INVALID };

    // Default constructor and construct with all values.
    Overlap() = default;
    
    Overlap(unsigned seq0, unsigned seq1, unsigned len0, unsigned len1, unsigned begin0, unsigned begin1,
            unsigned errors = INVALID) :
            seq0(seq0), seq1(seq1), len0(len0), len1(len1), begin0(begin0), begin1(begin1), errors(errors)
    {}

    // Returns a "flipped" alignment, i.e. seq0 and seq1 change roles.
    Overlap flip() const { return Overlap(seq1, seq0, len1, len0, begin1, begin0, errors); }

    // Returns a normalized overlap, i.e. if begin0 != begin1 then seq0 < seq1.  If begin0 == begin1 then ties are
    // broken by length (such that the container comes first), in case of stacks, ties are broken by id.
    Overlap normalize() const
    {
        if (begin0 != begin1)
        {
            return (begin0 < begin1) ? *this : flip();
        }
        else  // begin0 == begin1
        {
            if (len0 != len1)
                return (len0 > len1) ? *this : flip();
            else  // begin0 == begin1 && len0 == len1
                return (seq0 < seq1) ? *this : flip();
        }
    }

    // Returns the length of the overlap.
    int length() const
    {
        auto o = normalize();
        return (o.len0 - o.begin1);
    }

    // Returns "key", i.e. (seq0, seq1).
    std::pair<unsigned, unsigned> key() const { return std::make_pair(seq0, seq1); }

    // Compares two Overlap objects lexicographically.
    bool operator<(Overlap const & other) const { return (key() < other.key()); }
};

inline std::ostream & operator<<(std::ostream & out, Overlap const & ovl)
{
    return out << "Overlap(seq0=" << ovl.seq0 << ", seq1=" << ovl.seq1 << ", len0=" << ovl.len0
               << ", len1=" << ovl.len1 << ", begin0=" << ovl.begin0 << ", begin1=" << ovl.begin1
               << ", errors=" << ovl.errors << ")";
}

// ----------------------------------------------------------------------------
// Class OverlapStore
// ----------------------------------------------------------------------------

// Store for the overlaps.
//
// After calling insert(), you have to call refresh() before you can call any query function to resort the alignments by
// (seq0, seq1) and remove possible duplicates.  After calling eraseIf(), the store does not require refreshing.

class OverlapStore
{
public:
    typedef std::vector<Overlap>::const_iterator TIterator;

    // Constructor and destructor.
    OverlapStore();
    ~OverlapStore();  // for pimpl

    // Add the given alignment to the store, note that the alignment will be normalized.
    void insert(Overlap const & ovl, TFragments const & alignment);

    // Refresh the store such that the query functions work.
    void refresh();

    // Filter out alignments if they do not fit the predicate.
    void eraseIf(std::function<bool(Overlap)> pred);

    // returns number of overlaps.
    size_t size() const;

    // Retrieve the raw overlaps.
    std::vector<Overlap> const & overlaps() const;
    
    // Retrieve iterator to first / last element.
    TIterator begin() const;
    TIterator end() const;

    // Find alignment between the ordered (!) sequence pair.
    TIterator find(unsigned seq0, unsigned seq1) const;

    // Retrieve alignment fragments for the given iterator, must be valid.
    TFragments const & alignment(TIterator it) const;

    // Print store to stream.
    void print(std::ostream & out) const;

private:
    std::unique_ptr<OverlapStoreImpl> impl;
};

// ----------------------------------------------------------------------------
// Class OverlapperOptions
// ----------------------------------------------------------------------------

// Overlapper configuration.
struct OverlapperOptions
{
    double overlapErrorRate { 0.05 };
    int overlapMinLength { 40 };
    bool logging { false };

    OverlapperOptions() = default;
};

// ----------------------------------------------------------------------------
// Class Overlapper
// ----------------------------------------------------------------------------

class Overlapper
{
public:
    Overlapper(OverlapperOptions options = OverlapperOptions());
    ~Overlapper();  // for pimpl

    // Compute overlap from candidate and store it in overlap and alignment on success.  Return false if the overlap did
    // not exist.
    bool computeOverlap(Overlap & overlap,
                        TFragments & alignment,
                        seqan::Dna5String const & seqH,
                        seqan::Dna5String const & seqV,
                        OverlapCandidate const & candidate) const;    

    // Computes the overlap and stores it together with the alignment in the output parameters.  Returns false if the
    // overlap did not exist.
    bool computeOverlap(Overlap & overlap,
                        TFragments & alignment,
                        seqan::Dna5String const & seqH,
                        seqan::Dna5String const & seqV,
                        unsigned seq0,
                        unsigned seq1,
                        int diag,
                        int band) const
    {
        return computeOverlap(overlap, alignment, seqH, seqV, { seq0, seq1, diag - band, diag + band });
    }

    // Computes the overlap and directly stores it in the OverlapStore.  Returns false if the overlap did not exist.
    bool computeOverlap(OverlapStore & store,
                        seqan::Dna5String const & seqH,
                        seqan::Dna5String const & seqV,
                        OverlapCandidate const & candidate) const;

    // Computes the overlap and directly stores it in the OverlapStore.  Returns false if the overlap did not exist.
    bool computeOverlap(OverlapStore & store,
                        seqan::Dna5String const & seqH,
                        seqan::Dna5String const & seqV,
                        unsigned seq0,
                        unsigned seq1,
                        int diag,
                        int band) const
    {
        return computeOverlap(store, seqH, seqV, { seq0, seq1, diag - band, diag + band });
    }

    // Returns an overlap of the two sequences from candidate info.
    Overlap computeOverlap(seqan::Dna5String const & seqH,
                           seqan::Dna5String const & seqV,
                           OverlapCandidate const & candidate) const;


    // Returns an overlap of the two sequences at the given diagonal using the given bandwidth for overlapping.  Maximal
    // error rate and minimal overlap length are taken from the configuration.  Returns a default constructed Overlap
    // (fields set to INVALID or -1) if no such overlap could be found.
    Overlap computeOverlap(seqan::Dna5String const & seqH,
                           seqan::Dna5String const & seqV,
                           unsigned seq0,
                           unsigned seq1,
                           int diag,
                           int band)
    {
        return computeOverlap(seqH, seqV, { seq0, seq1, diag - band, diag + band });
    }

    // Test whether an overlap exists between seqH and seqV at the given diagonal with the given bandwidth.  Maximal
    // error rate and minimal overlap length are taken from configuration.
    bool overlapExists(seqan::Dna5String const & seqH,
                       seqan::Dna5String const & seqV,
                       int diagonal,
                       int bandwidth)
    {
        Overlap ovl;
        TFragments aln;
        return computeOverlap(ovl, aln, seqH, seqV, 0, 1, diagonal, bandwidth);
    }

private:
    std::unique_ptr<OverlapperImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace assembler

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASSEMBLER_OVERLAPPER_H_
