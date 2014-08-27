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

#include "overlapper.h"

#include <seqan/index.h>
#include <seqan/index/find_pigeonhole.h>

namespace {  // anonymous namespace
}  // anonymous namespace

namespace assembler {

// --------------------------------------------------------------------------
// Class OverlapCandidateGeneratorImpl
// --------------------------------------------------------------------------

class OverlapCandidateGeneratorImpl
{
public:

    OverlapCandidateGeneratorImpl(OverlapCandidateGeneratorOptions options = OverlapCandidateGeneratorOptions()) :
            options(options)
    {}

    void run(std::function<void(OverlapCandidate)> func,
             seqan::StringSet<seqan::Dna5String> const & seqs) const;

private:
    OverlapCandidateGeneratorOptions options;
};

void OverlapCandidateGeneratorImpl::run(std::function<void(OverlapCandidate)> func,
                                        seqan::StringSet<seqan::Dna5String> const & seqs) const
{
    typedef seqan::StringSet<seqan::Dna5String>                          TStringSet;
    typedef seqan::Shape<seqan::Dna5, seqan::OneGappedShape>             TShape;
    typedef seqan::IndexQGram<TShape, seqan::OpenAddressing>             TIndexSpec;
    typedef seqan::Index<TStringSet const, TIndexSpec>                   TIndex;
    typedef seqan::Pattern<TIndex, seqan::Pigeonhole<>>                  TFilterPattern;
    typedef seqan::Finder<seqan::Dna5String const, seqan::Pigeonhole<>>  TFilterFinder;

    // TODO(holtgrew): Replace by just using q-grams?

    // Convert q-gram size to max error rate.
    double maxErrorRate = 1.0 / options.k;
    
    // Build q-gram index.
    TIndex index(seqs);
    TFilterPattern filterPattern(index);
    _patternInit(filterPattern, maxErrorRate);

    // Perform the pigeonhole-based search.
    for (auto it = begin(seqs, seqan::Rooted()); !atEnd(it); ++it)
    {
        unsigned seq0 = position(it);
        TFilterFinder filterFinder(*it);
        while (find(filterFinder, filterPattern, maxErrorRate))
        {
            if (length(countOccurrencesMultiple(index, filterPattern.shape)) > options.kMerMaxOcc)
                continue;  // Ignore, too many matching sequences.

            unsigned seq1 = filterFinder.curHit->ndlSeqNo;
            // TODO(holtgrew): Care about dupes, i.e. use >=.
            if (seq1 == seq0)
                continue;  // Skip hits with self.
            __int64 lDiag = filterFinder.curHit->hstkPos;
            __int64 uDiag = filterFinder.curHit->hstkPos + filterFinder.curHit->bucketWidth - length(seqs[seq1]);
            SEQAN_ASSERT_GEQ(uDiag, lDiag);

            OverlapCandidate cand(positionToId(seqs, seq0), positionToId(seqs, seq1), lDiag, uDiag);
            func(cand);
            if (options.logging)
                std::cerr << "CANDIDATE\t" << cand << "\t" << seqs[seq0] << "\t" << seqs[seq1] << "\n"
                          << "    COUNTS\t" << length(countOccurrencesMultiple(index, filterPattern.shape)) << "\n";
        }
    }

}

// --------------------------------------------------------------------------
// Class OverlapCandidateGenerator
// --------------------------------------------------------------------------

OverlapCandidateGenerator::~OverlapCandidateGenerator()
{}

OverlapCandidateGenerator::OverlapCandidateGenerator(OverlapCandidateGeneratorOptions options) :
        impl(new OverlapCandidateGeneratorImpl(options))
{}

void OverlapCandidateGenerator::run(std::function<void(OverlapCandidate)> func,
                                    seqan::StringSet<seqan::Dna5String> const & seqs) const
{
    return impl->run(func, seqs);
}

}  // namespace assembler
