// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef APP_YARA_MAPPER_ALIGNER_H_
#define APP_YARA_MAPPER_ALIGNER_H_

using namespace seqan;

template <typename TGaps>
inline typename Size<TGaps>::Type
countLeadingGaps(TGaps const & gaps)
{
    return toViewPosition(gaps, 0);
}

template <typename TGaps>
inline typename Size<TGaps>::Type
countTrailingGaps(TGaps const & gaps)
{
    return length(gaps) - toViewPosition(gaps, length(source(gaps)) - 1) - 1;
}

template <typename TGlobalGaps, typename TLocalGaps>
inline void clipSemiGlobal(TGlobalGaps & global, TLocalGaps & local)
{
    typedef typename Size<TLocalGaps>::Type  TGapsSize;

    TGapsSize leadingGaps = countLeadingGaps(local);
    TGapsSize trailingGaps = countTrailingGaps(local);
    TGapsSize globalLenght = length(global);
    TGapsSize localLength = length(local);

    setClippedBeginPosition(global, leadingGaps);
    setClippedBeginPosition(local, leadingGaps);
    setClippedEndPosition(global, globalLenght - trailingGaps);
    setClippedEndPosition(local, localLength - trailingGaps);
}

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class MatchesAligner
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
struct MatchesAligner
{
    typedef typename Traits::TContigSeqs       TContigSeqs;
    typedef typename Traits::TReadSeqs         TReadSeqs;
    typedef typename Traits::TMatches          TMatches;
    typedef typename Traits::TCigar            TCigar;
    typedef typename Traits::TCigarSet         TCigarSet;
    typedef typename Traits::TCigarLimits      TCigarLimits;

    typedef String<GapAnchor<int> >            TGapAnchors;
    typedef AnchorGaps<TGapAnchors>            TAnchorGaps;

    // Thread-private data.
    TGapAnchors contigAnchors;
    TGapAnchors readAnchors;
    TCigar      cigar;
//    CharString  md;

    // Shared-memory read-write data.
    TCigarSet &             cigarSet;
    TCigarLimits &          cigarLimits;
    TMatches &              matches;

    // Shared-memory read-only data.
    TContigSeqs const &     contigSeqs;
    TReadSeqs const &       readSeqs;
    Options const &         options;

    MatchesAligner(TCigarSet & cigarSet,
                   TCigarLimits & cigarLimits,
                   TMatches & matches,
                   TContigSeqs const & contigSeqs,
                   TReadSeqs const & readSeqs,
                   Options const & options) :
        cigarSet(cigarSet),
        cigarLimits(cigarLimits),
        matches(matches),
        contigSeqs(contigSeqs),
        readSeqs(readSeqs),
        options(options)
    {
        _alignMatches(*this);
    }

    template <typename TMatch>
    void operator() (TMatch & match)
    {
        _alignMatchImpl(*this, match);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _alignMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _alignMatches(MatchesAligner<TSpec, Traits> & me)
{
    typedef typename Traits::TCigarSet                  TCigarSet;
    typedef typename StringSetLimits<TCigarSet>::Type   TCigarLimits;
    typedef typename Suffix<TCigarLimits>::Type         TCigarSuffix;

    TCigarLimits & cigarLimits = stringSetLimits(me.cigarSet);

    // Fill limits with cigar length estimates.
    resize(cigarLimits, length(me.matches) + 1, Exact());
    front(cigarLimits) = 0;
    TCigarSuffix const & cigarSuffix = suffix(cigarLimits, 1);
    transform(cigarSuffix, me.matches, [&](Match<void> const & match)
    {
        return 2 * getCigarLength<void>(match);
    },
    typename Traits::TThreading());

    // Bucket the cigar set.
    partialSum(cigarLimits, typename Traits::TThreading());
    assign(me.cigarSet.positions, prefix(cigarLimits, length(cigarLimits) - 1));
    resize(host(me.cigarSet), lengthSum(me.cigarSet), Exact());

    // Fill the cigars.
    resize(me.cigarLimits, length(me.matches) + 1, 0, Exact());
    iterate(me.matches, me, Standard(), typename Traits::TThreading());

    // Shrink the cigar limits.
    assign(cigarLimits, me.cigarLimits);
    partialSum(cigarLimits, typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _alignMatchImpl()
// ----------------------------------------------------------------------------
// Aligns one match.

template <typename TSpec, typename Traits, typename TMatchIt>
inline void _alignMatchImpl(MatchesAligner<TSpec, Traits> & me, TMatchIt & matchIt)
{
    typedef typename Value<TMatchIt>::Type          TMatch;
    typedef typename Traits::TContigSeqs            TContigSeqs;
    typedef typename Value<TContigSeqs const>::Type TContigSeq;
    typedef typename Infix<TContigSeq>::Type        TContigInfix;

    typedef typename Traits::TReadSeqs              TReadSeqs;
    typedef typename Value<TReadSeqs const>::Type   TReadSeq;

    typedef MatchesAligner<TSpec, Traits>           TMatchesAligner;
    typedef typename TMatchesAligner::TAnchorGaps   TAnchorGaps;
    typedef Gaps<TContigInfix, TAnchorGaps>         TContigGaps;
    typedef Gaps<TReadSeq, TAnchorGaps>             TReadGaps;

    TMatch & match = value(matchIt);

    if (isInvalid(match)) return;

    unsigned errors = getMember(match, Errors());
    TReadSeq const & readSeq = me.readSeqs[getReadSeqId(match, me.readSeqs)];
    TContigInfix const & contigInfix = infix(me.contigSeqs[getMember(match, ContigId())],
                                             getMember(match, ContigBegin()),
                                             getMember(match, ContigEnd()));

    clear(me.contigAnchors);
    clear(me.readAnchors);
    TContigGaps contigGaps(contigInfix, me.contigAnchors);
    TReadGaps readGaps(readSeq, me.readAnchors);

    int dpErrors = errors;

    // Do not align if the match contains no gaps.
    // TODO(esiragusa): reuse DP matrix.
    if (!(errors == 0 || (errors == 1 && length(contigInfix) == length(readSeq))))
    {
//        dpErrors = -globalAlignment(contigGaps, readGaps, Score<short, EditDistance>(), -(int)errors, (int)errors);

        dpErrors = globalAlignment(contigGaps, readGaps,
                                   Score<int>(0, -999, -1001, -1000),           // Match, mismatch, extend, open.
                                   AlignConfig<true, false, false, true>(),     // Top, left, right, bottom.
                                   -(int)errors, (int)errors, Gotoh()) / -999;

        if (me.options.verbose > 1 && countLeadingGaps(readGaps))
        {
            typedef Align<Dna5String, ArrayGaps>  TAlign;

            TAlign align;
            resize(rows(align), 2);
            assignSource(row(align, 0), contigInfix);
            assignSource(row(align, 1), readSeq);

            globalAlignment(align,
                            Score<int>(0, -999, -1001, -1000),           // Match, mismatch, extend, open.
                            AlignConfig<true, false, false, true>(),     // Top, left, right, bottom.
                            -(int)errors, (int)errors, Gotoh());

            std::cout << align << std::endl << std::endl;

            std::cout << "Leading gaps: " <<  countLeadingGaps(row(align, 0)) << " : " << countLeadingGaps(row(align, 1)) << std::endl;
            std::cout << "Trailing gaps: " <<  countTrailingGaps(row(align, 0)) << " : " << countTrailingGaps(row(align, 1)) << std::endl;

            clipSemiGlobal(row(align, 0), row(align, 1));
            std::cout << align << std::endl << std::endl;
        }

        match.contigBegin += clippedBeginPosition(contigGaps);
        match.contigEnd += clippedBeginPosition(contigGaps);
//        setContigPosition(match, getMember(match, ContigBegin()) + clippedBeginPosition(readGaps), getMember(match, ContigEnd()));
        clipSemiGlobal(contigGaps, readGaps);
    }

    SEQAN_ASSERT_EQ(dpErrors, (int)errors);
    ignoreUnusedVariableWarning(dpErrors);

    // Compute cigar.
    clear(me.cigar);
    getCigarString(me.cigar, contigGaps, readGaps);

    // Copy cigar to set.
    // TODO(esiragusa): use assign if possible.
//    me.cigarSet[getReadId(match)] = me.cigar;
    SEQAN_ASSERT_LEQ(length(me.cigar), length(me.cigarSet[getMember(match, ReadId())]));
    std::copy(begin(me.cigar, Standard()), end(me.cigar, Standard()), begin(me.cigarSet[getMember(match, ReadId())], Standard()));
    assignValue(me.cigarLimits, getMember(match, ReadId()) + 1, length(me.cigar));

//    clear(me.md);
//    getMDString(me.md, contigGaps, readGaps);
}

#endif  // #ifndef APP_YARA_MAPPER_ALIGNER_H_
