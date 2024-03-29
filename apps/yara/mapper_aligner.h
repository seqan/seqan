// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2024, Enrico Siragusa, FU Berlin
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

//#define YARA_PRINT_ALIGN

using namespace seqan2;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class MatchesAligner
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches>
struct MatchesAligner
{
    typedef typename Traits::TContigSeqs       TContigSeqs;
    typedef typename Traits::TReadSeqs         TReadSeqs;
    typedef typename Traits::TCigarString      TCigarString;
    typedef typename Traits::TCigars           TCigars;
    typedef typename Traits::TCigarLimits      TCigarLimits;

    typedef String<GapAnchor<int> >            TGapAnchors;
    typedef AnchorGaps<TGapAnchors>            TAnchorGaps;

    // Thread-private data.
    TGapAnchors contigAnchors;
    TGapAnchors readAnchors;
    TCigarString cigar;
//    CharString  md;

    // Shared-memory read-write data.
    TCigars &               cigarSet;
    TCigarLimits &          cigarLimits;
    TMatches &              matches;

    // Shared-memory read-only data.
    TContigSeqs const &     contigSeqs;
    TReadSeqs const &       readSeqs;
    Options const &         options;

    MatchesAligner(TCigars & cigarSet,
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

    template <typename TMatchIt>
    void operator() (TMatchIt & matchIt)
    {
        _alignMatchImpl(*this, matchIt);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _alignMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches>
inline void _alignMatches(MatchesAligner<TSpec, Traits, TMatches> & me)
{
    typedef typename Traits::TCigars                    TCigars;
    typedef typename StringSetLimits<TCigars>::Type     TCigarLimits;
    typedef typename Suffix<TCigarLimits>::Type         TCigarSuffix;
    typedef typename Traits::TMatchSpec                 TMatchSpec;

    TCigarLimits & cigarLimits = stringSetLimits(me.cigarSet);

    // Fill limits with cigar length estimates.
    resize(cigarLimits, length(me.matches) + 1, Exact());
    front(cigarLimits) = 0;
    TCigarSuffix const & cigarSuffix = suffix(cigarLimits, 1);
    transform(cigarSuffix, me.matches, getCigarLength<TMatchSpec>, typename Traits::TThreading());

    // Bucket the cigar set.
    partialSum(cigarLimits, typename Traits::TThreading());
    assign(me.cigarSet.positions, prefix(cigarLimits, length(cigarLimits) - 1));
    resize(host(me.cigarSet), lengthSum(me.cigarSet), Exact());

    // Fill the cigars.
    resize(me.cigarLimits, length(me.matches) + 1, 0, Exact());
    iterate(me.matches, me, Rooted(), typename Traits::TThreading());

    // Shrink the cigar limits.
    assign(cigarLimits, me.cigarLimits);
    partialSum(cigarLimits, typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _align(AffineGaps)
// ----------------------------------------------------------------------------

template <typename TContigGaps, typename TReadGaps, typename TErrors>
inline int _align(TContigGaps & contigGaps, TReadGaps & readGaps, TErrors errors, AffineGaps)
{
    return globalAlignment(contigGaps, readGaps,
                           Score<int>(0, -1000, -999, -1001),           // Match, mismatch, extend, open.
                           AlignConfig<true, false, false, true>(),     // Top, left, right, bottom.
                           -(int)errors, (int)errors, Gotoh()) / -999;
}

// ----------------------------------------------------------------------------
// Function _align(LinearGaps)
// ----------------------------------------------------------------------------

template <typename TContigGaps, typename TReadGaps, typename TErrors>
inline int _align(TContigGaps & contigGaps, TReadGaps & readGaps, TErrors errors, LinearGaps)
{
    return -globalAlignment(contigGaps, readGaps, Score<short, EditDistance>(), -(int)errors, (int)errors);
}

// ----------------------------------------------------------------------------
// Function _alignMatchImpl()
// ----------------------------------------------------------------------------
// Aligns one match.

template <typename TSpec, typename Traits, typename TMatches, typename TMatchIt>
inline void _alignMatchImpl(MatchesAligner<TSpec, Traits, TMatches> & me, TMatchIt & matchIt)
{
    typedef typename Value<TMatchIt>::Type          TMatch;
    typedef typename Traits::TContigSeqs            TContigSeqs;
    typedef typename Value<TContigSeqs const>::Type TContigSeq;
    typedef typename StringSetPosition<TContigSeqs const>::Type TContigPos;
    typedef typename Infix<TContigSeq>::Type        TContigInfix;
    typedef typename Size<TMatches>::Type           TSize;

    typedef typename Traits::TReadSeqs              TReadSeqs;
    typedef typename Value<TReadSeqs const>::Type   TReadSeq;

    typedef MatchesAligner<TSpec, Traits, TMatches> TMatchesAligner;
    typedef typename TMatchesAligner::TAnchorGaps   TAnchorGaps;
    typedef Gaps<TContigInfix, TAnchorGaps>         TContigGaps;
    typedef Gaps<TReadSeq, TAnchorGaps>             TReadGaps;

    TMatch & match = value(matchIt);
    TSize matchId = matchIt - begin(me.matches, Standard());


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

    // Do not align if the match contains no gaps.
    // TODO(esiragusa): reuse DP matrix.
    if (!(errors == 0 || (errors == 1 && length(contigInfix) == length(readSeq))))
    {
        int dpErrors = _align(contigGaps, readGaps, errors, TSpec());

        SEQAN_ASSERT_LEQ(dpErrors, (int)errors);
        ignoreUnusedVariableWarning(dpErrors);

        clipSemiGlobal(contigGaps, readGaps);

        // Shrink the match after realigning and clipping.
        TContigPos contigBegin(getMember(match, ContigId()), getMember(match, ContigBegin()));
        contigBegin = posAdd(contigBegin, beginPosition(contigGaps));
        TContigPos contigEnd = contigBegin;
        contigEnd = posAdd(contigEnd, endPosition(contigGaps));
        setContigPosition(match, contigBegin, contigEnd);
    }

#ifdef YARA_PRINT_ALIGN
        write(std::cerr, match);

        Align<Dna5String, TAnchorGaps> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), contigInfix);
        assignSource(row(align, 1), readSeq);
        _align(row(align, 0), row(align, 1), errors, TSpec());
        clipSemiGlobal(row(align, 0), row(align, 1));
        std::cerr << row(align, 0) << std::endl;
        std::cerr << row(align, 1) << std::endl;
#endif // YARA_PRINT_ALIGN

    // Compute cigar.
    clear(me.cigar);
    getCigarString(me.cigar, contigGaps, readGaps, length(contigInfix));
    SEQAN_CHECK(_getQueryLength(me.cigar) == length(readSeq), "CIGAR error.");
//    SEQAN_ASSERT_EQ(_getQueryLength(me.cigar), length(readSeq));

    // Copy cigar to set.
    // TODO(esiragusa): use assign if possible.
//    me.cigarSet[getReadId(match)] = me.cigar;
    SEQAN_CHECK(length(me.cigar) <= length(me.cigarSet[matchId]), "CIGAR error.");
//    SEQAN_ASSERT_LEQ(length(me.cigar), length(me.cigarSet[getMember(match, ReadId())]));
    std::copy(begin(me.cigar, Standard()), end(me.cigar, Standard()), begin(me.cigarSet[matchId], Standard()));
    assignValue(me.cigarLimits, matchId + 1, length(me.cigar));

//    clear(me.md);
//    getMDString(me.md, contigGaps, readGaps);
}

#endif  // #ifndef APP_YARA_MAPPER_ALIGNER_H_
