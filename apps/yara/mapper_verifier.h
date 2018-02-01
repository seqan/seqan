// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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

#ifndef APP_YARA_MAPPER_VERIFIER_H_
#define APP_YARA_MAPPER_VERIFIER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class AnchorsVerifier
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename Traits>
struct AnchorsVerifier
{
    typedef typename Traits::TContigSeqs       TContigSeqs;
    typedef typename Traits::TContigsPos       TContigsPos;
    typedef typename Traits::TReadSeqs         TReadSeqs;
    typedef typename Traits::TReadSeq          TReadSeq;
    typedef typename Traits::TReadsContext     TReadsContext;
    typedef typename Traits::TMatchesViewSet   TMatchesSet;
    typedef typename Traits::TMatchesAppender  TMatches;
    typedef typename Traits::TMatch            TMatch;

//    typedef Myers<>                                     TAlgorithm;
    typedef AffineGaps                                  TAlgorithm;
    typedef Verifier<TContigSeqs, TReadSeq, TAlgorithm> TVerifier;

    // Thread-private data.
    TVerifier           verifier;
    TMatch              prototype;

    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TMatches &          mates;

    // Shared-memory read-only data.
    TContigSeqs const & contigSeqs;
    TReadSeqs /*const*/ & readSeqs;
    TMatchesSet const & anchorsSet;
    unsigned            libraryLength;
    unsigned            libraryDev;
    Options const &     options;

    AnchorsVerifier(TReadsContext & ctx,
                    TMatches & mates,
                    TContigSeqs const & contigSeqs,
                    TReadSeqs /*const*/ & readSeqs,
                    TMatchesSet const & anchorsSet,
                    unsigned libraryLength,
                    unsigned libraryDev,
                    Options const & options) :
        verifier(contigSeqs),
        prototype(),
        ctx(ctx),
        mates(mates),
        contigSeqs(contigSeqs),
        readSeqs(readSeqs),
        anchorsSet(anchorsSet),
        libraryLength(libraryLength),
        libraryDev(libraryDev),
        options(options)
    {
        _verifyAnchorsImpl(*this);
    }

    template <typename TMatchesSetIt>
    void operator() (TMatchesSetIt const & anchorsSetIt)
    {
        _findMatesImpl(*this, anchorsSetIt);
    }

    template <typename TMatchPos, typename TMatchErrors>
    void operator() (TMatchPos matchBegin, TMatchPos matchEnd, TMatchErrors matchErrors)
    {
        _addMatchImpl(*this, matchBegin, matchEnd, matchErrors);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _verifyAnchorsImpl()
// ----------------------------------------------------------------------------
// Verifies all anchors.

template <typename TSpec, typename Traits>
inline void _verifyAnchorsImpl(AnchorsVerifier<TSpec, Traits> & me)
{
    // Guess the number of mates.
    reserve(me.mates, length(me.mates) + length(me.anchorsSet), Exact());

    // Iterate over all anchors.
    iterate(me.anchorsSet, me, Standard(), typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _findMatesImpl()
// ----------------------------------------------------------------------------
// Verifies all anchors.

template <typename TSpec, typename Traits, typename TMatchesSetIt>
inline void _findMatesImpl(AnchorsVerifier<TSpec, Traits> & me, TMatchesSetIt const & anchorsSetIt)
{
    auto const & anchors = *anchorsSetIt;
    auto readId = position(anchorsSetIt, me.anchorsSet);
    auto mateId = getMateId(me.readSeqs, readId);

    if (empty(anchors) || isMapped(me.ctx, mateId))
        return;

    SEQAN_ASSERT(empty(me.anchorsSet[mateId]));

//    for (auto const & anchor : anchors)
//        _findMateImpl(me, anchor);

    if (length(anchors) == 1)
        _findMateImpl(me, front(anchors));
}

// ----------------------------------------------------------------------------
// Function _findMateImpl()
// ----------------------------------------------------------------------------
// Verifies one anchor.

template <typename TSpec, typename Traits, typename TMatch>
inline void _findMateImpl(AnchorsVerifier<TSpec, Traits> & me, TMatch const & anchor)
{
    typedef typename Traits::TContigsPos                TContigsPos;

    typedef typename Traits::TReadSeqs                  TReadSeqs;
    typedef typename Size<TReadSeqs>::Type              TReadId;
    typedef typename Value<TReadSeqs>::Type             TReadSeq;
    typedef typename Size<TReadSeq>::Type               TErrors;

    // Get mate seq.
    TReadId mateSeqId = getMateSeqId(me.readSeqs, getReadSeqId(anchor, me.readSeqs));
    TReadSeq const & mateSeq = me.readSeqs[mateSeqId];

    TContigsPos contigBegin;
    TContigsPos contigEnd;

    if (isRevReadSeq(me.readSeqs, mateSeqId))
        _getMateContigPos(me, contigBegin, contigEnd, anchor, RightMate());
    else
        _getMateContigPos(me, contigBegin, contigEnd, anchor, LeftMate());

    // Fill readId.
    setReadId(me.prototype, me.readSeqs, mateSeqId);

    // Get absolute number of errors.
    TErrors maxErrors = getReadErrors<TMatch>(me.options, length(mateSeq));
    TErrors maxIndels = getReadIndels<TMatch>(me.options, length(mateSeq));

    verify(me.verifier, mateSeq, contigBegin, contigEnd, maxErrors, maxIndels, me);
}

// ----------------------------------------------------------------------------
// Function _addMatchImpl()
// ----------------------------------------------------------------------------
// Adds one mate.

template <typename TSpec, typename Traits, typename TMatchPos, typename TMatchErrors>
inline void _addMatchImpl(AnchorsVerifier<TSpec, Traits> & me,
                          TMatchPos matchBegin,
                          TMatchPos matchEnd,
                          TMatchErrors matchErrors)
{
    setContigPosition(me.prototype, matchBegin, matchEnd);
    me.prototype.errors = matchErrors;
    appendValue(me.mates, me.prototype, Generous(), typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _getMateContigPos()
// ----------------------------------------------------------------------------
// Computes the insert window.

// --> ... mate
template <typename TSpec, typename Traits, typename TContigPos, typename TMatch>
inline void _getMateContigPos(AnchorsVerifier<TSpec, Traits> const & me,
                              TContigPos & contigBegin,
                              TContigPos & contigEnd,
                              TMatch const & anchor,
                              RightMate)
{
    typedef typename Traits::TContig                TContig;
    typedef typename Size<TContig>::Type            TContigSize;

    TContigSize contigLength = length(me.contigSeqs[getMember(anchor, ContigId())]);

    setValueI1(contigBegin, getMember(anchor, ContigId()));
    setValueI1(contigEnd, getMember(anchor, ContigId()));

    contigBegin.i2 = 0;
    if (getMember(anchor, ContigBegin()) + me.libraryLength > 3 * me.libraryDev)
        contigBegin.i2 = getMember(anchor, ContigBegin()) + me.libraryLength - 3 * me.libraryDev;
    contigBegin.i2 = _min(contigBegin.i2, contigLength);

    contigEnd.i2 = _min(getMember(anchor, ContigBegin()) + me.libraryLength + 3 * me.libraryDev, contigLength);

    SEQAN_ASSERT_LEQ(getValueI2(contigBegin), getValueI2(contigEnd));
    SEQAN_ASSERT_LEQ(getValueI2(contigEnd) - getValueI2(contigBegin), 6 * me.libraryDev);
}

// mate ... <--
template <typename TSpec, typename Traits, typename TContigPos, typename TMatch>
inline void _getMateContigPos(AnchorsVerifier<TSpec, Traits> const & me,
                              TContigPos & contigBegin,
                              TContigPos & contigEnd,
                              TMatch const & anchor,
                              LeftMate)
{
    setValueI1(contigBegin, getMember(anchor, ContigId()));
    setValueI1(contigEnd, getMember(anchor, ContigId()));

    contigBegin.i2 = 0;
    if (getMember(anchor, ContigEnd()) > me.libraryLength + 3 * me.libraryDev)
        contigBegin.i2 = getMember(anchor, ContigEnd()) - me.libraryLength - 3 * me.libraryDev;

    contigEnd.i2 = 0;
    if (getMember(anchor, ContigEnd()) + me.libraryDev > me.libraryLength)
        contigEnd.i2 = getMember(anchor, ContigEnd()) - me.libraryLength + 3 * me.libraryDev;

    SEQAN_ASSERT_LEQ(getValueI2(contigBegin), getValueI2(contigEnd));
    SEQAN_ASSERT_LEQ(getValueI2(contigEnd) - getValueI2(contigBegin), 6 * me.libraryDev);
}

#endif  // #ifndef APP_YARA_MAPPER_VERIFIER_H_
