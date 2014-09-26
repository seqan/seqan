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
    typedef typename Traits::TMatchesSet       TMatchesSet;
    typedef typename Traits::TMatches          TMatches;
    typedef typename Traits::TMatch            TMatch;

    typedef Myers<>                                     TAlgorithm;
//    typedef Filter<MultipleShiftAnd>                    TAlgorithm;
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
    Options const &     options;

    AnchorsVerifier(TReadsContext & ctx,
                    TMatches & mates,
                    TContigSeqs const & contigSeqs,
                    TReadSeqs /*const*/ & readSeqs,
                    TMatchesSet const & anchorsSet,
                    Options const & options) :
        verifier(contigSeqs),
        prototype(),
        ctx(ctx),
        mates(mates),
        contigSeqs(contigSeqs),
        readSeqs(readSeqs),
        anchorsSet(anchorsSet),
        options(options)
    {
        _verifyAnchorsImpl(*this);
    }

    void operator() (TMatch const & anchor)
    {
        _findMateImpl(*this, anchor);
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
    // TODO(esiragusa): guess the number of pairs.
    reserve(me.mates, lengthSum(me.anchorsSet), Exact());

    // Iterate over all anchors.
    forEach(concat(me.anchorsSet), me, typename Traits::TThreading());
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

//    typedef Myers<>                                     TAlgorithm;
//    typedef Filter<MultipleShiftAnd>                    TAlgorithm;
//    typedef Verifier<TContigSeqs, TReadSeq, TAlgorithm> TVerifier;

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
    TErrors maxErrors = getReadErrors(me.options, length(mateSeq));

    verify(me.verifier, mateSeq, contigBegin, contigEnd, maxErrors, me);
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
    appendValue(me.mates, me.prototype, Insist(), typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _getMateContigPos()
// ----------------------------------------------------------------------------
// Computes the insert window.

template <typename TSpec, typename Traits, typename TContigPos, typename TMatch>
inline void _getMateContigPos(AnchorsVerifier<TSpec, Traits> const & me,
                              TContigPos & contigBegin,
                              TContigPos & contigEnd,
                              TMatch const & anchor,
                              RightMate)
{
    typedef typename Traits::TContig                TContig;
    typedef typename Size<TContig>::Type            TContigSize;

    TContigSize contigLength = length(me.contigSeqs[getContigId(anchor)]);

    setValueI1(contigBegin, getContigId(anchor));
    setValueI1(contigEnd, getContigId(anchor));

    contigBegin.i2 = 0;
    if (getMember(anchor, ContigBegin()) + me.options.libraryLength > me.options.libraryError)
        contigBegin.i2 = getMember(anchor, ContigBegin()) + me.options.libraryLength - me.options.libraryError;
    contigBegin.i2 = _min(contigBegin.i2, contigLength);

    contigEnd.i2 = _min(getMember(anchor, ContigBegin()) + me.options.libraryLength + me.options.libraryError, contigLength);

    SEQAN_ASSERT_LEQ(getValueI2(contigBegin), getValueI2(contigEnd));
    SEQAN_ASSERT_LEQ(getValueI2(contigEnd) - getValueI2(contigBegin), 2 * me.options.libraryError);
}

template <typename TSpec, typename Traits, typename TContigPos, typename TMatch>
inline void _getMateContigPos(AnchorsVerifier<TSpec, Traits> const & me,
                              TContigPos & contigBegin,
                              TContigPos & contigEnd,
                              TMatch const & anchor,
                              LeftMate)
{
    setValueI1(contigBegin, getContigId(anchor));
    setValueI1(contigEnd, getContigId(anchor));

    contigBegin.i2 = 0;
    if (getMember(anchor, ContigEnd()) > me.options.libraryLength + me.options.libraryError)
        contigBegin.i2 = getMember(anchor, ContigEnd()) - me.options.libraryLength - me.options.libraryError;

    contigEnd.i2 = 0;
    if (getMember(anchor, ContigEnd()) + me.options.libraryError > me.options.libraryLength)
        contigEnd.i2 = getMember(anchor, ContigEnd()) - me.options.libraryLength + me.options.libraryError;

    SEQAN_ASSERT_LEQ(getValueI2(contigBegin), getValueI2(contigEnd));
    SEQAN_ASSERT_LEQ(getValueI2(contigEnd) - getValueI2(contigBegin), 2 * me.options.libraryError);
}

#endif  // #ifndef APP_YARA_MAPPER_VERIFIER_H_
