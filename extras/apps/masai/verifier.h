// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_ESIRAGUSA_APPS_MASAI_VERIFIER_H_
#define SANDBOX_ESIRAGUSA_APPS_MASAI_VERIFIER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find.h>

#include "store.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TMatchesDelegate, typename TDistance = HammingDistance, typename TSpec = void>
struct Verifier
{
    TFragmentStore & store;
    TMatchesDelegate & matchesDelegate;
    String<TContigSeqSize>  contigSizes;

    TReadSeqSize            maxErrorsPerRead;
    TContigSeqSize          libraryLength;
    TContigSeqSize          libraryError;

    unsigned long           pairsCount;
    bool                    disabled;

    Verifier(TFragmentStore & store,
             TMatchesDelegate & matchesDelegate,
             TContigSeqSize libraryLength,
             TContigSeqSize libraryError,
             bool disabled = false) :
        store(store),
        matchesDelegate(matchesDelegate),
        maxErrorsPerRead(0),
        libraryLength(libraryLength),
        libraryError(libraryError),
        pairsCount(0),
        disabled(disabled)
    {
        _init(*this);
    }

};

template <typename TMatchesDelegate, typename TSpec>
struct Verifier<TMatchesDelegate, EditDistance, TSpec>:
    public Verifier<TMatchesDelegate>
{
    typedef Verifier<TMatchesDelegate>  TBase;
    typedef Pattern<TReadSeq, Myers<> > TPattern;

    TPattern    pattern;

    Verifier(TFragmentStore & store,
             TMatchesDelegate & matchesDelegate,
             TContigSeqSize libraryLength,
             TContigSeqSize libraryError,
             bool disabled = false) :
        TBase(store, matchesDelegate, libraryLength, libraryError, disabled)
    {
        _patternMatchNOfPattern(pattern, false);
        _patternMatchNOfFinder(pattern, false);
    }

};

//____________________________________________________________________________

template <typename THitsDelegate, typename TSpec>
struct Verifier<THitsDelegate, EditDistance, Filter<TSpec> >:
    public Verifier<THitsDelegate>
{
    typedef Verifier<THitsDelegate>  TBase;

    typedef Pattern<String<TReadSeq>, MultipleShiftAnd>   TPattern;

    unsigned seedLength;

    String<TReadSeq> seedsFwd;
    String<TReadSeq> seedsRev;

    TPattern patternFwd;
    TPattern patternRev;

    Verifier(TFragmentStore & store,
             THitsDelegate & hitsDelegate,
             TContigSeqSize libraryLength,
             TContigSeqSize libraryError,
             bool disabled = false) :
        TBase(store, hitsDelegate, libraryLength, libraryError, disabled),
        seedLength(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TMatchesDelegate, typename TDistance, typename TSpec>
inline void _init(Verifier<TMatchesDelegate, TDistance, TSpec> & verifier)
{
    reserve(verifier.contigSizes, length(verifier.store.contigStore), Exact());
    for (TContigStoreSize contigId = 0; contigId < length(verifier.store.contigStore); ++contigId)
        appendValue(verifier.contigSizes, length(verifier.store.contigStore[contigId].seq));
}

template <typename TMatchesDelegate, typename TDistance, typename TSpec, typename TReadId>
inline void preprocessMate(Verifier<TMatchesDelegate, TDistance, TSpec> &, TReadId)
{}

//____________________________________________________________________________

template <typename TMatchesDelegate, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Verifier<TMatchesDelegate, HammingDistance, TSpec> & verifier,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    typedef Iterator<TContigInfix, Standard>::Type          TContigInfixIterator;
    typedef Iterator<TReadSeq, Standard>::Type              TReadSeqIterator;

    if (verifier.disabled)
        return;

    // TODO Verifier must specify if readIds refer to left or right file
    // NOTE readId corresponds to left file
    readId *= 2;
    TReadId mateId = readId + 1;
    if (!reverseComplemented)
        mateId *= 2;

    TReadSeq & mate = verifier.store.readSeqStore[mateId];
    TReadSeqIterator mateBegin = begin(mate, Standard());
    TReadSeqIterator mateEnd = end(mate, Standard());
    TReadSeqSize mateLength = mateEnd - mateBegin;

    TContigSeq & contig = verifier.store.contigStore[contigId].seq;

    TContigSeqSize infixBeginPos;
    TContigSeqSize infixEndPos;

    if (reverseComplemented)
        _getContigInfix(verifier, contigId, beginPos, endPos, readId, infixBeginPos, infixEndPos, LeftMate());
    else
        _getContigInfix(verifier, contigId, beginPos, endPos, readId, infixBeginPos, infixEndPos, RightMate());

    TContigInfix contigInfix(contig, infixBeginPos, infixEndPos);
    TContigInfixIterator infixBegin = begin(contigInfix, Standard());
    TContigInfixIterator infixBorder = end(contigInfix, Standard()) - mateLength;

    for (TContigInfixIterator infixIt = infixBegin; infixIt != infixBorder; ++infixIt)
    {
        unsigned mateErrors = 0;

        for (TReadSeqIterator mateIt = mateBegin; mateIt != mateEnd; ++mateIt)
            if (*mateIt != *(infixIt + (mateIt - mateBegin)))
                if (++mateErrors > verifier.maxErrorsPerRead)
                    break;

        if (mateErrors <= verifier.maxErrorsPerRead)
        {
            verifier.pairsCount++;

//            if (reverseComplemented)
//            {
//                onMatch(verifier.matchesDelegate, contigId,
//                        (unsigned)(infixIt - infixBegin),
//                        (unsigned)((infixIt - infixBegin) + mateLength),
//                        mateId,
//                        mateErrors,
//                        beginPos, endPos, readId, errors);
//            }
//            else
//            {
//                onMatch(verifier.matchesDelegate, contigId,
//                        beginPos, endPos, readId, errors,
//                        (unsigned)(infixIt - infixBegin),
//                        (unsigned)((infixIt - infixBegin) + mateLength),
//                        mateId,
//                        mateErrors);
//            }
        }
    }
}

template <typename TMatchesDelegate, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Verifier<TMatchesDelegate, EditDistance, TSpec> & verifier,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    typedef Finder<TContigInfix>                    TFinder;

    if (verifier.disabled)
        return;

    // TODO Verifier must specify if readIds refer to left or right file
    // NOTE readId corresponds to left file
    readId *= 2;
    TReadId mateId = readId + 1;
    if (!reverseComplemented)
        mateId *= 2;

    TReadSeq & mate = verifier.store.readSeqStore[mateId];
    TContigSeq & contig = verifier.store.contigStore[contigId].seq;

    TContigSeqSize infixBegin;
    TContigSeqSize infixEnd;

    if (reverseComplemented)
        _getContigInfix(verifier, contigId, beginPos, endPos, readId, infixBegin, infixEnd, LeftMate());
    else
        _getContigInfix(verifier, contigId, beginPos, endPos, readId, infixBegin, infixEnd, RightMate());

    // DEBUG
//    std::cout << infixBegin << " " << infixEnd << std::endl;

    TContigInfix contigInfix(contig, infixBegin, infixEnd);

    TFinder finder(contigInfix);
    setHost(verifier.pattern, mate);

    bool paired = false;
    while (find(finder, verifier.pattern, -verifier.maxErrorsPerRead))
    {
//        findBegin(finder, verifier.pattern, getScore(verifier.pattern));

        // DEBUG
//        std::cout << "END " << position(finder) << std::endl;
//        std::cout << "BEGIN " << beginPosition(finder) << std::endl;

        paired = true;
//        verifier.pairsCount++;

//        if (reverseComplemented)
//        {
//            onMatch(verifier.matchesDelegate, contigId,
//                    (unsigned)(infixBegin + beginPosition(finder)),
//                    (unsigned)(infixBegin + endPosition(finder)),
//                    mateId,
//                    (unsigned char)-getScore(verifier.pattern),
//                    beginPos, endPos, readId, errors);
//        }
//        else
//        {
//            onMatch(verifier.matchesDelegate, contigId,
//                    beginPos, endPos, readId, errors,
//                    (unsigned)(infixBegin + beginPosition(finder)),
//                    (unsigned)(infixBegin + endPosition(finder)),
//                    mateId,
//                    (unsigned char)-getScore(verifier.pattern));
//        }
    }

    if (paired)
        verifier.pairsCount++;
}

//____________________________________________________________________________

template <typename THitsDelegate, typename TSpec, typename TReadId>
inline void preprocessMate(Verifier<THitsDelegate, EditDistance, Filter<TSpec> > & verifier, TReadId readId)
{
    // TODO Verifier must specify if readIds refer to left or right file
    // NOTE readId corresponds to left file
    readId *= 2;
    TReadId mateId = readId + 1;

    // Preprocess fwd mate
    _preprocessMate(verifier, verifier.patternFwd, verifier.seedsFwd, verifier.store.readSeqStore[mateId]);

    // Preprocess rev mate
    _preprocessMate(verifier, verifier.patternRev, verifier.seedsRev, verifier.store.readSeqStore[mateId * 2]);
}

template <typename THitsDelegate, typename TSpec, typename TPattern>
inline void _preprocessMate(Verifier<THitsDelegate, EditDistance, Filter<TSpec> > & verifier,
                            TPattern & pattern,
                            String<TReadSeq> & seeds,
                            TReadSeq & read)
{
    typedef Iterator<TReadSeq, Standard>::Type  TReadSeqIterator;
    typedef unsigned int                        TWord;

    TReadSeqSize seedsCount = verifier.maxErrorsPerRead + 1;
    TReadSeqSize readLength = length(read);
    verifier.seedLength = readLength / seedsCount;
//    verifier.seedLength = BitsPerValue<TWord>::VALUE / seedsCount;
//    verifier.seedLength = 64 / seedsCount;
    verifier.matchesDelegate.seedLength = verifier.seedLength;

    clear(seeds);
    reserve(seeds, seedsCount, Exact());

    for (TReadSeqSize seed = 0; seed < seedsCount; ++seed)
    {
        TReadSeqSize seedBegin = verifier.seedLength * seed;
        TReadSeqSize seedEnd = seedBegin + verifier.seedLength;

        appendValue(seeds, infix(read, seedBegin, seedEnd), Exact());
    }

    setHost(pattern, seeds);
}

template <typename THitsDelegate, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Verifier<THitsDelegate, EditDistance, Filter<TSpec> > & verifier,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    typedef Finder<TContigInfix>                    TFinder;

    if (verifier.disabled)
        return;

    // TODO Verifier must specify if readIds refer to left or right file
    // NOTE readId corresponds to left file
    readId *= 2;
    TReadId mateId = readId + 1;
    if (!reverseComplemented)
        mateId *= 2;

//    TReadSeq & mate = verifier.store.readSeqStore[mateId];
    TContigSeq & contig = verifier.store.contigStore[contigId].seq;

    TContigSeqSize infixBegin;
    TContigSeqSize infixEnd;

    if (reverseComplemented)
        _getContigInfix(verifier, contigId, beginPos, endPos, readId, infixBegin, infixEnd, LeftMate());
    else
        _getContigInfix(verifier, contigId, beginPos, endPos, readId, infixBegin, infixEnd, RightMate());

    // DEBUG
//    std::cout << infixBegin << " " << infixEnd << std::endl;

    TContigInfix contigInfix(contig, infixBegin, infixEnd);

    TFinder finder(contigInfix);

    // NOTE patternFwd and patternRev must be preprocessed with preprocessMate()
    unsigned hits = 0;
    if (reverseComplemented)
    {
        while (find(finder, verifier.patternFwd))
        {
            hits++;
//            verifier.matchesDelegate.minErrorsPerRead = position(verifier.patternFwd);
//            onSeedHit(verifier.matchesDelegate,
//                      contigId,
//                      (unsigned)(infixBegin + beginPosition(finder)),
//                      mateId,
//                      verifier.seedLength * position(verifier.patternFwd),
//                      0);
        }
    }
    else
    {
        while (find(finder, verifier.patternRev))
        {
            hits++;
//            verifier.matchesDelegate.minErrorsPerRead = position(verifier.patternRev);
//            onSeedHit(verifier.matchesDelegate,
//                      contigId,
//                      (unsigned)(infixBegin + beginPosition(finder)),
//                      mateId,
//                      verifier.seedLength * position(verifier.patternRev),
//                      0);
        }
    }
//    std::cout << "Hits: " << hits << std::endl;
    verifier.pairsCount += hits;
}

//____________________________________________________________________________

template <typename TMatchesDelegate, typename TDistance, typename TSpec, typename TReadId, typename TContigId, typename TContigPos>
inline void _getContigInfix(Verifier<TMatchesDelegate, TDistance, TSpec> & verifier,
                            TContigId contigId,
                            TContigPos beginPos,
                            TContigPos,
                            TReadId readId,
                            TContigSeqSize & infixBegin,
                            TContigSeqSize & infixEnd,
                            RightMate)
{
    TReadId mateId = readIdToPairId(readId, RightFile());
    TReadSeqSize mateLength = length(verifier.store.readSeqStore[mateId]);

    infixBegin = verifier.contigSizes[contigId];
    if (infixBegin > beginPos + verifier.libraryLength - verifier.libraryError - mateLength)
        infixBegin = beginPos + verifier.libraryLength - verifier.libraryError - mateLength;

    infixEnd = verifier.contigSizes[contigId];
    if (infixEnd > beginPos + verifier.libraryLength + verifier.libraryError)
        infixEnd = beginPos + verifier.libraryLength + verifier.libraryError;
}

template <typename TMatchesDelegate, typename TDistance, typename TSpec, typename TReadId, typename TContigId, typename TContigPos>
inline void _getContigInfix(Verifier<TMatchesDelegate, TDistance, TSpec> & verifier,
                            TContigId,
                            TContigPos,
                            TContigPos endPos,
                            TReadId readId,
                            TContigSeqSize & infixBegin,
                            TContigSeqSize & infixEnd,
                            LeftMate)
{
    TReadId mateId = readIdToPairId(readId, RightFile());
    TReadSeqSize mateLength = length(verifier.store.readSeqStore[mateId]);

    infixBegin = 0;
    if (endPos > verifier.libraryLength + verifier.libraryError)
        infixBegin = endPos - verifier.libraryLength - verifier.libraryError;

    infixEnd = 0;
    if (endPos > verifier.libraryLength - verifier.libraryError - mateLength)
        infixEnd = endPos - verifier.libraryLength + verifier.libraryError + mateLength;
}

#endif  // #ifndef SANDBOX_ESIRAGUSA_APPS_MASAI_MAPPER_H_
