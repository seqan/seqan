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
// This file contains the Pairer class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_PAIRER_H_
#define SEQAN_EXTRAS_MASAI_PAIRER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "tags.h"
#include "store.h"
#include "matches.h"
//#include "verifier.h"
//#include "extender.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Pairer
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec = void>
struct Pairer
{
    typedef MatchStore<>    TMatchStore;

    Holder<TReads>      reads;
    TMatchStore         _storeLeft;
    TMatchStore         _storeRight;
    TDelegate           & delegate;
    unsigned long       pairsCount;
    unsigned            libraryLength;
    unsigned            libraryError;

    Pairer(TDelegate & delegate, unsigned libraryLength, unsigned libraryError) :
        delegate(delegate),
        pairsCount(0),
        libraryLength(libraryLength),
        libraryError(libraryError)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readIdToPairId()
// ----------------------------------------------------------------------------

template <typename TReadId>
TReadId readIdToPairId(TReadId readId, LeftFile)
{
    return readId * 2;
}

template <typename TReadId>
TReadId readIdToPairId(TReadId readId, RightFile)
{
    return readId * 2 + 1;
}

// ----------------------------------------------------------------------------
// Function open()                                                     [Pairer]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TString>
bool open(Pairer<TReads, TDelegate, TSpec> & sorter, TString const & mappedReadsLeftFile, TString const & mappedReadsRightFile)
{
    return open(sorter._storeLeft, mappedReadsLeftFile) && open(sorter._storeRight, mappedReadsRightFile);
}

// ----------------------------------------------------------------------------
// Function close()                                                    [Pairer]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec>
bool close(Pairer<TReads, TDelegate, TSpec> & pairer)
{
    return close(pairer._storeLeft) && close(pairer._storeRight);
}

// ----------------------------------------------------------------------------
// Function setReads()                                                 [Pairer]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec>
void setReads(Pairer<TReads, TDelegate, TSpec> & pairer, TReads & reads)
{
    setValue(pairer.reads, reads);
}

// ----------------------------------------------------------------------------
// Function getReads()                                                 [Pairer]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec>
inline typename Reference<TReads>::Type
getReads(Pairer<TReads, TDelegate, TSpec> & pairer)
{
    return value(pairer.reads);
}

// ----------------------------------------------------------------------------
// Function pair()                                                     [Pairer]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec>
void pair(Pairer<TReads, TDelegate, TSpec> & pairer)
{
//    typedef Pairer<TReads, TDelegate, TSpec>        TPairer;
//    typedef typename TPairer::TMatchStore   TMatchStore;
//    typedef typename TMatchStore::TMatches  TMatches;
    typedef String<Match<> >                TMatches;

    TMatches matchesLeft;
    TMatches matchesRight;

    if (!getNext(pairer._storeLeft, matchesLeft) || !getNext(pairer._storeRight, matchesRight)) return;

    do
    {
        while (front(matchesLeft).readId > front(matchesRight).readId && getNext(pairer._storeRight, matchesRight)) ;

        if (front(matchesLeft).readId == front(matchesRight).readId)
        {
            removeDuplicateMatches(matchesLeft);
            removeDuplicateMatches(matchesRight);
            _matePair(pairer, matchesLeft, matchesRight);
        }
    }
    while (getNext(pairer._storeLeft, matchesLeft));
}

// ----------------------------------------------------------------------------
// Function _matePair()                                                [Pairer]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TRecordSpec, typename TStringSpec>
inline void _matePair(Pairer<TReads, TDelegate, TSpec> & pairer,
                      String<Match<TRecordSpec>, TStringSpec> const & matchesLeft,
                      String<Match<TRecordSpec>, TStringSpec> const & matchesRight)
{
    typedef Match<TRecordSpec>                                      TMatch;
    typedef MatchIterator<TMatch, TStringSpec>                      TMatchIterator;
    typedef typename MatchIterator<TMatch, TStringSpec>::TIterator  TIterator;

    TIterator matchesLeftFwdBegin;
    TIterator matchesLeftRevBegin;
    TIterator matchesLeftFwdEnd;
    TIterator matchesLeftRevEnd;

    TIterator matchesRightFwdBegin;
    TIterator matchesRightRevBegin;
    TIterator matchesRightFwdEnd;
    TIterator matchesRightRevEnd;

    TMatchIterator leftIt(matchesLeft);
    TMatchIterator rightIt(matchesRight);

    unsigned char leftContigId;
    unsigned char rightContigId;

    while (getNextContig(leftIt, matchesLeftFwdBegin, matchesLeftFwdEnd,
                         matchesLeftRevBegin, matchesLeftRevEnd))
    {
        if (matchesLeftFwdEnd != matchesLeftFwdBegin)
            leftContigId = (*matchesLeftFwdBegin).contigId;
        else
            leftContigId = (*matchesLeftRevBegin).contigId;

        while (getNextContig(rightIt, matchesRightFwdBegin, matchesRightFwdEnd,
                             matchesRightRevBegin, matchesRightRevEnd))
        {
            if (matchesRightFwdEnd != matchesRightFwdBegin)
                rightContigId = (*matchesRightFwdBegin).contigId;
            else
                rightContigId = (*matchesRightRevBegin).contigId;

            if (leftContigId == rightContigId)
            {
                if (matchesLeftFwdEnd != matchesLeftFwdBegin && matchesRightRevEnd != matchesRightRevBegin)
                {
//                    std::cout << "Left FWD vs Right REV" << std::endl;
//                    printMatches(matchesLeft, matchesLeftFwdBegin, matchesLeftFwdEnd);
//                    printMatches(matchesRight, matchesRightRevBegin, matchesRightRevEnd);
//
                    _matePair(pairer,
                              matchesLeftFwdBegin, matchesLeftFwdEnd,
                              matchesRightRevBegin, matchesRightRevEnd,
                              LeftFile(), RightFile());
//
//                    std::cout << "==================" << std::endl;
                }

                if (matchesRightFwdEnd != matchesRightFwdBegin && matchesLeftRevEnd != matchesLeftRevBegin)
                {
//                    std::cout << "Left REV vs Right FWD" << std::endl;
//                    printMatches(matchesLeft, matchesLeftRevBegin, matchesLeftRevEnd);
//                    printMatches(matchesRight, matchesRightFwdBegin, matchesRightFwdEnd);
//
                    _matePair(pairer,
                              matchesRightFwdBegin, matchesRightFwdEnd,
                              matchesLeftRevBegin, matchesLeftRevEnd,
                              RightFile(), LeftFile());
//
//                    std::cout << "==================" << std::endl;
                }
            }
            else if (leftContigId < rightContigId)
                break;
        }
    }
}

// ----------------------------------------------------------------------------
// Function _matePair()                                                [Pairer]
// ----------------------------------------------------------------------------

template <typename TReads, typename TDelegate, typename TSpec, typename TIterator, typename TMateFwd, typename TMateRev>
inline void _matePair(Pairer<TReads, TDelegate, TSpec> & pairer,
                      TIterator mateFwdBegin, TIterator mateFwdEnd,
                      TIterator mateRevBegin, TIterator mateRevEnd,
                      TMateFwd const & /* tag */,
                      TMateRev const & /* tag */)
{
    typedef Match<> TMatch;

    // NOTE mateFwd queue C= mateRev queue, i.e. mateFwdTail >= mateRevTail && mateFwdHead <= mateRevHead

    SEQAN_ASSERT_NEQ(mateFwdBegin, mateFwdEnd);
    unsigned mateFwdId = readIdToPairId((*mateFwdBegin).readId, TMateFwd());
    TReadSeqSize mateFwdLength = length(getSeqs(getReads(pairer))[mateFwdId]);

    TIterator mateFwdIt = mateFwdBegin;

    // Get next mateRev match.
    for (TIterator mateRevIt = mateRevBegin; mateRevIt != mateRevEnd; ++mateRevIt)
    {
        // Compute mateRev match queue.
        TContigSeqSize mateRevHead = endPos(*mateRevIt) - pairer.libraryLength + pairer.libraryError + mateFwdLength;
        TContigSeqSize mateRevTail = endPos(*mateRevIt) - pairer.libraryLength - pairer.libraryError;

        // Seek next mateFwd match after the tail of mateRev match queue.
        while (mateFwdIt != mateFwdEnd && (*mateFwdIt).beginPos < mateRevTail)
            ++mateFwdIt;

        TIterator mateFwdTailIt = mateFwdIt;
        if (mateFwdTailIt == mateFwdEnd)
            break;

        // Continue if mateFwd tail is beyond mateRev head.
        TContigSeqSize mateFwdTail = (*mateFwdTailIt).beginPos;
        if (mateFwdTail >= mateRevHead)
            continue;

        // Seek next mateFwd match after the head of mateRev match queue.
        while (mateFwdIt != mateFwdEnd && endPos(*mateFwdIt) <= mateRevHead)
            ++mateFwdIt;
        TIterator mateFwdHeadIt = mateFwdIt;

        // Couple all fwds mates vs rev mate.
        for (TIterator mateFwdQueueIt = mateFwdTailIt; mateFwdQueueIt != mateFwdHeadIt; ++mateFwdQueueIt)
        {
            pairer.pairsCount++;

            // Convert readIds to pairIds
            TMatch mateFwd;
            TMatch mateRev;
            assign(mateFwd, *mateFwdQueueIt);
            assign(mateRev, *mateRevIt);
            mateFwd.readId = readIdToPairId((*mateFwdQueueIt).readId, TMateFwd());
            mateRev.readId = readIdToPairId((*mateRevIt).readId, TMateRev());

            onMatch(pairer.delegate, mateFwd, mateRev);

//            onMatch(pairer.delegate, *mateFwdQueueIt, *mateRevIt);
        }

        // Empty mateFwd queue.
//        if (mateFwdTailIt == mateFwdHeadIt) continue;
    }
}

/*
template <typename TReads, typename TDelegate, typename TSpec, typename TString, typename TErrors, typename TDistance>
bool mateMappedReads(Pairer<TReads, TDelegate, TSpec> & pairer,
                     TString const & mappedReadsLeftFile,
                     TString const & mappedPairsFile,
                     TErrors errors,
                     TDistance, Raw)
{
    typedef String<Match<>, TStream>                        TWriterStream;
    typedef MatchWriter<TWriterStream, TDistance, Raw>      TMatchWriter;
    typedef Verifier<TMatchWriter, TDistance>               TVerifier;

//    typedef Extender<TMatchWriter, TDistance>               TExtender;
//    typedef Verifier<TExtender, TDistance, Filter<void> >   TVerifier;

    TWriterStream file;
    if (pairer.dumpResults)
        open(file, toCString(mappedPairsFile), OPEN_RDWR | OPEN_CREATE);

    TMatchWriter writer(file, pairer.store, pairer.readsCount, pairer.dumpResults);
    TVerifier verifier(pairer.store, writer, pairer.libraryLength, pairer.libraryError);

//    TExtender extender(pairer.store, writer, pairer.readsCount, 0, true);
//    extender.minErrorsPerRead = 0;
//    extender.maxErrorsPerRead = errors;
//
//    TVerifier verifier(pairer.store, extender, pairer.libraryLength, pairer.libraryError);
    verifier.maxErrorsPerRead = errors;

    mateMappedReads(pairer, mappedReadsLeftFile, verifier);

//    std::cout << "Pairs:\t\t\t\t" << verifier.pairsCount << std::endl;

    return true;
}

template <typename TReads, typename TDelegate, typename TSpec, typename TString, typename TErrors, typename TDistance>
bool mateMappedReads(Pairer<TReads, TDelegate, TSpec> & pairer,
                     TString const & mappedReadsLeftFile,
                     TString const & mappedPairsFile,
                     TErrors errors,
                     TDistance, Sam)
{
    typedef String<char, TStream>                           TWriterStream;
    typedef MatchWriter<TWriterStream, TDistance, Sam>      TMatchWriter;
    typedef Verifier<TMatchWriter, TDistance>               TVerifier;

    TWriterStream file;
    if (pairer.dumpResults)
        open(file, toCString(mappedPairsFile), OPEN_RDWR | OPEN_CREATE);

    TMatchWriter writer(file, pairer.store, pairer.readsCount, pairer.dumpResults);
    TVerifier verifier(pairer.store, writer, pairer.libraryLength, pairer.libraryError);
    verifier.maxErrorsPerRead = errors;

    // TODO(esiragusa):Remove writeCigar from members.
    writer.writeCigar = pairer.writeCigar;

    mateMappedReads(pairer, mappedReadsLeftFile, verifier);

//    std::cout << "Pairs:\t\t\t\t" << verifier.pairsCount << std::endl;

    return true;
}

template <typename TReads, typename TDelegate, typename TSpec, typename TString>
bool mateMappedReads(Pairer<TReads, TDelegate, TSpec> & pairer,
                     TString const & mappedReadsLeftFile,
                     TMatchesDelegate & matchesDelegate)
{
    typedef Match<>                  TMatch;
    typedef MatchStore<TMatch>       TMatchStore;
    typedef String<TMatch>           TMatches;

    TMatchStore storeLeft;

    if (!open(storeLeft, mappedReadsLeftFile))
        return false;

    TMatches matchesLeft;

    if (!getNext(storeLeft, matchesLeft))
        return false;

    do
    {
        removeDuplicateMatches(matchesLeft);
        _matePair(pairer, matchesLeft, matchesDelegate);
    }
    while (getNext(storeLeft, matchesLeft));

    close(storeLeft);

    return true;
}

template <typename TReads, typename TDelegate, typename TSpec, typename TRecordSpec, typename TStringSpec>
inline void _matePair(Pairer<TReads, TDelegate, TSpec> &,
                      String<Match<TRecordSpec>, TStringSpec> const & matchesLeft,
                      TMatchesDelegate & matchesDelegate)
{
    typedef Match<TRecordSpec>                              TMatch;
    typedef String<TMatch, TStringSpec>                     TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;

    TIterator matchesIt = begin(matchesLeft, Standard());
    TIterator matchesEnd = end(matchesLeft, Standard());

    preprocessMate(matchesDelegate, (*matchesIt).readId);

    for (; matchesIt != matchesEnd; ++matchesIt)
        onMatch(, *matchesIt);
}
*/

#endif  // #ifndef SEQAN_EXTRAS_MASAI_PAIRER_H_
