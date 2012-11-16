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
#include <seqan/stream.h>

#include "tags.h"
#include "store.h"
#include "indexer.h"
#include "matches.h"
//#include "verifier.h"
#include "extender.h"
#include "stream.h"
#include "writer.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Pairer
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Pairer
{
    typedef Indexer<Nothing>    TIndexer;

    TFragmentStore      store;
    TIndexer            indexer;

    unsigned            readsCount;
    unsigned long       pairsCount;

    bool                writeCigar;
    bool                dumpResults;

    unsigned            libraryLength;
    unsigned            libraryError;

    // TODO(esiragusa): Remove writeCigar from Pairer members.
    Pairer(unsigned libraryLength, unsigned libraryError, bool writeCigar = true, bool dumpResults = true) :
        indexer(store),
        readsCount(0),
        pairsCount(0),
        writeCigar(writeCigar),
        dumpResults(dumpResults),
        libraryLength(libraryLength),
        libraryError(libraryError)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function loadReads()                                                [Pairer]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString>
bool loadReads(Pairer<TSpec> & pairer, TString const & readsLeftFile, TString const & readsRightFile)
{
    // TODO(esiragusa): Use loadReads() from store.h
    if (!loadReads(pairer.store, readsLeftFile, readsRightFile))
        return false;

    pairer.readsCount = length(pairer.store.readSeqStore);

    _loadReadsRC(pairer);

    return true;
}

template <typename TSpec>
bool _loadReadsRC(Pairer<TSpec> & pairer)
{
    for (TReadSeqStoreSize readId = 0; readId < pairer.readsCount; ++readId)
    {
        TReadSeq & read = pairer.store.readSeqStore[readId];
        appendValue(pairer.store.readSeqStore, read);
        reverseComplement(back(pairer.store.readSeqStore));
    }

    return true;
}

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

/*
template <typename TSpec, typename TString, typename TErrors, typename TDistance>
bool mateMappedReads(Pairer<TSpec> & pairer,
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

template <typename TSpec, typename TString, typename TErrors, typename TDistance>
bool mateMappedReads(Pairer<TSpec> & pairer,
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

template <typename TSpec, typename TString, typename TMatchesDelegate>
bool mateMappedReads(Pairer<TSpec> & pairer,
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

template <typename TSpec, typename TRecordSpec, typename TStringSpec, typename TMatchesDelegate>
inline void _matePair(Pairer<TSpec> &,
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
        onMatch(matchesDelegate, *matchesIt);
}
*/

// ----------------------------------------------------------------------------
// Function mateMappedReads()                                          [Pairer]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString, typename TDistance>
bool mateMappedReads(Pairer<TSpec> & pairer,
                     TString const & mappedReadsLeftFile,
                     TString const & mappedReadsRightFile,
                     TString const & mappedPairsFile,
                     TDistance const & /*tag*/,
                     Raw const & /*tag*/)
{
    typedef External<ExternalConfigLarge<> >            TStream;
    typedef String<Match<>, TStream>                    TWriterStream;
//    typedef Stream<FileStream<Match<>, MMapWriter> >    TWriterStream;
    typedef MatchWriter<TWriterStream, TDistance, Raw>  TMatchWriter;

    TWriterStream file;
    
    if (pairer.dumpResults)
        if (!open(file, toCString(mappedPairsFile), OPEN_RDWR | OPEN_CREATE))
            return false;

    TMatchWriter writer(file, pairer.store, pairer.readsCount, pairer.dumpResults);

    _mateMappedReads(pairer, mappedReadsLeftFile, mappedReadsRightFile, writer);

    std::cout << "Pairs:\t\t\t\t" << pairer.pairsCount << std::endl;

    return true;
}

template <typename TSpec, typename TString, typename TDistance>
bool mateMappedReads(Pairer<TSpec> & pairer,
                     TString const & mappedReadsLeftFile,
                     TString const & mappedReadsRightFile,
                     TString const & mappedPairsFile,
                     TDistance const & /*tag*/,
                     Sam const & /*tag*/)
{
    typedef External<ExternalConfigLarge<> >            TStream;
    typedef String<char, TStream>                       TWriterStream;
//    typedef Stream<FileStream<char, MMapWriter> >       TWriterStream;
    typedef MatchWriter<TWriterStream, TDistance, Sam> TMatchWriter;

    TWriterStream file;
    
    if (pairer.dumpResults)
        if (!open(file, toCString(mappedPairsFile), OPEN_RDWR | OPEN_CREATE))
            return false;

    TMatchWriter writer(file, pairer.store, pairer.readsCount, pairer.dumpResults);

    // TODO(esiragusa): Remove writeCigar from members.
    writer.writeCigar = pairer.writeCigar;

    _mateMappedReads(pairer, mappedReadsLeftFile, mappedReadsRightFile, writer);

    std::cout << "Pairs:\t\t\t\t" << pairer.pairsCount << std::endl;

    return true;
}

// ============================================================================

template <typename TSpec, typename TString, typename TMatchesDelegate>
bool _mateMappedReads(Pairer<TSpec> & pairer,
                      TString const & mappedReadsLeftFile,
                      TString const & mappedReadsRightFile,
                      TMatchesDelegate & matchesDelegate)
{
    typedef Match<>                  TMatch;
    typedef MatchStore<TMatch>       TMatchStore;
    typedef String<TMatch>           TMatches;

    TMatchStore storeLeft;
    TMatchStore storeRight;

    if (!open(storeLeft, mappedReadsLeftFile))
        return false;

    if (!open(storeRight, mappedReadsRightFile))
    {
        close(storeLeft);
        return false;
    }

    TMatches matchesLeft;
    TMatches matchesRight;

    if (!getNext(storeLeft, matchesLeft) || !getNext(storeRight, matchesRight))
        return false;

    do
    {
        while (front(matchesLeft).readId > front(matchesRight).readId && getNext(storeRight, matchesRight)) ;

        if (front(matchesLeft).readId == front(matchesRight).readId)
        {
            removeDuplicateMatches(matchesLeft);
            removeDuplicateMatches(matchesRight);
            _matePair(pairer, matchesLeft, matchesRight, matchesDelegate);
        }
    }
    while (getNext(storeLeft, matchesLeft));

    close(storeRight);
    close(storeLeft);

    return true;
}

// ============================================================================

template <typename TSpec, typename TRecordSpec, typename TStringSpec, typename TMatchesDelegate>
inline void _matePair(Pairer<TSpec> & pairer,
                      String<Match<TRecordSpec>, TStringSpec> const & matchesLeft,
                      String<Match<TRecordSpec>, TStringSpec> const & matchesRight,
                      TMatchesDelegate & matchesDelegate)
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
                              LeftFile(), RightFile(),
                              matchesDelegate);
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
                              RightFile(), LeftFile(),
                              matchesDelegate);
//
//                    std::cout << "==================" << std::endl;
                }
            }
            else if (leftContigId < rightContigId)
                break;
        }
    }
}

// ============================================================================

template <typename TSpec, typename TIterator, typename TMateFwd, typename TMateRev, typename TMatchesDelegate>
inline void _matePair(Pairer<TSpec> & pairer,
                      TIterator mateFwdBegin, TIterator mateFwdEnd,
                      TIterator mateRevBegin, TIterator mateRevEnd,
                      TMateFwd const &,
                      TMateRev const &,
                      TMatchesDelegate & matchesDelegate)
{
    typedef Match<> TMatch;

    // NOTE mateFwd queue C= mateRev queue, i.e. mateFwdTail >= mateRevTail && mateFwdHead <= mateRevHead

    SEQAN_ASSERT_NEQ(mateFwdBegin, mateFwdEnd);
    unsigned mateFwdId = readIdToPairId((*mateFwdBegin).readId, TMateFwd());
    TReadSeqSize mateFwdLength = length(pairer.store.readSeqStore[mateFwdId]);

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

            onMatch(matchesDelegate, mateFwd, mateRev);

//            onMatch(matchesDelegate, *mateFwdQueueIt, *mateRevIt);
        }

        // Empty mateFwd queue.
//        if (mateFwdTailIt == mateFwdHeadIt) continue;
    }
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_PAIRER_H_
