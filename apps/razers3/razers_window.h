/*==========================================================================
 RazerS - Fast Read Mapping with Controlled Loss Rate
 http://www.seqan.de/projects/razers.html

 ============================================================================
 Copyright (C) 2010

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your options) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_RAZERS_WINDOW_H
#define SEQAN_HEADER_RAZERS_WINDOW_H

#include <iostream>

namespace seqan {

//////////////////////////////////////////////////////////////////////////////
// Find read matches in a single genome sequence
//
// Creates finder on the contig given by the ID and a verifier.
// Searches through the contig using the findWindowNext() function.
// The results are dumped in the (aligned) store.
template <
    typename TFragmentStore,
    typename TReadIndex,
    typename TSwiftSpec,
    typename TPreprocessing,
    typename TCounts,
    typename TRazerSOptions,
    typename TRazerSMode>
void _mapSingleReadsToContigWindow(
    TFragmentStore & store,
    unsigned                                  contigId,                     // ... and its sequence number
    Pattern<TReadIndex, Swift<TSwiftSpec> > & swiftPattern,
    TPreprocessing & preprocessing,
    TCounts & cnts,
    char                                      orientation,                      // q-gram index of reads
    TRazerSOptions & options,
    TRazerSMode                       const & mode)
{
    _mapSingleReadsToContigWindow(store, store, contigId, swiftPattern, swiftPattern, preprocessing, cnts, orientation, options, mode);
}

//////////////////////////////////////////////////////////////////////////////
// Find read matches in a single genome sequence
//
// Creates finder on the contig given by the ID and a verifier.
// Searches through the contig using the findWindowNext() function.
// The results are dumped in the (aligned) store.
// Specialized version: Contigs are taken from the main store but the results
//   are written to the block store. Used in the parallel reads over whole
//   genome version.
template <
    typename TFragmentStore,
    typename TReadIndex,
    typename TSwiftSpec,
    typename TSwiftPatternHandler,
    typename TPreprocessing,
    typename TCounts,
    typename TRazerSOptions,
    typename TRazerSMode>
void _mapSingleReadsToContigWindow(
    TFragmentStore & mainStore,
    TFragmentStore & blockStore,
    unsigned                                  contigId,                     // ... and its sequence number
    Pattern<TReadIndex, Swift<TSwiftSpec> > & swiftPattern,
    TSwiftPatternHandler & swiftPatternHandler,
    TPreprocessing & preprocessing,
    TCounts & cnts,
    char                                      orientation,                      // q-gram index of reads
    TRazerSOptions & options,
    TRazerSMode                      const & mode)
{
    // FILTRATION
    typedef typename TFragmentStore::TContigSeq             TContigSeq;
    typedef Finder<TContigSeq, Swift<TSwiftSpec> >          TSwiftFinder;
    typedef Pattern<TReadIndex, Swift<TSwiftSpec> >         TSwiftPattern;

    // VERIFICATION
    typedef MatchVerifier<
        TFragmentStore,
        TRazerSOptions,
        TRazerSMode,
        TPreprocessing,
        TSwiftPattern,
        TCounts>                                           TVerifier;
    typedef typename Fibre<TReadIndex, FibreText>::Type TReadSet;

    // HITS
    typedef typename TSwiftFinder::THitString               THitString;
    typedef typename Value<THitString>::Type                TSwiftHit;
    typedef typename Size<THitString>::Type                 THitStringSize;

    typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;

    // output what is done if verbous
    if (options._debugLevel >= 1)
    {
        std::cerr << std::endl << "Process genome seq #" << contigId;
        if (orientation == 'F')
            std::cerr << "[fwd]";
        else
            std::cerr << "[rev]";
    }
    // lock contig
    lockContig(mainStore, contigId);
    TContigSeq & contigSeq = mainStore.contigStore[contigId].seq;
    if (orientation == 'R')
        reverseComplement(contigSeq);

    // Create finder and verifier
    TSwiftFinder    swiftFinder(contigSeq, options.repeatLength, 1);
    TVerifier       verifier(blockStore, options, preprocessing, swiftPattern, cnts);

    // initialize verifier
    verifier.onReverseComplement = (orientation == 'R');
    verifier.genomeLength = length(contigSeq);
    verifier.m.contigId = contigId;

    // if the pattern can be initialized and there is a non-repeat region in the contig that fits a qgram.
    if (windowFindBegin(swiftFinder, swiftPattern, options.errorRate))
    {

        _proFloat myTime = sysTime();

        bool sequenceLeft = true;

        // while there is more contig sequence to search through
        while (sequenceLeft)
        {
            sequenceLeft = windowFindNext(swiftFinder, swiftPattern, 1000000);

            printf("filter: %f sec\n", sysTime() - myTime);
            myTime = sysTime();

            // get the found hits from the finder
            THitString hits = getSwiftHits(swiftFinder);
            // verifiy them
            for (THitStringSize h = 0; h < length(hits); ++h)
            {
                verifier.m.readId = hits[h].ndlSeqNo;             //array oder jedesmal berechnen
                matchVerify(verifier, swiftInfix(hits[h], contigSeq), hits[h].ndlSeqNo, host(host(swiftPattern))[hits[h].ndlSeqNo], mode);
                ++options.countFiltration;
            }

            printf("verify: %f sec\n", sysTime() - myTime);
            myTime = sysTime();

            // compact matches if neccessary
            if (length(blockStore.alignedReadStore) > options.compactThresh)
            {
                TAlignedReadStoreSize oldSize = length(blockStore.alignedReadStore);
                if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE || options.threshold == 0)
                    maskDuplicates(blockStore, mode);       // overlapping parallelograms cause duplicates

                compactMatches(blockStore, cnts, options, mode, swiftPatternHandler, COMPACT);
                if (options._debugLevel >= 2)
                    std::cerr << '(' << oldSize - length(blockStore.alignedReadStore) << " matches removed)";
            }

        }

        // clear finders
        windowFindEnd(swiftFinder, swiftPattern);

    }

    if (!unlockAndFreeContig(mainStore, contigId))                              // if the contig is still used
        if (orientation == 'R')
            reverseComplement(contigSeq);
    // we have to restore original orientation
}

}

#endif
