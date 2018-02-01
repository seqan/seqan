// ==========================================================================
//                              casbar
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Sabrina Krakau <sabrina.krakau@fu-berlin.de>
// ==========================================================================

//#define SEQAN_ENABLE_PARALLELISM

//#define CALL_PROFILE
//#define SNPSTORE_DEBUG

#include <seqan/platform.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/parallel.h>

//#ifdef STDLIB_VS
//#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
//#else
//#define SEQAN_DEFAULT_TMPDIR "./"
//#endif

#include <cmath>
#include <iostream>

#include <map>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>
#include <seqan/bed_io.h>

#include <boost/detail/workaround.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp>

#include "casbar_alphabets.h"
#include "casbar_util.h"
#include "casbar_realignment.h"
#include "casbar.h"
#include "casbar_calling.h"

//using namespace std;
using namespace seqan;

// load entire genome into memory
template <typename TGenomeSet, typename TGenomeNames>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList, ::std::map<CharString,unsigned> &gIdStringToIdNumMap, TGenomeNames & genomeNames)
{
    clear(genomes);
    clear(genomeNames);

    seqan::CharString id, seq;
    for (unsigned filecount = 0; filecount < length(fileNameList); ++filecount)
    {
        seqan::SeqFileIn seqFileIn;
        if (!open(seqFileIn, toCString(fileNameList[filecount])))
            return false;

        unsigned seqCount = 0;
        for (; !atEnd(seqFileIn); ++seqCount)
        {
            readRecord(id, seq, seqFileIn);
            cropAfterFirst(id, IsBlank());

            gIdStringToIdNumMap.insert(std::make_pair(id, length(genomes)));

            appendValue(genomes, seq);
            appendValue(genomeNames, id);
        }
    }

    return !empty(genomes);
}

// transform global cooridnates to coordinates relative to chromosomal segment
template<typename TFragmentStore, typename TContigPos, typename TOptions>
void
transformCoordinates(TFragmentStore &fragmentStore, TContigPos startCoord, TOptions&)
{
    typedef typename TFragmentStore::TAlignedReadStore          TMatches;
    //typedef typename Value<TMatches>::Type                      TMatch;
    typedef typename Iterator<TMatches,Standard>::Type          TMatchIt;

    TMatchIt mIt        = begin(fragmentStore.alignedReadStore,Standard());
    TMatchIt mItEnd     = end(fragmentStore.alignedReadStore,Standard());

    while(mIt != mItEnd)
    {
        (*mIt).endPos -= startCoord;
        (*mIt).beginPos -= startCoord;
        ++mIt;
    }

}

// copy matches relevant for next window
template<typename TFragmentStore, typename TSetContigAnchorGaps, typename TSize, typename TContigPos, typename TOptions>
void
copyNextWindowMatchesAndReads(TFragmentStore &fragmentStore,
                              TSetContigAnchorGaps &setContigAnchorGaps,
                              typename TFragmentStore::TReadSeqStore &tmpReads,
                              typename TFragmentStore::TReadStore &tmpRs,
                              typename TFragmentStore::TAlignedReadStore &tmpMatches,
                              typename TFragmentStore::TAlignQualityStore &tmpQualities,
                              TSetContigAnchorGaps &tmpSetContigAnchorGaps,
                              TSize ,
                              TContigPos currentWindowEnd,
                              TOptions &options)
{

    typedef typename TFragmentStore::TAlignedReadStore          TMatches;
    typedef typename Value<TMatches>::Type                      TMatch;
    typedef typename Iterator<TMatches,Standard>::Type          TMatchIt;
    typedef typename Iterator<TSetContigAnchorGaps,Standard>::Type          TSetContigGapsIter;
    typedef typename Id<TFragmentStore>::Type                   TId;
    //typedef typename Value<TReadClips>::Type                    TPair;

    SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore),length(fragmentStore.alignQualityStore));

    ::std::sort(begin(fragmentStore.alignedReadStore, Standard()), end(fragmentStore.alignedReadStore, Standard()), LessGPos<TMatch>());

    if(options._debugLevel > 1 )::std::cout << "Copying matches overlapping more than one window ... \n";

    TMatchIt mIt        = end(fragmentStore.alignedReadStore,Standard());
    TMatchIt mItBegin   = begin(fragmentStore.alignedReadStore,Standard());
    --mIt;

    TSetContigGapsIter itG = end(setContigAnchorGaps);
    //TSetContigGapsIter itGBegin = begin(setContigAnchorGaps);
    --itG;

    // We will use minCoord/maxCoord to store the temporarily minimal and maximal coordinates in the window.
    int minCoord = std::numeric_limits<int>::max();
    int maxCoord = std::numeric_limits<int>::min();
    //CharString str = "discBef";
    //_dumpMatches(fragmentStore, str);

    //std::cout << " max hit length = " << options.maxHitLength << std::endl;
    // look for matches that are inside our window of interest, copy corresponding matches,reads,qualities
    while(mIt >= mItBegin && _min((*mIt).beginPos,(*mIt).endPos) + (TContigPos)options.maxHitLength + (TContigPos)options.windowBuff >= currentWindowEnd )
    {
        if( _max((*mIt).beginPos,(*mIt).endPos) + (TContigPos)options.windowBuff > currentWindowEnd )
        {
            TId id = length(tmpMatches);
            appendValue(tmpMatches,*mIt);
            tmpMatches[id].id = id;
            tmpMatches[id].readId = id;
            appendValue(tmpReads,fragmentStore.readSeqStore[(*mIt).readId]);
            appendValue(tmpRs,fragmentStore.readStore[(*mIt).readId]);
            appendValue(tmpQualities,fragmentStore.alignQualityStore[(*mIt).readId]);
            appendValue(tmpSetContigAnchorGaps, *itG);
            maxCoord = std::max(maxCoord, (int)std::max(mIt->beginPos, mIt->endPos));
            minCoord = std::min(minCoord, (int)std::min(mIt->beginPos, mIt->endPos));
        }
        --mIt;
        --itG;
    }

    // Write minimal and maximal coordinate from reads in this window to options.minCoord/options.maxCoord.
    if (minCoord != std::numeric_limits<int>::max())
        options.minCoord = minCoord;
    if (maxCoord != std::numeric_limits<int>::min())
        options.maxCoord = maxCoord;

    if(options._debugLevel > 1)
        std::cout << length(tmpMatches)<<" matches left over from previous window.\n";
}


// little helper
template<typename TMatch>
char
orientation(TMatch & match)
{
    if(match.endPos > match.beginPos)
        return 'F';
    else
        return 'R';
}

// discard reads stacking up, give preference to high quality reads
template<typename TFragmentStore, typename TSize, typename TOptions>
void
applyPileupCorrection(TFragmentStore    &fragmentStore,
                      TSize                         arrayBeginPos,
                      TSize                         arrayEndPos,
                      TOptions                      &options)
{
    typedef StringSet<String<Dna5Q>, Owner<ConcatDirect<> > >    TFalseType;
    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    //typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename TFragmentStore::TReadSeqStore      TReads;
    //typedef typename Value<TReads>::Type                TRead;
    //typedef typename TFragmentStore::TContigStore       TGenomeSet;
    //typedef typename Value<TGenomeSet>::Type            TGenome;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;

    if(IsSameType<TReads, TFalseType >::VALUE)
        std::cout << "Hier stimmt was nciht. strinsetspec concat direct\n";

    if(options._debugLevel > 0) std::cout << arrayEndPos - arrayBeginPos  << " matches subject to pile up correction." << std::endl;

    if(options.orientationAware)
        ::std::sort(iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard()),
                    iter(fragmentStore.alignedReadStore, arrayEndPos, Standard()),
                    LessGStackOaMQ<TMatches,TMatchQualities>(fragmentStore.alignQualityStore));
    else
        ::std::sort(iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard()),
                    iter(fragmentStore.alignedReadStore, arrayEndPos, Standard()),
                    LessGStackMQ<TMatches,TMatchQualities>(fragmentStore.alignQualityStore));


    TMatchIterator matchIt          = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
    TMatchIterator matchRangeEnd    = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
    TMatchIterator matchItKeep      = matchIt;

    while(matchIt != matchRangeEnd)
    {
        TContigPos currentBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
        TContigPos currentEnd = _max((*matchIt).beginPos,(*matchIt).endPos);
        unsigned currentSeqno = (*matchIt).contigId;
        char currentOrientation = orientation(*matchIt);
        unsigned currPile = 0;
        while(matchIt != matchRangeEnd
              && (*matchIt).contigId == currentSeqno
              && _min((*matchIt).beginPos,(*matchIt).endPos) == currentBegin
              && _max((*matchIt).beginPos,(*matchIt).endPos) == currentEnd
              && (!options.orientationAware || orientation(*matchIt) == currentOrientation)
              && currPile < options.maxPile)
        {
            *matchItKeep = *matchIt;
            ++currPile;
            ++matchIt;
            ++matchItKeep;
        }
        //if(matchRangeEnd > matchItEnd) ::std::cerr <<"neeeeeeee\n";
        while(matchIt != matchRangeEnd
              && (*matchIt).contigId == currentSeqno
              && _min((*matchIt).beginPos,(*matchIt).endPos) == currentBegin
              && _max((*matchIt).beginPos,(*matchIt).endPos) == currentEnd
              && (!options.orientationAware || orientation(*matchIt) == currentOrientation))
            ++matchIt;

    }

    if(options._debugLevel > 0) std::cout << matchIt - matchItKeep << " matches discarded." << std::endl;
    resize(fragmentStore.alignedReadStore,matchItKeep - begin(fragmentStore.alignedReadStore,Standard()));

    //  ::std::cout << "numMatches = " << length(fragmentStore.alignedReadStore) << ::std::endl;
    SEQAN_ASSERT_LEQ(length(fragmentStore.alignedReadStore), length(fragmentStore.alignQualityStore));
    SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.alignQualityStore));

    //  str="pileAft";
    //  _dumpMatches(fragmentStore,str);
}

// average quality of read is kept as extra info for each match
// used for prioritization in pile up correction
template<typename TFragmentStore, typename TSize, typename TOptions>
void
addReadQualityToMatches(TFragmentStore  &fragmentStore,
                        TSize                           arrayBeginPos,
                        TSize                           arrayEndPos,
                        TOptions &)
{
    typedef typename TFragmentStore::TAlignedReadStore      TMatches;
    typedef typename TFragmentStore::TReadSeqStore          TReads;
    typedef typename Value<TReads>::Type                    TRead;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;

    TIterator it = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
    TIterator itEnd = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());

    int avgRQ;
    for (; it != itEnd; ++it)
    {
        TRead const &read = fragmentStore.readSeqStore[(*it).readId];
        avgRQ = 0;
        for(unsigned i = 0; i < length(read); ++i)
            avgRQ += (int) getQualityValue(read[i]);
        // watch out, this is new: use mapping quality if given
        //if((fragmentStore.alignQualityStore[(*it).id]).score == 0 || (char)(avgRQ/length(read))<(fragmentStore.alignQualityStore[(*it).id]).score)
        if((char)(avgRQ/length(read))<(fragmentStore.alignQualityStore[(*it).id]).score)

            (fragmentStore.alignQualityStore[(*it).id]).score = (char)(avgRQ/length(read));
    }
}

// checks whether an alignment has indels
template<typename TValue, typename TSpec>
bool alignmentHasIndel(Align<TValue,TSpec> &align)
{
    typedef Align<TValue,TSpec>  TAlign;
    typedef typename Row<TAlign>::Type      TAlignRow;
    typedef typename Iterator<TAlignRow>::Type  TAlignIterator;

    bool hasIndel = false;
    for(unsigned i = 0; !hasIndel && i < length(rows(align)); ++i)
    {
        TAlignIterator it = begin(row(align,i));
        TAlignIterator itEnd = end(row(align,i));
        while (it != itEnd)
        {
            if(isGap(it))
            {
                hasIndel = true;
                break;
            }
            ++it;
        }
    }
    return hasIndel;
}

// perform read clipping
template<typename TFragmentStore, typename TReadClips, typename TSize, typename TOptions>
void
clipReads(TFragmentStore    &fragmentStore,
          TReadClips    &readClips,
          TSize     arrayBeginPos,
          TSize     arrayEndPos,
          TOptions  &options)
{
    typedef typename TFragmentStore::TAlignedReadStore      TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;     // TMatchQualities
    typedef typename Value<TAlignQualityStore>::Type        TAlignQuality;
    typedef typename TFragmentStore::TReadSeqStore          TReads;
    typedef typename Value<TReads>::Type                TRead;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;             // TGenome
    typedef typename Position<TContigSeq>::Type         TContigPos;

    TIterator it = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
    TIterator itEnd = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
    Align<TRead, ArrayGaps> align;
    resize(rows(align), 2);
#ifdef SNPSTORE_DEBUG
    bool extraV = true;
#endif

    Score<int> scoreType = Score<int>(0, -999, -1001, -1000);
    if(length(readClips) < (arrayEndPos-arrayBeginPos)) ::std::cout << length(readClips) << " readclips but " << (arrayEndPos-arrayBeginPos) << " many reads.\n";
    TContigSeq &currGenome = fragmentStore.contigStore[0].seq;
    int kickout = 0;
    for (; it != itEnd; ++it)
    {
        TMatch &m = *it;
        TRead &read = fragmentStore.readSeqStore[m.readId];
        int clipLeft = readClips[m.readId].i1;
        int clipRight = readClips[m.readId].i2;
        TContigPos beginPos = (m.beginPos < m.endPos ) ? m.beginPos : m.endPos;
        TContigPos endPos = (m.beginPos > m.endPos ) ? m.beginPos : m.endPos;
        TAlignQuality &aliQ = fragmentStore.alignQualityStore[m.id];

#ifdef SNPSTORE_DEBUG
        TContigPos beginBefore = beginPos;
#endif
        if(m.id != m.readId) ::std::cout << "match id != readId \n";
        if(clipLeft+clipRight > (int)length(read) || clipLeft > (int)length(read) || clipRight > (int)length(read))
        {
            if(options._debugLevel > 0)::std::cout << "WARNING: clipLeft+clipRight > readLen!!\n";
#ifdef SNPSTORE_DEBUG
            ::std::cout << "readlength = "<<length(read)<< " \n";
            ::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
            ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
            ::std::cout << "read=" << read << std::endl;
            ::std::cout << "beginPos=" << beginPos << std::endl;
            ::std::cout << "endPos=" << endPos << std::endl;
#endif
            clipLeft = length(read);
            clipRight = 0;
        }
#ifdef SNPSTORE_DEBUG
            ::std::cout << "readlength = "<<length(read)<< " \n";
            ::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
            ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
            ::std::cout << "read=" << read << std::endl;
            ::std::cout << "beginPos=" << beginPos << std::endl;
            ::std::cout << "endPos=" << endPos << std::endl;
#endif
        if(clipLeft > 0 || clipRight > 0)
        {
            //  if(extraV) ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << std::endl;
            if((int)length(read)-clipLeft-clipRight < options.minClippedLength)
            {
                if(options._debugLevel > 1 )
                    ::std::cout << "Discarded: "<<read<<" at position "<< beginPos <<"\n";
                m.endPos = m.beginPos = 0; // "mask" read
                read = "";
                ++kickout;
                continue;
            }
            // adjust read sequence
            read = infix(read,clipLeft,length(read)-clipRight);

            // upate avg read quality
            int avgRQ = 0;
            for(unsigned i = 0; i < length(read); ++i)
                avgRQ += (int) getQualityValue(read[i]);
            aliQ.score = (char)(avgRQ/length(read));
            //      if(extraV) ::std::cout << "aliQ.score = " << (int)aliQ.score << ::std::endl;


            //do semi-global alignment
            assignSource(row(align, 0), read);
            assignSource(row(align, 1), infix(currGenome, beginPos, endPos));
            if ((*it).endPos < (*it).beginPos)
                reverseComplement(source(row(align, 1)));

            int score = globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());
            aliQ.errors = (unsigned) round((float)-score/1000);

#ifdef SNPSTORE_DEBUG
            if(extraV) ::std::cout << align << std::endl;
            if(extraV) ::std::cout << "aliQ.errors = " << (int) aliQ.errors << std::endl;
#endif

            // transform first and last read character to genomic positions
            if(aliQ.pairScore == 1)
            {
                unsigned viewPosReadFirst = toViewPosition(row(align, 0), 0);
                unsigned viewPosReadLast  = toViewPosition(row(align, 0), length(read) - 1);
                unsigned genomePosReadFirst = toSourcePosition(row(align, 1), viewPosReadFirst);
                unsigned genomePosReadLast  = toSourcePosition(row(align, 1), viewPosReadLast);
                //              if(isGap(row(align,1),viewPosReadFirst))
                //              {
                //                  std::cout << "bgein gap --> do nothing " << std::endl;
                //
                //              }

                if(isGap(row(align,1),viewPosReadLast))
                {
                    genomePosReadLast--;
                }
#ifdef SNPSTORE_DEBUG
                if(extraV)
                {   ::std::cout << "viewPosReadFirst = " << viewPosReadFirst << ::std::endl;
                    ::std::cout << "viewPosReadLast = " << viewPosReadLast << ::std::endl;
                    ::std::cout << "genomePosReadFirst = " << genomePosReadFirst << ::std::endl;
                }
#endif
                if(m.beginPos < m.endPos) // forward
                {
                    endPos = beginPos + (genomePosReadLast + 1);
                    beginPos += genomePosReadFirst;
                }
                else
                {
                    beginPos = endPos - genomePosReadLast - 1;
                    endPos = endPos - genomePosReadFirst;
                }

                if(!alignmentHasIndel(align)) aliQ.pairScore = 0;
            }
            else
            {
                if(m.beginPos < m.endPos) // forward
                {
                    endPos -= clipRight;
                    beginPos += clipLeft;
                }
                else
                {
                    endPos -= clipLeft;
                    beginPos += clipRight;
                }
            }

            // set genomic positions
            if(m.beginPos < m.endPos) // forward
            {
                m.endPos = endPos;
                m.beginPos = beginPos;
            }
            else
            {
                m.endPos = beginPos;
                m.beginPos = endPos;
            }
            if(beginPos > 300000000 || endPos > 300000000)
            {
                ::std::cout << "bgein groesser 300mio neu beginPos = "<< beginPos << " endpos=" << endPos << ::std::endl;
#ifdef SNPSTORE_DEBUG
                ::std::cout << "WARNING: clipLeft+clipRight > readLen!!??\n";
                ::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
                ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
                ::std::cout << "read=" << read << std::endl;
                ::std::cout << "beginPos before=" << beginBefore << std::endl;
                ::std::cout << "beginPos=" << beginPos << std::endl;
                ::std::cout << "endPos=" << endPos << std::endl;
#endif
            }


        }
    }
    if(options._debugLevel > 0)
        ::std::cout << kickout <<" reads too short after clipping, discarding!\n";

}

template<typename TTmpReads, typename TTmpMatches, typename TTmpQualities, typename TStr>
void
_dumpMatches(TTmpReads & reads, TTmpMatches & matches, TTmpQualities & qualities, TStr str)
{

    std::cout << "Length of matches = " << length(matches)  << "\n";
    std::cout << "Length of reads   = " << length(reads)  << "\n";
    std::cout << "Length of matchqs = " << length(qualities)  << "\n";

    for(unsigned i = 0 ; i < length(matches); ++i)
    {
        char ori = (matches[i].beginPos < matches[i].endPos) ? 'F' : 'R';
        std::cout << "--"<<str<<"Match number " << i << ":\n";
        std::cout << "--"<<str<<"MatchId  = " << matches[i].id << "\n";
        std::cout << "--"<<str<<"ReadId   = " << matches[i].readId << "\n";
        std::cout << "--"<<str<<"ContigId = " << matches[i].contigId << std::flush << "\n";
        std::cout << "--"<<str<<"gBegin   = " << _min(matches[i].beginPos, matches[i].endPos) << "\n";
        std::cout << "--"<<str<<"gEnd     = " << _max(matches[i].beginPos, matches[i].endPos) << "\n";
        std::cout << "--"<<str<<"orient   = " << ori << std::flush << std::endl;
        if(length(qualities) > matches[i].id)
        {
            std::cout << "--"<<str<<"EditDist = " << (int) qualities[matches[i].id].errors << "\n";
            std::cout << "--"<<str<<"AvgQ     = " << (int) qualities[matches[i].id].score << "\n";
        }
        std::cout << "--"<<str<<"Readseq  = " << reads[matches[i].readId] << std::flush << "\n";

    }
}

template<typename TContigIntervals, typename TFragmentStore, typename TMethOptions>
inline void
assignIntervalsToContigs(TContigIntervals &contigIntervals, TFragmentStore &fragStore, TMethOptions &methOptions)
{
    typedef typename Value<TContigIntervals>::Type      TIntervals;
    typedef typename Value<TIntervals>::Type            TInterval;

    resize(contigIntervals, length(fragStore.contigStore));
    unsigned i = 0;
    for (unsigned contigId = 0; contigId < length(fragStore.contigStore); ++contigId)
    {
        while (i < length(methOptions.intervals) && methOptions.intervals[i].contigName == fragStore.contigNameStore[contigId])
        {
            TInterval interval;
            interval.i1 = methOptions.intervals[i].startPos;
            interval.i2 = methOptions.intervals[i].endPos;
            appendValue(contigIntervals[contigId], interval);
            ++i;
        }
    }
}

template <typename TSpec, typename TContigId, typename TBamFileIns, typename TRecords, typename TContigIntervals, typename TOptions, typename TMethOptions>
inline bool
detectSNPsForContig(seqan::VcfFileOut & vcfFileOut,
                    seqan::BedFileOut & bedFileOut,
                    FragmentStore<TSpec> &fragmentStore1,
                    TContigId &currContigId,
                    TBamFileIns & bamFileIns,
                    TRecords &records,
                    TContigIntervals &contigIntervals,
                    TOptions &options,
                    TMethOptions &methOptions)
{
    typedef          FragmentStore<TSpec>               TFragmentStore;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TAlignedReadStore  TAlignedReadStore;      // TMatches
    typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;     // TMatchQualities
    typedef typename TFragmentStore::TReadStore         TReadStore;             // TReadSet
    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;              // TReadSet
    typedef typename TFragmentStore::TContigStore       TContigStore;           // TGenomeSet
    typedef typename Value<TContigStore>::Type          TContig;

    typedef String<String<typename TFragmentStore::TContigGapAnchor> >      TSetContigAnchorGaps;

    if (!empty(methOptions.intervals) && empty(contigIntervals[currContigId])) return 0;
    unsigned currInterval = 0;

    // parse matches batch by batch
    TContigPos currentWindowBegin = 0;
    if(options._debugLevel > 0) ::std::cout << "Scanning genome #" << currContigId << " ..." << ::std::endl;

    // containers for those matches that overlap two windows
    TAlignedReadStore tmpMatches;
    TAlignQualityStore tmpQualities;
    TReadStore tmpRs;
    TReadSeqStore tmpReads;     // Something went wrong when keeping all reads of contig, so keep it like that for the moment
    TSetContigAnchorGaps tmpSetContigAnchorGaps;
    options.minCoord = std::numeric_limits<unsigned>::max();
    options.maxCoord = 0;

    // snp calling is done for all positions between windowBegin and windowEnd
    // bs_change
    if (!empty(methOptions.intervals))
    {
        SEQAN_ASSERT_LT((unsigned)contigIntervals[currContigId][currInterval].i2, (unsigned)length(fragmentStore1.contigStore[currContigId].seq));
        currentWindowBegin = contigIntervals[currContigId][currInterval].i1;
    }
    while(currentWindowBegin <  (TContigPos)length(fragmentStore1.contigStore[currContigId].seq))
    {
        TContigPos currentWindowEnd = currentWindowBegin + options.windowSize;
        if(currentWindowEnd > (TContigPos)length(fragmentStore1.contigStore[currContigId].seq)) currentWindowEnd = (TContigPos)length(fragmentStore1.contigStore[currContigId].seq);
        // bs_change
        if (!empty(methOptions.intervals) && currentWindowEnd > contigIntervals[currContigId][currInterval].i2)
        {
            currentWindowEnd = contigIntervals[currContigId][currInterval].i2;
            //std::cout << "currInterval startPos: " << contigIntervals[currContigId][currInterval].i1 << "  endPos: " << contigIntervals[currContigId][currInterval].i2 << std::endl;
        }
        //std::cout << "currentWindowBegin: " << currentWindowBegin << "  currentWindowEnd: " << currentWindowEnd << std::endl;

        if(options._debugLevel > 0)
            std::cout << "Sequence number " << currContigId << " window " << currentWindowBegin << ".." << currentWindowEnd << "\n";

        TFragmentStore fragmentStoreTmp;  // Use temp. fragmentStore to store only current reads and to use only infix of contig seq for realigning etc.
        TSetContigAnchorGaps setContigAnchorGaps;

        // add the matches that were overlapping with this and the last window (copied in order to avoid 2 x makeGlobal)
        if(!empty(tmpMatches))
        {
            //sumreads -=  length(tmpReads);  // count these reads only once
            resize(fragmentStoreTmp.alignQualityStore,length(tmpMatches));
            resize(fragmentStoreTmp.alignedReadStore,length(tmpMatches));
            resize(fragmentStoreTmp.readSeqStore,length(tmpMatches));
            resize(fragmentStoreTmp.readStore,length(tmpMatches));
            resize(setContigAnchorGaps, length(tmpMatches));

            arrayMoveForward(begin(tmpQualities,Standard()), end(tmpQualities,Standard()), begin(fragmentStoreTmp.alignQualityStore,Standard()));
            arrayMoveForward(begin(tmpMatches,Standard()), end(tmpMatches,Standard()), begin(fragmentStoreTmp.alignedReadStore,Standard()));
            arrayMoveForward(begin(tmpReads,Standard()), end(tmpReads,Standard()), begin(fragmentStoreTmp.readSeqStore,Standard()));
            arrayMoveForward(begin(tmpRs,Standard()), end(tmpRs,Standard()), begin(fragmentStoreTmp.readStore,Standard()));
            arrayMoveForward(begin(tmpSetContigAnchorGaps,Standard()), end(tmpSetContigAnchorGaps,Standard()), begin(setContigAnchorGaps,Standard()));
#ifdef SNPSTORE_DEBUG
            CharString strstr = "afterCopyInFrag";
            _dumpMatches(fragmentStoreTmp,strstr);
#endif
        }
        TContig contig;
        contig.seq = fragmentStore1.contigStore[currContigId].seq;
        appendValue(fragmentStoreTmp.contigStore, contig);
        appendValue(fragmentStoreTmp.contigNameStore, fragmentStore1.contigNameStore[currContigId]);

        // parse matches for current window
        if(options._debugLevel > 0) std::cout << "Parsing reads up to position " << currentWindowEnd << "...\n";
        for(unsigned j = 0; j < length(options.readFNames); ++j)
        {
            unsigned sizeBefore = length(fragmentStoreTmp.alignedReadStore);

            int result = readMatchesFromSamBam(setContigAnchorGaps, *bamFileIns[j], records[j], fragmentStoreTmp, fragmentStore1,
                                               currContigId, currentWindowBegin, currentWindowEnd, options);

            if(result == CALLSNPS_GFF_FAILED)
            {
                std::cerr << "Failed to open read file " << options.readFNames[j] << "or reads are not sorted correctly. " << ::std::endl;
                return result;
            }
            if(result > 0) return result;
            if(options._debugLevel > 0) std::cout << "parsed reads of file " << j << "\n";

            // store average quality of each read if mapq is 0 or if average quality is < mapq
            // addReadQualityToMatches(fragmentStoreTmp,sizeBefore,(unsigned)length(fragmentStoreTmp.alignedReadStore),options); // bs_change: we will deal with real mapqs
            // do pile up correction if lane-based. lane-specific pileup correction not really used
            if(options.maxPile != 0 && options.laneSpecificMaxPile)
                applyPileupCorrection(fragmentStoreTmp,sizeBefore,(unsigned)length(fragmentStoreTmp.alignedReadStore),options);
        }
        SEQAN_ASSERT_EQ(length(fragmentStoreTmp.alignedReadStore), length(setContigAnchorGaps));

        if (options._debugLevel > 1)  // number of reads currently in memory
            std::cerr << lengthSum(fragmentStoreTmp.readSeqStore) << " bps of " << length(fragmentStoreTmp.readSeqStore) << " reads in memory." << ::std::endl;
        //sumreads +=  length(fragmentStoreTmp.readSeqStore);  // total count of reads

        // do merged pile up correction
        if(options.maxPile != 0 && !options.laneSpecificMaxPile)
            applyPileupCorrection(fragmentStoreTmp,(unsigned)0,(unsigned)length(fragmentStoreTmp.alignedReadStore),options);

        // these were set while parsing matches, first and last position of parsed matches
        TContigPos startCoord = _max((int)options.minCoord-options.realignAddBorder,0);// can be < currentWindowBegin
        TContigPos endCoord = _min(options.maxCoord+options.realignAddBorder,length(fragmentStore1.contigStore[currContigId].seq)); // can be > currentWindoEnd

        if(!empty(fragmentStoreTmp.alignedReadStore))
        {
            //initial values of min and max coords for next round are set here
            if(currentWindowEnd != (TContigPos)length(fragmentStore1.contigStore[currContigId].seq))
            {
                clear(tmpMatches);
                clear(tmpQualities);
                clear(tmpRs);
                clear(tmpReads);
                clear(tmpSetContigAnchorGaps);
                copyNextWindowMatchesAndReads(fragmentStoreTmp, setContigAnchorGaps, tmpReads, tmpRs, tmpMatches, tmpQualities, tmpSetContigAnchorGaps, currContigId, currentWindowEnd, options);
            }
#ifdef SNPSTORE_DEBUG
            std::cout << "Min = " << options.minCoord << " Max = " << options.maxCoord << std::endl;
            std::cout << "startCoord = " << startCoord << " endCoord = " << endCoord << std::endl;
            std::cout << "currentWindowBegin = " << currentWindowBegin << " currentWindowEnd = " << currentWindowEnd << std::endl;
#endif
            // Aligned read coordinates are relative to current chromosomal window (segment)
            transformCoordinates(fragmentStoreTmp,startCoord,options);
            // set the current contig segment as contig sequence
            fragmentStoreTmp.contigStore[0].seq = infix(fragmentStore1.contigStore[currContigId].seq,startCoord,endCoord);

            if(options.realign)
            {
                doCheckRealignCall(fragmentStoreTmp, startCoord, currentWindowBegin, currentWindowEnd, setContigAnchorGaps, vcfFileOut, bedFileOut, methOptions, options);
            }
            else
            {
                // Convert gaps of pairwise matches to msa
                //std::cout << "convertPairWiseToGlobalAlignment..." << std::endl;
#ifdef CALL_PROFILE
                double timeStamp = sysTime();
#endif
                convertPairWiseToGlobalAlignment(fragmentStoreTmp, setContigAnchorGaps);
#ifdef CALL_PROFILE
                Times::instance().time_convertPWToGlobal += (sysTime() - timeStamp);
#endif
                doSnpAndMethCalling(fragmentStoreTmp,  startCoord, currentWindowBegin, currentWindowEnd, false, vcfFileOut, bedFileOut, methOptions, options);    //bs
            }
        }
        if (!empty(methOptions.intervals) && currentWindowEnd == contigIntervals[currContigId][currInterval].i2)   // Reached end off current interval
        {
            if (currInterval < length(contigIntervals[currContigId])-1)   // If existent, go to next interval within the current contig
            {
                currentWindowBegin = contigIntervals[currContigId][currInterval+1].i1;
                clear(tmpReads);
                clear(tmpRs);
                clear(tmpMatches);
                clear(tmpQualities);
                clear(tmpSetContigAnchorGaps);
            }
            else currentWindowBegin = length(fragmentStore1.contigStore[currContigId].seq);    // No interval int this contig to analyze anymore -> jump to end of contig
            ++currInterval;
        }
        else currentWindowBegin = currentWindowEnd;
    }
    return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TMethOptions>
int detectSNPs(SNPCallingOptions &options, TMethOptions &methOptions)
{
    typedef FragmentStore<SnpStoreSpec_>                                        TFragmentStore;
    typedef typename TFragmentStore::TContigPos                                 TContigPos;

    SEQAN_PROTIMESTART(load_time);

    // Load genomes in FragmentStore
    TFragmentStore      fragmentStore1;
    StringSet<CharString> genomeFileNameList; // possible for later to load multiple files
    appendValue(genomeFileNameList,options.genomeFName);
    loadContigs(fragmentStore1, genomeFileNameList);
    refresh(fragmentStore1.contigNameStoreCache);


    // Assign intervals to analyze to string which can be accessed by contigId
    String<String<Interval<TContigPos> > > contigIntervals;
    if(!empty(methOptions.intervals)) assignIntervalsToContigs(contigIntervals, fragmentStore1, methOptions);
    // Prepare genotype priors
    computeGenotypePriors(methOptions, options);
    // Store fileName for temp file for each contig, needed later to merge in correct order
    String<CharString> contigTempFileNamesVcf;
    String<CharString> contigTempFileNamesBed;
    resize(contigTempFileNamesVcf, length(fragmentStore1.contigStore));
    resize(contigTempFileNamesBed, length(fragmentStore1.contigStore));

    // TODO(holtgrew): The new's below are not paired with delete's, thus there are leaks!

#if !defined(SEQAN_ENABLE_PARALLELISM)
    // For each contig we need own record readers (for each given sam file)
    std::vector<seqan::BamFileIn *> bamFileIns(length(options.readFNames));
    std::vector<seqan::BamAlignmentRecord> records(bamFileIns.size());
    for (unsigned i = 0; i < length(options.readFNames); ++i)
    {
        bamFileIns[i] = new seqan::BamFileIn;
        if (!open(*bamFileIns[i], toCString(options.readFNames[i])))
        {
            ::std::cerr << "Failed to open read file " << options.readFNames[i] << ::std::endl;
            return CALLSNPS_GFF_FAILED;
        }

        // Read header.
        seqan::BamHeader header;
        readHeader(header, *bamFileIns[i]);
    }
#endif  // #if !defined(SEQAN_ENABLE_PARALLELISM)
    bool abort = false;
#if defined(SEQAN_ENABLE_PARALLELISM)
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) // TODO Check if guided is faster
#endif  // #if defined(SEQAN_ENABLE_PARALLELISM)

    for (int currContigId = 0; currContigId < (int)length(fragmentStore1.contigStore); ++currContigId)
    {
#if defined(SEQAN_ENABLE_PARALLELISM)
        // For each contig we need our own BamFileIn (for each given sam file)
        std::vector<seqan::BamFileIn *> bamFileIns(length(options.readFNames));
        std::vector<seqan::BamAlignmentRecord> records(bamFileIns.size());

        for (unsigned i = 0; i < length(options.readFNames); ++i)
        {
            bamFileIns[i] = new seqan::BamFileIn;
            if (!open(*bamFileIns[i], toCString(options.readFNames[i])))
            {
                ::std::cerr << "Failed to open read file " << options.readFNames[i] << ::std::endl;
                abort = true;
            }
            else
            {
                // Read header.
                seqan::BamHeader header;
                readHeader(header, *bamFileIns[i]);
            }
        }
#endif  // #if defined(SEQAN_ENABLE_PARALLELISM)

        if (!abort)
        {
        // Temp. contig output
            CharString tempFileNameVcf;
            CharString tempFileNameBed;
            append(tempFileNameVcf, SEQAN_TEMP_FILENAME());
            append(tempFileNameBed, SEQAN_TEMP_FILENAME());
            //std::cout << "temp file Name: " << tempFileNameVcf << std::endl;
            stringstream ss;
            ss << currContigId;
            append(tempFileNameVcf, ss.str());
            append(tempFileNameBed, ss.str());
            append(tempFileNameVcf, ".vcf");
            append(tempFileNameBed, ".bed");
            contigTempFileNamesVcf[currContigId] = tempFileNameVcf;
            contigTempFileNamesBed[currContigId] = tempFileNameBed;

            seqan::VcfFileOut vcfFileOut(toCString(tempFileNameVcf));
            seqan::VcfHeader vcfHeader;

            appendName(contigNamesCache(context(vcfFileOut)), fragmentStore1.contigNameStore[currContigId]);
            appendName(sampleNamesCache(context(vcfFileOut)), "NA00001");

            appendValue(vcfHeader, VcfHeaderRecord("fileformat", "VCFv4.1"));
            appendValue(vcfHeader, VcfHeaderRecord("fileDate", "20090805"));
            appendValue(vcfHeader, VcfHeaderRecord("source", "temporary file snp_meth_store"));
            appendValue(vcfHeader, VcfHeaderRecord("reference", "file:///seq/references/1000GenomesPilot-NCBI36.fasta"));
            appendValue(vcfHeader, VcfHeaderRecord("contig", "<ID=21,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"));
            appendValue(vcfHeader, VcfHeaderRecord("phasing", "partial"));
            appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"));
            appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"));
            appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"));
            appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"));
            appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"));
            appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"));
            appendValue(vcfHeader, VcfHeaderRecord("FILTER", "<ID=q10,Description=\"Quality below 10\">"));
            appendValue(vcfHeader, VcfHeaderRecord("FILTER", "<ID=s50,Description=\"Less than 50% of samples have data\">"));
            appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
            appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
            appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"));
            appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"));
            writeHeader(vcfFileOut, vcfHeader);

            seqan::BedFileOut bedFileOut(toCString(tempFileNameBed));
            //XXX addSequenceName(tempBedStream, fragmentStore1.contigNameStore[currContigId]);
            detectSNPsForContig(vcfFileOut, bedFileOut, fragmentStore1, currContigId, bamFileIns, records, contigIntervals, options, methOptions);
        }
    }

    if (abort)
        return CALLSNPS_GFF_FAILED;

    // Prepare VCF output
    seqan::VcfFileOut vcfFileOut(toCString(options.vcfOut));
    seqan::VcfHeader vcfHeader;

    appendName(sampleNamesCache(context(vcfFileOut)), "NA00001");
    appendValue(vcfHeader, VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(vcfHeader, VcfHeaderRecord("fileDate", "20090805"));
    appendValue(vcfHeader, VcfHeaderRecord("source", "snp_meth_store"));
    appendValue(vcfHeader, VcfHeaderRecord("reference", "file:///seq/references/1000GenomesPilot-NCBI36.fasta"));
    appendValue(vcfHeader, VcfHeaderRecord("contig", "<ID=21,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"));
    appendValue(vcfHeader, VcfHeaderRecord("phasing", "partial"));
    appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"));
    appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"));
    appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"));
    appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"));
    appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"));
    appendValue(vcfHeader, VcfHeaderRecord("INFO", "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"));
    appendValue(vcfHeader, VcfHeaderRecord("FILTER", "<ID=q10,Description=\"Quality below 10\">"));
    appendValue(vcfHeader, VcfHeaderRecord("FILTER", "<ID=s50,Description=\"Less than 50% of samples have data\">"));
    appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
    appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"));
    appendValue(vcfHeader, VcfHeaderRecord("ID", "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"));
    writeHeader(vcfFileOut, vcfHeader);

    // Prepare BED output
    seqan::BedFileOut bedFileOut;
    if (!open(bedFileOut, toCString(options.bedOut)))
    {
        std::cerr << "Could not open output BED file.\n";
        return 1;
    }

   // Append content of temp files to final output in contig order
    for (unsigned i = 0; i < length(fragmentStore1.contigStore); ++i)
    {
        // VCF
        seqan::VcfFileIn vcfFileIn;
        if (!open(vcfFileIn, toCString(contigTempFileNamesVcf[i])))
        {
            std::cerr << "ERROR: Could not open temporary vcf file: " << contigTempFileNamesVcf[i] << "\n";
            return 1;
        }

        seqan::VcfHeader vcfHeader;
        readHeader(vcfHeader, vcfFileIn);

        seqan::VcfRecord vcfRecord;
        while (!atEnd(vcfFileIn))
        {
            readRecord(vcfRecord, vcfFileIn);
            if (vcfRecord.rID >= (int)length(contigNames(context(vcfFileOut))))
                appendName(contigNamesCache(context(vcfFileOut)), fragmentStore1.contigNameStore[i]);
            writeRecord(vcfFileOut, vcfRecord);
        }
        ClassTest::_deleteTempFile(ClassTest::_stripFileName(toCString(contigTempFileNamesVcf[i])));
        remove(toCString(contigTempFileNamesVcf[i]));  // Delete temp. files

        // BED
        seqan::BedFileIn bedFileIn;
        if (!open(bedFileIn, toCString(contigTempFileNamesBed[i])))
        {
            std::cerr << "ERROR: Could not open temporary bed file: " << contigTempFileNamesBed[i] << "\n";
            return 1;
        }

        seqan::BedRecord<seqan::Bed6> bedRecord;
        while (!atEnd(bedFileIn))
        {
            readRecord(bedRecord, bedFileIn);
            writeRecord(bedFileOut, bedRecord);
        }
        ClassTest::_deleteTempFile(ClassTest::_stripFileName(toCString(contigTempFileNamesBed[i])));
        remove(toCString(contigTempFileNamesBed[i]));
    }

    methOptions.statsCGMethylated = methOptions.statsCGMethylated/methOptions.countCG;
    methOptions.statsCHGMethylated = methOptions.statsCHGMethylated/methOptions.countCHG;
    methOptions.statsCHHMethylated = methOptions.statsCHHMethylated/methOptions.countCHH;
    std::cout << "Average methylation rate in CG context: " << methOptions.statsCGMethylated << std::endl;
    std::cout << "Average methylation rate in CHG context: " << methOptions.statsCHGMethylated << std::endl;
    std::cout << "Average methylation rate in CHH context: " << methOptions.statsCHHMethylated << std::endl;
    return 0;
}


template <typename TMethOptions>
ArgumentParser::ParseResult
parseCommandLine(SNPCallingOptions & options, TMethOptions &methOptions, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("casbar");
    // Set short description, version, and date.
    setShortDescription(parser, "SNP and methylation level calling");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "BS-Seq Analysis");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIALIGNMENT FILE\\fP> -o <\\fISNP FILE\\fP> -b <\\fIMETH-LEVEL FILE\\fP>");
    addDescription(parser, "SNP and methylation level calling in mapped bisulfite read data.");

    // We require two arguments.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "GENOME"));
    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A reference genome file.");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "ALIGNMENTS"));
    setValidValues(parser, 1, BamFileIn::getFileExtensions());
    setHelpText(parser, 1, "SAM input file containing four-letter read alignments (must be sorted by coordinates).");

    addSection(parser, "Options");
    addOption(parser, ArgParseOption("o", "output", "Output file for SNPs.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "output", VcfFileOut::getFileExtensions());
    setRequired(parser, "output", true);
    addOption(parser, ArgParseOption("b", "bed", "Bed output file for methylation level calls.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "bed", BedFileOut::getFileExtensions());
    setRequired(parser, "bed", true);

    addOption(parser, ArgParseOption("mu", "multi", "Keep non-unique reads."));
    addOption(parser, ArgParseOption("sqo", "solexa-qual-offset", "Base qualities are encoded as char value - 64 (instead of char - 33)."));
    hideOption(parser, "sqo");
    addOption(parser, ArgParseOption("mp", "max-pile", "Maximal number of matches allowed to pile up at the same genome position.", ArgParseArgument::INTEGER));
    setMinValue(parser, "max-pile", "0");
    setDefaultValue(parser, "max-pile", options.maxPile);
    addOption(parser, ArgParseOption("mmp", "merged-max-pile", "Do pile up correction on merged lanes."));
    addOption(parser, ArgParseOption("mc", "min-coverage", "Minimal required number of reads covering a candidate position.", ArgParseArgument::INTEGER));
    setMinValue(parser, "min-coverage", "1");
    setDefaultValue(parser, "min-coverage", options.minCoverage);
    addOption(parser, ArgParseOption("eb", "exclude-border", "Exclude read positions within eb base pairs of read borders for SNP calling.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "exclude-border", options.excludeBorderPos);
    addOption(parser, ArgParseOption("su", "suboptimal", "Keep suboptimal reads."));
    addOption(parser, ArgParseOption("re", "realign", "Realign reads around indel candidates."));
    hideOption(parser, "re");
    addOption(parser, ArgParseOption("cq", "corrected-quality", "New quality calibration factor.", ArgParseArgument::DOUBLE));
    hideOption(parser, "cq");
    addOption(parser, ArgParseOption("pws", "parse-window-size", "Genomic window size for parsing reads (concerns memory consumption, choose smaller windows for higher coverage).", ArgParseArgument::INTEGER));
    setMinValue(parser, "parse-window-size", "1");
    setMaxValue(parser, "parse-window-size", "100000");
    setDefaultValue(parser, "parse-window-size", options.windowSize);
    addOption(parser, ArgParseOption("reb", "realign-border", "Realign border.", ArgParseArgument::INTEGER));
    setMinValue(parser, "realign-border", "0");
    setMaxValue(parser, "realign-border", "10");
    setDefaultValue(parser, "realign-border", options.realignAddBorder);
    hideOption(parser, "reb");
    addOption(parser, ArgParseOption("it", "indel-threshold", "Minimal number of indel-supporting reads required for realignment.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "indel-threshold", options.indelCountThreshold);
    hideOption(parser, "it");
    addOption(parser, ArgParseOption("I", "intervals", "Genomic intervals to analyze. E.g. 21:1000-2000.",  ArgParseArgument::STRING, "TEXT"));

    addSection(parser, "Calling options");
    addOption(parser, ArgParseOption("bcr", "bs-conv-rate", "Bisulfite conversion rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "bs-conv-rate", "0");
    setMaxValue(parser, "bs-conv-rate", "1");
    setDefaultValue(parser, "bs-conv-rate", methOptions.convRate);
    addOption(parser, ArgParseOption("mm", "min-mutations", "Minimal number of deviating bases for calling.", ArgParseArgument::INTEGER));
    setMinValue(parser, "min-mutations", "1");
    setDefaultValue(parser, "min-mutations", options.minMutT);
    addOption(parser, ArgParseOption("mq", "min-quality", "Minimal average quality for calling.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "min-quality", "0");
    setDefaultValue(parser, "min-quality", options.avgQualT);
    addOption(parser, ArgParseOption("mmq", "min-map-quality", "Minimum base call quality for a match to be considered.", ArgParseArgument::INTEGER));
    setMinValue(parser, "min-map-quality", "0");
    setDefaultValue(parser, "min-map-quality", options.minMapQual);
    addOption(parser, ArgParseOption("hes", "prob-hetero-snp", "Heterozygous SNP probability to compute genotype prior probabilities.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "prob-hetero-snp", "0");
    setMaxValue(parser, "prob-hetero-snp", "1");
    setDefaultValue(parser, "prob-hetero-snp", options.pHetSnp);
    addOption(parser, ArgParseOption("hos", "prob-homo-snp", "Homozygous SNP probability to compute genotype prior probabilities.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "prob-homo-snp", "0");
    setMaxValue(parser, "prob-homo-snp", "1");
    setDefaultValue(parser, "prob-homo-snp", options.pHomoSnp);
    addOption(parser, ArgParseOption("spl", "beta-sampling", "Sample beta values (instead of using MLE and newton method)."));
    hideOption(parser, "spl");
    addOption(parser, ArgParseOption("msc", "min-score", "Minimum score to call.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "min-score", methOptions.minScoreToCallSnp);
    addOption(parser, ArgParseOption("mpc", "min-prob", "Minimum genotype probability to call.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "min-prob", "0");
    setMaxValue(parser, "min-prob", "1");
    setDefaultValue(parser, "min-prob", methOptions.minProbToCallSnp);
    addOption(parser, ArgParseOption("umq", "use-mapq", "Use mapqs as weights for SNP/BS calling."));
    hideOption(parser, "umq");
    addOption(parser, ArgParseOption("gp", "genotype-priors", "Use non-uniform genotype prior probabilities."));
    hideOption(parser, "gp");   // simulate data and check
    addOption(parser, ArgParseOption("nec", "ns-errors-calling", "Use empirical error frequencies of Illumina sequencing data to compute likelihoods in bayesian model (corresponding to Dohm et al. 2008)."));


    // Realignment
    addOption(parser, ArgParseOption("nse", "ns-subst-errors", "Use non-uniform substitution error frequencies for realigning."));
    hideOption(parser, "nse");
    addOption(parser, ArgParseOption("nsi", "ns-ins-errors", "Use non-uniform insertion error frequencies for realigning."));
    hideOption(parser, "nsi");
    addOption(parser, ArgParseOption("nsd", "ns-del-errors", "Use non-uniform deletion error frequencies for realigning."));
    hideOption(parser, "nsd");
    addOption(parser, ArgParseOption("dr", "del-rate", "Genomic deletion rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "del-rate", "0");
    setMaxValue(parser, "del-rate", "1");
    hideOption(parser, "dr");
    addOption(parser, ArgParseOption("der", "del-error-rate", "Deletion error rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "del-error-rate", "0");
    setMaxValue(parser, "del-error-rate", "1");
    hideOption(parser, "der");
    addOption(parser, ArgParseOption("ier", "ins-error-rate", "Insertion error rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "ins-error-rate", "0");
    setMaxValue(parser, "ins-error-rate", "1");
    hideOption(parser, "ier");
    addOption(parser, ArgParseOption("egs", "end-gap-score", "Simple score for end gaps (must be in the range of internal gaps) to avoid introducing of gaps.", ArgParseArgument::DOUBLE));
    hideOption(parser, "egs");
    addOption(parser, ArgParseOption("sl", "score-limit", "Score limit to avoid to high negative scores in profile realignment.", ArgParseArgument::DOUBLE));
    hideOption(parser, "sl");

    addSection(parser, "Other options");
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBcasbar\\fP \\fB-nec\\fP \\fB-o\\fP \\fBsnps.vcf\\fP \\fB-b\\fP \\fBmeth_levels.bed\\fP \\fBgenome.fa\\fP \\fBmapped_bisulfite_reads.sam\\fP",
                "SNP and methylation level calling by taking nonuniform sequencing error probabilities into account.");
    addListItem(parser, "\\fBcasbar\\fP \\fB-nec\\fP \\fB-mc\\fP \\fB6 \\fB-msc\\fP \\fB5\\fP \\fB-o\\fP \\fBsnps.vcf\\fP \\fB-b\\fP \\fBmeth_levels.bed\\fP \\fBgenome.fa\\fP \\fBmapped_bisulfite_reads.sam\\fP",
                "SNP and methylation level calling by taking nonuniform sequencing error probabilities into account for genomic positions with a coverage >= 6. The minimum score a SNP or methylation level to be called is 5.");


    // ${CASBAR} -nec -mc 6 -msc 5 -hes 0.005 -o snps_se_0.vcf -b meths_se_0.bed  hg18_chr21_3000.fa reads_se_N6000_2.CT_GA.verified.pos_so.sam

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    // Options:
    getOptionValue(options.vcfOut, parser, "output");
    getOptionValue(options.bedOut, parser, "bed");
    options.keepMultiReads = isSet(parser, "multi");
    if (isSet(parser, "solexa-qual-offset"))
        options.asciiQualOffset = 64;
    getOptionValue(options.maxPile, parser, "max-pile");
    options.laneSpecificMaxPile = !isSet(parser, "merged-max-pile");
    getOptionValue(options.minCoverage, parser, "min-coverage");
    getOptionValue(options.excludeBorderPos, parser, "exclude-border");
    options.keepSuboptimalReads = isSet(parser, "suboptimal");
    options.realign = isSet(parser, "realign");
    if (options.realign)
        options.windowSize = 1000;
    hideOption(parser, "re");   // hide for the beginning !

    getOptionValue(options.windowSize, parser, "parse-window-size");
    getOptionValue(options.realignAddBorder, parser, "realign-border");
    // SNP Calling Options:
    getOptionValue(methOptions.convRate, parser, "bs-conv-rate");
    getOptionValue(options.minMutT, parser, "min-mutations");
    getOptionValue(options.avgQualT, parser, "min-quality");
    getOptionValue(options.minMapQual, parser, "min-map-quality");
    getOptionValue(options.pHetSnp, parser, "prob-hetero-snp");
    getOptionValue(options.pHomoSnp, parser, "prob-homo-snp");

    getOptionValue(options.indelCountThreshold, parser, "indel-threshold");
    // Bs Options:
    clear(methOptions.intervals);
    if (isSet(parser, "intervals"))
    {
        resize(methOptions.intervals, getOptionValueCount(parser, "intervals"), Exact());
        for (unsigned i = 0; i < getOptionValueCount(parser, "intervals"); ++i)
        {
            CharString startPos;
            CharString endPos;
            CharString value;
            getOptionValue(value, parser, "intervals", i);
            unsigned j;
            for (j = 0; value[j] != ':' && j < length(value); ++j)
            {
                appendValue(methOptions.intervals[i].contigName, value[j], Generous());
            }
            ++j;
            for (; value[j] != '-' && j < length(value); ++j)
            {
                appendValue(startPos, value[j], Generous());
            }
            ++j;
            methOptions.intervals[i].startPos = lexicalCast<unsigned>(startPos);
            for (; j < length(value); ++j)
            {
                appendValue(endPos, value[j], Generous());
            }
            methOptions.intervals[i].endPos = lexicalCast<unsigned>(endPos);

            std::cout << " contigName: " << methOptions.intervals[i].contigName << "   startPos: " << methOptions.intervals[i].startPos << "  endPos: " << methOptions.intervals[i].endPos;
        }
    }
    if (isSet(parser, "beta-sampling"))
        methOptions.betaSampling = true;
    else
        methOptions.betaSampling = false;

    getOptionValue(methOptions.minScoreToCallSnp, parser, "min-score");
    getOptionValue(methOptions.minProbToCallSnp, parser, "min-prob");

    if (isSet(parser, "use-mapq")) methOptions.useMapq = true;
    else methOptions.useMapq = false;

    if (isSet(parser, "genotype-priors")) methOptions.uniformGenPriors = false;
    if (isSet(parser, "ns-errors-calling")) methOptions.uniformSeqErrorsCalling = false;
    if (isSet(parser, "ns-subst-errors")) methOptions.nonSimpleSubstErrors = true;
    if (isSet(parser, "ns-ins-errors")) methOptions.nonSimpleInsErrors = true;
    if (isSet(parser, "ns-del-errors")) methOptions.nonSimpleDelErrors = true;

    getOptionValue(methOptions.delRate, parser, "del-rate");
    getOptionValue(methOptions.delErrorRate, parser, "del-error-rate");
    getOptionValue(methOptions.insErrorRate, parser, "ins-error-rate");
    getOptionValue(methOptions.endGapScore, parser, "end-gap-score");
    getOptionValue(methOptions.scoreLimit, parser, "score-limit");

    // Other Options:
    //getOptionValue(options.outputLog, parser, "log-file");
    if (isSet(parser, "verbose"))
        options._debugLevel = max(options._debugLevel, 1);
    if (isSet(parser, "very-verbose"))
        options._debugLevel = max(options._debugLevel, 2);

    getArgumentValue(options.genomeFName, parser, 0);
    unsigned countFiles = getArgumentValueCount(parser, 1);
    resize(options.readFNames, countFiles);
    for (unsigned i = 0; i < countFiles; ++i)
        getArgumentValue(options.readFNames[i], parser, 1, i);

    if(options.windowSize > 1000000)
        options.windowSize = 10000;

    return ArgumentParser::PARSE_OK;
}


int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    SNPCallingOptions      options;
    MethCallingOptions      methOptions;
    ArgumentParser::ParseResult res = parseCommandLine(options, methOptions, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    for(int arg = 0; arg < argc; ++arg) {
        options.programCall << argv[arg] << " ";
    }

    //////////////////////////////////////////////////////////////////////////////
    // check for variants

#ifdef CALL_PROFILE
    double timeStamp = sysTime();
#endif
    int result = detectSNPs(options, methOptions);
    if (result > 0)
    {
        cerr << "ERROR: Something went wrong. Try 'snpStore --help' for more information." << endl << endl;
        return 0;
    }
#ifdef CALL_PROFILE
    Times::instance().time_all = sysTime() - timeStamp;
    std::cout << "  Time needed for all: " << Times::instance().time_all/60.0 << "min" << std::endl;
    std::cout << "  Time needed for doBsCalling: " << Times::instance().time_doBsCalling/60.0 << "min" << std::endl;
    std::cout << "  Time needed for optimization: " << Times::instance().time_optimization/60.0 << "min" << std::endl;
    std::cout << "  Time needed for IO: " << Times::instance().time_IO/60.0 << "min" << std::endl;
    std::cout << "  Time needed for convertPairwiseToGlobal: " << Times::instance().time_convertPWToGlobal/60.0 << "min" << std::endl;
#endif

    if (options._debugLevel > 0)
    {
        std::cout << "CounteBViolated: " << methOptions.counteBViolated << std::endl;
        std::cout << "CountPlanB (Newton, error bound violated or f'' > 0: sampling applied): " << methOptions.countPlanB << std::endl;
        std::cout << "CountNoPlanB: " << methOptions.countNoPlanB << std::endl;
        std::cout << "CountCovTooLow: " << methOptions.countCovTooLow << std::endl;
        std::cout << "CountScoreTooLow: " << methOptions.countScoreTooLow << std::endl;
    }

    return 0;
}

