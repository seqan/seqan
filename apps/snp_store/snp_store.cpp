/*==========================================================================
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


// switches on corrected heterozygote variance in allele ratios (higher sensitivity for hets, at the cost of more false positive snp calls)
#define CORRECTED_HET  // requires boost library


// for longer reads such as 454 data (will use different scoring scheme)
//#define READS_454


///// switch on extreme debug output
// #define SNPSTORE_DEBUG
//#define SNPSTORE_DEBUG_CANDPOS


#include <seqan/platform.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>


#ifdef STDLIB_VS
#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
#define SEQAN_DEFAULT_TMPDIR "./"
#endif


//#include "../../../apps/rep_sep/utils.h"
//#include "../../../apps/rep_sep/assembly_parser.h"
//#include "../../../apps/rep_sep/column_scanner.h"
//#include "../../../apps/rep_sep/rgraph.h"
//#include "../../../apps/rep_sep/rep_sep.h"

#include "snp_store.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <seqan/arg_parse.h>

using namespace std;
using namespace seqan;



// load entire genome into memory
template <typename TGenomeSet, typename TGenomeSetSize, typename TGenomeNames>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList, ::std::map<CharString, TGenomeSetSize> &gIdStringToIdNumMap, TGenomeNames & genomeNames)
{
    unsigned gSeqNo = 0;
    unsigned filecount = 0;
    CharString temp;
    clear(genomeNames);
    while(filecount < length(fileNameList))
    {
        clear(temp);
        SeqFileIn seqFileIn;
        if (!open(seqFileIn, toCString(fileNameList[filecount])))
            return false;

        unsigned seqCount = 0;
        Dna5String seq;
        for (; !atEnd(seqFileIn); ++seqCount)
        {
            readRecord(temp, seq, seqFileIn);
            // Trim sequence ID and insert into mapping.
            cropAfterFirst(temp, NotFunctor<IsGraph>());
            gIdStringToIdNumMap.insert(::std::make_pair(temp, length(genomes))); // keeps the whole fasta ID including white spaces
            appendValue(genomeNames, temp);
            // Register genome sequence.
            appendValue(genomes, seq);
        }

        gSeqNo += seqCount;
        ++filecount;
    }
    resize(genomes, gSeqNo);
    return (gSeqNo > 0);
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
template<typename TFragmentStore, typename TReadCigars, typename TReadCounts, typename TReadClips, typename TSize, typename TContigPos, typename TOptions>
void
copyNextWindowMatchesAndReads(TFragmentStore &fragmentStore,
                              TReadCounts &readCounts,
                              TReadCigars &readCigars,
                              TReadClips  &readClips,
                              TReadCounts &tmpReadCounts,
                              typename TFragmentStore::TReadSeqStore &tmpReads,
                              typename TFragmentStore::TReadStore &tmpRs,
                              typename TFragmentStore::TAlignedReadStore &tmpMatches,
                              typename TFragmentStore::TAlignQualityStore &tmpQualities,
                              TReadClips &tmpReadClips,
                              TReadCigars &tmpReadCigars,
                              TSize ,
                              TContigPos currentWindowEnd,
                              TOptions &options)
{

    typedef typename TFragmentStore::TAlignedReadStore          TMatches;
    typedef typename Value<TMatches>::Type                      TMatch;
    typedef typename Iterator<TMatches,Standard>::Type          TMatchIt;
    typedef typename Id<TFragmentStore>::Type                   TId;
    //typedef typename Value<TReadClips>::Type                    TPair;

    SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore),length(fragmentStore.alignQualityStore));

    ::std::sort(begin(fragmentStore.alignedReadStore, Standard()), end(fragmentStore.alignedReadStore, Standard()), LessGPos<TMatch>());

    if(options._debugLevel > 1 )::std::cout << "Copying matches overlapping more than one window ... \n";

    TMatchIt mIt        = end(fragmentStore.alignedReadStore,Standard());
    TMatchIt mItBegin   = begin(fragmentStore.alignedReadStore,Standard());
    --mIt;

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
            if(!empty(readCounts))appendValue(tmpReadCounts,readCounts[(*mIt).readId]);
            appendValue(tmpQualities,fragmentStore.alignQualityStore[(*mIt).readId]);
            appendValue(tmpReadCigars,readCigars[(*mIt).readId]);
//            appendValue(tmpReadClips,TPair(0,0));
            appendValue(tmpReadClips,readClips[(*mIt).readId]);
            maxCoord = std::max(maxCoord, (int)std::max(mIt->beginPos, mIt->endPos));
            minCoord = std::min(minCoord, (int)std::min(mIt->beginPos, mIt->endPos));
        }
        --mIt;
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

    //CharString str = "pileBef";
    //_dumpMatches(fragmentStore, str);

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
    //typedef typename Value<TMatches>::Type                  TMatch;
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
            avgRQ += (int) getQualityValue(read[i]) ;
        // watch out, this is new: use mapping quality if given
        if((fragmentStore.alignQualityStore[(*it).id]).score == 0 || (char)(avgRQ/length(read))<(fragmentStore.alignQualityStore[(*it).id]).score)
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
                avgRQ += (int) getQualityValue(read[i]) ;
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


template<typename TTable>
void
printHetTable(TTable & hetTable)
{
    cout << "Printing het table:" << std::endl;
    for (int n1=0; n1<256; ++n1)
    {
        for (int n2=0; n2<256; ++n2)
        {
            cout << hetTable[n1<<8|n2] << "\t";

        }
        cout << std::endl;
    }
    cout << std::endl;
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



//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int detectSNPs(SNPCallingOptions<TSpec> &options)
{

    typedef FragmentStore<SnpStoreSpec_>            TFragmentStore;
    //typedef typename TFragmentStore::TReadSeq       TReadSeq;               // TRead
    typedef typename TFragmentStore::TContigSeq     TContigSeq;             // TGenome
    //typedef typename Position<TReadSeq>::Type     TReadPos;               // TPos
    //typedef typename TFragmentStore::TReadPos       TReadPos;               // TPos
    //typedef typename Position<TContigSeq>::Type       TContigPos;             // TContigPos
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TAlignedReadStore  TAlignedReadStore;      // TMatches
    typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;     // TMatchQualities
    typedef typename TFragmentStore::TReadStore     TReadStore;             // TReadSet
    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;              // TReadSet
    typedef typename TFragmentStore::TContigStore       TContigStore;           // TGenomeSet
    typedef typename Value<TContigStore>::Type      TContig;
    typedef              TContigSeq                    TGenome;
    typedef StringSet<TGenome>                         TGenomeSet;

    typedef String<String<TContigPos > >            TPositions;
    typedef typename Iterator<String<TContigPos > >::Type TPosIterator;
    typedef typename Size<TGenomeSet>::Type         TGenomeSetSize;
    typedef ::std::map<CharString, TGenomeSetSize>  TGenomeMap;
    //typedef typename TGenomeMap::iterator           TMapIter;
    typedef String<unsigned>                TReadCounts;
    typedef String<Pair<int,int> >              TReadClips;
    typedef StringSet<String<Pair<char,int> > >     TReadCigars;

    TGenomeSet              genomes;
    StringSet<CharString>           genomeFileNameList; // filenamen
    StringSet<CharString>           genomeNames;                // todo: raus
    TGenomeMap              gIdStringToIdNumMap;    // name to id
    TPositions              positions; // list of positions to inspect


    // dump configuration in verbose mode
    if (options._debugLevel >= 1)
    {
        ::std::cerr << "___SETTINGS____________" << ::std::endl;
        ::std::cerr << "Genome file:                             \t" << options.genomeFName << ::std::endl;
        ::std::cerr << "Read files:                              \t" << options.readFNames[0] << ::std::endl;
        for(unsigned i = 1; i < length(options.readFNames); ++i)
            ::std::cerr << "                                         \t" << options.readFNames[i] << ::std::endl;
/*      if(options.inputFormat == 1)
        {
            ::std::cerr << "Quality files:                           \t" << qualityFNames[0] << ::std::endl;
            for(unsigned i = 1; i < length(qualityFNames); ++i)
                ::std::cerr << "                                         \t" << qualityFNames[i] << ::std::endl;
        }*/::std::cerr << "MaxPile:                                 \t" << options.maxPile << ::std::endl;
        if(options.laneSpecificMaxPile)::std::cerr << "Lane specific:                           \tYES" << ::std::endl;
        else ::std::cerr << "Lane specific:                           \tNO" << ::std::endl;
        ::std::cerr << "MinCoverage:                             \t" << options.minCoverage << ::std::endl;
        if(options.method == 0)
        {
            ::std::cerr << "MinMutThreshold:                         \t" << options.minMutT << ::std::endl;
            ::std::cerr << "MinPercentageThreshold:                  \t" << options.percentageT << ::std::endl;
            ::std::cerr << "MinQualityThreshold:                     \t" << options.avgQualT << ::std::endl;
        }
        else
        {
            ::std::cerr << "MinMappingQuality:                       \t" << options.minMapQual << ::std::endl;
        }
        if(options.outputIndel != "")
        {
            ::std::cerr << "IndelCountThreshold:                     \t" << options.indelCountThreshold << ::std::endl;
            ::std::cerr << "IndelPercentageThreshold:                \t" << options.indelPercentageT << ::std::endl;
            ::std::cerr << "IndelWindow:                             \t" << options.indelWindow << ::std::endl;
        }
        ::std::cerr << ::std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Step 1: Determine genome file type and load genomes
    SEQAN_PROTIMESTART(load_time);

    int result = getGenomeFileNameList(genomeFileNameList, options);
    if(result == CALLSNPS_GENOME_FAILED || !loadGenomes(genomes, genomeFileNameList, gIdStringToIdNumMap, genomeNames))
    {
        ::std::cerr << "Failed to open genome file " << options.genomeFName << ::std::endl;
        return result;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Step 2: Load fragmentStore.readSeqStore and fragmentStore.alignedReadStore
    // open read files and  store open file pointers
    String<int> highestChrId;
    resize(highestChrId,length(options.readFNames),0);
    vector< ::std::fstream* > readFileStreams;
    readFileStreams.resize(length(options.readFNames));
    for (unsigned i = 0; options.inputFormat == 0 && i < length(options.readFNames); ++i)
    {
        readFileStreams[i] = new fstream(toCString(options.readFNames[i]), ios_base::in | ios::binary);
        if(!(*(readFileStreams[i])).is_open())
        {
            ::std::cerr << "Failed to open read file " << options.readFNames[i] << ::std::endl;
            return CALLSNPS_GFF_FAILED;
        }
    }

    // Open all input SAM and BAM files.
    std::vector<BamFileIn *> bamFileIns;
    if (options.inputFormat >= 1)  // 1 -- SAM, 2 -- BAM
        for (unsigned i = 0; i < length(options.readFNames); ++i)
        {
            bamFileIns.push_back(new BamFileIn());
            //std::cout << "Opening SAM/BAM file\n";
            if (!open(*bamFileIns.back(), toCString(options.readFNames[i])))
            {
                std::cerr << "[ERROR] Could not open BAM file" << options.readFNames[i] << "\n";
                return 1;
            }
        }
    // Allocate one BamAlignmentRecord for each input BAM file.
    std::vector<BamAlignmentRecord> records(length(options.readFNames));

    /////////////////////////////////////////////////////////////////////
    // open out file streams and store open file pointers
    ::std::ofstream snpFileStream;
    if (options.outputSNP != "")
    {

        // prepare lookup tables for maq statistics
        if (options.method == 1 )
        {
            computeCnks(options.cnks,options.fks,options);
#ifdef CORRECTED_HET
            if(options.correctedHetTable) options.priorHetQ = computeHetTable(options.hetTable,options); // corrected het table computed with normal distribution
            else
#endif
                options.priorHetQ = computeHetTable(options.hetTable,options,MaqMethod()); // original Maq method
            //printHetTable(options.hetTable);
            //printHetTable(options.hetTable2);

        }

        snpFileStream.open(toCString(options.outputSNP),::std::ios_base::out);
        if(!snpFileStream.is_open())
            return CALLSNPS_OUT_FAILED;
        snpFileStream << "#" << (options.programCall).str() << "\n";
        if(options.outputFormat < 2)
        {
            if(options.orientationAware)
                snpFileStream << "#chr\tpos\tref\t[A+]\t[C+]\t[G+]\t[T+]\t[A-]\t[C-]\t[G-]\t[T-]\tcov\tcall";
            else
                snpFileStream << "#chr\tpos\tref\tA\tC\tG\tT\tcov\tcall";
            //if(options.method==1)
            snpFileStream << "\tquality\tsnpQ\n";
            //else
            //  file <<"\n";
        }
    }
    ::std::ofstream indelFileStream;
    if (options.outputIndel != "")
    {
        indelFileStream.open(toCString(options.outputIndel),::std::ios_base::out);
        if(!indelFileStream.is_open())
            return CALLSNPS_OUT_FAILED;
    }
    //  ::std::ofstream cnvFileStream;
    //  if (*options.outputCNV != 0)
    //  {
    //      cnvFileStream.open(options.outputCNV,::std::ios_base::out);
    //      if(!cnvFileStream.is_open())
    //          return CALLSNPS_OUT_FAILED;
    //  }

    ::std::ofstream posFileStream;
    if(options.inputPositionFile != "")
    {
        resize(positions,length(genomeNames));
        result = loadPositions(positions,gIdStringToIdNumMap, toCString(options.inputPositionFile),options);
        if(result != 0)
        {
            ::std::cerr << "Failed to read position file " << options.inputPositionFile << ::std::endl;
            return result;
        }
        posFileStream.open(toCString(options.outputPosition),::std::ios_base::out);
        if(!posFileStream.is_open())
            return CALLSNPS_OUT_FAILED;
        if(options.orientationAware)
            posFileStream << "#chr\tpos\t[A+]\t[C+]\t[G+]\t[T+]\t[gap+]\t[A-]\t[C-]\t[G-]\t[T-]\t[gap-]\n";
        else
            posFileStream << "#chr\tpos\tA\tC\tG\tT\tgap\n";

    }
    /////////////////////////////////////////////////////////////////////////////
    // helper variables
    Pair<int,int> zeroPair(0,0);
    int sumreads = 0;
    int sumwindows = 0;

    bool positionStatsOnly = (options.outputSNP == "" && options.outputPosition != "");
    TPosIterator inspectPosIt, inspectPosItEnd;
    // The iteration of the positions entries is convoluted with the rest of the code.  We create a
    // fake positions entry of the correct size (length(genomes)) but for positionStatsOnly to false
    // if no positions are given.
    if (empty(positions))
    {
        resize(positions, length(genomes));
        positionStatsOnly = false;  // ignor einspectPosIt below
    }

    bool firstCall = true;

    /////////////////////////////////////////////////////////////////////////////
    // Start scanning for SNPs/indels
    // for each chromosome
    for(unsigned i=0; i < length(genomes); ++i)
    {
        //std::cout << genomeNames[i] << "\n";
        inspectPosIt = begin(positions[i], Standard());
        inspectPosItEnd = end(positions[i], Standard());
        if (positionStatsOnly && inspectPosIt == inspectPosItEnd)
            continue;
        // parse matches batch by batch
        TContigPos currentWindowBegin = 0;
        if(options._debugLevel > 0) ::std::cout << "Scanning genome #" << i << " ..." << ::std::endl;

        // containers for those matches that overlap two windows
        TAlignedReadStore tmpMatches;
        TAlignQualityStore tmpQualities;
        TReadStore tmpRs;
        TReadSeqStore tmpReads;
        TReadCounts tmpReadCounts;
        TReadClips tmpReadClips;
        TReadCigars tmpReadCigars;
        options.minCoord = std::numeric_limits<unsigned>::max();
        options.maxCoord = 0;

        // snp calling is done for all positions between windowBegin and windowEnd
        while(currentWindowBegin < (TContigPos)length(genomes[i]))
        {
            TContigPos currentWindowEnd = currentWindowBegin + options.windowSize;
            if(currentWindowEnd > (TContigPos)length(genomes[i])) currentWindowEnd = (TContigPos)length(genomes[i]);

            if(positionStatsOnly && !(currentWindowBegin <= *inspectPosIt && *inspectPosIt< currentWindowEnd))
            {
                SEQAN_ASSERT_LEQ(currentWindowBegin,*inspectPosIt);
                clear(tmpReads); clear(tmpRs); clear(tmpMatches); clear(tmpQualities);
                clear(tmpReadClips); clear(tmpReadCounts); clear(tmpReadCigars);
                currentWindowBegin = currentWindowEnd;
                ++sumwindows;
                continue;
            }
            if(options._debugLevel > 0)
                ::std::cout << "Sequence number " << i << " window " << currentWindowBegin << ".." << currentWindowEnd << "\n";

            TFragmentStore fragmentStore;
            TReadCounts readCounts;  // Count number of reads that are identical to the given one. Useful for micro RNA data where there were millions of identical reads. Must be in GFF input, not supported for SAM input.
            TReadClips readClips;  // Soft clipping information and/or clipping information from GFF/SAM tag. Clipping is postponed after pileup correction.
            TReadCigars readCigars; // Currently only stored for split-mapped reads. Split-mapped reads need special handling, especially for realignment.

            // add the matches that were overlapping with this and the last window (copied in order to avoid 2 x makeGlobal)
            if(!empty(tmpMatches))
            {
                sumreads -=  length(tmpReads);  // count these reads only once
                resize(fragmentStore.alignQualityStore,length(tmpMatches));
                resize(fragmentStore.alignedReadStore,length(tmpMatches));
                resize(fragmentStore.readSeqStore,length(tmpMatches));
                resize(fragmentStore.readStore,length(tmpMatches));
                if(!empty(tmpReadClips))resize(readClips,length(tmpMatches));
                if(!empty(tmpReadCounts)) resize(readCounts,length(tmpMatches));
                if(!empty(tmpReadCigars))resize(readCigars,length(tmpMatches));

                arrayMoveForward(begin(tmpQualities,Standard()), end(tmpQualities,Standard()), begin(fragmentStore.alignQualityStore,Standard()));
                arrayMoveForward(begin(tmpMatches,Standard()), end(tmpMatches,Standard()), begin(fragmentStore.alignedReadStore,Standard()));
                arrayMoveForward(begin(tmpReads,Standard()), end(tmpReads,Standard()), begin(fragmentStore.readSeqStore,Standard()));
                arrayMoveForward(begin(tmpRs,Standard()), end(tmpRs,Standard()), begin(fragmentStore.readStore,Standard()));
                if(!empty(tmpReadCounts)) arrayMoveForward(begin(tmpReadCounts,Standard()), end(tmpReadCounts,Standard()), begin(readCounts,Standard()));
                if(!empty(tmpReadCigars)) arrayMoveForward(begin(tmpReadCigars,Standard()), end(tmpReadCigars,Standard()), begin(readCigars,Standard()));
                if(!empty(tmpReadClips)) arrayMoveForward(begin(tmpReadClips,Standard()), end(tmpReadClips,Standard()), begin(readClips,Standard()));
#ifdef SNPSTORE_DEBUG
                CharString strstr = "afterCopyInFrag";
                _dumpMatches(fragmentStore,strstr);
#endif

            }

            // parse matches for current window
            if(options._debugLevel > 0)
                ::std::cout << "Parsing reads up to position " << currentWindowEnd << "...\n";
            for(unsigned j = 0; j < length(options.readFNames); ++j)
            {
                unsigned sizeBefore = length(fragmentStore.alignedReadStore);


                // currently only gff supported
                if (options.inputFormat == 0) // GFF
                    result = readMatchesFromGFF_Batch(readFileStreams[j], fragmentStore, readCounts, readClips,
                                                      readCigars, genomes[i], gIdStringToIdNumMap,
                                                      i, currentWindowBegin, currentWindowEnd, highestChrId[j], options);
                if (options.inputFormat >= 1) // SAM
                    result = readMatchesFromSamBam_Batch(*bamFileIns[j], records[j], fragmentStore, readCounts, readClips,
                                                         readCigars, genomes[i], gIdStringToIdNumMap,
                                                         i, currentWindowBegin, currentWindowEnd, highestChrId[j], options, firstCall);
                firstCall = false;

                if (result == CALLSNPS_GFF_FAILED)
                {
                    std::cerr << "Failed to open read file " << options.readFNames[j] << ::std::endl;
                    std::cerr << "or reads are not sorted correctly. " << ::std::endl;
                    return result;
                }
                if (result > 0)
                    return result;

                if (options._debugLevel > 0)
                    std::cout << "parsed reads of file " << j << "\n";

                // store average quality of each read
                addReadQualityToMatches(fragmentStore,sizeBefore,(unsigned)length(fragmentStore.alignedReadStore),options);

                // do pile up correction if lane-based. lane-specific pileup correction not really used
                if(options.maxPile != 0 && options.laneSpecificMaxPile)
                    applyPileupCorrection(fragmentStore,sizeBefore,(unsigned)length(fragmentStore.alignedReadStore),options);

            }

            if (options._debugLevel > 1)  // number of reads currently in memory
                ::std::cerr << lengthSum(fragmentStore.readSeqStore) << " bps of " << length(fragmentStore.readSeqStore) << " reads in memory." << ::std::endl;
            sumreads +=  length(fragmentStore.readSeqStore);  // total count of reads

            // do merged pile up correction
            if(options.maxPile != 0 && !options.laneSpecificMaxPile)
                applyPileupCorrection(fragmentStore,(unsigned)0,(unsigned)length(fragmentStore.alignedReadStore),options);

            // these were set while parsing matches, first and last position of parsed matches
            //          TContigPos startCoord = options.minCoord;// can be < currentWindowBegin
            //          TContigPos endCoord = options.maxCoord; // can be > currentWindoEnd

            // check
            TContigPos startCoord = _max((int)options.minCoord-options.realignAddBorder,0);// can be < currentWindowBegin
            TContigPos endCoord = _min(options.maxCoord+options.realignAddBorder,length(genomes[i])); // can be > currentWindoEnd


            if(!empty(fragmentStore.alignedReadStore))
            {
                //initial values of min and max coords for next round are set here
                if(currentWindowEnd != (TContigPos)length(genomes[i]))
                {
                    clear(tmpMatches);
                    clear(tmpQualities);
                    clear(tmpRs);
                    clear(tmpReads);
                    clear(tmpReadCounts);
                    clear(tmpReadClips);
                    clear(tmpReadCigars);
                    copyNextWindowMatchesAndReads(fragmentStore,readCounts,readCigars,readClips,tmpReadCounts,tmpReads,tmpRs,tmpMatches,tmpQualities,tmpReadClips,tmpReadCigars,i,currentWindowEnd,options);
#ifdef SNPSTORE_DEBUG
                    CharString strstr = "afterCopyInTmp";
                    _dumpMatches(tmpReads,tmpMatches,tmpQualities,strstr);
#endif
                }

                //  ::std::cout << "Min = " << options.minCoord << " Max = " << options.maxCoord << std::endl;
                //  ::std::cout << "Min = " << startCoord << " Max = " << endCoord << std::endl;

                // coordinates are relative to current chromosomal window (segment)
                transformCoordinates(fragmentStore,startCoord,options);

                // set the current chromosomal segment as contig sequence
                TContig conti;
                conti.seq = infix(genomes[i],startCoord,endCoord);
                appendValue(fragmentStore.contigStore, conti, Generous() );
                appendValue(fragmentStore.contigNameStore, genomeNames[i], Generous() );// internal id is always 0

                // clip Reads if clipping is switched on and there were clip tags in the gff file
                if((!options.dontClip && options.clipTagsInFile) || options.softClipTagsInFile)
                {
                    options.useBaseQuality = false; // activate "average read quality"-mode for snp calling, low quality bases should be clipped anyway
                    clipReads(fragmentStore,readClips,(unsigned)0,(unsigned)length(fragmentStore.alignedReadStore),options);
                }

                // check for indels
                if (options.outputIndel != "")
                {
                    if(options._debugLevel > 1) ::std::cout << "Check for indels..." << std::endl;
                    if(!options.realign) dumpShortIndelPolymorphismsBatch(fragmentStore, readCigars, fragmentStore.contigStore[0].seq, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, indelFileStream, options);
                }

                // // check for CNVs
                //              if (*options.outputCNV != 0)
                //                  dumpCopyNumberPolymorphismsBatch(fragmentStore, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, cnvFileStream, options);

#ifdef SNPSTORE_DEBUG
                CharString strstr = "test";
                //              _dumpMatches(fragmentStore, strstr );
#endif
                if (options.outputSNP != "")
                {
                    if(options._debugLevel > 1) ::std::cout << "Check for SNPs..." << std::endl;
                    if(options.realign)
                        dumpVariantsRealignBatchWrap(fragmentStore, readCigars, readCounts, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, snpFileStream,indelFileStream,options);
                    else
                        dumpSNPsBatch(fragmentStore, readCigars, readCounts, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, snpFileStream,options);
                }
                if(positionStatsOnly)
                {
                    if(options._debugLevel > 1) ::std::cout << "Dumping info for query positions..." << std::endl;
                    if(options.realign)
                        dumpPositionsRealignBatchWrap(fragmentStore, inspectPosIt, inspectPosItEnd, readCigars, readCounts, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, posFileStream, options);
                    else dumpPosBatch(fragmentStore, inspectPosIt, inspectPosItEnd, readCigars, readCounts, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, posFileStream, options);
                }
            }
//            else
//            {
                if(positionStatsOnly && !empty(positions))
                {
                    while(inspectPosIt != inspectPosItEnd && *inspectPosIt < currentWindowEnd)
                    {
                        if(options.orientationAware)
                            posFileStream << genomeNames[i] << '\t' << *inspectPosIt + options.positionFormat << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" << std::endl;
                        else
                            posFileStream << genomeNames[i] << '\t' << *inspectPosIt + options.positionFormat << "\t0\t0\t0\t0\t0\t0" << std::endl;
                        ++inspectPosIt;
                    }
                }
//            }
            currentWindowBegin = currentWindowEnd;
            ++sumwindows;
        }

    }
    if (options.outputSNP != "")
        snpFileStream.close();

    if (options.outputIndel != "")
        indelFileStream.close();

    if (options.outputPosition != "")
        posFileStream.close();

    //  if (options.outputCNV != "")
    //      cnvFileStream.close();

    return 0;
}




// log file to keep track of happenings
template <typename TSpec>
int writeLogFile(int argc, const char *argv[], SNPCallingOptions<TSpec> &options)
{
    ::std::ofstream logfile;
    logfile.open(toCString(options.outputLog), ::std::ios_base::out | ::std::ios_base::trunc);
    if (!logfile.is_open())
    {
        ::std::cerr << "Failed to open log file" << ::std::endl;
        return 1;
    }
    logfile << "#call" << std::endl;
    for (int i = 0; i < argc; ++i)
        logfile << argv[i] << " ";
    logfile << std::endl;

    logfile << "#files" << ::std::endl;
    logfile << "Genome=\"" << options.genomeFName << "\""<<::std::endl;
    logfile << "Reads=\"" << options.readFNames[0];
    for(unsigned i = 1; i < length(options.readFNames); ++i)
        logfile << " " << options.readFNames[i] << ::std::endl;
    logfile << "\"" << std::endl;
    if(options.outputSNP != "")
    {
        logfile << "OutputSnp=\"" << CharString(options.outputSNP) << "\"" << ::std::endl;
    }
    if(options.outputIndel != "")
    {
        logfile << "OutputIndel=\"" << CharString(options.outputIndel) << "\""<< ::std::endl;
    }
    logfile << "#settings" << ::std::endl;
    logfile << "MaxPile=" << options.maxPile << ::std::endl;
    logfile << "MinCov=" << options.minCoverage << ::std::endl;
    logfile << "Method=" << options.method << ::std::endl;
    if(options.method == 0)
    {
        logfile << "MinMutT=" << options.minMutT << ::std::endl;
        logfile << "MinPercT=" << options.percentageT << ::std::endl;
        logfile << "MinQualT=" << options.avgQualT << ::std::endl;
    }
    else
    {
        logfile << "MinMapQual=" << options.minMapQual << ::std::endl;
    }
    if(options.outputIndel != "")
    {
        logfile << "MinIndel=" << options.indelCountThreshold << ::std::endl;
        logfile << "MinPercIndelT=" << options.indelPercentageT << ::std::endl;
    }
    logfile.close();
    return 0;
}


template <typename TSpec>
ArgumentParser::ParseResult
parseCommandLine(SNPCallingOptions<TSpec> & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("snp_store");
    // Set short description, version, and date.
    setShortDescription(parser, "SnpStore");
    setCategory(parser, "Variant Detection");
	setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIALIGNMENT FILE\\fP> [<\\fIALIGNMENT FILE\\fP> ...]");

    // We require two mandatory arguments: genome file and read file(s)
    addArgument(parser, ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "GENOME"));
    setValidValues(parser, 0, ".fa .fasta");
    setHelpText(parser, 0, "A reference genome file.");

    std::vector<std::string> alignmentFormats(BamFileIn::getFileExtensions());
    alignmentFormats.push_back(".gff");
    addArgument(parser, ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "ALIGNMENTS", true));
    setValidValues(parser, 1, alignmentFormats);
    setHelpText(parser, 1, "Read alignment file(s) sorted by genomic position.");

    addDescription(parser, "SNP and Indel Calling in Mapped Read Data.");
    addSection(parser, "Main Options");
    addOption(parser, ArgParseOption("o", "output", "SNP output file (mandatory).", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "output", ".vcf");
    setRequired(parser, "output");
    addOption(parser, ArgParseOption("osc", "only-successful-candidates", "Output only successfully called SNP candidates. Default: Output all candidates."));
    addOption(parser, ArgParseOption("dc", "dont-clip", "Ignore clip tags in gff. Default: off."));

    addOption(parser, ArgParseOption("mu", "multi", "Keep non-unique fragmentStore.alignedReadStore. Default: off."));
    addOption(parser, ArgParseOption("hq", "hide-qualities", "Only show coverage (no qualities) in SNP output file. Default: off."));
    addOption(parser, ArgParseOption("sqo", "solexa-qual-offset", "Base qualities are encoded as char value - 64 (instead of char - 33)."));
    addOption(parser, ArgParseOption("id", "indel-file", "Output file for called indels in gff format. Default: off.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "indel-file", ".gff");
    addOption(parser, ArgParseOption("m", "method", "Set method used for SNP calling either threshold based or Maq method.", ArgParseArgument::STRING));
    setValidValues(parser, "method", "thresh maq");
    setDefaultValue(parser, "method", "maq");
    addOption(parser, ArgParseOption("mp", "max-pile", "Maximal number of matches allowed to pile up at the same genome position.", ArgParseArgument::INTEGER));
    setMinValue(parser, "max-pile", "1");
    setDefaultValue(parser, "max-pile", options.maxPile);
    addOption(parser, ArgParseOption("mmp", "merged-max-pile", "Do pile up correction on merged lanes. Default: off."));
    addOption(parser, ArgParseOption("mc", "min-coverage", "Minimal required number of reads covering a candidate position.", ArgParseArgument::INTEGER));
    setMinValue(parser, "min-coverage", "1");
    setDefaultValue(parser, "min-coverage", options.minCoverage);
    addOption(parser, ArgParseOption("fc", "force-call", "Always call base if count is >= fc, ignore other parameters. Default: off.", ArgParseArgument::INTEGER));
    setMinValue(parser, "force-call", "1");
    setDefaultValue(parser, "force-call", options.forceCallCount);
    addOption(parser, ArgParseOption("oa", "orientation-aware", "Distinguish between forward and reverse reads. Default: off."));
//    addOption(parser, ArgParseOption("cnv", "output-cnv", "Name of CNV result file.", ArgParseArgument::OUTPUT_FILE));
//    hideOption(parser, "cnv");
    addOption(parser, ArgParseOption("op", "output-positions", "Name of positions output file.", ArgParseArgument::STRING));
    hideOption(parser, "op");
    addOption(parser, ArgParseOption("ip", "input-positions", "Name of positions input file.", ArgParseArgument::STRING));
    hideOption(parser, "ip");
    addOption(parser, ArgParseOption("mpr", "max-polymer-run", "Discard indels in homopolymer runs longer than mpr.", ArgParseArgument::INTEGER));
    setMinValue(parser, "max-polymer-run", "0");
    setDefaultValue(parser, "max-polymer-run", options.maxPolymerRun);
    addOption(parser, ArgParseOption("dp", "diff-pos", "Minimal number of different read positions supporting the mutation.", ArgParseArgument::INTEGER));
    setMinValue(parser, "diff-pos", "0");
    setDefaultValue(parser, "diff-pos", options.minDifferentReadPos);
    addOption(parser, ArgParseOption("eb", "exclude-border", "Exclude read positions within eb base pairs of read borders for SNV calling. Default: off.", ArgParseArgument::INTEGER));
    setMinValue(parser, "exclude-border", "0");
    setDefaultValue(parser, "exclude-border", options.excludeBorderPos);
    addOption(parser, ArgParseOption("su", "suboptimal", "Keep suboptimal reads. Default: off"));
    addOption(parser, ArgParseOption("re", "realign", "Realign reads around indel candidates. Default: off"));
    addOption(parser, ArgParseOption("cq", "corrected-quality", "New quality calibration factor.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "corrected-quality", "0");
    setDefaultValue(parser, "corrected-quality", options.newQualityCalibrationFactor);
    hideOption(parser, "cq");
    addOption(parser, ArgParseOption("pws", "parse-window-size", "Genomic window size for parsing reads (concerns memory consumption, choose smaller windows for higher coverage).", ArgParseArgument::INTEGER));
    setMinValue(parser, "parse-window-size", "1");
    setDefaultValue(parser, "parse-window-size", options.windowSize);
    addOption(parser, ArgParseOption("reb", "realign-border", "Realign border.", ArgParseArgument::INTEGER));
    setMinValue(parser, "realign-border", "0");
    setMaxValue(parser, "realign-border", "10");
    setDefaultValue(parser, "realign-border", options.realignAddBorder);
    hideOption(parser, "reb");

    addSection(parser, "SNP calling options");
    addSection(parser, " Threshold method related");
    addOption(parser, ArgParseOption("mm", "min-mutations", "Minimal number of observed mutations for mutation to be called.", ArgParseArgument::INTEGER));
    setMinValue(parser, "min-mutations", "1");
    setDefaultValue(parser, "min-mutations", options.minMutT);
    addOption(parser, ArgParseOption("pt", "perc-threshold", "Minimal percentage of mutational base for mutation to be called.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "perc-threshold", "0");
    setDefaultValue(parser, "perc-threshold", options.percentageT);
    addOption(parser, ArgParseOption("mq", "min-quality", "Minimal average quality of mutational base for mutation to be called.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "min-quality", "0");
    setDefaultValue(parser, "min-quality", options.avgQualT);
    addSection(parser, " Maq method related");
    addOption(parser, ArgParseOption("th", "theta", "Dependency coefficient.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "theta", "0");
    setDefaultValue(parser, "theta", options.theta);
    addOption(parser, ArgParseOption("hr", "hetero-rate", "Heterozygote rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "hetero-rate", "0");
    setMaxValue(parser, "hetero-rate", "1");
    setDefaultValue(parser, "hetero-rate", options.hetRate);
    addOption(parser, ArgParseOption("mmq", "min-map-quality", "Minimum base call (mapping) quality for a match to be considered.", ArgParseArgument::INTEGER));
    setMinValue(parser, "min-map-quality", "0");
    setDefaultValue(parser, "min-map-quality", options.minMapQual);
    addOption(parser, ArgParseOption("ch", "corrected-het", "Use amplification bias corrected distribution for heterozygotes. Default: off."));
    addOption(parser, ArgParseOption("maf", "mean-alleleFreq", "Mean ref allele frequency in heterozygotes.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "mean-alleleFreq", "0");
    setDefaultValue(parser, "mean-alleleFreq", options.meanAlleleFrequency);
    addOption(parser, ArgParseOption("ac", "amp-cycles", "Number of amplification cycles.", ArgParseArgument::INTEGER));
    setMinValue(parser, "amp-cycles", "0");
    setDefaultValue(parser, "amp-cycles", options.amplificationCycles);
    addOption(parser, ArgParseOption("ae", "amp-efficiency", "Polymerase efficiency, probability of amplification.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "amp-efficiency", "0");
    setMaxValue(parser, "amp-efficiency", "1");
    setDefaultValue(parser, "amp-efficiency", options.amplificationEfficiency);
    addOption(parser, ArgParseOption("in", "initial-N", "Initial allele population size.", ArgParseArgument::INTEGER));
    setMinValue(parser, "initial-N", "0");
    setDefaultValue(parser, "initial-N", options.initialN);
    addOption(parser, ArgParseOption("pht", "print-hetTable", "Print het table. Default: off."));
    hideOption(parser, "pht");
    addOption(parser, ArgParseOption("mec", "min-explained-column", "Minimum fraction of alignment column reads explained by genotype call.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "min-explained-column", "0");
    setMaxValue(parser, "min-explained-column", "1");
    setDefaultValue(parser, "min-explained-column", options.minExplainedColumn);

    addSection(parser, "Indel calling options");
    addOption(parser, ArgParseOption("it", "indel-threshold", "Minimal number of indel-supporting reads required for indel calling.", ArgParseArgument::INTEGER));
    setMinValue(parser, "indel-threshold", "1");
    setDefaultValue(parser, "indel-threshold", options.indelCountThreshold);
    addOption(parser, ArgParseOption("ipt", "indel-perc-threshold", "Minimal ratio of indel-supporting/covering reads for indel to be called.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "indel-perc-threshold", "0");
    setMaxValue(parser, "indel-perc-threshold", "1");
    setDefaultValue(parser, "indel-perc-threshold", options.indelPercentageT);
    addOption(parser, ArgParseOption("iqt", "indel-quality-thresh", "Minimal average quality of inserted base/deletion-neighboring bases for indel to be called.", ArgParseArgument::INTEGER));
    setMinValue(parser, "indel-quality-thresh", "0");
    setDefaultValue(parser, "indel-quality-thresh", options.indelQualityThreshold);
    addOption(parser, ArgParseOption("bsi", "both-strands-indel", "Both strands need to be observed for indel to be called. Default: off."));
    addOption(parser, ArgParseOption("iw", "indel-window", "Overlap window used for indel calling.", ArgParseArgument::INTEGER));
    setMinValue(parser, "indel-window", "0");
    setDefaultValue(parser, "indel-window", options.indelWindow);
    hideOption(parser, "iw");
    addOption(parser, ArgParseOption("ebi", "exclude-border-indel", "Same as option -eb but for indel candidates.", ArgParseArgument::INTEGER));
    setMinValue(parser, "exclude-border-indel", "0");
    setDefaultValue(parser, "exclude-border-indel", options.indelDepthMinOverlap);
    addOption(parser, ArgParseOption("cws", "cnv-window-size", "CNV window size.", ArgParseArgument::INTEGER));
    setMinValue(parser, "cnv-window-size", "1");
    setMaxValue(parser, "cnv-window-size", "10000");
    setDefaultValue(parser, "cnv-window-size", options.cnvWindowSize);
    hideOption(parser, "cws");

    addSection(parser, "Other options");
    addOption(parser, ArgParseOption("lf", "log-file", "Write log to FILE.", ArgParseArgument::STRING));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));

    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBsnp_store\\fP \\fB-mc\\fP \\fB2\\fP \\fB-it\\fP \\fB2\\fP \\fBexampleGenome.fa\\fP \\fBexampleReads.gff\\fP \\fB-o\\fP \\fBexampleSNPs.vcf\\fP \\fB-id\\fP \\fBexampleIndels.gff\\fP",
                "Call SNPs and indels of a low-coverage example (minimum coverage and indel threshold were reduced to 2).");
    addListItem(parser,
                "\\fBsnp_store\\fP \\fB-re\\fP \\fB-mc\\fP \\fB2\\fP \\fB-it\\fP \\fB2\\fP \\fBexampleGenome.fa\\fP \\fBexampleReads.gff\\fP \\fB-o\\fP \\fBexampleSNPs.vcf\\fP \\fB-id\\fP \\fBexampleIndels.gff\\fP",
                "Computes a realignment before variant calling. Now, the two 1bp insertions should have been merged into one 2bp insertion.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    // Options:
    getOptionValue(options.outputSNP, parser, "output");

    if (isSet(parser, "only-successful-candidates"))
        options.outputFormat = 1;
    options.dontClip = isSet(parser, "dont-clip");
    options.keepMultiReads = isSet(parser, "multi");
    options.showQualityStrings = !isSet(parser, "hide-qualities");
    if (isSet(parser, "solexa-qual-offset"))
        options.asciiQualOffset = 64;
    getOptionValue(options.outputIndel, parser, "indel-file");
    std::string tmp;
    getOptionValue(tmp, parser, "method");
    options.method = (tmp == "maq")? 1: 0;
    getOptionValue(options.maxPile, parser, "max-pile");
    options.laneSpecificMaxPile = !isSet(parser, "merged-max-pile");
    getOptionValue(options.minCoverage, parser, "min-coverage");
    getOptionValue(options.forceCallCount, parser, "force-call");
    options.orientationAware = isSet(parser, "orientation-aware");
//    getOptionValue(options.outputCNV, parser, "output-cnv");
    getOptionValue(options.outputPosition, parser, "output-positions");
    getOptionValue(options.inputPositionFile, parser, "input-positions");
    getOptionValue(options.maxPolymerRun, parser, "max-polymer-run");
    getOptionValue(options.minDifferentReadPos, parser, "diff-pos");
    getOptionValue(options.excludeBorderPos, parser, "exclude-border");
    options.keepSuboptimalReads = isSet(parser, "suboptimal");
    options.realign = isSet(parser, "realign");
    getOptionValue(options.newQualityCalibrationFactor, parser, "corrected-quality");
    getOptionValue(options.windowSize, parser, "parse-window-size");
    getOptionValue(options.realignAddBorder, parser, "realign-border");
    // SNP Calling Options:
    getOptionValue(options.minMutT, parser, "min-mutations");
    getOptionValue(options.percentageT, parser, "perc-threshold");
    getOptionValue(options.avgQualT, parser, "min-quality");
    getOptionValue(options.theta, parser, "theta");
    getOptionValue(options.hetRate, parser, "hetero-rate");
    getOptionValue(options.minMapQual, parser, "min-map-quality");
    options.correctedHetTable = isSet(parser, "corrected-het");
    getOptionValue(options.meanAlleleFrequency, parser, "mean-alleleFreq");
    getOptionValue(options.amplificationCycles, parser, "amp-cycles");
    getOptionValue(options.amplificationEfficiency, parser, "amp-efficiency");
    getOptionValue(options.initialN, parser, "initial-N");
    options.printHetTable = isSet(parser, "print-hetTable");
    getOptionValue(options.minExplainedColumn, parser, "min-explained-column");
    // Indel Calling Options:
    getOptionValue(options.indelCountThreshold, parser, "indel-threshold");
    getOptionValue(options.indelPercentageT, parser, "indel-perc-threshold");
    getOptionValue(options.indelQualityThreshold, parser, "indel-quality-thresh");
    options.bothIndelStrands = isSet(parser, "both-strands-indel");
    getOptionValue(options.indelWindow, parser, "indel-window");
    getOptionValue(options.indelDepthMinOverlap, parser, "exclude-border-indel");
    getOptionValue(options.cnvWindowSize, parser, "cnv-window-size");
    // Other Options:
    getOptionValue(options.outputLog, parser, "log-file");
    if (isSet(parser, "verbose"))
        options._debugLevel = max(options._debugLevel, 1);
    if (isSet(parser, "very-verbose"))
        options._debugLevel = max(options._debugLevel, 2);

    getArgumentValue(options.genomeFName, parser, 0);
    unsigned countFiles = getArgumentValueCount(parser, 1);

    if (countFiles == 0)
    {
        cerr << "No mapping files specified." << endl;
        return ArgumentParser::PARSE_ERROR;
    }
    resize(options.readFNames, countFiles);
    for (unsigned i = 0; i < countFiles; ++i)
    {
        getArgumentValue(options.readFNames[i], parser, 1, i);

        // Get lower case of the output file name.  File endings are accepted in both upper and lower case.
        CharString tmp = options.readFNames[i];
        toLower(tmp);
        unsigned format = 0;
        if (endsWith(tmp, ".gff"))
            format = 0;
        else if (endsWith(tmp, ".sam"))
            format = 1;
        else if (endsWith(tmp, ".bam"))
            format = 2;

        if (i == 0)
        {
            options.inputFormat = format;
        }
        else
        {
            if (options.inputFormat != format)
            {
                cerr << "All mapping files must have the same format." << endl;
                return ArgumentParser::PARSE_ERROR;
            }
        }
    }

    // some additional option checking:

    //if(options.inputFormat == 1 && (!options.qualityFile || (length(qualityFNames)!=length(options.readFNames))))
    //{
    //  cerr << "If mapped read file is in Eland format, a .qual or .fastq file containing read qualities needs to be specified." << endl << endl;
    //  return 0;
    //}
    if(options.inputPositionFile == "" && options.outputPosition != "")
    {
        cerr << "Position analysis output specified, but no position file given." << endl << endl;
        return ArgumentParser::PARSE_ERROR;
    }

    if((options.realign && options.windowSize > 50000) || options.windowSize > 1000000)
        options.windowSize = 10000;

    if (options.outputLog != "")
        writeLogFile(argc, argv, options);

    if(options.runID == "")
    {
        ::std::string tempStr = toCString(options.readFNames[0]);
        size_t lastPos = tempStr.find_last_of("/\\");
        if (lastPos == tempStr.npos)
            lastPos = 0;
        else
            ++lastPos;
        options.runID = tempStr.substr(lastPos);
    }

    return ArgumentParser::PARSE_OK;
}



int main(int argc, const char *argv[])
{

    ArgumentParser parser;
    SNPCallingOptions<>     options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    for(int arg = 0; arg < argc; ++arg) {
        options.programCall << argv[arg] << " ";
    }

    //////////////////////////////////////////////////////////////////////////////
    // check for variants
    int result = detectSNPs(options);
    if (result > 0)
    {
        cerr << "ERROR: Something went wrong. Try 'snpStore --help' for more information." << endl << endl;
        return 0;
    }
    return result;
}
