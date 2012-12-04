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


#ifdef PLATFORM_WINDOWS
#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
#define SEQAN_DEFAULT_TMPDIR "./"
#endif


//#include "../../../extras/apps/rep_sep/utils.h"
//#include "../../../extras/apps/rep_sep/assembly_parser.h"
//#include "../../../extras/apps/rep_sep/column_scanner.h"
//#include "../../../extras/apps/rep_sep/rgraph.h"
//#include "../../../extras/apps/rep_sep/rep_sep.h"

#include "snp_store.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;
using namespace seqan;



// load entire genome into memory
template <typename TGenomeSet, typename TGenomeNames>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList, ::std::map<CharString,unsigned> &gIdStringToIdNumMap, TGenomeNames & genomeNames)
{
    unsigned gSeqNo = 0;
    unsigned filecount = 0;
    CharString temp;
    clear(genomeNames);
    while(filecount < length(fileNameList))
    {
        clear(temp);
        MultiFasta multiFasta;
        if (!open(multiFasta.concat, toCString(fileNameList[filecount]), OPEN_RDONLY)) return false;
        split(multiFasta, Fasta());
        
        unsigned seqCount = length(multiFasta);
        if(length(genomes) < gSeqNo+seqCount) 
            resize(genomes,gSeqNo+seqCount);
        for(unsigned i = 0; i < seqCount; ++i)
        {
            assignSeq(genomes[gSeqNo+i], multiFasta[i], Fasta());       // read Genome sequence
            assignSeqId(temp, multiFasta[i], Fasta());
            for (unsigned pos = 0; pos < length(temp); ++pos)
            {
                if(temp[pos]=='\t' || temp[pos]=='\b' || temp[pos]==' ')
                {
                    resize(temp,pos);
                    break;
                }
            }
            gIdStringToIdNumMap.insert(::std::make_pair(temp,gSeqNo+i)); // keeps the whole fasta ID including white spaces
            appendValue(genomeNames,temp);
        }
        gSeqNo += seqCount;
        ++filecount;
    }
    resize(genomes,gSeqNo);
    return (gSeqNo > 0);
}





// transform global cooridnates to coordinates relative to chromosomal segment
template<typename TFragmentStore, typename TContigPos, typename TOptions>
void 
transformCoordinates(TFragmentStore &fragmentStore, TContigPos startCoord, TOptions&)
{
    typedef typename TFragmentStore::TAlignedReadStore          TMatches;
    typedef typename Value<TMatches>::Type                      TMatch;
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
    typedef typename Value<TReadClips>::Type                    TPair;

    SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore),length(fragmentStore.alignQualityStore));

    ::std::sort(begin(fragmentStore.alignedReadStore, Standard()), end(fragmentStore.alignedReadStore, Standard()), LessGPos<TMatch>());    

    if(options._debugLevel > 1 )::std::cout << "Copying matches overlapping more than one window ... \n";
    
    TMatchIt mIt        = end(fragmentStore.alignedReadStore,Standard());
    TMatchIt mItBegin   = begin(fragmentStore.alignedReadStore,Standard());
    --mIt;

    // We will use minCoord/maxCoord to store the temporarily minimal and maximal coordinates in the window.
    int minCoord = maxValue<int>();
    int maxCoord = minValue<int>();
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
    if (minCoord != maxValue<int>())
        options.minCoord = minCoord;
    if (maxCoord != minValue<int>())
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
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename TFragmentStore::TReadSeqStore      TReads;
    typedef typename Value<TReads>::Type                TRead;
    typedef typename TFragmentStore::TContigStore       TGenomeSet;
    typedef typename Value<TGenomeSet>::Type            TGenome;
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
    typedef typename Value<TMatches>::Type                  TMatch;
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
int detectSNPs(
               const char *genomeFileName,
               String<CharString> & readFNames,
               String<CharString> &,
               SNPCallingOptions<TSpec> &options)
{
    
    typedef FragmentStore<SnpStoreSpec_>            TFragmentStore;
    typedef typename TFragmentStore::TReadSeq       TReadSeq;               // TRead
    typedef typename TFragmentStore::TContigSeq     TContigSeq;             // TGenome
    //typedef typename Position<TReadSeq>::Type     TReadPos;               // TPos
    typedef typename TFragmentStore::TReadPos       TReadPos;               // TPos
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
    typedef ::std::map<CharString,unsigned>         TGenomeMap;
    typedef typename TGenomeMap::iterator           TMapIter;
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
        ::std::cerr << "Genome file:                             \t" << genomeFileName << ::std::endl;
        ::std::cerr << "Read files:                              \t" << readFNames[0] << ::std::endl;
        for(unsigned i = 1; i < length(readFNames); ++i)
            ::std::cerr << "                                         \t" << readFNames[i] << ::std::endl;
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
        if(*options.outputIndel != 0)
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
    
    int result = getGenomeFileNameList(genomeFileName, genomeFileNameList, options);
    if(result == CALLSNPS_GENOME_FAILED || !loadGenomes(genomes, genomeFileNameList,gIdStringToIdNumMap,genomeNames))
    {
        ::std::cerr << "Failed to open genome file " << genomeFileName << ::std::endl;
        return result;
    }
    
    //////////////////////////////////////////////////////////////////////////////
    // Step 2: Load fragmentStore.readSeqStore and fragmentStore.alignedReadStore
    // open read files and  store open file pointers
    String<int> highestChrId;
    resize(highestChrId,length(readFNames),0);
    vector< ::std::fstream* > readFileStreams;
    readFileStreams.resize(length(readFNames));
    for(unsigned i = 0; options.inputFormat != 2 && i < length(readFNames); ++i)
    {
        readFileStreams[i] = new fstream(toCString(readFNames[i]), ios_base::in | ios::binary);
        if(!(*(readFileStreams[i])).is_open())
        {
            ::std::cerr << "Failed to open read file " << readFNames[i] << ::std::endl;
            return CALLSNPS_GFF_FAILED;
        }
    }
    String<RecordReader<std::fstream, SinglePass< > >* > recordReaders;

    String<Stream<Bgzf>* > bgzfStreams;
    if(options.inputFormat == 2)
    {
        resize(bgzfStreams,length(readFNames));
        for(unsigned i = 0; i < length(readFNames); ++i)
        {
            bgzfStreams[i] = new Stream<Bgzf>();
            std::cout <<"Opening bam file" << std::endl;
            if (!open(*bgzfStreams[i], toCString(readFNames[i]), "r"))
            {
                std::cerr << "[ERROR] Could not open BAM file" << readFNames[i] << std::endl;
              return 1;
            }
        }
    }

    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;


    TNameStoreCache refNameStoreCache(genomeNames);
    String<BamIOContext<TNameStore> > contexts;
    String<BamAlignmentRecord> records;
    if(options.inputFormat > 0)
    {
        if(options.inputFormat == 1) resize(recordReaders, length(readFNames));
        resize(records, length(readFNames));
        //resize(contexts, length(readFNames));
        for (unsigned i = 0; i < length(readFNames); ++i)
        {
            clear(records[i].qName);
            if(options.inputFormat == 1) recordReaders[i] = new RecordReader<std::fstream,SinglePass< > >(*readFileStreams[i]);
            //appendValue(contexts,BamIOContext<TNameStore>(refNameStore, refNameStoreCache));
            appendValue(contexts,BamIOContext<TNameStore>(genomeNames, refNameStoreCache));
        }
    }
    /////////////////////////////////////////////////////////////////////
    // open out file streams and store open file pointers
    ::std::ofstream snpFileStream; 
    if (*options.outputSNP != 0)
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
        
        snpFileStream.open(options.outputSNP,::std::ios_base::out);
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
    if (*options.outputIndel != 0)
    {
        indelFileStream.open(options.outputIndel,::std::ios_base::out);
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
    if(*options.inputPositionFile != 0)
    {
        resize(positions,length(genomeNames));
        result = loadPositions(positions,gIdStringToIdNumMap,options.inputPositionFile,options);
        if(result != 0)
        {
            ::std::cerr << "Failed to read position file " << options.inputPositionFile << ::std::endl;
            return result;
        }
        posFileStream.open(options.outputPosition,::std::ios_base::out);
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
    
    bool positionStatsOnly = (*options.outputSNP == 0 && *options.outputPosition != 0);
    TPosIterator inspectPosIt, inspectPosItEnd;

    bool firstCall = true;

    /////////////////////////////////////////////////////////////////////////////
    // Start scanning for SNPs/indels
    // for each chromosome
    for(unsigned i=0; i < length(genomes); ++i)
    {
        //std::cout << genomeNames[i] << "\n";
        if(!empty(positions))
        {
            inspectPosIt = begin(positions[i],Standard());
            inspectPosItEnd = end(positions[i],Standard());
            if(inspectPosIt == inspectPosItEnd )
                if(positionStatsOnly)
                    continue;
        }
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
        options.minCoord = MaxValue<unsigned>::VALUE;
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
            for(unsigned j = 0; j < length(readFNames); ++j)
            {
                unsigned sizeBefore = length(fragmentStore.alignedReadStore);
                
                
                // currently only gff supported
                if(options.inputFormat == 0) // GFF
                    result = readMatchesFromGFF_Batch(readFileStreams[j], fragmentStore, readCounts, readClips,
                                              readCigars, genomes[i], gIdStringToIdNumMap, 
                                              i, currentWindowBegin, currentWindowEnd, highestChrId[j], options);
                if(options.inputFormat == 1) // SAM
                    result = readMatchesFromSamBam_Batch(*recordReaders[j], contexts[j], records[j], fragmentStore, readCounts, readClips,
                                              readCigars, genomes[i], gIdStringToIdNumMap, 
                                              i, currentWindowBegin, currentWindowEnd, highestChrId[j], options, Sam(),firstCall);
               if(options.inputFormat == 2) // BAM
                    result = readMatchesFromSamBam_Batch(*bgzfStreams[j], contexts[j], records[j], fragmentStore, readCounts, readClips,
//                    result = readMatchesFromSamBam_Batch(*recordReaders[j], contexts[j], records[j], fragmentStore, readCounts, readClips,
                                              readCigars, genomes[i], gIdStringToIdNumMap, 
                                              i, currentWindowBegin, currentWindowEnd, highestChrId[j], options, Bam(),firstCall);

               firstCall = false;
                if(result == CALLSNPS_GFF_FAILED)
                {
                    ::std::cerr << "Failed to open read file " << readFNames[j] << ::std::endl;
                    ::std::cerr << "or reads are not sorted correctly. " << ::std::endl;
                    return result;
                }
                if(result > 0)
                    return result;
                
                if(options._debugLevel > 0)
                    ::std::cout << "parsed reads of file " << j << "\n";
                
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
                if (*options.outputIndel != 0)
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
                if (*options.outputSNP != 0)
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
    if (*options.outputSNP != 0)
        snpFileStream.close();
    
    if (*options.outputIndel != 0)
        indelFileStream.close();
    
    if (*options.outputPosition != 0)
        posFileStream.close();

    //  if (*options.outputCNV != 0)
    //      cnvFileStream.close();
    
    return 0;
}




// log file to keep track of happenings
template <typename TSpec>
int writeLogFile(
                 int argc, const char *argv[],
                 const char *genomeFileName,
                 String<CharString> & readFNames,
                 String<CharString> & ,
                 SNPCallingOptions<TSpec> &options)
{
    
    ::std::ofstream logfile;
    logfile.open(options.outputLog, ::std::ios_base::out | ::std::ios_base::trunc);
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
    logfile << "Genome=\"" << genomeFileName << "\""<<::std::endl;
    logfile << "Reads=\"" << readFNames[0];
    for(unsigned i = 1; i < length(readFNames); ++i)
        logfile << " " << readFNames[i] << ::std::endl;
    logfile << "\"" << std::endl;
    if(*options.outputSNP != 0)
    {
        logfile << "OutputSnp=\"" << CharString(options.outputSNP) << "\"" << ::std::endl;
    }
    if(*options.outputIndel != 0)
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
    if(*options.outputIndel != 0)
    {
        logfile << "MinIndel=" << options.indelCountThreshold << ::std::endl;
        logfile << "MinPercIndelT=" << options.indelPercentageT << ::std::endl;
    }
    logfile.close();
    return 0;
}



template <typename TOptions>
void printHelp(int, const char *[], TOptions &options, bool longHelp = false)
{
    
    cerr << "Usage: snpStore [OPTION]... <GENOME FILE> <MAPPED READ FILE>" << endl;
    if (longHelp) {
        cerr << endl << "Options:" << endl;
        cerr << "  -o,  --output FILE               \t" << "change output filename (default <READ FILE>.snp)" << endl;
        cerr << "  -if, --input-format NUM          \t" << "input format:" << endl;
        cerr << "                                   \t" << "0 = GFF sorted according to genome pos (default)" << endl;
        cerr << "                                   \t" << "1 = SAM sorted according to genome pos" << endl;
        cerr << "  -of, --output-format NUM         \t" << "output format:" << endl;
        cerr << "                                   \t" << "0 = output all candidate snps (default)" << endl;
        cerr << "                                   \t" << "1 = output succesful candidate snps only" << endl;
        cerr << "  -dc, --dont-clip                 \t" << "ignore clip tags in gff (off)" << endl;
        cerr << "  -mu, --multi                     \t" << "keep non-unique fragmentStore.alignedReadStore (off)" << endl;
        cerr << "  -hq, --hide-qualities            \t" << "only show coverage (no qualities) in SNP output file (off)" << endl;
        cerr << "  -sqo,--solexa-qual-offset        \t" << "base qualities are encoded as Ascii value - 64 (instead of Ascii - 33)" << endl;
        cerr << "  -id, --indel-file FILE           \t" << "output file for called indels in gff format (off)" << endl;
        cerr << "  -m,  --method NUM                \t" << "set method used for SNP calling" << endl;
        cerr << "                                   \t" << "0 = threshold method" << endl;
        cerr << "                                   \t" << "1 = maq (default)" << endl;
        //      cerr << "                                   \t" << "(default = "<<options.method << ")" << endl;
        cerr << "  -mp, --max-pile NUM              \t" << "maximal number of matches allowed to pile up at the same genome position ("<<options.maxPile<<")" << endl;
        cerr << "  -mmp,--merged-max-pile           \t" << "do pile up correction on merged lanes (off)" << endl;
        cerr << "  -mc, --min-coverage NUM          \t" << "minimal required number of reads covering a candidate position ("<< options.minCoverage<<")" << endl;
        cerr << "  -fc, --force-call NUM            \t" << "always call base if count is >= fc, ignore other parameters (off)" << endl;
        cerr << "  -oa, --orientation-aware         \t" << "distinguish between forward and reverse reads (off)" << endl;
        cerr << "  -fl, --force-length NUM          \t" << "read length to be used (ignores suffix of read) (off)" << endl;
        cerr << endl;
        cerr << "SNP calling options: " << endl;
        cerr << "  Threshold method related: " << endl;
        cerr << "  -mm, --min-mutations NUM         \t" << "minimal number of observed mutations for mutation to be called ("<<options.minMutT<<")" << endl;
        cerr << "  -pt, --perc-threshold NUM        \t" << "minimal percentage of mutational base for mutation to be called (" << options.percentageT << ")" << endl;
        cerr << "  -mq, --min-quality NUM           \t" << "minimal average quality of mutational base for mutation to be called ("<<options.avgQualT <<")" << endl;
        cerr << "  Maq method related: " << endl;
        cerr << "  -th, --theta NUM                 \t" << "dependency coefficient ("<< options.theta <<")" << endl;
        cerr << "  -hr, --hetero-rate NUM           \t" << "heterozygote rate ("<< options.hetRate <<")" <<  endl;
        cerr << "  -mmq,--min-map-quality NUM       \t" << "minimum base call (mapping) quality for a match to be considered ("<< options.minMapQual <<")" <<  endl;
        cerr << "  -ch, --corrected-het             \t" << "use amplification bias corrected distribution for heterozygotes (off)" << endl;
        cerr << "  -maf,--mean-alleleFreq NUM       \t" << "mean ref allele frequency in heterozygotes (" << options.meanAlleleFrequency << ")"  << endl;
        cerr << "  -ac, --amp-cycles NUM            \t" << "number of amplification cycles (" << options.amplificationCycles << ")"  << endl;
        cerr << "  -ae, --amp-efficiency NUM        \t" << "polymerase efficiency, probability of amplification (" << options.amplificationEfficiency << ")"  << endl;
        cerr << "  -in, --initialN NUM              \t" << "initial allele population size (" << options.initialN << ")" << endl;
        cerr << "  -pht,--print-hetTable            \t" << "print het table (off)" << endl;
//        cerr << "  -cq, --corrected-quality NUM     \t" << "do simple quality recalibration with x=NUM (" << options.newQualityCalibrationFactor << ")" << endl;
        cerr << "  -mec,--min-explained-col         \t" << "minimum fraction of alignment column reads explained by genotype call (" << options.minExplainedColumn << ")"  << endl;
        cerr << "Indel calling options: " << endl;
        cerr << "  -it, --indel-threshold NUM       \t" << "minimal number of indel-supporting reads required for indel calling (" << options.indelCountThreshold<<")"<< endl;
        cerr << "  -ipt,--indel-perc-threshold NUM  \t" << "minimal ratio of indel-supporting/covering reads for indel to be called (" << options.indelPercentageT<<")" << endl;
        cerr << "  -iqt,--indel-quality-thresh NUM  \t" << "minimal average quality of inserted base/deletion-neighboring bases for indel to be called (" << options.indelQualityThreshold<<")" << endl;
        cerr << "  -bsi,--both-strands-indel        \t" << "both strands need to be observed for indel to be called (off)" << endl;
//      cerr << "  -iw, --indel-window              \t" << "overlap window used for indel calling (" << options.indelWindow<<")"<< endl;
        
        cerr << endl<< "Other options: " << endl;
        cerr << "  -lf, --log-file FILE             \t" << "write log file to FILE" << endl;
        cerr << "  -v,  --verbose                   \t" << "verbose mode" << endl;
        cerr << "  -vv, --very-verbose              \t" << "very verbose mode" << endl;
        cerr << "  -h,  --help                      \t" << "print this help" << endl << endl;
    }
    else {
        cerr << "Try 'snpStore --help' for more information." << endl <<endl;
    }
}



int main(int argc, const char *argv[]) 
{
    //////////////////////////////////////////////////////////////////////////////
    // Parse command line
    
    SNPCallingOptions<>     options;
    
    unsigned                fnameCount = 0;
    const char              *genomeFName = "";
    String<CharString>      readFNames;
    String<CharString>      qualityFNames;
    
    for(int arg = 0; arg < argc; ++arg) {
        options.programCall << argv[arg] << " ";
    }
    
    /*  std::cout << "lgamma(1) = " << lgamma(1) << std::endl;
     std::cout << "lgamma(2) = " << lgamma(2) << std::endl;
     std::cout << "lgamma(3) = " << lgamma(3) << std::endl;
     std::cout << "lgamma(4) = " << lgamma(4) << std::endl;
     std::cout << "lgamma(25) = " << lgamma(25) << std::endl;
     std::cout << "lgamma(105) = " << lgamma(105) << std::endl;
     std::cout << "lgamma(255) = " << lgamma(255) << std::endl;*/
    for(int arg = 1; arg < argc; ++arg) {
        if (argv[arg][0] == '-') {
            // parse options
            if (strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--method") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.method;
                    if (!istr.fail())
                    {
                        if (options.method > 1)
                            cerr << "Invalid method option." << endl << endl;
                        else
                            continue;
                    }
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-mq") == 0 || strcmp(argv[arg], "--min-quality") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.avgQualT;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-mm") == 0 || strcmp(argv[arg], "--min-mutations") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.minMutT;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-fc") == 0 || strcmp(argv[arg], "--force-call") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.forceCallCount;
                    if (!istr.fail())
                    {
                        if(options.forceCallCount < 1)
                            cerr << "--force-call expects a positive integer." << endl; 
                        else continue;
                    }
                    
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-mpr") == 0 || strcmp(argv[arg], "--max-polymer-run") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.maxPolymerRun;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-bsi") == 0 || strcmp(argv[arg], "--both-strands-indel") == 0) {
                options.bothIndelStrands = true;
                continue;
            }
            if (strcmp(argv[arg], "-iqt") == 0 || strcmp(argv[arg], "--indel-quality-thresh") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.indelQualityThreshold;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-mec") == 0 || strcmp(argv[arg], "--min-explained-column") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.minExplainedColumn;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-ipt") == 0 || strcmp(argv[arg], "--indel-perc-threshold") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.indelPercentageT;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-pt") == 0 || strcmp(argv[arg], "--perc-threshold") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.percentageT;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-mp") == 0 || strcmp(argv[arg], "--max-pile") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.maxPile;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-dp") == 0 || strcmp(argv[arg], "--diff-pos") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.minDifferentReadPos;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-eb") == 0 || strcmp(argv[arg], "--exclude-border") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.excludeBorderPos;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-ebi") == 0 || strcmp(argv[arg], "--exclude-border-indel") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.indelDepthMinOverlap;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-mc") == 0 || strcmp(argv[arg], "--min-coverage") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.minCoverage;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-dc") == 0 || strcmp(argv[arg], "--dont-clip") == 0) {
                options.dontClip = true;
                continue;
            }
            if (strcmp(argv[arg], "-su") == 0 || strcmp(argv[arg], "--suboptimal") == 0) {
                options.keepSuboptimalReads = true;
                continue;
            }
            if (strcmp(argv[arg], "-ch") == 0 || strcmp(argv[arg], "--corrected-het") == 0) {
                options.correctedHetTable = true;
                continue;
            }
            if (strcmp(argv[arg], "-mu") == 0 || strcmp(argv[arg], "--multi") == 0) {
                options.keepMultiReads = true;
                continue;
            }
            if (strcmp(argv[arg], "-re") == 0 || strcmp(argv[arg], "--realign") == 0) {
                options.realign = true;
                continue;
            }
            if (strcmp(argv[arg], "-hq") == 0 || strcmp(argv[arg], "--hide-qualities") == 0) {
                options.showQualityStrings = false;
                continue;
            }
            if (strcmp(argv[arg], "-mmp") == 0 || strcmp(argv[arg], "--merged-max-pile") == 0) {
                options.laneSpecificMaxPile = false;
                continue;
            }
            if (strcmp(argv[arg], "-oa") == 0 || strcmp(argv[arg], "--orientation-aware") == 0) {
                options.orientationAware = true;
                continue;
            }
            if (strcmp(argv[arg], "-pht") == 0 || strcmp(argv[arg], "--print-hetTable") == 0) {
                options.printHetTable = true;
                continue;
            }
 
            if (strcmp(argv[arg], "-ae") == 0 || strcmp(argv[arg], "--amp-efficiency") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.amplificationEfficiency;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-cq") == 0 || strcmp(argv[arg], "--corrected-quality") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.newQualityCalibrationFactor;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }

            if (strcmp(argv[arg], "-maf") == 0 || strcmp(argv[arg], "--mean-alleleFreq") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.meanAlleleFrequency;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
    
            if (strcmp(argv[arg], "-ac") == 0 || strcmp(argv[arg], "--amp-cycles") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.amplificationCycles;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-in") == 0 || strcmp(argv[arg], "--initial-N") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.initialN;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }

            if (strcmp(argv[arg], "-iw") == 0 || strcmp(argv[arg], "--indel-window") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.indelWindow;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-it") == 0 || strcmp(argv[arg], "--indel-threshold") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.indelCountThreshold;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-mmq") == 0 || strcmp(argv[arg], "--min-map-quality") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.minMapQual;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-th") == 0 || strcmp(argv[arg], "--theta") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.theta;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-pws") == 0 || strcmp(argv[arg], "--parse-window-size") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.windowSize;
                    if (!istr.fail() && options.windowSize > 0)
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-sqo") == 0 || strcmp(argv[arg], "--solexa-qual-offset") == 0) {
                options.asciiQualOffset = 64;
                continue;
            }
            if (strcmp(argv[arg], "-reb") == 0 || strcmp(argv[arg], "--realign-border") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.realignAddBorder;
                    if (!istr.fail() && options.realignAddBorder >= 0 && options.realignAddBorder <= 100)
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-cws") == 0 || strcmp(argv[arg], "--cnv-window-size") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.cnvWindowSize;
                    if (!istr.fail() && options.cnvWindowSize > 0)
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-hr") == 0 || strcmp(argv[arg], "--hetero-rate") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.hetRate;
                    if (!istr.fail())
                        continue;
            }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-if") == 0 || strcmp(argv[arg], "--input-format") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.inputFormat;
                    if (!istr.fail())
                    {
                        if(options.inputFormat > 2 )
                            cerr << "--input-format expects 0 or 1." << endl;   
                        else continue;
                    }
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-of") == 0 || strcmp(argv[arg], "--output-format") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    istringstream istr(argv[arg]);
                    istr >> options.outputFormat;
                    if (!istr.fail())
                        continue;
                }
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-lf") == 0 || strcmp(argv[arg], "--log-file") == 0) {
                if (arg + 1 == argc) {
                    printHelp(argc, argv, options, true);
                    return 0;
                }
                ++arg;
                options.outputLog = argv[arg];
                continue;
            }
            if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
                if (arg + 1 == argc) {
                    printHelp(argc, argv, options, true);
                    return 0;
                }
                ++arg;
                options.outputSNP = argv[arg];
                continue;
            }
            if (strcmp(argv[arg], "-id") == 0 || strcmp(argv[arg], "--indel-file") == 0) {
                if (arg + 1 == argc) {
                    printHelp(argc, argv, options, true);
                    return 0;
                }
                ++arg;
                options.outputIndel = argv[arg];
                continue;
            }
            if (strcmp(argv[arg], "-cnv") == 0 || strcmp(argv[arg], "--output-cnv") == 0) {
                if (arg + 1 == argc) {
                    printHelp(argc, argv, options, true);
                    return 0;
                }
                ++arg;
                options.outputCNV = argv[arg];
                continue;
            }
            if (strcmp(argv[arg], "-op") == 0 || strcmp(argv[arg], "--output-positions") == 0) {
                if (arg + 1 == argc) {
                    printHelp(argc, argv, options, true);
                    return 0;
                }
                ++arg;
                options.outputPosition = argv[arg];
                continue;
            }
            if (strcmp(argv[arg], "-ip") == 0 || strcmp(argv[arg], "--input-positions") == 0) {
                if (arg + 1 == argc) {
                    printHelp(argc, argv, options, true);
                    return 0;
                }
                ++arg;
                options.inputPositionFile = argv[arg];
                continue;
            }
            if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
                // print help
                printHelp(argc, argv, options, true);
                return 0;
            }
            if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
                options._debugLevel = max(options._debugLevel, 1);
                continue;
            }
            if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--very-verbose") == 0) {
                options._debugLevel = max(options._debugLevel, 2);
                continue;
            }
            cerr << "Unknown option: " << argv[arg] << endl << endl;
            printHelp(argc, argv, options);
            return 0;
        } else {
            // parse file name
            if (fnameCount == 0)
                genomeFName = argv[arg];
            if (fnameCount == 1)
            {
                if(argv[arg][0]=='[')
                {
                    String<char> tempStr = argv[arg];
                    appendValue(readFNames,suffix(tempStr,1),Generous());
                    ++arg;
                    while(arg < argc && argv[arg][0] != '-')
                    {
                        tempStr = argv[arg];
                        appendValue(readFNames,tempStr,Generous());
                        ++arg;
                    }
                    if(readFNames[length(readFNames)-1][length(readFNames[length(readFNames)-1])-1] != ']' )
                        cerr << "Something wrong with read file list?" << endl;
                    resize(readFNames[length(readFNames)-1],length(readFNames[length(readFNames)-1])-1);
                    --arg;
                }
                else
                {
                    //split by whitesapce and append each read file
                    ::std::string tempStr(argv[arg]);
                    size_t firstPos = 0;
                    size_t lastPos = tempStr.find(' ');
                    ::std::string tempFile = tempStr.substr(firstPos,lastPos);
                    appendValue(readFNames,String<char>(tempFile),Generous());
                    while (lastPos != 0 && lastPos != tempStr.npos)
                    {
                        while(tempStr[lastPos]==' ')
                            ++lastPos;
                        firstPos = lastPos; 
                        lastPos = tempStr.find(' ',firstPos);
                        if (lastPos != tempStr.npos) tempFile = tempStr.substr(firstPos,lastPos-firstPos);
                        else tempFile = tempStr.substr(firstPos,length(tempStr));
                        appendValue(readFNames,String<char>(tempFile),Generous());
                    }
                }
            }
            if (fnameCount == 2) {
                cerr << "More than 2 input files specified." <<endl;
                cerr << "If more than 2 mapped read files are to be parsed, use quotation marks directly before first file name and directly after last file name (e.g. \"lane1.gff lane2.gff\")." << endl << endl;
                printHelp(argc, argv, options);
                return 0;
            }
            ++fnameCount;
        }
    }
    
    // some option checking
    if (fnameCount != 2) {
        if (argc > 1 && !options.printVersion)
            cerr << "Exactly 2 input files need to be specified." << endl << endl;
        printHelp(argc, argv, options);
        return 0;
    }
    //if(options.inputFormat == 1 && (!options.qualityFile || (length(qualityFNames)!=length(readFNames))))
    //{
    //  cerr << "If mapped read file is in Eland format, a .qual or .fastq file containing read qualities needs to be specified." << endl << endl;
    //  return 0;
    //}
    if(*options.inputPositionFile == 0 && *options.outputPosition != 0)
    {
        cerr << "Position analysis output specified, but no position file given." << endl << endl;
        return 0;
    }
    
    if((options.realign && options.windowSize > 50000) || options.windowSize > 1000000) 
        options.windowSize = 10000;
    
    if (*options.outputLog != 0)
        writeLogFile(argc, argv, genomeFName, readFNames, qualityFNames, options);
    
    if(options.runID == "")
    {
        ::std::string tempStr = toCString(readFNames[0]);
        size_t lastPos = tempStr.find_last_of('/') + 1;
        if (lastPos == tempStr.npos) lastPos = tempStr.find_last_of('\\') + 1;
        if (lastPos == tempStr.npos) lastPos = 0;
        options.runID = tempStr.substr(lastPos);
    }
    
    
    
    //////////////////////////////////////////////////////////////////////////////
    // check for variants
    int result = detectSNPs(genomeFName, readFNames, qualityFNames, options);
    if (result > 0)
    {
        printHelp(argc, argv, options);
        return 0;
    }
    return result;
}
