 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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

#ifndef SEQAN_HEADER_RAZERS_MATEPAIRS_H
#define SEQAN_HEADER_RAZERS_MATEPAIRS_H

#include <seqan/misc/misc_dequeue.h>

namespace SEQAN_NAMESPACE_MAIN
{

// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.


//////////////////////////////////////////////////////////////////////////////
// Definitions

typedef StringSet<TRead const, Dependent<> >	TMPReadSet;

#ifdef RAZERS_MEMOPT

	template <typename TShape>
	struct SAValue< Index<TMPReadSet, IndexQGram<TShape> > > {
		typedef Pair<
			unsigned,				
			unsigned,
			BitPacked<24, 8>	// max. 16M reads of length < 256
		> Type;
	};
	
#else

	template <typename TShape>
	struct SAValue< Index<TMPReadSet, IndexQGram<TShape> > > {
		typedef Pair<
			unsigned,			// many reads
			unsigned,			// of arbitrary length
			Pack
		> Type;
	};

#endif

	
template <typename TShape, typename TSpec>
struct Cargo< Index<TMPReadSet, IndexQGram<TShape, TSpec> > > {
	typedef struct {
		double		abundanceCut;
		int			_debugLevel;
	} Type;
};

#ifdef RAZERS_PRUNE_QGRAM_INDEX

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TShape>
inline bool _qgramDisableBuckets(Index<TMPReadSet, IndexQGram<TShape> > &index) 
{
	typedef Index<TMPReadSet, IndexQGram<TShape>	>	TReadIndex;
	typedef typename Fibre<TReadIndex, QGramDir>::Type	TDir;
	typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
	typedef typename Value<TDir>::Type					TSize;
	
	TDir &dir    = indexDir(index);
	bool result  = false;
	unsigned counter = 0;
	TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
	if (thresh < 100) thresh = 100;
	
	TDirIterator it = begin(dir, Standard());
	TDirIterator itEnd = end(dir, Standard());
	for (; it != itEnd; ++it)
		if (*it > thresh) 
		{
			*it = (TSize)-1;
			result = true;
			++counter;
		}
	
	if (counter > 0 && cargo(index)._debugLevel >= 1)
		::std::cerr << "Removed " << counter << " k-mers" << ::std::endl;
	
	return result;
}

#endif

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TFSSpec, typename TFSConfig, typename TRazerSOptions>
bool loadReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	const char							* fileNameL,		// left mates file
	const char							* fileNameR,		// right mates file
	TRazerSOptions						& options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);

	MultiFasta leftMates;
	MultiFasta rightMates;

	if (!open(leftMates.concat, fileNameL, OPEN_RDONLY)) return false;
	if (!open(rightMates.concat, fileNameR, OPEN_RDONLY)) return false;

	AutoSeqFormat formatL;
	guessFormat(leftMates.concat, formatL);
	split(leftMates, formatL);

	AutoSeqFormat formatR;
	guessFormat(rightMates.concat, formatR);
	split(rightMates, formatR);

	unsigned seqCount = length(leftMates);
	if (seqCount != length(rightMates))
	if (options._debugLevel > 1) 
	{
		::std::cerr << "Numbers of mates differ: " << seqCount << "(left) != " << length(rightMates) << "(right).\n";
		return false;
	}

	String<Dna5Q>	seq[2];
	CharString		qual[2];
	CharString		id[2];
	
	unsigned kickoutcount = 0;
	unsigned maxReadLength = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0 || options.readNaming == 3)
		{
			assignSeqId(id[0], leftMates[i], formatL);              // read left Fasta id
			assignSeqId(id[1], rightMates[i], formatR);             // read right Fasta id
			if (options.readNaming == 0)
			{
				append(id[0], "/L");
				append(id[1], "/R");
			}
		}
		
		assignSeq(seq[0], leftMates[i], formatL);                   // read left Read sequence
		assignSeq(seq[1], rightMates[i], formatR);                  // read right Read sequence
		assignQual(qual[0], leftMates[i], formatL);                 // read left ascii quality values  
		assignQual(qual[1], rightMates[i], formatR);                // read right ascii quality values  
		
		if (countN)
		{
			for (int j = 0; j < 2; ++j)
			{
				int maxBase = (int)(0.8 * length(seq[j]));
				int allowed[5] = 
					{ maxBase, maxBase, maxBase, maxBase, (int)(options.errorRate * length(seq[j]))};
				for (unsigned k = 0; k < length(seq[j]); ++k)
					if (--allowed[ordValue(getValue(seq[j], k))] == 0)
					{
//						std::cout << "Ignoring mate-pair: " << seq[0] << " " << seq[1] << std::endl;
						clear(seq[0]);
						clear(seq[1]);
						clear(id[0]);
						clear(id[1]);
						++kickoutcount;
						break;
					}
			}
		}
		
		for (int j = 0; j < 2; ++j)
		{
			// store dna and quality together
			assignQualities(seq[j], qual[j]);
			
			if (options.trimLength > 0 && length(seq[j]) > (unsigned)options.trimLength)
				resize(seq[j], options.trimLength);
		}
		appendMatePair(store, seq[0], seq[1], id[0], id[1]);
		if (maxReadLength < length(seq[0]))
			maxReadLength = length(seq[0]);
		if (maxReadLength < length(seq[1]))
			maxReadLength = length(seq[1]);
	}
	// memory optimization
	reserve(store.readSeqStore.concat, length(store.readSeqStore.concat), Exact());
//	reserve(store.readNameStore.concat, length(store.readNameStore.concat), Exact());

	typedef Shape<Dna, SimpleShape> TShape;
	typedef typename SAValue< Index<TMPReadSet, IndexQGram<TShape, OpenAddressing> > >::Type TSAValue;
	TSAValue sa(0, 0);
	sa.i1 = ~sa.i1;
	sa.i2 = ~sa.i2;
	
	if ((unsigned)sa.i1 < length(store.readSeqStore) - 1)
	{
		::std::cerr << "Maximal read number of " << (unsigned)sa.i1 + 1 << " exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}
	if ((unsigned)sa.i2 < maxReadLength - 1)
	{
		::std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality mate-pairs.\n";
	return (seqCount > 0);
}

	template <typename TFragmentStore>
	struct LessPairScore : 
		public ::std::binary_function <
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type, 
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type,
			bool >
	{
		TFragmentStore &store;
		
		LessPairScore(TFragmentStore &_store):
			store(_store) {}
		
		inline bool operator() (
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type const &a, 
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type const &b) const 
		{
			typedef typename TFragmentStore::TReadStore			TReadStore;
			typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;
			typedef typename TFragmentStore::TAlignQualityStore	TAlignQualityStore;
			typedef typename Value<TReadStore>::Type			TRead;
			typedef typename Value<TAlignedReadStore>::Type		TAlignedRead;
			typedef typename Value<TAlignQualityStore>::Type	TQual;
			typedef typename Id<TRead>::Type					TId;

			// pair number
			TRead const &ra = store.readStore[a.readId];
			TRead const &rb = store.readStore[b.readId];
			if (ra.matePairId < rb.matePairId) return true;
			if (ra.matePairId > rb.matePairId) return false;

			// quality
			SEQAN_ASSERT_NEQ(a.id, TAlignedRead::INVALID_ID);
			SEQAN_ASSERT_NEQ(b.id, TAlignedRead::INVALID_ID);
			TQual const &qa = store.alignQualityStore[a.id];
			TQual const &qb = store.alignQualityStore[b.id];
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			
			return a.pairMatchId < b.pairMatchId;
		}
	};
	

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TFragmentStore, typename TSpec, typename TSwiftL, typename TSwiftR >
void compactPairMatches(
	TFragmentStore			& store,
	RazerSOptions<TSpec>	& options,
	TSwiftL					& swiftL, 
	TSwiftR					& swiftR)
{
	typedef typename TFragmentStore::TReadStore						TReadStore;
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	typedef typename Value<TReadStore>::Type						TRead;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Value<TAlignQualityStore>::Type				TQual;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;

	SEQAN_ASSERT_EQ(length(store.alignedReadStore) % 2, 0u);
	
	unsigned matePairId = -2;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreDistCutOff = MinValue<int>::VALUE;
    int scoreRangeBest = (options.distanceRange == 0u)? MinValue<int>::VALUE: -(int)options.distanceRange;

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	// sort 
	sortAlignedReads(store.alignedReadStore, LessPairScore<TFragmentStore>(store));

	for (; it != itEnd; ++it) 
	{
		SEQAN_ASSERT_NEQ((*it).id, TAlignedRead::INVALID_ID);
		SEQAN_ASSERT_NEQ((*it).readId, TAlignedRead::INVALID_ID);
        SEQAN_ASSERT_EQ(it->pairMatchId, (it + 1)->pairMatchId);
		
		// ignore pair alignments if one of the mates is marked as deleted (<=> contigId is invalid)
		if (it->contigId == TAlignedRead::INVALID_ID || (it + 1)->contigId == TAlignedRead::INVALID_ID)
		{
			++it;
			continue;
		}
		TRead &r = store.readStore[(*it).readId];
		TQual &q = store.alignQualityStore[(*it).id];
		if (matePairId == r.matePairId)
		{ 
			if (q.pairScore <= scoreDistCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = -1 - q.pairScore;
					if (options.purgeAmbiguous && q.pairScore > scoreRangeBest)
						maxErrors = -1;

					setMaxErrors(swiftL, matePairId, maxErrors);
					setMaxErrors(swiftR, matePairId, maxErrors);

					if (maxErrors == -1 && options._debugLevel >= 2)
						::std::cerr << "(pair #" << matePairId << " disabled)";

					if (options.purgeAmbiguous)
					{
						if (q.pairScore > scoreRangeBest || IsSameType<TSwiftL, Nothing>::VALUE)
							dit = ditBeg;
						else {
							*dit = *it;	++dit; ++it;
							*dit = *it; ++dit;
							continue;
						}
					}
				}
#endif
				++it;
				continue;
			}
		}
		else
		{
			matePairId = r.matePairId;
			hitCount = 0;
			if (options.distanceRange > 0)
				scoreDistCutOff = q.pairScore - options.distanceRange;
			ditBeg = dit;
		}
		*dit = *it;	++dit; ++it;
		*dit = *it;	++dit;
	}
	resize(store.alignedReadStore, dit - begin(store.alignedReadStore, Standard()));
	compactAlignedReads(store);
}



#ifndef RAZERS_PARALLEL
//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TPreprocessing,
	typename TRazerSOptions >
void mapMatePairReads(
	FragmentStore<TFSSpec, TFSConfig>		& store,
	TGenome									& genome,				// genome ...
	unsigned								  contigId,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPatternL,
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPatternR,
	TPreprocessing							& preprocessingL,
	TPreprocessing							& preprocessingR,
	char									  orientation,			// q-gram index of reads
	TRazerSOptions							& options)
{
	typedef FragmentStore<TFSSpec, TFSConfig>				TFragmentStore;
	typedef typename TFragmentStore::TMatePairStore			TMatePairStore;
	typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
	typedef typename Value<TMatePairStore>::Type			TMatePair;
	typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
	typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	typedef typename Id<TAlignedRead>::Type					TId;

	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename MakeSigned_<TGPos>::Type				TSignedGPos;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >			TSwiftFinderR;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// MATE-PAIR FILTRATION
	typedef Triple<__int64, TAlignedRead, TAlignQuality>	TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;

	// VERIFICATION
	typedef MatchVerifier <
		TFragmentStore, 
		TRazerSOptions, 
		TPreprocessing, 
		TSwiftPattern >										TVerifier;

	const unsigned NOT_VERIFIED = 1u << (8*sizeof(unsigned)-1);

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	TReadSet	&readSetL = host(host(swiftPatternL));
	TReadSet	&readSetR = host(host(swiftPatternR));
	TVerifier	verifierL(store, options, preprocessingL, swiftPatternL);
	TVerifier	verifierR(store, options, preprocessingR, swiftPatternR);

	verifierL.oneMatchPerBucket = true;
	verifierR.oneMatchPerBucket = true;
	verifierL.m.contigId = contigId;
	verifierR.m.contigId = contigId;

	if (empty(readSetL))
		return;

	// distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
	// distance >= libLen - libErr - 2*parWidth + shapeLen
	TSize readLength = length(readSetL[0]);
	TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(swiftPatternL)));
	TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(swiftPatternL)));
	TGPos scanShift = (minDistance < 0)? 0: minDistance;

	// exit if contig is shorter than library size
	if (length(genome) <= scanShift)
		return;

	TGenomeInf genomeInf = infix(genome, scanShift, length(genome));
	TSwiftFinderL swiftFinderL(genome, options.repeatLength, 1);
	TSwiftFinderR swiftFinderR(genomeInf, options.repeatLength, 1);

	TDequeue fifo;						// stores left-mate potential matches
	String<__int64> lastPotMatchNo;		// last number of a left-mate potential
	__int64 lastNo = 0;					// last number over all left-mate pot. matches in the queue
	__int64 firstNo = 0;				// first number over all left-mate pot. match in the queue
	Pair<TGPos> gPair;

	resize(lastPotMatchNo, length(host(swiftPatternL)), (__int64)-2, Exact());

	TSize gLength = length(genome);

	TAlignedRead mR;
	TAlignQuality qR;
	TDequeueValue fL(-2, mR, qR);	// to supress uninitialized warnings
	
	// iterate all verification regions returned by SWIFT
	while (find(swiftFinderR, swiftPatternR, options.errorRate)) 
	{
		++options.countFiltration;
		unsigned matePairId = swiftPatternR.curSeqNo;
		TGPos rEndPos = endPosition(swiftFinderR) + scanShift;
		TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;

		// (1) Remove out-of-window left mates from fifo.
		while (!empty(fifo) && (TSignedGPos)front(fifo).i2.endPos + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
		{
///            std::cerr << "  -Left [" << front(fifo).i2.beginPos << "\t" << front(fifo).i2.endPos << ')' << std::endl;
			popFront(fifo);
			++firstNo;
		}

        // (2) Add within-window left mates to fifo.
		while (empty(fifo) || (TSignedGPos)back(fifo).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
		{
			if (find(swiftFinderL, swiftPatternL, options.errorRate))
			{
				++options.countFiltration;
				gPair = positionRange(swiftFinderL);
				if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
				{
					// link in
					fL.i1 = lastPotMatchNo[swiftPatternL.curSeqNo];
					lastPotMatchNo[swiftPatternL.curSeqNo] = lastNo++;

					fL.i2.readId = store.matePairStore[swiftPatternL.curSeqNo].readId[0] | NOT_VERIFIED;
					fL.i2.beginPos = gPair.i1;
					fL.i2.endPos = gPair.i2;
					
///            std::cerr << "  +Left \t" << firstNo + length(fifo) << ":\t[" << fL.i2.beginPos << "\t" << fL.i2.endPos << ')' << std::endl;
					pushBack(fifo, fL);
				}
			} else {
				break;
			}
		}

		int	bestLeftScore = MinValue<int>::VALUE;
		int bestLibSizeError = MaxValue<int>::VALUE;
		TDequeueIterator bestLeft = TDequeueIterator();

		bool rightVerified = false;
		TDequeueIterator it;
		unsigned leftReadId = store.matePairStore[matePairId].readId[0];
		__int64 last = (__int64)-1;
		__int64 lastValid = (__int64)-1;
		__int64 i;
		for (i = lastPotMatchNo[matePairId]; firstNo <= i; last = i, i = (*it).i1)
		{
///            std::cout<< "\t[" << i << "]" << "\t" << fifo[3].i1 << std::endl;
			it = &value(fifo, i - firstNo);

			// search left mate
			{
				// verify left mate (equal seqNo), if not done already
				if ((*it).i2.readId & NOT_VERIFIED)
				{
                    if ((TSignedGPos)(*it).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
					{
                        ++options.countVerification;
///                        if (i==0)
///                        std::cout<<"here"<<std::endl;
						if (matchVerify(verifierL, infix(genome, (TSignedGPos)(*it).i2.beginPos, (TSignedGPos)(*it).i2.endPos), 
								matePairId, readSetL, TSwiftSpec()))
						{
///                            std::cerr << "  Left+ " << verifierL.m.endPos << std::endl;
							verifierL.m.readId = (*it).i2.readId & ~NOT_VERIFIED;		// has been verified positively
							(*it).i2 = verifierL.m;
							(*it).i3 = verifierL.q;
						} else {
							(*it).i2.readId = ~NOT_VERIFIED;							// has been verified negatively 
							continue;													// we intentionally do not set lastPositive to i 
						}																// to remove i from linked list 
					} else { 
						lastValid = i; 
						continue;														// left pot. hit is out of tolerance window 
					} 
				} //else {}																// left match is verified already 
 
				// short-cut negative matches 
				if (last != lastValid) 
				{ 
					SEQAN_ASSERT_NEQ(lastValid, i); 
					if (lastValid == (__int64)-1) 
						lastPotMatchNo[matePairId] = i; 
					else 
						value(fifo, lastValid - firstNo).i1 = i;
				}
				lastValid = i;
				
                if (!rightVerified)														// here a verfied left match is available
                {
                    ++options.countVerification;
					if (matchVerify(verifierR, infix(swiftFinderR, genomeInf), 
							matePairId, readSetR, TSwiftSpec()))
                    {
                        rightVerified = true;
                        mR = verifierR.m;
						qR = verifierR.q;
                    } else {
                        // Break out of lastPotMatch loop, rest of find(right SWIFT results loop will not
                        // be executed since bestLeftScore remains untouched.
                        i = (*it).i1;
                        break;
                    }
                }

				if ((*it).i2.readId == leftReadId)
				{
					int score = (*it).i3.score;
					if (bestLeftScore <= score)
					{
                        // distance between left mate beginning and right mate end
                        __int64 dist = (__int64)verifierR.m.endPos - (__int64)(*it).i2.beginPos;
						int libSizeError = options.libraryLength - dist;
                        
///                        if (orientation == 'F')
///                            std::cout << (__int64)(*it).i2.beginPos << "\t" << (__int64)verifierR.m.beginPos;
///                        else
///                            std::cout << (__int64)(*it).i2.endPos << "\t" << (__int64)verifierR.m.endPos;
///                        std::cout << '\t' << dist << '\t' << libSizeError << std::endl;

						if (libSizeError < 0)
							libSizeError = -libSizeError;
                        if (libSizeError > options.libraryError)
                            continue;
						if (bestLeftScore == score)
						{
							if (bestLibSizeError > libSizeError)
							{
								bestLibSizeError = libSizeError;
								bestLeft = it;
							}
						}
						else
						{
							bestLeftScore = score;
							bestLibSizeError = libSizeError;
							bestLeft = it;
//							if (bestLeftScore == 0) break;	// TODO: replace if we have real qualities
						}
					}
				}
			}
		}

        // (3) Short-cut negative matches.
        if (last != lastValid)
        {
            SEQAN_ASSERT_NEQ(lastValid, i);
            if (lastValid == (__int64)-1)
                lastPotMatchNo[matePairId] = i;
            else
                value(fifo, lastValid - firstNo).i1 = i;
        }
		
		// verify right mate, if left mate matches
		if (bestLeftScore != MinValue<int>::VALUE)
		{
			fL.i2 = (*bestLeft).i2;
			fL.i3 = (*bestLeft).i3;

			// transform mate readNo to global readNo
			TMatePair &mp = store.matePairStore[matePairId];
			fL.i2.readId = mp.readId[0];
			mR.readId    = mp.readId[1];

			// transform coordinates to the forward strand
			if (orientation == 'F')
			{
				TSize temp = mR.beginPos;
				mR.beginPos = mR.endPos;
				mR.endPos = temp;
			} else 
			{
				fL.i2.beginPos = gLength - fL.i2.beginPos;
				fL.i2.endPos = gLength - fL.i2.endPos;
				TSize temp = mR.beginPos;
				mR.beginPos = gLength - mR.endPos;
				mR.endPos = gLength - temp;
//					dist = -dist;
			}
			
			// set a unique pair id
			fL.i2.pairMatchId = mR.pairMatchId = options.nextPairMatchId;
			if (++options.nextPairMatchId == TAlignedRead::INVALID_ID)
				options.nextPairMatchId = 0;

			// score the whole match pair
			fL.i3.pairScore = qR.pairScore = fL.i3.score + qR.score;

			// both mates match with correct library size
			if (!options.spec.DONT_DUMP_RESULTS)
			{
				fL.i2.id = length(store.alignedReadStore);
				appendValue(store.alignedReadStore, fL.i2, Generous());
				appendValue(store.alignQualityStore, fL.i3, Generous());
				mR.id = length(store.alignedReadStore);
				appendValue(store.alignedReadStore, mR, Generous());
				appendValue(store.alignQualityStore, qR, Generous());

				if (length(store.alignedReadStore) > options.compactThresh)
				{
					typename Size<TAlignedReadStore>::Type oldSize = length(store.alignedReadStore);
					maskDuplicates(store);	// overlapping parallelograms cause duplicates
					compactPairMatches(store, options, swiftPatternL, swiftPatternR);
					
					if (length(store.alignedReadStore) * 4 > oldSize)			// the threshold should not be raised
						options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed
					
					if (options._debugLevel >= 2)
						::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
				}
			}
		} 
	}
}
#endif


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TGNoToFile, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapMatePairReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	StringSet<CharString>				& genomeFileNameList,
	String<TGNoToFile>					& gnoToFileMap,
	RazerSOptions<TSpec>				& options,
	TShape const						& shape,
	Swift<TSwiftSpec> const)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	
	typedef typename Value<TReadSeqStore>::Type			TRead;
	typedef StringSet<TRead>							TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

//	std::cout << "SA-TYPE:" <<sizeof(typename SAValue<TIndex>::Type)<<std::endl;

	// split mate-pairs over two indices
	TReadSet readSetL, readSetR;
	unsigned pairCount = length(store.matePairStore);
	resize(readSetL, pairCount, Exact());
	resize(readSetR, pairCount, Exact());

	for (unsigned i = 0; i < pairCount; ++i)
	{
		assign(readSetL[i], store.readSeqStore[store.matePairStore[i].readId[0]]);
		assign(readSetR[i], store.readSeqStore[store.matePairStore[i].readId[1]]);
	}
	reverseComplement(readSetR);

	// configure q-gram index
	TIndex swiftIndexL(readSetL, shape);
	TIndex swiftIndexR(readSetR, shape);
	reverse(indexShape(swiftIndexR));		// right mate qualities are reversed -> reverse right shape
	
	cargo(swiftIndexL).abundanceCut = options.abundanceCut;
	cargo(swiftIndexR).abundanceCut = options.abundanceCut;
	cargo(swiftIndexL)._debugLevel = options._debugLevel;
	cargo(swiftIndexR)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPatternL(swiftIndexL);
	TSwiftPattern swiftPatternR(swiftIndexR);
	swiftPatternL.params.minThreshold = options.threshold;
	swiftPatternR.params.minThreshold = options.threshold;
	swiftPatternL.params.tabooLength = options.tabooLength;
	swiftPatternR.params.tabooLength = options.tabooLength;
	swiftPatternL.params.printDots = false;
	swiftPatternR.params.printDots = options._debugLevel > 0;

	// init edit distance verifiers
	String<TMyersPattern> forwardPatternsL;
	String<TMyersPattern> forwardPatternsR;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		resize(forwardPatternsL, pairCount, Exact());
		resize(forwardPatternsR, pairCount, Exact());
		for(unsigned i = 0; i < pairCount; ++i)
		{
#ifdef RAZERS_NOOUTERREADGAPS
			if (!empty(readSetL[i]) && !empty(readSetR[i])) {
				setHost(forwardPatternsL[i], prefix(readSetL[i], length(readSetL[i]) - 1));
				setHost(forwardPatternsR[i], prefix(readSetR[i], length(readSetR[i]) - 1));
			}
#else
			setHost(forwardPatternsL[i], readSetL[i]);
			setHost(forwardPatternsR[i], readSetR[i]);
#endif
			_patternMatchNOfPattern(forwardPatternsL[i], options.matchN);
			_patternMatchNOfPattern(forwardPatternsR[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsL[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsR[i], options.matchN);
		}
	}

	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	unsigned filecount = 0;
	unsigned numFiles = length(genomeFileNameList);
	unsigned contigId = 0;

	// open genome files, one by one	
	while (filecount < numFiles)
	{
		// open genome file	
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;

		// remove the directory prefix of current genome file
		::std::string genomeFile(toCString(genomeFileNameList[filecount]));
		size_t lastPos = genomeFile.find_last_of('/') + 1;
		if (lastPos == genomeFile.npos) lastPos = genomeFile.find_last_of('\\') + 1;
		if (lastPos == genomeFile.npos) lastPos = 0;
		::std::string genomeName = genomeFile.substr(lastPos);
		
		CharString	id;
		Dna5String	genome;
		unsigned gseqNoWithinFile = 0;
		// iterate over genome sequences
		SEQAN_PROTIMESTART(find_time);
		for(; !_streamEOF(file); ++contigId)
		{
			if (options.genomeNaming == 0)
			{
				//readID(file, id, Fasta());			// read Fasta id
				readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
				appendValue(store.contigNameStore, id, Generous());
			}
			read(file, genome, Fasta());			// read Fasta sequence
			
			appendValue(gnoToFileMap, TGNoToFile(genomeName, gseqNoWithinFile));
			
			if (options.forward)
				mapMatePairReads(store, genome, contigId, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, 'F', options);

			if (options.reverse)
			{
				reverseComplement(genome);
				mapMatePairReads(store, genome, contigId, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
		file.close();
		++filecount;
	}

	reverseComplement(readSetR);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
	}
	return 0;
}


}

#endif
