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

namespace seqan
{

// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.



template <typename TValue, typename TSpec = Alloc<> >
class Dequeue
{
public:
	typedef String<TValue, TSpec>						TString;
	typedef typename Iterator<TString, Standard>::Type	TIter;

	String<TValue, TSpec> data_string;

	TIter data_begin;	// string beginning
	TIter data_end;		// string end

	TIter data_front;	// front fifo character
	TIter data_back;	// back fifo character
	bool data_empty;	// fifo is empty

//____________________________________________________________________________

public:
	inline Dequeue()
	{
		clear(*this);
	}
};

//////////////////////////////////////////////////////////////////////////////
// Iterators
//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Standard> 
{
	typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Standard> 
{
	typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Rooted> 
{
	typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Rooted> 
{
	typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
empty(Dequeue<TValue, TSpec> const &me)
{
	return me.data_empty;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(Dequeue<TValue, TSpec> &me)
{
	clear(me.data_string);
	me.data_begin = begin(me.data_string, Standard());
	me.data_end = end(me.data_string, Standard());

	me.data_front = me.data_back = me.data_begin;
	me.data_empty = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPos>
inline TValue &
value(Dequeue<TValue, TSpec> &me, TPos pos)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
	TSize wrap = length(me.data_string) - (me.data_front - me.data_begin);
	
	if ((TSize)pos < wrap)
		return value(me.data_front + pos);
	else
		return value(me.data_begin + (pos - wrap));
}

template <typename TValue, typename TSpec, typename TPos>
inline TValue const &
value(Dequeue<TValue, TSpec> const &me, TPos pos)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
	TSize wrap = length(me.data_string) - (me.data_front - me.data_begin);
	
	if ((TSize)pos < wrap)
		return value(me.data_front + pos);
	else
		return value(me.data_begin + (pos - wrap));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
front(Dequeue<TValue, TSpec> &me)
{
	return *me.data_front;
}

template <typename TValue, typename TSpec>
inline TValue const &
front(Dequeue<TValue, TSpec> const &me)
{
	return *me.data_front;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
back(Dequeue<TValue, TSpec> &me)
{
	return *me.data_back;
}

template <typename TValue, typename TSpec>
inline TValue const &
back(Dequeue<TValue, TSpec> const &me)
{
	return *me.data_back;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
popFront(Dequeue<TValue, TSpec> &me)
{
	if (me.data_empty) return false;

	if (me.data_front == me.data_back)
		me.data_empty = true;
	else
	{
		if (++me.data_front == me.data_end)
			me.data_front = me.data_begin;
	}

	return true;
}

template <typename TValue, typename TSpec>
inline bool
popBack(Dequeue<TValue, TSpec> &me)
{
	if (me.data_empty) return false;

	if (me.data_front == me.data_back)
		me.data_empty = true;
	else
	{
		if (me.data_back == me.data_begin)
			me.data_back = me.data_end;
		--me.data_back;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
pushFront(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	typedef typename Dequeue<TValue, TSpec>::TIter TIter;

	if (me.data_empty) 
	{
		if (me.data_begin == me.data_end)
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
		me.data_empty = false;
	}
	else 
	{
		TIter new_front = me.data_front;
		if (new_front == me.data_begin)
			new_front = me.data_end;
		--new_front;

		if (new_front == me.data_back)
		{
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());

			if (me.data_front == me.data_begin)
				me.data_front = me.data_end;
			--me.data_front;
		} else
			me.data_front = new_front;
	}
	assign(*me.data_front, _value);
}

template <typename TValue, typename TSpec>
inline void
pushBack(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
	typedef typename Dequeue<TValue, TSpec>::TIter TIter;

	if (me.data_empty) 
	{
		if (me.data_begin == me.data_end)
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
		me.data_empty = false;
	}
	else 
	{
		TIter new_back = me.data_back;
		if (++new_back == me.data_end)
			new_back = me.data_begin;

		if (new_back == me.data_front)
		{
			reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
			// in this case reserve adds new space behind data_back
			++me.data_back;
		} else
			me.data_back = new_back;
	}
	assign(*me.data_back, _value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<Dequeue<TValue, TSpec> >::Type
length(Dequeue<TValue, TSpec> const &me)
{
	if (empty(me)) return 0;

	if (me.data_front <= me.data_back)
		return (me.data_back - me.data_front) + 1;
	else
		return (me.data_end - me.data_begin) - (me.data_front - me.data_back) + 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size<Dequeue<TValue, TSpec> >::Type
reserve(Dequeue<TValue, TSpec> &me, TSize_ new_capacity, Tag<TExpand> tag)
{
	typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
//	::std::cout << "resize to "<<new_capacity<<::std::endl;
	TSize len = length(me);
	if (len < new_capacity && length(me.data_string) != new_capacity)
	{
		TSize pos_front = me.data_front - me.data_begin;
		TSize pos_back  = me.data_back  - me.data_begin;
		TSize new_freeSpace = new_capacity - len;

		if (pos_front <= pos_back)
		{
			// |empty|data|empty|
			// 0
			TSize freeSpace = length(me.data_string) - len;
			if (new_freeSpace > freeSpace)
				resize(me.data_string, new_capacity, tag);
			else
			{
				freeSpace -= new_freeSpace;	// reduce the free space by <freeSpace>
				if (pos_front >= freeSpace)
				{
					resizeSpace(me.data_string, pos_front - freeSpace, (TSize)0, pos_front, tag);
					pos_back -= freeSpace;
					pos_front -= freeSpace;
				}
				else
				{
					freeSpace -= pos_front;
					resizeSpace(me.data_string, length(me.data_string) - freeSpace, pos_back + 1, length(me.data_string), tag);
					resizeSpace(me.data_string, (TSize)0, (TSize)0, pos_front, tag);
					pos_back -= pos_front;
					pos_front = 0;
				}
			}
		}
		else
		{
			// |data|empty|data|
			// 0
			resizeSpace(me.data_string, new_freeSpace, pos_back + 1, pos_front, tag);
			pos_front += new_freeSpace;
		}

		me.data_begin = begin(me.data_string, Standard());
		me.data_end = end(me.data_string, Standard());
		me.data_front = me.data_begin + pos_front;
		me.data_back = me.data_begin + pos_back;
	}
	return length(me.data_string);
}

//////////////////////////////////////////////////////////////////////////////
// Definitions

typedef StringSet<TRead const, Dependent<> >	TMPReadSet;

#ifdef RAZERS_MEMOPT

	template <typename TShape>
	struct SAValue< Index<TMPReadSet, TShape> > {
		typedef Pair<
			unsigned,				
			unsigned,
			BitPacked<24, 8>	// max. 16M reads of length < 256
		> Type;
	};
	
#else

	template <typename TShape>
	struct SAValue< Index<TMPReadSet, TShape> > {
		typedef Pair<
			unsigned,			// many reads
			unsigned,			// of arbitrary length
			Pack
		> Type;
	};

#endif

	
template <typename TShape>
struct Cargo< Index<TMPReadSet, TShape> > {
	struct Type
	{
		double		abundanceCut;
		int			_debugLevel;

		Type() : abundanceCut(0), _debugLevel(0)
        {}
	};
};

#ifdef RAZERS_PRUNE_QGRAM_INDEX

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TShape>
inline bool _qgramDisableBuckets(Index<TMPReadSet, IndexQGram<TShape> > &index) 
{
	typedef Index<TReadSet, IndexQGram<TShape>	>		TReadIndex;
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
template <typename TReadSet, typename TNameSet, typename TRazerSOptions>
bool loadReads(
	TReadSet &reads,			// resulting mate sequences
	TNameSet &fastaIDs,			// resulting mate ids
	const char *fileNameL,		// left mates file
	const char *fileNameR,		// right mates file
	TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);

    SeqFileIn leftMates, rightMates;

    bool success;
    if (!isEqual(fileNameL, "-"))
        success = open(leftMates, fileNameL);
    else
        success = open(leftMates, std::cin);
    if (!success)
        return false;

    if (!isEqual(fileNameR, "-"))
        success = open(rightMates, fileNameR);
    else
        success = open(rightMates, std::cin);
    if (!success)
        return false;

    CharString fastaId[2];
    String<Dna5Q> seq[2];
    CharString qual[2];

	unsigned kickoutcount = 0;
	while (!atEnd(leftMates) && !atEnd(rightMates))
	{
        readRecord(fastaId[0], seq[0], qual[0], leftMates);         // read Fasta id, sequence and qualities
        readRecord(fastaId[1], seq[1], qual[1], rightMates);        // read Fasta id, sequence and qualities

		if (options.readNaming == 0)
		{
            append(fastaId[0], "/L");
            append(fastaId[1], "/R");
			appendValue(fastaIDs, fastaId[0]);
			appendValue(fastaIDs, fastaId[1]);
		}
		
		reverseComplement(seq[1]);
		reverse(qual[1]);
		
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
						++kickoutcount;
						break;
					}
			}
		}
		
		for (int j = 0; j < 2; ++j)
		{
			// store dna and quality together
			for (unsigned p = 0; p < length(qual[j]) && p < length(seq[j]); ++p)
				assignQualityValue(seq[j][p], (int)(ordValue(qual[j][p]) - 33));
			
			if (options.trimLength > 0 && length(seq[j]) > (unsigned)options.trimLength)
				resize(seq[j], options.trimLength);
		}
		appendValue(reads, seq[0], Generous());
		appendValue(reads, seq[1], Generous());
	}

	if (atEnd(leftMates) != atEnd(rightMates) && options._debugLevel > 1)
	{
		::std::cerr << "Numbers of mates differ.\n";
		return false;
	}

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality mate-pairs.\n";
	return !empty(reads);
}

	template <typename TReadMatch>
	struct LessPairErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if ((a.rseqNo >> 1) < (b.rseqNo >> 1)) return true;
			if ((a.rseqNo >> 1) > (b.rseqNo >> 1)) return false;

			// quality
#ifdef RAZERS_MATEPAIRS
			if (a.pairScore > b.pairScore) return true;
			if (a.pairScore < b.pairScore) return false;
			return a.pairId < b.pairId;
#else
			return a.editDist < b.editDist;
#endif
		}
	};
	


//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TCounts, typename TSpec, typename TSwiftL, typename TSwiftR >
void compactPairMatches(TMatches &matches, TCounts & /*cnts*/, RazerSOptions<TSpec> &options, TSwiftL &swiftL, TSwiftR &swiftR)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreDistCutOff = std::numeric_limits<int>::min();

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	// sort 
	::std::sort(it, itEnd, LessPairErrors<TMatch>());

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == ((*it).rseqNo >> 1))
		{ 
			if ((int)(*it).pairScore <= scoreDistCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = -1 - (*it).pairScore;
					if (options.purgeAmbiguous)
						maxErrors = -1;

					setMaxErrors(swiftL, readNo, maxErrors);
					setMaxErrors(swiftR, readNo, maxErrors);

					if (maxErrors == -1 && options._debugLevel >= 2)
						::std::cerr << "(read #" << readNo << " disabled)";

					if (options.purgeAmbiguous)
						dit = ditBeg;
				}
#endif
				continue;
			}
		}
		else
		{
			readNo = (*it).rseqNo >> 1;
			hitCount = 0;
			if (options.distanceRange > 0)
				scoreDistCutOff = (*it).pairScore - options.distanceRange;
			ditBeg = dit;
		}
		*dit = *it;	++dit; ++it;
		*dit = *it;	++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
}



#ifndef RAZERS_PARALLEL
//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TCounts,
	typename TSpec >
void mapMatePairReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternL,
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPatternR,
	TVerifier &forwardPatternsL,
	TVerifier &forwardPatternsR,
	TCounts & cnts,
	char orientation,				// q-gram index of reads
	RazerSOptions<TSpec> &options)
{
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename MakeSigned_<TGPos>::Type				TSignedGPos;
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >			TSwiftFinderR;
	//typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// MATE-PAIR FILTRATION
	typedef Pair<int64_t,TMatch>							TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;

	const unsigned NOT_VERIFIED = 1u << (8*sizeof(unsigned)-1);

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	TReadSet &readSetL = host(host(swiftPatternL));
	TReadSet &readSetR = host(host(swiftPatternR));

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
	String<int64_t> lastPotMatchNo;		// last number of a left-mate potential
	int64_t lastNo = 0;					// last number over all left-mate pot. matches in the queue
	int64_t firstNo = 0;				// first number over all left-mate pot. match in the queue
	Pair<TGPos> gPair;

	resize(lastPotMatchNo, length(host(swiftPatternL)), (int64_t)-1, Exact());

	TSize gLength = length(genome);
	TMatch mR = {	// to supress uninitialized warnings
		0, 0, 0, 0,
#ifdef RAZERS_MATEPAIRS
		0, 0, 0,
#endif
		0,
#ifdef RAZERS_EXTENDED_MATCH
		0, 0,
#endif
#ifdef RAZERS_SPLICED
		0, 0,
#endif
		0
	};
	TDequeueValue fL(-1, mR);	// to supress uninitialized warnings
	fL.i2.gseqNo = gseqNo;
	mR.gseqNo = gseqNo;
	fL.i2.orientation = orientation;
	mR.orientation = (orientation == 'F')? 'R': 'F';
	
//	unsigned const preFetchMatches = 2048;

	// iterate all verification regions returned by SWIFT
	while (find(swiftFinderR, swiftPatternR, options.errorRate)) 
	{
		unsigned rseqNo = swiftPatternR.curSeqNo;
		TGPos rEndPos = endPosition(swiftFinderR) + scanShift;
		TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;

		// remove out-of-window left mates from fifo
		while (!empty(fifo) && front(fifo).i2.gEnd + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
		{
			popFront(fifo);
			++firstNo;
		}
/*		
		if (empty(fifo) || back(fifo).gEnd + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
			for (unsigned i = 0; i < preFetchMatches; ++i)
				if (find(swiftFinderL, swiftPatternL, options.errorRate, false))
					pushBack(fifo, mL);
				else
					break;
*/
		// add within-window left mates to fifo
		while (empty(fifo) || back(fifo).i2.gEnd + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
		{
			if (find(swiftFinderL, swiftPatternL, options.errorRate))
			{
				gPair = positionRange(swiftFinderL);
				if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
				{
					// link in
					fL.i1 = lastPotMatchNo[swiftPatternL.curSeqNo];
					lastPotMatchNo[swiftPatternL.curSeqNo] = lastNo++;

					fL.i2.rseqNo = swiftPatternL.curSeqNo | NOT_VERIFIED;
					fL.i2.gBegin = gPair.i1;
					fL.i2.gEnd = gPair.i2;
					
					pushBack(fifo, fL);
				}
			} else
				break;
		}

		int	bestLeftErrors = std::numeric_limits<int>::max();
		int bestLibSizeError = std::numeric_limits<int>::max();
		TDequeueIterator bestLeft = TDequeueIterator();

		TDequeueIterator it;
		int64_t lastPositive = (int64_t)-1;
		for (int64_t i = lastPotMatchNo[rseqNo]; firstNo <= i; i = (*it).i1)
		{
			it = &value(fifo, i - firstNo);
			
			// search left mate
//			if (((*it).i2.rseqNo & ~NOT_VERIFIED) == rseqNo)
			{
				// verify left mate (equal seqNo), if not done already
				if ((*it).i2.rseqNo & NOT_VERIFIED)
				{
					if (matchVerify(
							(*it).i2, infix(genome, (*it).i2.gBegin, (*it).i2.gEnd), 
							rseqNo, readSetL, forwardPatternsL, 
							options, TSwiftSpec()))
					{
						(*it).i2.rseqNo &= ~NOT_VERIFIED;		// has been verified positively
						
						// short-cut negative matches
						if (lastPositive == (int64_t)-1)
							lastPotMatchNo[rseqNo] = i;
						else
							value(fifo, lastPositive - firstNo).i1 = i;
						lastPositive = i;
					} else
						(*it).i2.rseqNo = ~NOT_VERIFIED;		// has been verified negatively
				}
/*
				if ((*it).i2.rseqNo == rseqNo)
				{
					bestLeft = it;
					bestLeftErrors = (*it).i2.editDist;
					break;
				}
*/
				if ((*it).i2.rseqNo == rseqNo)
					if (bestLeftErrors >= (*it).i2.editDist)
					{
						int libSizeError = options.libraryLength - (int)((int64_t)mR.gEnd - (int64_t)(*it).i2.gBegin);
						if (libSizeError < 0) libSizeError = -libSizeError;
						if (bestLeftErrors == (*it).i2.editDist)
						{
							if (bestLibSizeError > libSizeError)
							{
								bestLibSizeError = libSizeError;
								bestLeft = it;
							}
						}
						else
						{
							bestLeftErrors = (*it).i2.editDist;
							bestLibSizeError = libSizeError;
							bestLeft = it;
							if (bestLeftErrors == 0) break;
						}
					}
			}
//			else
//				std::cout << "HUH?" << std::endl;
		}

		// short-cut negative matches
		if (lastPositive == (int64_t)-1)
			lastPotMatchNo[rseqNo] = (int64_t)-1;
		else
			value(fifo, lastPositive - firstNo).i1 = (int64_t)-1;
		
		// verify right mate, if left mate matches
		if (bestLeftErrors != std::numeric_limits<int>::max())
		{
			if (matchVerify(
					mR, infix(swiftFinderR),
					rseqNo, readSetR, forwardPatternsR,
					options, TSwiftSpec()))
			{
				// distance between left mate beginning and right mate end
				int64_t dist = (int64_t)mR.gEnd - (int64_t)(*bestLeft).i2.gBegin;
				if (dist <= options.libraryLength + options.libraryError &&
					options.libraryLength <= dist + options.libraryError)
				{
					fL.i2 = (*bestLeft).i2;

					// transform mate readNo to global readNo
					fL.i2.rseqNo = 2*rseqNo;
					mR.rseqNo   = 2*rseqNo + 1;

					// transform coordinates to the forward strand
					if (orientation == 'R') 
					{
						TSize temp = fL.i2.gBegin;
						fL.i2.gBegin = gLength - fL.i2.gEnd;
						fL.i2.gEnd = gLength - temp;
						temp = mR.gBegin;
						mR.gBegin = gLength - mR.gEnd;
						mR.gEnd = gLength - temp;
						dist = -dist;
					}

					// set a unique pair id
					fL.i2.pairId = mR.pairId = options.nextMatePairId;
					if (++options.nextMatePairId == 0)
						options.nextMatePairId = 1;

					// score the whole match pair
					fL.i2.pairScore = mR.pairScore = 0 - fL.i2.editDist - mR.editDist;

					// relative positions
					fL.i2.mateDelta = dist;
					mR.mateDelta = -dist;

					// both mates match with correct library size
/*								std::cout << "found " << rseqNo << " on " << orientation << gseqNo;
					std::cout << " dist:" << dist;
					if (orientation=='F')
						std::cout << " \t_" << fL.i2.gBegin+1 << "_" << mR.gEnd;
					else
						std::cout << " \t_" << mR.gBegin+1 << "_" << mL.gEnd;
//							std::cout << " L_" << (*bestLeft).gBegin << "_" << (*bestLeft).gEnd << "_" << (*bestLeft).editDist;
//							std::cout << " R_" << mR.gBegin << "_" << mR.gEnd << "_" << mR.editDist;
					std::cout << std::endl;
*/
					if (!options.spec.DONT_DUMP_RESULTS)
					{
						appendValue(matches, fL.i2, Generous());
						appendValue(matches, mR, Generous());
						if (length(matches) > options.compactThresh)
						{
							typename Size<TMatches>::Type oldSize = length(matches);
//									maskDuplicates(matches);	// overlapping parallelograms cause duplicates
							compactPairMatches(matches, cnts, options, swiftPatternL, swiftPatternR);
							options.compactThresh += (options.compactThresh >> 1);
							if (options._debugLevel >= 2)
								::std::cerr << '(' << oldSize - length(matches) << " matches removed)";
						}
					}
					++options.TP;
				}
				else
					++options.FP;
			}
		} 
	}
}
#endif


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TMatches, 
	typename TReadSet_, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapMatePairReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet_ const &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_CONCATREADS
	typedef TReadSet_									TReadSet;
#else
	typedef typename Value<TReadSet_ const>::Type		TRead;
	typedef StringSet<TRead, Dependent<> >				TReadSet;
#endif
	typedef Index<TReadSet, IndexQGram<TShape> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

//	std::cout << "SA-TYPE:" <<sizeof(typename SAValue<TIndex>::Type)<<std::endl;

	// split mate-pairs over two indices
	TReadSet readSetL, readSetR;
	unsigned readCount = length(readSet) / 2;
	reserve(readSetL, readCount, Exact());
	reserve(readSetR, readCount, Exact());

	for (unsigned i = 0; i < readCount; ++i)
	{
#ifdef RAZERS_CONCATREADS
		appendValue(readSetL, readSet[2*i], Generous());
		appendValue(readSetR, readSet[2*i+1], Generous());
#else
		assign(readSetL[i], readSet[2*i], Exact());
		assign(readSetR[i], readSet[2*i+1], Exact());
#endif
	}

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
		resize(forwardPatternsL, readCount, Exact());
		resize(forwardPatternsR, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatternsL[i], readSetL[i]);
			setHost(forwardPatternsR[i], readSetR[i]);
			_patternMatchNOfPattern(forwardPatternsL[i], options.matchN);
			_patternMatchNOfPattern(forwardPatternsR[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsL[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsR[i], options.matchN);
		}
	}

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	unsigned numFiles = length(genomeFileNameList);
	unsigned gseqNo = 0;

	// open genome files, one by one	
	for (unsigned filecount = 0; filecount < numFiles; ++filecount)
	{
		// open genome file	
		SeqFileIn file;
		if (!open(file, toCString(genomeFileNameList[filecount])))
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
		for(; !atEnd(file); ++gseqNo)
		{
            readRecord(id, genome, file);			// read Fasta id and sequence

			if (options.genomeNaming == 0)
			{
                cropAfterFirst(id, IsWhitespace());     // crop id after the first whitespace
				appendValue(genomeNames, id, Generous());
			}
			
			gnoToFileMap.insert(::std::make_pair(gseqNo,::std::make_pair(genomeName,gseqNoWithinFile)));
			
			if (options.forward)
				mapMatePairReads(matches, genome, gseqNo, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'F', options);

			if (options.reverse)
			{
				reverseComplement(genome);
				mapMatePairReads(matches, genome, gseqNo, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
	}

	compactPairMatches(matches, cnts, options, swiftPatternL, swiftPatternR);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Swift FP: " << options.FP << ::std::endl;
		::std::cerr << "Swift TP: " << options.TP << ::std::endl;
	}
	return 0;
}


}

#endif
