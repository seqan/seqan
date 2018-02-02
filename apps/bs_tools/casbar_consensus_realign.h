#ifndef __APPS_BS_TOOLS_CASBAR_CONSENSUS_REALIGN_H__
#define __APPS_BS_TOOLS_CASBAR_CONSENSUS_REALIGN_H__

#include <seqan/misc/svg.h>

using namespace std;
using namespace seqan;

// For read bases
template<typename TProfileChar>
inline void
addToProfileChar(TProfileChar &profileChar, bool &top, Dna5 const &readBase, double &qual, double &mapq)
{
    typedef typename ValueSize<TProfileChar>::Type TSize;

    double e = pow(10, -qual/10.0);
    double pm = 1.0 - pow(10, -mapq/10.0);
    if (mapq > 254.5) pm = 1.0;
    if (top)
    {
        if (readBase == 'N')
            for (TSize i = 0; i < 4; ++i)
                profileChar.count[i] += (1.0/4.0)*pm;   // Not equal distr. !
        else
        {
            profileChar.count[ordValue(readBase)] += (1.0-e)*pm;
            for (TSize i = 0; i < 4; ++i)
                if (i != static_cast<TSize>(ordValue(readBase)))
                    profileChar.count[i] += (e/3.0)*pm;
        }
    }
    else
    {
        if (readBase == 'N')
            for (TSize i = 4; i < 8; ++i)
                profileChar.count[i] += (1.0/4.0)*pm;
        else
        {
            profileChar.count[ordValue(readBase)+4] += (1.0-e)*pm;
            for (TSize i = 4; i < 8; ++i)
                if (i != static_cast<TSize>(ordValue(readBase)+4))
                    profileChar.count[i] += (e/3.0)*pm;
        }
    }
    ++profileChar.count[9]; // Count (non-gap) bases
}

// For reference bases
template<typename TProfileChar>
inline void
addToProfileChar(TProfileChar &profileChar, Dna5 const &refBase)
{
    profileChar.count[8] = ordValue(refBase);
    ++profileChar.count[9]; // Count (non-gap) bases
 }

// For read bases
template<typename TProfileChar>
inline void
removeFromProfileChar(TProfileChar &profileChar, bool &top, Dna5 const &readBase, double &qual, double &mapq)
{
    typedef typename ValueSize<TProfileChar>::Type TSize;

    double e = pow(10, -qual/10.0);
    double pm = 1.0 - pow(10, -mapq/10.0);
    if (mapq > 254.5) pm = 1.0;

    if (top)
    {
        if (readBase == 'N')
            for (TSize i = 0; i < 4; ++i)
                profileChar.count[i] -= (1.0/4.0)*pm;
        else
        {
            profileChar.count[ordValue(readBase)] -= (1.0-e)*pm;
            for (TSize i = 0; i < 4; ++i)
                if (i != static_cast<TSize>(ordValue(readBase)))
                    profileChar.count[i] -= (e/3.0)*pm;
        }
    }
    else
    {
        if (readBase == 'N')
            for (TSize i = 4; i < 8; ++i)
                profileChar.count[i] -= (1.0/4.0)*pm;
        else
        {
            profileChar.count[ordValue(readBase)+4] -= (1.0-e)*pm;
            for (TSize i = 4; i < 8; ++i)
                if (i != static_cast<TSize>(ordValue(readBase)+4))
                    profileChar.count[i] -= (e/3.0)*pm;
        }
    }
    --profileChar.count[9];     // Count (non-gap) bases

}

// For ref bases
template<typename TProfileChar>
inline void
removeFromProfileChar(TProfileChar &profileChar, Dna5 const &/*readBase*/)
{
    profileChar.count[8] = -1;  // Empty ref base entry in profile
    --profileChar.count[9];     // Count (non-gap) bases
}

template <typename TAlignedReads, typename TSpec, typename TFragmentStore, typename TGapPos>
inline void
insertGap(String<TAlignedReads, TSpec>& alignedReadStoreTmp,
          int &numGapsTop,
          int &numGapsBottom,
          TFragmentStore& fragmentStoreOrig,   // Needed to get original strand to separate gap counts
		  TGapPos const gapPos)
{
	typedef String<TAlignedReads, TSpec> TAlignedReadStore;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;

    numGapsTop = 0;
    numGapsBottom = 0;
	TAlignIter alignIt = begin(alignedReadStoreTmp, Standard());
	TAlignIter alignItEnd = end(alignedReadStoreTmp, Standard());
	for(;alignIt != alignItEnd; ++alignIt)
    {
        if (alignIt->id == length(fragmentStoreOrig.alignedReadStore))    // Ref
            insertGap(*alignIt, gapPos);
        else if (fragmentStoreOrig.alignedReadStore[alignIt->id].beginPos <  fragmentStoreOrig.alignedReadStore[alignIt->id].endPos) // Top
            numGapsTop += (int)(insertGap(*alignIt, gapPos) * (1.0 -pow(10, -fragmentStoreOrig.alignQualityStore[alignIt->id].score/10.0)));
        else                                                // Bottom
            numGapsBottom += (int)(insertGap(*alignIt, gapPos) * (1.0 -pow(10, -fragmentStoreOrig.alignQualityStore[alignIt->id].score/10.0)));
    }
}

// Perform one realignment round.
// TODO(holtgrew): Rename to reflect this more clearly.
// TODO(holtgrew): TConsensus/consensus are profiles, really.
template<typename TFragSpec, typename TConfig, typename TAlignedRead, typename TSpec, typename TConsensus, typename TBandwidth, typename TOptions, typename TModel>
void
reAlign(double &profileScore,
        FragmentStore<TFragSpec, TConfig>& fragStore,
		String<TAlignedRead, TSpec>& contigReads,
		TConsensus& consensus,
		TBandwidth const bandwidth,
		bool includeReference,
		TOptions &options,
        double & timeBeforeAlign,
        double & timeAlign,
        double & timeAfterAlign,
        TModel const &
        )
{
    //std::cout << " Do reAlign() ... " << std::endl;
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef String<TAlignedRead, TSpec> TAlignedReadStore;
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TReadSeq TReadSeq;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadGapAnchor TGapAnchor;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
	typedef typename Iterator<TConsensus, Standard>::Type TConsIter;

    profileScore = 0;

	// Initialization
	typedef typename Value<TConsensus>::Type TProfileChar;

	// Remove each fragment and realign it to the profile.
	TAlignedReadIter alignIt = begin(contigReads, Standard());
	TAlignedReadIter alignItEnd = end(contigReads, Standard());
	if (includeReference)
        --alignItEnd;
	TConsensus bandConsensus;
	TConsensus myRead;
	TConsensus newConsensus;
    int i = 0;
    int help = 0;   // TODO rm
	for (; alignIt != alignItEnd; ++alignIt) {
        double tBegin = sysTime();
        if (i++ > 1000 || i == 1) {
            if (i != 1)
                i = 0;
            fflush(stdout);
        }

        bool isRef = false;
	    if (alignIt->id == length(fragStore.alignedReadStore)) isRef = true;
        bool top = (fragStore.alignedReadStore[alignIt->id].beginPos < fragStore.alignedReadStore[alignIt->id].endPos);
        double mapE = 0.0; // ?
        if (!isRef) mapE = pow(10, -fragStore.alignQualityStore[alignIt->id].score/10.0);
        double mapq = 0;
        if (!isRef) mapq = fragStore.alignQualityStore[alignIt->id].score;  // TODO Get rid of it

        ++help;   // TODO rm
		TSize itConsPos = 0;
		TConsIter itCons = begin(consensus, Standard());
		TConsIter itConsEnd = end(consensus, Standard());

		// Initialize the consensus of the band. -> part of the whole profile where current read mapped
		clear(myRead);
		resize(myRead, length(fragStore.readSeqStore[alignIt->readId]), TProfileChar());
		resize(bandConsensus, 2 * bandwidth + (alignIt->endPos - alignIt->beginPos), Generous());
		TConsIter bandConsIt = begin(bandConsensus);
		TConsIter bandConsItEnd = end(bandConsensus);
		TConsIter myReadIt = begin(myRead);
		TReadPos bandOffset = 0;
		if (bandwidth < (TBandwidth) alignIt->beginPos) {
			bandOffset = alignIt->beginPos - bandwidth;
			itCons += bandOffset; itConsPos += bandOffset;
			SEQAN_ASSERT_LEQ(itCons, itConsEnd);
		}
		int leftDiag = (alignIt->beginPos - bandOffset) - bandwidth;
		int rightDiag = leftDiag + 2 * bandwidth;
		//int increaseBand = 0;
		int increaseBandLeft = 0;
		int increaseBandRight = 0;
		int removedBeginPos = 0;
		int removedEndPos = 0;
		for (TReadPos iPos = bandOffset; iPos < alignIt->beginPos && itCons != itConsEnd && bandConsIt != bandConsItEnd; ++itCons, ++bandConsIt, ++itConsPos, ++iPos)
			*bandConsIt = *itCons; // fill in positions left of readbegin
		TSize itConsPosBegin = itConsPos;  // start position of read basically, right? if(itConsPosBegin != alignIt->beginPos) std::cout <<"nicht unbedingt gleich\n";
		alignIt->beginPos = alignIt->endPos = 0; // So this read is discarded in all gap operations

		// Remove sequence from profile (and add to the consensus??)  // TODO(holtgrew): Add to consensus part right?
		typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
		TReadIter itRead = begin(fragStore.readSeqStore[alignIt->readId], Standard());
		TReadIter itReadEnd = end(fragStore.readSeqStore[alignIt->readId], Standard());
		typedef typename Iterator<String<TGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(alignIt->gaps, Standard());
		TReadGapsIter itGapsEnd = end(alignIt->gaps, Standard());
		TReadPos old = 0;
		int diff = 0;
		TReadPos clippedBeginPos = 0;
		TReadPos clippedEndPos = 0;
		SEQAN_ASSERT_LT(itRead, itReadEnd);
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			old = itGaps->seqPos;
			clippedBeginPos = old; // gaps at beginning? or really clipped?
			itRead += old;
			diff -= old;
			++itGaps;
			SEQAN_ASSERT_LT(itRead, itReadEnd);
		}
		for (; itGaps != itGapsEnd && itCons != itConsEnd; ++itGaps) {
			// limit should never be larger than read length
			TReadPos limit = itGaps->seqPos;
			SEQAN_ASSERT_LT(itGaps->seqPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
			int newDiff = (itGaps->gapPos - limit);
			SEQAN_ASSERT_LT(itGaps->gapPos, (TReadPos)length(consensus));
			if (diff > newDiff) {
				clippedEndPos = diff - newDiff;
				limit -= clippedEndPos;
			}
			for (; old < limit && itCons != itConsEnd && itRead != itReadEnd && bandConsIt != bandConsItEnd; ++old, ++itRead) {
                // Recalculate profile value for current read to subtract it
                if (!isRef)
                {
                    double qual = getQualityValue(*itRead);
                    double mapq = fragStore.alignQualityStore[alignIt->id].score;
                    bool top = (fragStore.alignedReadStore[alignIt->id].beginPos < fragStore.alignedReadStore[alignIt->id].endPos);
                    removeFromProfileChar(*itCons, top, (Dna5)(*itRead), qual, mapq);
                }
                else
                {
                    removeFromProfileChar(*itCons, (Dna5)(*itRead));
                }
				if (!empty(*itCons)) {
					*bandConsIt = *itCons;
					++bandConsIt;
					++itConsPos;
					removedEndPos = 0;
				}
				else
				{
					if (itConsPosBegin != itConsPos)
					{
                        ++increaseBandLeft; // insertion --> increaseBandLeft, read has character here, consensus doesnt
						++removedEndPos;
					}
					else ++removedBeginPos; // begin gaps
					removeGap(contigReads, itConsPos);
				}
				(*myReadIt).count[0] = ordValue(*itRead);
				if (!isRef) (*myReadIt).count[1] = getQualityValue(*itRead);            // Me: store quality value
				++myReadIt;
				++itCons;
			}
			for (; diff < newDiff && itCons != itConsEnd && bandConsIt != bandConsItEnd; ++diff) {
                ++increaseBandRight; // deletion --> increaseBandRight, read has gaps here, consensus doesnt
				if (isRef) (*itCons).count[8] = -2;
                else
                {
                    if(top) (*itCons).count[10] -= (1.0-mapE);
                    else (*itCons).count[11] -= (1.0-mapE);
                }

				if (!empty(*itCons)) {
					*bandConsIt = *itCons;
					++bandConsIt;
					++itConsPos;
				}
				else
				    removeGap(contigReads, itConsPos);  //++increaseBandRight;}
				++itCons;
			}
		}
		if (!clippedEndPos) {
			for (; itRead!=itReadEnd && itCons != itConsEnd && bandConsIt != bandConsItEnd; ++itRead) {
				//SEQAN_ASSERT_LT(itCons, itConsEnd);
                //SEQAN_ASSERT_LT(itRead, itReadEnd);
                // Recalculate profile value for current read to subtract it
                if (!isRef)
                {
                    double qual = getQualityValue(*itRead);
                    double mapq = fragStore.alignQualityStore[alignIt->id].score;
                    bool top = (fragStore.alignedReadStore[alignIt->id].beginPos < fragStore.alignedReadStore[alignIt->id].endPos);
                    removeFromProfileChar(*itCons, top, (Dna5)(*itRead), qual, mapq);
                }
                else
                    removeFromProfileChar(*itCons, (Dna5)(*itRead));
				if (!empty(*itCons))
				{
					*bandConsIt = *itCons;
					++bandConsIt;
					++itConsPos;
					removedEndPos = 0;
				}
				else
				{  // only gaps left in this column after removing myRead
					if (itConsPosBegin != itConsPos)
					{
                        ++increaseBandLeft; // insertion --> increaseBandLeft, read is longer than consensus here
                        ++removedEndPos;
					}
					else ++removedBeginPos;
					removeGap(contigReads, itConsPos);
				}
				(*myReadIt).count[0] = ordValue(*itRead);
				if (!isRef) (*myReadIt).count[1] = getQualityValue(*itRead);
				++myReadIt;
				++itCons;
			}
		}
		bool singleton = (itConsPosBegin == itConsPos);
		increaseBandLeft -= removedEndPos;
        //increaseBand = increaseBandLeft + increaseBandRight;

		// Go further up to the bandwidth
		for (TReadPos iPos = 0; ((itCons != itConsEnd) && (iPos < (TReadPos) bandwidth)) && bandConsIt != bandConsItEnd; ++itCons, ++iPos, ++bandConsIt)
            *bandConsIt = *itCons;
		resize(bandConsensus, bandConsIt - begin(bandConsensus, Standard()), Generous());
		resize(myRead, myReadIt - begin(myRead, Standard()), Generous());

		// Realign the consensus with the sequence.
		typedef StringSet<TConsensus, Dependent<> > TStringSet;
		TStringSet pairSet;
		appendValue(pairSet, bandConsensus);
		appendValue(pairSet, myRead);

		typedef String<Fragment<> > TFragmentString;
		TFragmentString matches;

        double tBegAlign = sysTime();
		leftDiag -= removedBeginPos;
		rightDiag -= removedBeginPos;

        // TODO precompute als this....
        double const *seqErrorFreqs;
        double const *delErrorFreqs;
        double const *insErrorFreqs;
        if (options.nonSimpleSubstErrors) seqErrorFreqs = SeqErrorFreqs<double, BsNonSimple>::getData();
        else seqErrorFreqs = SeqErrorFreqs<double, BsSimple>::getData();
        if (options.nonSimpleInsErrors) insErrorFreqs = InsErrorFreqs<double, BsNonSimple>::getData();
        else insErrorFreqs = InsErrorFreqs<double, BsSimple>::getData();
        double scalingFactorDelErrors;
        if (options.nonSimpleDelErrors)
        {
            delErrorFreqs = DelErrorFreqs<double, BsNonSimple>::getData();
            scalingFactorDelErrors = 1.0; //options.scalingFactorDelErrorsNonSimple;
        }
        else
        {
            delErrorFreqs = DelErrorFreqs<double, BsSimple>::getData();
            scalingFactorDelErrors = options.scalingFactorDelErrorsSimple;
        }
		if (!singleton) {
            if (isRef)
            {
                Score<double, BsTagList<BsProfileScoreRef, TModel, InnerCell> > scoringScheme(options, seqErrorFreqs, insErrorFreqs, delErrorFreqs, scalingFactorDelErrors);
                assignTargetFreqs(scoringScheme, bandConsensus, options);
                profileScore += globalAlignment(matches, pairSet, scoringScheme, AlignConfig<false,false,false,false>(), _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBandRight, (int) length(pairSet[0])), Gotoh());
            }
            else if (top && (fragStore.readStore[alignIt->readId].matePairId == TAlignedRead::INVALID_ID || fragStore.readStore[alignIt->readId].matePairId  == 1)) // top, original
            {
                Score<double, BsTagList<BsProfileScoreCT, TModel, InnerCell> > scoringScheme(options, seqErrorFreqs, insErrorFreqs, delErrorFreqs, scalingFactorDelErrors);
                assignTargetFreqs(scoringScheme, bandConsensus, options, BsTop());
                profileScore += globalAlignment(matches, pairSet, scoringScheme, AlignConfig<false,false,false,false>(), _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBandRight, (int) length(pairSet[0])), Gotoh());
            }
            else if (top && fragStore.readStore[alignIt->readId].matePairId  == 2)    // top, reverse mapped
            {
                Score<double, BsTagList<BsProfileScoreCTRight, TModel, InnerCell> > scoringScheme(options, seqErrorFreqs, insErrorFreqs, delErrorFreqs, scalingFactorDelErrors);
                assignTargetFreqs(scoringScheme, bandConsensus, options, BsTop());
                profileScore += globalAlignment(matches, pairSet, scoringScheme, AlignConfig<false,false,false,false>(), _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBandRight, (int) length(pairSet[0])), Gotoh());
            }
            else if (!top && (fragStore.readStore[alignIt->readId].matePairId  == TAlignedRead::INVALID_ID || fragStore.readStore[alignIt->readId].matePairId  == 1))    // bottom, original
            {
                Score<double, BsTagList<BsProfileScoreGA, TModel, InnerCell> > scoringScheme(options, seqErrorFreqs, insErrorFreqs, delErrorFreqs, scalingFactorDelErrors);
                assignTargetFreqs(scoringScheme, bandConsensus, options, BsBottom());
                profileScore += globalAlignment(matches, pairSet, scoringScheme, AlignConfig<false,false,false,false>(), _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBandRight, (int) length(pairSet[0])), Gotoh());
            }
            else if (!top && fragStore.readStore[alignIt->readId].matePairId  == 2)   // bottom, forward mapped
            {
                Score<double, BsTagList<BsProfileScoreGARight, TModel, InnerCell> > scoringScheme(options, seqErrorFreqs, insErrorFreqs, delErrorFreqs, scalingFactorDelErrors);
                assignTargetFreqs(scoringScheme, bandConsensus, options, BsBottom());
                profileScore += globalAlignment(matches, pairSet, scoringScheme, AlignConfig<false,false,false,false>(), _max(leftDiag - increaseBandLeft, -1 * (int) length(pairSet[1])), _min(rightDiag + increaseBandRight, (int) length(pairSet[0])), Gotoh());
            }
            //std::cout << "in between profileScore: " << profileScore << std::endl;
   		}
        double tEndAlign = sysTime();

		// Add the read back to the consensus and build the new consensus.
		resize(newConsensus, length(bandConsensus) + length(myRead), Generous());
		TConsIter newConsIt = begin(newConsensus, Standard());
		TConsIter bandIt = begin(bandConsensus, Standard());
		TConsIter bandItEnd = end(bandConsensus, Standard());
		typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
		TFragIter fragIt = end(matches, Standard());
		TFragIter fragItEnd = begin(matches, Standard());
		TReadPos consPos = 0;
		TReadPos readPos = 0;
		TReadPos alignPos = 0;
		clear(alignIt->gaps);
		diff = 0;
		if (clippedBeginPos) {
			appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos, 0), Generous() );
			diff -= clippedBeginPos;
		}
		bool firstMatch = true;
		if (fragIt != fragItEnd) { // walk through segment matches that represent read-msa alignment
			do {
				--fragIt;
				int gapLen = fragIt->begin1 - consPos;  // Number of gaps in profile before next fragment starts
				if (firstMatch) gapLen = 0;
				// While we havn't reached the next fragment in the current profile yet (-> gaps in read)
				while (consPos < (TReadPos)fragIt->begin1) {
					SEQAN_ASSERT_LT(bandIt, bandItEnd);
					SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
					if (isRef) (*bandIt).count[8] = -1;
                    else
                    {
                        if (!firstMatch && top) (*bandIt).count[10] += (1.0-mapE);  // Update gap counts for current gap
                        else if (!firstMatch) (*bandIt).count[11] += (1.0-mapE);
                    }
					*newConsIt = *bandIt;   // Fill new profile with updated old profile columns at read gap poitions
					++newConsIt;
					++bandIt;
					++consPos;
					++alignPos;
				}
				// We are already at the next fragment begin position in the profile, but there read bases left
                // While we haven't reached the next fragment begin position in the read with readPos (-> gaps in profile)
				while (readPos < (TReadPos)fragIt->begin2) {
					SEQAN_ASSERT_LT(readPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
					SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
                    // Append new gap anchor to read gaps presenting given length
                    if (gapLen) {
						diff += gapLen; // add gap of length gaplen to readGaps
						appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff), Generous() );
						gapLen = 0; // do this only once
					}
					int numGapsTop;
					int numGapsBottom;
					// Is this done for all reads? shouldn't we skip the current read???
			        insertGap(contigReads, numGapsTop, numGapsBottom, fragStore, bandOffset + alignPos);    // insert gap into each read of contigReads at corresponding position and return numbers/freqs
					TProfileChar tmpChar;   // insert new column into profile, containing only current read base
					if (!isRef)
                    {
                        double qual = myRead[readPos].count[1];
						addToProfileChar(tmpChar, top, (unsigned)myRead[readPos].count[0], qual, mapq);
                        tmpChar.count[8] = -1;
                    }
                    else
                    {
						addToProfileChar(tmpChar, (unsigned)myRead[readPos].count[0]);
                    }
					tmpChar.count[10] = numGapsTop;
                    tmpChar.count[11] = numGapsBottom;
					*newConsIt = tmpChar; ++newConsIt;
					++readPos; ++alignPos;
				}
                // For each position in the current fragment
				for (TSize i = 0; i<fragIt->len; ++i, ++bandIt, ++consPos, ++readPos, ++alignPos, ++newConsIt) {
				    SEQAN_ASSERT_LT(bandIt, bandItEnd);
					SEQAN_ASSERT_LT(readPos, (TReadPos)length(fragStore.readSeqStore[alignIt->readId]));
					SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
					if (firstMatch) {
						firstMatch = false;
						alignIt->beginPos = bandOffset + consPos;
					} else if (gapLen) {
						diff += gapLen; // Why is this done here?
						appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos, clippedBeginPos + readPos + diff), Generous() );
						gapLen = 0;
					}
					SEQAN_ASSERT_LT(bandIt, bandItEnd);
					// Add current base of read to profile
					if (!isRef)
                    {
                        double qual = myRead[readPos].count[1];
						addToProfileChar(*bandIt, top, (unsigned)myRead[readPos].count[0], qual, mapq);
                    }
                    else
                    {
						addToProfileChar(*bandIt, (unsigned)myRead[readPos].count[0]);
                    }
					*newConsIt = *bandIt;
				}
			} while (fragIt != fragItEnd);
		}
		
		for (; readPos < (TReadPos)length(myRead); ++readPos) { // For rest not matched read positions
            int numGapsTop;
            int numGapsBottom;
            insertGap(contigReads, numGapsTop, numGapsBottom, fragStore, bandOffset + alignPos);    // insert gap in each read of contigReads at corresponding position
            TProfileChar tmpChar;  // insert new column into profile, containing only current read base
			if (!isRef)
            {
                double qual = myRead[readPos].count[1];
				addToProfileChar(tmpChar, top, (unsigned)myRead[readPos].count[0], qual, mapq);
                tmpChar.count[8] = -1;
            }
            else
            {
				addToProfileChar(tmpChar, (unsigned)myRead[readPos].count[0]);
            }
	    	tmpChar.count[10] = numGapsTop;
            tmpChar.count[11] = numGapsBottom;
			SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
			*newConsIt = tmpChar; ++newConsIt;
			++alignPos;
		}
		// Adjust entries in temporary aligned read store 'contigReads' by clipped positions
		if (singleton) alignIt->beginPos = bandOffset;
		alignIt->endPos = alignIt->beginPos + clippedBeginPos + readPos + diff;
		if (clippedEndPos) {
			diff -= clippedEndPos;
			appendValue(alignIt->gaps, TGapAnchor(clippedBeginPos + readPos + clippedEndPos, clippedBeginPos + readPos + clippedEndPos + diff), Generous() );
		}
		for (; bandIt != bandItEnd; ++bandIt, ++newConsIt)
		{
			SEQAN_ASSERT_LT(newConsIt, end(newConsensus,Standard()));
			*newConsIt = *bandIt;
		}
		resize(newConsensus, newConsIt - begin(newConsensus, Standard()), Generous());

		replace(consensus, bandOffset, itCons - begin(consensus), newConsensus);

        double tEnd = sysTime();

        timeBeforeAlign += tBegAlign - tBegin;
        timeAlign += tEndAlign - tBegAlign;
        timeAfterAlign += tEnd -tEndAlign;

 	}
}

//////////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): realignmentMethod should not be optional or moved to the end of the list.
// TODO(holtgrew): The method should be selected with an enum instead of an int.
template<typename TSpec, typename TConfig, typename TId, typename TBandwidth, typename TOptions, typename TModel>
void
reAlign(FragmentStore<TSpec, TConfig> & fragStore,
		TId const contigId,
		TBandwidth const bandwidth,
		bool includeReference,
		TOptions &options,
		TModel const &)
{
	typedef FragmentStore<TSpec, TConfig>                   TFragmentStore;
	typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
	typedef typename TFragmentStore::TReadPos               TReadPos;
	
	typedef typename TFragmentStore::TContigStore		                    TContigStore;
	typedef typename Value<TContigStore>::Type		                        TContig;
	typedef typename TFragmentStore::TContigPos 		                    TContigPos;
	typedef typename TFragmentStore::TContigSeq 		                    TContigSeq;
	typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
	
	typedef typename TFragmentStore::TReadSeq                       TReadSeq;
	typedef typename TFragmentStore::TReadGapAnchor                 TGapAnchor;
	typedef typename Value<TAlignedReadStore>::Type                 TAlignedElement;

	typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

	// Sort the reads according to their begin position.
	sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
	alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());


    if(true)
    {
        std::cout << "Before realigning " << std::endl;
        TContigGaps contigGaps2(fragStore.contigStore[0].seq, fragStore.contigStore[0].gaps);
        TContigPos maxPos2 = positionSeqToGap(contigGaps2,length(fragStore.contigStore[0].seq)-1)+1;
        maxPos2 = _max(maxPos2,(TContigPos)length(fragStore.contigStore[0].seq));
        std::cout << "maxPos visual = " << maxPos2 << std::endl;
        //std::cout << " genomeLen = " << genomeLen << std::endl;
        AlignedReadLayout layout;
        layoutAlignment(layout, fragStore);
        printAlignment(std::cout, layout, fragStore, 0, (TContigPos)0, (TContigPos)200, 0, 300);

        if (maxPos2 == 855 || maxPos2 == 424 || maxPos2 == 696)
        {
            CharString name = "/home/takifugu/sabrina7/Data/plotData/before_re";
            stringstream ss;
            ss << maxPos2;
            ss << "_";
            ss << length(fragStore.alignedReadStore);
            string str = ss.str();
            append(name, str);
            append(name, ".svg");
            SVGFile svg(toCString(name));
            printAlignment(svg, layout, fragStore, 0, (TContigPos)0, (TContigPos)200, 0, 300);
        }
    }

	// Copy all reads belonging to this contig and reverse complement them if necessary.
	TAlignedReadStore contigReads;  // TODO(holtgrew): Rather contigAlignedReads?
	TReadPos maxPos = 0;
	TReadPos minPos = std::numeric_limits<TReadPos>::max();
	for (; alignIt != alignItEnd; ++alignIt) {
		if (alignIt->beginPos > alignIt->endPos) {
			reverseComplement(fragStore.readSeqStore[alignIt->readId]);
			TAlignedElement alignedEl = *alignIt;
			alignedEl.id = position(alignIt, fragStore.alignedReadStore);   // to get original beginPos from fragStore later
			TReadPos tmp = alignedEl.beginPos;
			alignedEl.beginPos = alignedEl.endPos;
			alignedEl.endPos = tmp;
			if (alignedEl.beginPos < minPos)
                minPos = alignedEl.beginPos;
			if (alignedEl.endPos > maxPos)
                maxPos = alignedEl.endPos;
			appendValue(contigReads, alignedEl, Generous() );
		} else {
			if (alignIt->beginPos < minPos)
                minPos = alignIt->beginPos;
			if (alignIt->endPos > maxPos)
                maxPos = alignIt->endPos;
			TAlignedElement alignedEl = *alignIt;
            alignedEl.id = position(alignIt, fragStore.alignedReadStore);   // to get original beginPos from fragStore later
			appendValue(contigReads, alignedEl, Generous() );
		}
	}
    // Append reference sequence to aligned reads for contigs if requested to do so.
	if (includeReference) {
		TId dummyReadId = length(fragStore.readSeqStore);
		TId dummyMatchId = length(fragStore.alignedReadStore);
		appendRead(fragStore, fragStore.contigStore[contigId].seq);
		appendValue(fragStore.readNameStore, fragStore.contigNameStore[contigId], Generous());

		TAlignedElement el;
		el.id = dummyMatchId;
		el.readId = dummyReadId;
		el.contigId = contigId;
		minPos = el.beginPos = 0;
		TContigGaps contigGaps(fragStore.contigStore[contigId].seq, fragStore.contigStore[contigId].gaps);
		maxPos = el.endPos = _max(maxPos, (TReadPos)positionSeqToGap(contigGaps,length(fragStore.contigStore[contigId].seq)-1)+1);
		maxPos = el.endPos = _max(maxPos, (TReadPos)length(fragStore.contigStore[contigId].seq));
		el.gaps = fragStore.contigStore[contigId].gaps;
		appendValue(contigReads, el, Generous());
	}

	// Create the consensus sequence
	typedef ProfileChar<DnaMR, double>                              TProfileChar;   // A, C, G, T, C (G from reverse strand), T (A from reverse strand), N, R (ref. ord value)
    typedef typename ValueSize<TProfileChar>::Type                  TSizeP;
	typedef String<TProfileChar>                                    TProfileString;
	typedef typename Iterator<TProfileString, Standard>::Type       TConsIter;
	//TSizeP gapPos = ValueSize<DnaMR>::VALUE;
	TProfileString consensus;
	TProfileChar profChar;
	for (TSizeP i = 0; i < ValueSize<TProfileChar>::VALUE; ++i)
	    profChar.count[i] = 0.0;
	profChar.count[8] = -2.0;
	resize(consensus, maxPos - minPos, profChar);

	TConsIter itCons = begin(consensus, Standard() );
	TConsIter itConsEnd = end(consensus, Standard());
	TAlignIter contigReadsIt = begin(contigReads, Standard() );
	TAlignIter contigReadsItEnd = end(contigReads, Standard() );
	for(;contigReadsIt != contigReadsItEnd; ++contigReadsIt) {
        bool isRef;
	    if (contigReadsIt->id == length(fragStore.alignedReadStore)) isRef = true;
        else isRef = false;
        bool top = true;
        if (!isRef) top = (fragStore.alignedReadStore[contigReadsIt->id].beginPos < fragStore.alignedReadStore[contigReadsIt->id].endPos);
        double mapE = 0.0; // ?
        if (!isRef) mapE = pow(10, -fragStore.alignQualityStore[contigReadsIt->id].score/10.0);

		contigReadsIt->beginPos -= minPos;  // Me: adjust read positions to min observed postion
		contigReadsIt->endPos -= minPos;
		itCons = begin(consensus, Standard() );
		itCons += contigReadsIt->beginPos;  // Me: jump to first position of curr. read in consensus/profile

		typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
		TReadIter itRead = begin(fragStore.readSeqStore[contigReadsIt->readId], Standard() );
		TReadIter itReadEnd = end(fragStore.readSeqStore[contigReadsIt->readId], Standard() );
		typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(contigReadsIt->gaps, Standard() );
		TReadGapsIter itGapsEnd = end(contigReadsIt->gaps, Standard() );

		TReadPos old = 0;
		int diff = 0;
		bool clippedEnd = false;

		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {   // Me: skip first gap before (clipping at the beginning?)
			old = itGaps->seqPos;   // Me: negative?
			itRead += old;
			diff -= old;
			++itGaps;
		}
		for(;itGaps != itGapsEnd; ++itGaps) {   // Me: for each gap anchor
			TReadPos limit = itGaps->seqPos;
			int newDiff = (itGaps->gapPos - limit);     // Me: sum of gaps up to this position
			SEQAN_ASSERT_LEQ(itGaps->gapPos, (int)length(consensus));
			// changed LT -> LEQ: in case of gap at end of contig, last gapAnchor gapPos is same as length of gaps (could this cause problems somewhere else?)
			if (diff > newDiff) {
				limit -= (diff - newDiff);
				clippedEnd = true;
			}
            // Me: for each seq. position before next gap, increase count for corresponding base
			for(;old < limit && itRead != itReadEnd && itCons != itConsEnd; ++old, ++itRead)
			{
				SEQAN_ASSERT_LT(itRead, itReadEnd);
				if (!isRef)
                {
                    double qual = getQualityValue(*itRead);
                    double mapq = fragStore.alignQualityStore[contigReadsIt->id].score;
                    addToProfileChar(value(itCons++), top, (Dna5)(*itRead), qual, mapq);
                }
                else
                {
                    addToProfileChar(value(itCons++), (Dna5)(*itRead));
                }
			}
			for(;diff < newDiff; ++diff)    // Me: increase count for each gap
            {
                if (isRef) value(itCons).count[8] = -1;
                else
                {
			        if (top) (value(itCons)).count[10] += (1.0-mapE);
                    else (value(itCons)).count[11] += (1.0-mapE);
                }
                ++itCons;
            }
		}
		if (!clippedEnd)
		{
			for( ; itRead!=itReadEnd && itCons != itConsEnd;++itRead) // Me: count bases behind last gap
            {
 				if (!isRef)
                {
                    double qual = getQualityValue(*itRead);
                    double mapq = fragStore.alignQualityStore[contigReadsIt->id].score;
                    addToProfileChar(value(itCons++), top, (Dna5)(*itRead), qual, mapq);
                }
                else
                {
                        addToProfileChar(value(itCons++), (Dna5)(*itRead));
                }
            }
		}
	}
    /*
    std::cout << "Before realigning " << std::endl;
    for (unsigned i = 0; i < 200; ++i)
    {
        std::cout << "Pos: " << i << "   ";
        for (unsigned j = 0; j < 4; ++j)
        {
            std::cout << consensus[i].count[j] << " ";
        }
        std::cout << "   ";
        for (unsigned j = 4; j < 8; ++j)
        {
            std::cout << consensus[i].count[j] << " ";
        }
        std::cout << "   ";
        for (unsigned j = 8; j < 12; ++j)
        {
            std::cout << consensus[i].count[j] << " ";
        }
        std::cout << std::endl;
    }
    */
    double tBefore = 0, tAlign = 0, tAfter = 0;
	double profileScore;
	reAlign(profileScore, fragStore, contigReads, consensus, bandwidth, includeReference, options, tBefore, tAlign, tAfter, TModel());
    //std::cout << "profileScore: " << profileScore << std::endl;
	double oldProfileScore = profileScore;
    ++profileScore;
    unsigned limit = 0;
	while(profileScore > oldProfileScore && limit < 3) {     // Me: use sum of single alignment scores as profileScore to check for improvement
		oldProfileScore = profileScore;
        double tBefore = 0, tAlign = 0, tAfter = 0;
		reAlign(profileScore, fragStore, contigReads, consensus, bandwidth, includeReference, options, tBefore, tAlign, tAfter, TModel());
		//profileScore = 0; // TODO rm, for test only
        //std::cout << "new profileScore : " << profileScore << std::endl;
        if (profileScore < oldProfileScore - 5)
        {
            ++limit;                            // avoid running forever
            profileScore = oldProfileScore + 1; // if lower, run again
        }
	}

	// Update all the aligned reads and the new consensus
	alignIt = begin(fragStore.alignedReadStore);
	TAlignIter contigReadIt = begin(contigReads, Standard());
	for (; alignIt != alignItEnd; ++alignIt) {
		if (alignIt->beginPos > alignIt->endPos) {
			reverseComplement(fragStore.readSeqStore[alignIt->readId]);
			alignIt->beginPos = contigReadIt->endPos;
			alignIt->endPos = contigReadIt->beginPos;
		} else {
			alignIt->beginPos = contigReadIt->beginPos;
			alignIt->endPos = contigReadIt->endPos;
		}
		// Remove empty gap anchors
		clear(alignIt->gaps);
		typedef typename Iterator<TGapAnchor, Standard>::Type TGapIter;
		TGapIter gapIt = begin(contigReadIt->gaps, Standard());
		TGapIter gapItEnd = end(contigReadIt->gaps, Standard());
		int diff = 0;
		for(;gapIt != gapItEnd; ++gapIt) {
			if ((int) gapIt->gapPos - (int) gapIt->seqPos != diff) {
				diff = (int) gapIt->gapPos - (int) gapIt->seqPos;
				appendValue(alignIt->gaps, *gapIt, Generous() );
			}
		}
		++contigReadIt;
	}
	// Update contig gaps
	if (includeReference) // causes problems ?
    {
		clear(fragStore.contigStore[0].gaps);
		typedef typename Iterator<TGapAnchor, Standard>::Type TGapIter;
		TGapIter gapIt = begin(contigReads[length(contigReads) - 1].gaps, Standard());
		TGapIter gapItEnd = end(contigReads[length(contigReads) - 1].gaps, Standard());
		int diff = 0;
		for(;gapIt != gapItEnd; ++gapIt) {
			if ((int) gapIt->gapPos - (int) gapIt->seqPos != diff) { // Only if gap not empty
				diff = (int) gapIt->gapPos - (int) gapIt->seqPos;
				appendValue(fragStore.contigStore[0].gaps, *gapIt, Generous() );
			}
		}
        fragStore.contigStore[0].gaps = contigReads[length(contigReads) - 1].gaps;
    }
	/*
    std::cout << "After realigning " << std::endl;
    for (unsigned i = 0; i < 200; ++i)
    {
        std::cout << "pos: " << i << "   ";
        for (unsigned j = 0; j < 4; ++j)
        {
            std::cout << consensus[i].count[j] << " ";
        }
        std::cout << "   ";
        for (unsigned j = 4; j < 8; ++j)
        {
            std::cout << consensus[i].count[j] << " ";
        }
        std::cout << "   ";
        for (unsigned j = 8; j < 12; ++j)
        {
            std::cout << consensus[i].count[j] << " ";
        }
        std::cout << std::endl;
    }
    */
    if(true)
    {
        std::cout << "After realigning " << std::endl;
        TContigGaps contigGaps3(fragStore.contigStore[0].seq, fragStore.contigStore[0].gaps);
        TContigPos maxPos3 = positionSeqToGap(contigGaps3,length(fragStore.contigStore[0].seq)-1)+1;
        maxPos = _max(maxPos3,(TContigPos)length(fragStore.contigStore[0].seq));
        std::cout << "maxPos visual = " << maxPos3 << std::endl;
        //std::cout << " genomeLen = " << genomeLen << std::endl;
        AlignedReadLayout layout;
        layoutAlignment(layout, fragStore);
        printAlignment(std::cout, layout, fragStore, 0, (TContigPos)0, (TContigPos)200, 0, 300);

        if (maxPos3 == 851 || maxPos3 == 422|| maxPos3 == 692)
        {
            CharString name = "/home/takifugu/sabrina7/Data/plotData/after_re";
            stringstream ss;
            ss << maxPos3;
            //ss << "_";
            //ss << length(fragStore.alignedReadStore);
            string str = ss.str();
            append(name, str);
            append(name, ".svg");
            SVGFile svg(toCString(name));
            printAlignment(svg, layout, fragStore, 0, (TContigPos)0, (TContigPos)200, 0, 300);
        }
    }

	if (includeReference)
		appendValue(fragStore.alignedReadStore, contigReads[length(contigReads) - 1]); // Because we need the beginPos for case of gaps at the beginning
}





#endif
