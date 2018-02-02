 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by Anne-Katrin Emde

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

#ifndef SEQAN_HEADER_READ_SIMULATOR_H
#define SEQAN_HEADER_READ_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>

#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/modifier.h>


namespace seqan
{

template <typename TOperation, typename TAlphabet>
inline TAlphabet
sample(std::mt19937 & rng, TOperation m, TAlphabet base)
{
    std::uniform_int_distribution<int> dist(0, 2);

    int ret = Dna(dist(rng));
    if (m == SEQAN_MISMATCH && ret >= ordValue(base))
        ret += 1;
    return Dna(ret);
}

template < typename TGenome >
void simulateGenome(std::mt19937 & rng,
                    TGenome &genome, int size)
{
//	mtRandInit(); 
	resize(genome, size);
	for(int i = 0; i < size; ++i)
		genome[i] = sample(rng, 0, (Dna)0);
}


template<typename TPosString>
void
fillupStartpos(
    std::mt19937 & rng,
    TPosString & sortedStartPos, 
	int numReads,
	int readLength,
	int maxErrors,			// how many errors they may have
	int libSize,			// library size, 0 disables mate-pair simulation
	int libError,
	int seqLength,			// maximal library size deviation
	double forwardProb)
{	
	
	const int REVCOMP = 1 << (sizeof(int)*8-1);
	resize(sortedStartPos, numReads);

	int fragmentSize = readLength;
	if (libSize > 0)
		fragmentSize = libSize + libError;


	// sample positions
	std::uniform_real_distribution<double> dist(0, 1);
	for (int i = 0; i < numReads; ++i)
		sortedStartPos[i] = (int)((seqLength - fragmentSize - maxErrors + 1.0) * dist(rng));

	std::sort(begin(sortedStartPos),end(sortedStartPos));

	// sample orientations
	for (int i = 0; i < numReads; ++i)
	    if (dist(rng) >= forwardProb)
			sortedStartPos[i] |= REVCOMP;

	if (libSize > 0)
	{
		resize(sortedStartPos,2*numReads);
		// sample mate-pair positions and inverse orientations
		for(int i=0;i<numReads;++i)
		{
			int leftPos = sortedStartPos[i] & ~REVCOMP;
			int lSize = fragmentSize - (int)((2.0 * libError + 1.0) * dist(rng));
			int rightPos = leftPos + lSize - readLength;
			if ((sortedStartPos[i] & REVCOMP) == 0)
			{
				sortedStartPos[i+numReads] = rightPos | REVCOMP;
			}
			else
			{
				sortedStartPos[i] = rightPos | REVCOMP;
				sortedStartPos[i+numReads] = leftPos;
			}
		}
		numReads*=2;
	}
}



//////////
// Simulates a set of reads from a set of haplotypes with a certain error distribution
//////////
template < 
	typename TReadSet,
	typename TReadIDs,
	typename TGenomeSet,
	typename TDistr >
void simulateReads(
            std::mt19937 & rng,
	TReadSet &readSet,		// generated read sequences
	TReadIDs &readIDs,		// corresponding Fasta ids
	TGenomeSet &genomeSet,	// source genome sequences
	int numReads,			// how many reads should be generated
	int maxErrors,			// how many errors they may have
	TDistr &errorDist,		// error probability distribution
	int libSize,			// library size, 0 disables mate-pair simulation
	int libError,			// maximal library size deviation
	double forwardProb,
	bool verbose = false)
{
	//typedef typename Value<TReadSet>::Type				TRead;
	//typedef typename Value<TGenomeSet>::Type			TGenome;
	//typedef typename Infix<TGenome>::Type				TGenomeInfix;
	//typedef ModifiedString<TGenomeInfix, ModReverse>	TGenomeInfixRev;
	//typedef ModifiedString<TRead, ModReverse>			TReadRev;

	//typedef Finder<TGenomeInfix>						TMyersFinder;
	//typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;

//	typedef Finder<TGenomeInfix>						TMyersFinderRev;
	//typedef Finder<TGenomeInfixRev>						TMyersFinderRev;
	//typedef Pattern<TReadRev, MyersUkkonenGlobal>		TMyersPatternRev;

//	mtRandInit();
//	typedef TGenome TRevComp;

	int readLength = length(errorDist)/4;
	const int REVCOMP = 1 << (sizeof(int)*8-1);
	//int KJ = 2*maxErrors;

	String<int> bucketCounter;
	resize(bucketCounter,maxErrors,0);
	
	String<int> kickOutCount;
	resize(kickOutCount,maxErrors,0);

	if (verbose)
		std::cout << "\nSimulating...";

	String<int> modificationPattern;
	reserve(modificationPattern, readLength + maxErrors);
	int inValidModPat = 0;

	//at the moment: nur eine source sequenz --> nimm immer genomeSet[0]
	TGenome& currentSource = genomeSet[0];
	int seqLength = length(currentSource);

	// int fragmentSize = readLength;
	// if (libSize > 0)
	// 	fragmentSize = libSize + libError;

	String<int> sortedStartPos;
/*				# Pick a library size
				currentLibrary = sample(1:(length(librarySizes)), 1)
				lSize = round(rnorm(1, mean=librarySizes[currentLibrary], sd=librarySd[currentLibrary]))
				if (start < end) {
					if (start + libSize <= seqLength) {
						startMatePair = start + libSize - readLength + 1 
						endMatePair = startMatePair + readLength - 1
						readMate = sourceSeq[startMatePair:endMatePair]
						readMate = reverseComplement(readMate)
						tmp = startMatePair 
						startMatePair = endMatePair 
						endMatePair = tmp
						invalidMatePair = 0
					}
				} else {
					if (start - libSize >= 1) {
						startMatePair = start - libSize  
						endMatePair = startMatePair + readLength - 1
						readMate = sourceSeq[startMatePair:endMatePair]
						invalidMatePair = 0
					}
				}*/

	int realNumReads = numReads;
	unsigned int samplePosCounter = 0;
	int readCounter = 0;
	while (readCounter < numReads) {
		clear(modificationPattern);
		//# Pick a haplotype
		//currentHaplotype = sample(1:(length(haplotypes)), 1)

		// Sample a read
		if(samplePosCounter >= length(sortedStartPos))//get new start positions
		{
			fillupStartpos(rng, sortedStartPos, numReads, readLength, maxErrors,	libSize, libError, seqLength, forwardProb);
			samplePosCounter = 0;
		}
		int  startPos = sortedStartPos[samplePosCounter] & ~REVCOMP;
		bool revComp  = (sortedStartPos[samplePosCounter] & REVCOMP) != 0;
		int  maxEnd   = startPos + readLength + maxErrors - 1;
	
		TGenome read;
		resize(read,readLength);
		
		TGenome readTemplate;// = infix(currentSource,startPos,maxEnd);
		resize(readTemplate,maxEnd-startPos);
		arrayCopy(iter(currentSource,startPos), iter(currentSource,maxEnd), begin(readTemplate)); //infix(currentSource,startPos,maxEnd);
		
		if(revComp) reverseComplement(readTemplate);

		// int lastOp = 0;
		int currOp = 0;
		// Sequence the reads
		int countErrors = 0;
		int pos = 0;
		int trueLength = 0;
		bool successful = false;

		std::uniform_real_distribution<double> dist(0, 1);

		while(pos < maxEnd) {
			// lastOp = currOp;
			double prob = dist(rng);
			//	std::cout << "prob = " << prob << "\t";
			int m;
			for(m = 0; m < 4; ++m)
			{
				double modProb = _transformBack(errorDist[m*readLength + trueLength]);
				if (prob < modProb)
				{
					currOp = m;
					break;
				}
				prob -= modProb;
			}
			if (m == 4) std::cout << "HUH?";
		//	std::cout << "operation = " << operation << "\t";
			if(pos==0 &&  currOp == SEQAN_INSERT)// (currOp==SEQAN_DELETE || currOp == SEQAN_INSERT))
			{
				currOp = 0;
				continue;
			}
			appendValue(modificationPattern,currOp,Generous());

			// ignore reads with Ns
			if (currOp != SEQAN_INSERT && readTemplate[pos] == 'N')
				countErrors = maxErrors + 1;

			// Insert Delete is the same as Delete Insert and both are the same as Mismatch (or match)
			if(currOp == SEQAN_MATCH) read[trueLength] = readTemplate[pos];
			else
			{
				++countErrors;
				if(currOp != SEQAN_INSERT)
					read[trueLength] = sample(rng, currOp,readTemplate[pos]);
			}
			if(currOp != SEQAN_INSERT) ++trueLength; //if read nucleotide is not deleted
			if(currOp != SEQAN_DELETE) ++pos; //if read nucleotide is not an insert
//			if((lastOp==SEQAN_DELETE && currOp==SEQAN_INSERT) || (currOp==SEQAN_DELETE && lastOp==SEQAN_INSERT))
//			{
//				--countErrors;
//				currOp = SEQAN_MISMATCH;
//			//	std::cout << "ID=DI=M\n";
//			}
		//	std::cout << "true len = " << trueLength << std::endl;
			if(trueLength == readLength || countErrors >= maxErrors)
			{
				if(countErrors < maxErrors)// && currOp != SEQAN_INSERT)
					successful = true;
				break;
			}
		}
	
	/*
		if(successful)
		{
			int patLen = length(modificationPattern);
			int err = 0, del = 0;
			for (int j = 0; j < KJ; ++j)
			{
				switch (modificationPattern[j]) {
					case SEQAN_MATCH:
						++del;
						break;
	
					case SEQAN_DELETE:
						++del;
						SEQAN_FALLTHROUGH
		
					case SEQAN_INSERT:
						++err;
						break;
	
					default:;
				}
				if (del > 0 && del <= err)
				{
					successful = false;
					++inValidModPat;
					break;
				}
			}
			if(successful)
			{
				err = del = 0;
				for (int j = patLen - 1; j >= patLen - KJ; --j)
				{
					switch (modificationPattern[j]) {
						case SEQAN_MATCH:
							++del;
							break;
		
						case SEQAN_DELETE:
							++del;
							SEQAN_FALLTHROUGH
			
						case SEQAN_INSERT:
							++err;
							break;
					
						default:;
					}
					if (del > 0 && del <= err)
					{
						successful = false;
						++inValidModPat;
						break;
					}
				}
			}
		}
		*/
		
/*		countMateErrors = 0
		if (simulateMatePairs == 1) {
			for(pos in 1:length(readMate)) {
				if (runif(1) <= _transformBack(errorDist[pos])) {
					readMate[pos] = sample(alphabet[alphabet!= readMate[pos]], 1)
					countMateErrors = countMateErrors + 1
				}
			}
		}*/
	
		if(successful)
		{
			//verify that number of errors is correct
			bool kickOut = false;
/*			int start1 = startPos;
			int maxEnd1 = maxEnd;
			while(start1 > 0 && (startPos - start1) < countErrors) --start1;
			while(maxEnd1 > 0 && (maxEnd1 - maxEnd) < countErrors) ++maxEnd1;
			TGenomeInfix genomeInfix(currentSource,start1,maxEnd1);
			TMyersFinder myersFinder(genomeInfix);

			// init forward verifiers
			if(revComp) reverseComplement(read);
			TMyersPattern forwardPattern(read);
			TMyersPattern &myersPattern = forwardPattern;
			
			// find end of best semi-global alignment
			int maxScore = std::numeric_limits<int>::min();
			int minScore = -(int)countErrors;
			TMyersFinder maxPos;
			while (find(myersFinder, myersPattern, minScore))
				if (maxScore < getScore(myersPattern)) 
				{
					maxScore = getScore(myersPattern);
					maxPos = myersFinder;
				}
			
			if (maxScore >= minScore) 
			{
				TGenomeInfixRev		infRev(infix(currentSource, start1, start1+position(maxPos)+1));
				TReadRev		readRev(read);
				TMyersFinderRev		myersFinderRev(infRev);
				TMyersPatternRev	myersPatternRev(readRev);
				// find beginning of best semi-global alignment
				if (find(myersFinderRev, myersPatternRev, maxScore))
					start1 = start1 + position(maxPos) - (position(myersFinderRev) + 1);
				else {
					// this case should never occur
					if(revComp) std::cerr <<"reverse\n";
					std::cerr << "posMaxpos = " << position(maxPos) << std::endl;
					std::cerr << "startPos = " << startPos << std::endl;
					std::cerr << "maxEnd   = " << maxEnd << std::endl;
					std::cerr << "HUH?\n" << std::endl;
					std::cerr << "fGENOME: " << host(myersFinder) << std::endl;
					std::cerr << "fREAD:   " << read << std::endl;
					std::cerr << "iGENOME: " << infix(currentSource, start1,start1+position(maxPos)+1) << std::endl;
					std::cerr << "rGENOME: " << infRev << std::endl;
					std::cerr << "rREAD:   " << readRev << std::endl;
				}
			
			} 
			if(revComp) reverseComplement(read);
			SEQAN_ASSERT(maxScore >= -(int)countErrors);
			if(maxScore != -(int)countErrors)
				kickOut = true;*/
			if(!kickOut)
			{
				std::stringstream id;
				resize(read,trueLength);
				++bucketCounter[countErrors];
				++readCounter;
				if(verbose && readCounter%10000 == 0)std::cout << readCounter<<"..." << std::flush;
				//Add read to readSet
				if(!revComp) id << startPos << ',' << startPos+pos;
				else id << maxEnd << ',' << maxEnd-pos;
				id << "[id=" << readCounter << ",fragId=" << readCounter % realNumReads;
				id << ",repeatId=" << 0 <<",errors=" << countErrors;
				if (revComp) id << ",orientation=R]";
				else id << ",orientation=F]";

				appendValue(readIDs, id.str(),Generous());
				appendValue(readSet, read, Generous());
			}
			else ++kickOutCount[countErrors];

		}
		++samplePosCounter;
//		else std::cout << "Not successful\n";
	}
//	if (simulateMatePairs == 1) {
//		print(bucketCounter / (2* numOfReads))
//	} else {
	if (verbose)
	{
		std::cout << "\n\nBucket frequencies:\n";
		for(unsigned i = 0; i < length(bucketCounter); ++i)
			std::cout << (double) bucketCounter[i] / numReads << std::endl;
		std::cout << std::endl;
		std::cout << "\nBucket kickout count:\n";
		for(unsigned i = 0; i < length(kickOutCount); ++i)
		{
			if((kickOutCount[i] + bucketCounter[i]) > 0) std::cout << (double) kickOutCount[i] / (kickOutCount[i] + bucketCounter[i]) << std::endl;
			else std::cout << "0\n";
		}
		std::cout << std::endl;
		
		std::cout << "\nInvalid modification pattern count: "<<inValidModPat<<std::endl;
	}

//	if (simulateMatePairs == 1) {
//		write(scan(tmpPath, what = 'character'), file=readPath, sep=std::endl, append = TRUE)
//		unlink(tmpPath)
//	}
}

}

#endif
