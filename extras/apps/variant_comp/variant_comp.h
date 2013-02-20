 /*==========================================================================
  SNP Calling routine of RazerS - Fast Read Mapping with Controlled Loss Rate
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

#ifndef SEQAN_HEADER_COMPAREVAR_H
#define SEQAN_HEADER_COMPAREVAR_H


#include <iostream>
#include <fstream>
#include <cmath>

#define EPS_PREC 0.00000001
#define MAX_DUPLICATION_SIZE 400000


namespace SEQAN_NAMESPACE_MAIN
{

    
    //////////////////////////////////////////////////////////////////////////////
// Default options


//____________________________________________________________________________
// Global Parameters



struct IndelCheck{};

template<typename TSpec = IndelCheck>
struct IndelCompareOptions
{
	CharString output; 	    // output file for statistics and shared indels
	CharString outputFN; 	    // output file for unmatched reference indels
	CharString outputSnp; 	        // output file for statistics and shared SNPs
	CharString outputSnpFN; 	    // output file for unmatched reference SNPs
	CharString inputReference;	// reference indels+SNPs in gff format
	CharString inputPredicted;	// predicted indels+SNPs in gff format

	int positionTolerance;		// max. deviation of predicted indel position from reference indel position
	double sizeTolerance;		    // max. deviation of predicted indel size from reference indel size (in percent)
	bool sequenceContext;
	int annotateRepeats;
    bool genotypeAware;
	
	const char *attachTag;          // tag to attach to overlap output
	int _debugLevel;
	String<Pair<int,int> > ranges;  // seperate statistics for each range
	
	IndelCompareOptions()
	{
		output = "";		
		outputFN = "";		
		outputSnp = "";		
		outputSnpFN = "";		
		inputReference = "";		
		inputPredicted = "";		
		
		positionTolerance  = 10;  
		sizeTolerance = 0.0;
		sequenceContext = false;
        genotypeAware = false;
		
		annotateRepeats = 50;	// 0 -> dont add "percentRepeat"-tag in result file
								// otherwise take +- x flanking sequence
		attachTag = "";
		_debugLevel = 0;
		appendValue(ranges,Pair<int,int>(-50,-6));
		appendValue(ranges,Pair<int,int>(-6,0));
		appendValue(ranges,Pair<int,int>(1,7));
		appendValue(ranges,Pair<int,int>(7,51));
		appendValue(ranges,Pair<int,int>(51,501));
		appendValue(ranges,Pair<int,int>(501,5001));
		
	}
	
};

enum SVTypes {
		SNP = 0,
		INSERTION = 1,
		DELETION = 2,
		INVERSION = 3,
        TRANSLOCATION = 4
	};

struct IndelInfo{
	
	unsigned genomeId;
	unsigned originalPos;
	unsigned secondGenomeId;
	unsigned secondOriginalPos;
	unsigned simPos;
	int indelSize;
    int type;
	double quality; // should be int (and real quality value)
	bool duplication;
    CharString genotype; // het or hom
	Dna5String insertionSeq;
	CharString idStr;
	CharString ninethCol;
	CharString field2;
	CharString field3;
	
	// int percentRepeat;
};


struct SnpInfo{
	
	unsigned genomeId;
	unsigned originalPos;
	unsigned simPos;
	int quality;
    CharString genotype; // het or hom
    CharString idStr;
	CharString ninethCol;
	CharString field2;
	CharString field3;
	
	// int percentRepeat;
};

/*
inline void
dumpSnpInfo(SnpInfo & snp)
{
    std::cout << "genomeId" << snp.genomeId << std::endl;
	std::cout << "originalPos" << snp.originalPos << std::endl;
	std::cout << "simPos" << snp.simPos << std::endl;
	std::cout << "quality" << snp.quality << std::endl;
    std::cout << "genotype" << snp.genotype << std::endl; // het or hom
    std::cout << "idStr" << snp.idStr << std::endl;
	std::cout << "ninethCol" << snp.ninethCol << std::endl;
	std::cout << "field2" << snp.field2 << std::endl;
	std::cout << "field3" << snp.field3 << std::endl;
	

}
*/


//////////////////////////// Write Functions ///////////////////////////////////

// write one gff line
template <typename TFile, typename TOptions>
bool write(TFile &file, IndelInfo &refIndel, CharString &genomeID, CharString &tagAppend, TOptions & options)
{
	if(!file.is_open()) return false;
	file << genomeID << '\t';
	
	if(!empty(refIndel.field2)) file << refIndel.field2 << '\t';
	else file << "variantCmp\t";

	if(!empty(refIndel.field3)) file << refIndel.field3 << '\t';
	else
	{
		if(refIndel.indelSize < 0) file << "insertion\t";
		if(refIndel.indelSize == 0) file << "baseexchange\t";
		if(refIndel.indelSize > 0) file << "deletion\t";
	}

	if(refIndel.indelSize <= 0) file << refIndel.originalPos + 1 << '\t' << refIndel.originalPos + 1 ;
	else file << refIndel.originalPos + 1 << '\t' << refIndel.originalPos + refIndel.indelSize ;
	
	file << "\t" << refIndel.quality << "\t+\t.\t";
	if(!empty(refIndel.ninethCol)) file << refIndel.ninethCol;
	else
	{
		file << "ID=" << refIndel.idStr << ";size=" << refIndel.indelSize;
        if(options.genotypeAware && refIndel.genotype==0) file << ";geno=hom";
        else if(options.genotypeAware && refIndel.genotype==1) file << ";geno=het";

//		if(refIndel.duplication) file << ";duplication=1";
	}
	if(refIndel.duplication) file << ";duplication=1";
	
	file << toCString(tagAppend) << std::endl;
	
	return true;
	
}


// write one gff line
template <typename TFile, typename TOptions>
bool write(TFile &file, SnpInfo &refSnp, CharString &genomeID, CharString &tagAppend, TOptions & options)
{
	if(!file.is_open()) return false;
	file << genomeID << '\t';
	
	if(!empty(refSnp.field2)) file << refSnp.field2 << '\t';
	else file << "variantCmp\t";

	if(!empty(refSnp.field3)) file << refSnp.field3 << '\t';
	else file << "snp\t";
		
	file << refSnp.originalPos + 1 << '\t' << refSnp.originalPos + 1 ;
    file << "\t" << refSnp.quality << "\t+\t.\t";
	if(!empty(refSnp.ninethCol)) file << refSnp.ninethCol;
	else
	{
		file << "ID=" << refSnp.idStr;
        if(options.genotypeAware && refSnp.genotype==0) file << ";geno=hom";
        else if(options.genotypeAware && refSnp.genotype==1) file << ";geno=het";
	}
	
	file << toCString(tagAppend) << std::endl;
	
	return true;
	
}

///////////////////////// Compare Functions ////////////////////////////////////////    
    
//____________________________________________________________________________

template <typename TIndel>
struct LessGPosSize : public ::std::binary_function < TIndel, TIndel, bool >
{
	inline bool operator() (TIndel const &a, TIndel const &b) const
	{
		if (a.genomeId < b.genomeId ) return true;
		if (a.genomeId > b.genomeId ) return false;

		if (a.originalPos < b.originalPos ) return true;
		if (a.originalPos > b.originalPos ) return false;

		return a.indelSize < b.indelSize;
	}
};

//____________________________________________________________________________


template <typename TVariant>
struct LessGPos : public ::std::binary_function < TVariant, TVariant, bool >
{
	inline bool operator() (TVariant const &a, TVariant const &b) const 
	{
		// genome sequence
		if (a.genomeID < b.genomeID) return true;
		if (a.genomeID > b.genomeID) return false;

        // genome position
        return (a.originalPos < b.originalPos);
	}
};

//____________________________________________________________________________

template <typename TIndel, typename TOptions>
bool compareIndelPair(TIndel &refIndel, TIndel &predIndel, TOptions &options)
{
    if(refIndel.type != predIndel.type) return false;
	
    if(options.genotypeAware && predIndel.genotype != refIndel.genotype)
        return false;

	int sizeTol = int((double)abs(refIndel.indelSize) * options.sizeTolerance);
	if(!(refIndel.indelSize - sizeTol <= predIndel.indelSize && predIndel.indelSize <= refIndel.indelSize + sizeTol))
		return false; // --> doesnt match
	
	// check if begin position matches
	if((int)refIndel.originalPos - options.positionTolerance <= (int) predIndel.originalPos
			&& predIndel.originalPos <= refIndel.originalPos + options.positionTolerance)
		return true; // matches
	
	// if it is a deletion or inversion
	if(refIndel.indelSize > 0) // check if end position matches
	{
		if (((int)refIndel.originalPos+refIndel.indelSize-options.positionTolerance 
				<= (int)predIndel.originalPos+predIndel.indelSize)
			&& ((int) predIndel.originalPos+predIndel.indelSize
				<= (int)refIndel.originalPos+refIndel.indelSize+options.positionTolerance))
			return true;
	}
	else // if it is an insertion
	{
		if(refIndel.duplication) // dirty
		{
	//		if (abs((int)refIndel.originalPos - (int)predIndel.originalPos) < -refIndel.indelSize + options.positionTolerance)
		//		return true;
	
		}
	}
	return false;
	
}

template <typename TIndel, typename TPos, typename TOptions>
bool compareIndelPair(TIndel &refIndel, TIndel &predIndel, TPos beginPoint, TPos endPoint, TOptions &options)
{
    if(refIndel.type != predIndel.type) return false;

    if(options.genotypeAware && predIndel.genotype != refIndel.genotype)
        return false;

    int sizeTol = int((double)abs(refIndel.indelSize) * options.sizeTolerance);
	if(!(refIndel.indelSize - sizeTol <= predIndel.indelSize && predIndel.indelSize <= refIndel.indelSize + sizeTol))
		return false; // --> doesnt match
	
	// check if begin position matches
	if(refIndel.originalPos >= (unsigned)beginPoint
		&& refIndel.originalPos <= (unsigned)endPoint)
		return true; // matches
	return false;
	
}




template<typename TIndel, typename TGenome, typename TPos>
void 
computeEir(TIndel &indel, TGenome & genome, TPos & beginPoint,TPos & endPoint)
{

	SEQAN_ASSERT(indel.indelSize > 0 || (int)length(indel.insertionSeq) == -indel.indelSize || indel.duplication == true);

	if( indel.type == DELETION  // deletion
        || (indel.type == INSERTION && indel.duplication == true && empty(indel.insertionSeq) && indel.secondOriginalPos != indel.originalPos)) // duplication, sequence not explicitly given
	{
		unsigned iDel = indel.originalPos;
		unsigned iRef = indel.originalPos + indel.indelSize;
		while(iRef < length(genome) && genome[iRef]==genome[iDel])
		{
			++iRef; ++iDel;
		}
		endPoint = iRef;
		iDel = indel.originalPos + indel.indelSize - 1;
		iRef = _max(0,(int)indel.originalPos - 1);
		while(iRef > 0 && genome[iRef] == genome[iDel])
		{
			--iRef; --iDel;
		}
		beginPoint = iRef;
	}
	else if(indel.type == INSERTION && !empty(indel.insertionSeq))
	{
        SEQAN_ASSERT_EQ((int)length(indel.insertionSeq),-indel.indelSize);
		unsigned iIns = 0; // relative to insertion sequence
		unsigned iRef = indel.originalPos;
//		std::cout << "right: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
		while(iRef < length(genome) && genome[iRef]==indel.insertionSeq[iIns])
		{
//			std::cout << "right: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
			++iRef; ++iIns;
			if(iIns >= length(indel.insertionSeq)) iIns = 0;
		}
		endPoint = iRef;
		iIns = -indel.indelSize - 1;
		iRef = _max(0,(int)indel.originalPos - 1);
//		std::cout << "left: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
		while(iRef > 0 && genome[iRef] == indel.insertionSeq[iIns])
		{
//			std::cout << "left: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
			if(iIns == 0) iIns = -indel.indelSize;
			--iRef; --iIns;
		}
		beginPoint = iRef;
	}
    else if (indel.type == INVERSION)
    {
   		unsigned iFor = _max(0,(int)indel.originalPos - 1);
		unsigned iRev = indel.originalPos + indel.indelSize;
        Dna5String r(genome[iRev]);
        reverseComplement(r);
        while(iRev < length(genome) && iFor > 0 && r[0]==genome[iFor])
		{
			++iRev; --iFor;
            r[0] = genome[iRev];
            reverseComplement(r);
        }
		endPoint = iRev;
        beginPoint = iFor;
    }
//	std::cout <<"begin=" << beginPoint << " end=" << endPoint << std::endl;
	

}

//////////////////////////////////////////////////////////////////////////////
// 1) build interval tree of reference indels
// 2) go through predicted indels one by one
template <
	typename TIndelSet,
	typename TGenome,
	typename TGenomeIDs,
	typename TOptions
>
int compareIndels(
	TIndelSet			&refIndels,		  // reference indels
	TIndelSet			&predIndels, // predicted indels
	TGenome				&genomes,
	TGenomeIDs			&genomeIDs,
	TOptions 			&options)	  	  // options
{
	
	typedef typename Value<TIndelSet>::Type    TIndel;
	typedef typename Iterator<TIndelSet>::Type TIndelIt;

	::std::ostringstream fileName;
	if (!empty(options.output))
		fileName << options.output;
	else
		fileName << options.inputPredicted << ".overlap.gff";
	
	::std::ofstream file;
	file.open(fileName.str().c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) {
		::std::cerr << "\nFailed to open output file" << ::std::endl;
		return 1;
	}


	// a set to keep track of which of the reference indels were matched
	::std::set<int> refFoundSet;
	typename ::std::set<int>::iterator foundSetIt;
	
	::std::sort(begin(predIndels),end(predIndels),LessGPosSize<TIndel>());
	::std::sort(begin(refIndels),end(refIndels),LessGPosSize<TIndel>());

	TIndelIt refIndelIt = begin(refIndels);
	TIndelIt refIndelsEnd = end(refIndels);
	TIndelIt predIndelIt = begin(predIndels);
	TIndelIt predIndelsEnd = end(predIndels);
	TIndelIt currPredGenomeEnd,currPredGenomeBegin,currRefGenomeEnd,currRefGenomeBegin;
	
	int TP = 0;

	// rememer stats for each range
	String<unsigned> rangeCountRef,rangeCountPred,rangeTP,rangeFound;
    unsigned predInversions = 0, refInversions = 0, inversionsTP = 0, inversionsFound = 0;
	resize(rangeCountRef,length(options.ranges),0);
	resize(rangeCountPred,length(options.ranges),0);
	resize(rangeTP,length(options.ranges),0);         // to calculate FP
	resize(rangeFound,length(options.ranges),0);  // to calculate FN
	for(unsigned j = 0; j < length(options.ranges); ++j)
	{
		for(TIndelIt it = begin(refIndels); it != end(refIndels); ++it)
        {
    		if((*it).type != INVERSION  && options.ranges[j].i1 <= (*it).indelSize && (*it).indelSize < options.ranges[j].i2)
	    		++rangeCountRef[j];
        }
		for(TIndelIt it = begin(predIndels); it != end(predIndels); ++it)
        {
	    	if((*it).type != INVERSION  && options.ranges[j].i1 <= (*it).indelSize && (*it).indelSize < options.ranges[j].i2)
   				++rangeCountPred[j];
        }
	}
    for(TIndelIt it = begin(refIndels); it != end(refIndels); ++it)
        if((*it).type == INVERSION) ++refInversions;
	for(TIndelIt it = begin(predIndels); it != end(predIndels); ++it)
        if((*it).type == INVERSION) ++predInversions;
 
	CharString tagAttach = ";";
	if(*options.attachTag != 0) append(tagAttach,options.attachTag);
	else append(tagAttach,"matchID");
	append(tagAttach,"=");

	// foreach genomeID
	//   build interval tree from ref indels
	//   while predictedIndels on current genomeID
	//     check indel interval +- positionTolerance
	//     foreach relevant reference indel
	//       check size +- sizeTolerance
	//       if(match) output refIndel-predictedIndel-combi
	//                 and increase TP counter
	
	if(options._debugLevel > 0) std::cout << "Starting to compare indels..." <<std::endl;
	for(unsigned i = 0; i < length(genomeIDs); ++i)
	{
		if(options._debugLevel > 1 ) std::cout << "." <<std::flush;

		//skip ahead if necessary
		while(refIndelIt != refIndelsEnd && (*refIndelIt).genomeId < i)
			++refIndelIt;
		while(predIndelIt != predIndelsEnd && (*predIndelIt).genomeId < i)
			++predIndelIt;
		
		// get range of ref variants that are on current chromosome
		TIndelIt currRefGenomeBegin = refIndelIt;
		while(refIndelIt != refIndelsEnd && (*refIndelIt).genomeId == i)
			++refIndelIt;
		TIndelIt currRefGenomeEnd = refIndelIt;
//		if(currRefGenomeEnd - currRefGenomeBegin == 0)
//			continue;
			
		// get range of predicted variants that are on current chromosome
		TIndelIt currPredGenomeBegin = predIndelIt;
		while(predIndelIt != predIndelsEnd && (*predIndelIt).genomeId == i)
			++predIndelIt;
		TIndelIt currPredGenomeEnd = predIndelIt;
		if(currPredGenomeEnd - currPredGenomeBegin == 0)
			continue;

			
		// build interval tree of all variant intervals
		typedef int TCargo;
		typedef int TValue;
		String<IntervalAndCargo<TValue,TCargo> > intervals;
		int ii = currRefGenomeBegin - begin(refIndels);
		for(refIndelIt = currRefGenomeBegin; refIndelIt != currRefGenomeEnd; ++refIndelIt)
		{
			if((*refIndelIt).indelSize > 0) 
				appendValue(intervals, IntervalAndCargo<TValue,TCargo>((*refIndelIt).originalPos,(*refIndelIt).originalPos + (*refIndelIt).indelSize, ii ));
			else
			{
				if((*refIndelIt).duplication)
					appendValue(intervals, IntervalAndCargo<TValue,TCargo>((*refIndelIt).originalPos + (*refIndelIt).indelSize,(*refIndelIt).originalPos + 1 - (*refIndelIt).indelSize, ii ));
				else
					appendValue(intervals, IntervalAndCargo<TValue,TCargo>((*refIndelIt).originalPos,(*refIndelIt).originalPos + 1, ii ));
			}
			++ii;
		}
		IntervalTree<TValue,TCargo> itree(intervals);
		if(options._debugLevel > 1 ) std::cout << "Interval tree done.\n" <<std::flush;
		
		// for each predicted variant check if it overlaps with a reference variant
		for(predIndelIt = currPredGenomeBegin; predIndelIt != currPredGenomeEnd; ++predIndelIt)
		{
			TIndel &predIndel = *predIndelIt;
			if(options._debugLevel > 1 ) std::cout << "." <<std::flush;

			// hack! needs to be fixed in snpStore, insertion position is one position to the left!
			// if(predIndel.indelSize < 0) predIndel.originalPos += 1;
			
			// "core" interval
			TValue beginPoint = predIndel.originalPos;
			TValue endPoint = predIndel.originalPos;
			if(predIndel.indelSize > 0) endPoint += predIndel.indelSize;
			else endPoint += 1;
			
			// add tolerance in genomic position
			if(options.sequenceContext)
			{
				//compute eir, equivalent indel region, Krawitz et. al
             //   if(predIndel.type == DELETION || (predIndel.type == INSERTION && !empty(predIndel.insertionSeq)) )
    				computeEir(predIndel,genomes[i],beginPoint,endPoint);
				beginPoint -= options.positionTolerance;
				endPoint += options.positionTolerance;
				SEQAN_ASSERT_LT(beginPoint,endPoint);
			}
			else
			{
				beginPoint -= options.positionTolerance;
				endPoint += options.positionTolerance;
			}
			
			// retrieve overlapping intervals
			String<TCargo> intersectingIntervals;
			findIntervals(itree,beginPoint,endPoint,intersectingIntervals);

			if(options._debugLevel > 1) std::cout <<"begin=" << beginPoint << " end=" << endPoint << std::endl;
			
			// check whether one of the overlapping intervals is a variant match
			bool predFound = false;
			CharString tagAppend = "";
			for(unsigned j = 0; j < length(intersectingIntervals); ++j)
			{
				if(options._debugLevel > 1 ) std::cout << "+" <<std::flush;
				TIndel &refIndel = refIndels[intersectingIntervals[j]];
				bool res = false;
				res = compareIndelPair(refIndel,predIndel,options);
				if(!res && options.sequenceContext) res = compareIndelPair(refIndel,predIndel,beginPoint,endPoint,options);
			//	if(!options.sequenceContext) res = compareIndelPair(refIndel,predIndel,options);
			//	else res = compareIndelPair(refIndel,predIndel,beginPoint,endPoint,options);
				if(res)
				{
					if(options.annotateRepeats>0)
					{
						int countN = 0;
						for(TValue k = (TValue)_max((int)0,(int)refIndel.originalPos - (int)options.annotateRepeats); k < (TValue)refIndel.originalPos && k < (TValue)length(genomes[i]); k++) // count 'N's in upstream flanking sequence
							if(genomes[i][k] == 'N') ++countN;
						for(TValue k = refIndel.originalPos; k < _min((TValue)length(genomes[i]),(TValue)refIndel.originalPos+options.annotateRepeats); k++) // count 'N's in upstream flanking sequence
							if(genomes[i][k] == 'N') ++countN;
						int total = refIndel.originalPos - _max((int)0,(int)refIndel.originalPos - (int)options.annotateRepeats);
						total += _min((TValue)length(genomes[i]),(TValue)refIndel.originalPos+options.annotateRepeats) - refIndel.originalPos;
						std::stringstream tagAppendStr;
						tagAppendStr << ";percRep="  << int(100*countN/total);
						append(tagAppend,tagAppendStr.str());
							
					}
//					CharString tagAppend = ";matchID=";
					if(!predFound) append(tagAppend,tagAttach);
					else append(tagAppend,",");
					append(tagAppend,refIndel.idStr);
					predFound = true;
					refFoundSet.insert(intersectingIntervals[j]);
				}
				
			}
			if(predFound)
			{
				++TP;
                if(predIndel.type == INVERSION)
                    ++inversionsTP;
                else
    			{
                    for(unsigned j = 0; j < length(options.ranges); ++j)
	    			{
                        if(options.ranges[j].i1 <= predIndel.indelSize && predIndel.indelSize < options.ranges[j].i2)
		    				++rangeTP[j];
                    }
                }
			}
			else
			{
				if(options._debugLevel>1)std::cout << "FP:" << predIndel.ninethCol << " begin=" << beginPoint << " end=" << endPoint << std::endl;
			}
			write(file,predIndel,genomeIDs[i],tagAppend,options);
			
		}
	
		if(options._debugLevel > 1 ) std::cout << std::endl;
	}
        // check number of false negatives

	for(foundSetIt = refFoundSet.begin(); foundSetIt != refFoundSet.end(); ++foundSetIt)
	{
        if(refIndels[*foundSetIt].type == INVERSION)
            ++inversionsFound;
        else   
            for(unsigned j = 0; j < length(options.ranges); ++j)
	    		if(options.ranges[j].i1 <= refIndels[*foundSetIt].indelSize && refIndels[*foundSetIt].indelSize < options.ranges[j].i2)
		    		++rangeFound[j];
                
	}

	if(!empty(options.outputFN))
    {
        ::std::ofstream fileFN;
	    fileFN.open(toCString(options.outputFN), ::std::ios_base::out | ::std::ios_base::trunc);
    	if (!fileFN.is_open()) {
	    	::std::cerr << "\nFailed to open FN output file" << ::std::endl;
	    	return 1;
    	}

		CharString tagAppend = ";notMatched";
        int currentId = 0;
        for(foundSetIt = refFoundSet.begin(); foundSetIt != refFoundSet.end(); ++foundSetIt)
    	{
            while (*foundSetIt > currentId)
            {
			    write(fileFN,refIndels[currentId],genomeIDs[refIndels[currentId].genomeId],tagAppend,options);
                ++currentId;
            }
            ++currentId;
        }            
        while (currentId < (int)length(refIndels))
        {
		    write(fileFN,refIndels[currentId],genomeIDs[refIndels[currentId].genomeId],tagAppend,options);
            ++currentId;
        }
        fileFN.close(); 
	}


	// optionally append non-overlapped reference indels to output file
	file << "###################################################" << ::std::endl;
	file << "# Total Stats: " << ::std::endl;
	file << "# Number of reference indels: " << length(refIndels) << ::std::endl;
	file << "# Number of predicted indels: " << length(predIndels) << ::std::endl;
	file << "# Number of TP predictions  : " << TP << ::std::endl;
	file << "# Number of FP predictions  : " << length(predIndels)-TP << ::std::endl;
	file << "# Number of FN predictions  : " << length(refIndels)-refFoundSet.size() << ::std::endl;

	for(unsigned j = 0; j < length(options.ranges); ++j)
	{
		file << "###################################################" << ::std::endl;
		file << "# Bucket Stats [" << options.ranges[j].i1 << "," << options.ranges[j].i2 << "): " << ::std::endl;
		file << "# Number of reference indels: " << rangeCountRef[j] << ::std::endl;
		file << "# Number of predicted indels: " << rangeCountPred[j] << ::std::endl;
		file << "# Number of TP predictions  : " << rangeTP[j] << ::std::endl;
		file << "# Number of FP predictions  : " << rangeCountPred[j]-rangeTP[j] << ::std::endl;
		file << "# Number of FN predictions  : " << rangeCountRef[j]-rangeFound[j] << ::std::endl;
	}

    file << "###################################################" << ::std::endl;
	file << "# Inversion Stats: " << ::std::endl;
	file << "# Number of reference inversions: " << refInversions << ::std::endl;
	file << "# Number of predicted inversions: " << predInversions << ::std::endl;
	file << "# Number of TP predictions  : " << inversionsTP << ::std::endl;
	file << "# Number of FP predictions  : " << predInversions-inversionsTP << ::std::endl;
	file << "# Number of FN predictions  : " << refInversions-inversionsFound << ::std::endl;

	file.close();

	if(options._debugLevel > 0)
	{
		std::cout << "# Total Stats: " << ::std::endl;
		std::cout << "# Number of reference indels: " << length(refIndels) << ::std::endl;
		std::cout << "# Number of predicted indels: " << length(predIndels) << ::std::endl;
		std::cout << "# Number of TP predictions  : " << TP << ::std::endl;
		std::cout << "# Number of FP predictions  : " << length(predIndels)-TP << ::std::endl;
		std::cout << "# Number of FN predictions  : " << length(refIndels)-refFoundSet.size() << ::std::endl;
	}
	// FN = length(referenceIndels) - TP
	// FP = length(predictedIndels) - TP
	// output statistics at end of file

	
	return 0;
	
}




template<typename TOptions>
bool
compareSnpPair(SnpInfo & predSnp, SnpInfo & refSnp, TOptions & options)
{
    if (predSnp.genomeId != refSnp.genomeId || predSnp.originalPos != refSnp.originalPos)
        return false;
    if(options.genotypeAware && predSnp.genotype != refSnp.genotype)
        return false;
    return true;
}


//////////////////////////////////////////////////////////////////////////////
// 1) sort SNP sets according to genome and position
// 2) scan over them and record shared/missing entries
template <
	typename TSnpSet,
	typename TGenome,
	typename TGenomeIDs,
	typename TOptions
>
int compareSnps(
	TSnpSet			&refSnps,		  // reference SNPs
	TSnpSet			&predSnps, // predicted SNPs
	TGenome				&,
	TGenomeIDs			&genomeIDs,
	TOptions 			&options)	  	  // options
{
	
	//typedef typename Value<TSnpSet>::Type    TSnp;
	typedef typename Iterator<TSnpSet>::Type TSnpIt;

	::std::ostringstream fileName;
	if (!empty(options.outputSnp))
		fileName << options.outputSnp;
	else
		fileName << options.inputPredicted << ".snp_overlap.gff";
	
	::std::ofstream file;
	file.open(fileName.str().c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) {
		::std::cerr << "\nFailed to open SNP output file" << ::std::endl;
		return 1;
	}
    ::std::ofstream fileFN;
   	if(!empty(options.outputSnpFN))
    {
	    fileFN.open(toCString(options.outputSnpFN), ::std::ios_base::out | ::std::ios_base::trunc);
    	if (!fileFN.is_open()) {
	    	::std::cerr << "\nFailed to open SNP FN output file" << ::std::endl;
	    	return 1;
    	}
    }

	int predTP = 0; // true positive predicted SNPs
    int refTP = 0;  // recovered refset SNPs  // actually for SNPs we should have predTP==refTP

	TSnpIt refSnpIt = begin(refSnps);
	TSnpIt refSnpsEnd = end(refSnps);
	TSnpIt predSnpIt = begin(predSnps);
	TSnpIt predSnpsEnd = end(predSnps);
	TSnpIt currPredGenomeEnd,currPredGenomeBegin,currRefGenomeEnd,currRefGenomeBegin;

	CharString tagAppendEmpty = "";

	if(options._debugLevel > 0) std::cout << "Starting to compare SNPs..." <<std::endl;
	for(unsigned i = 0; i < length(genomeIDs); ++i)
	{
		if(options._debugLevel > 1 ) std::cout << "." <<std::flush;

		//skip ahead if necessary
		while(refSnpIt != refSnpsEnd && (*refSnpIt).genomeId < i)
			++refSnpIt;
		while(predSnpIt != predSnpsEnd && (*predSnpIt).genomeId < i)
			++predSnpIt;
		
		// get range of ref variants that are on current chromosome
		TSnpIt currRefGenomeBegin = refSnpIt;
		while(refSnpIt != refSnpsEnd && (*refSnpIt).genomeId == i)
			++refSnpIt;
		TSnpIt currRefGenomeEnd = refSnpIt;
//		if(currRefGenomeEnd - currRefGenomeBegin == 0)
//			continue;
			
		// get range of predicted variants that are on current chromosome
		TSnpIt currPredGenomeBegin = predSnpIt;
		while(predSnpIt != predSnpsEnd && (*predSnpIt).genomeId == i)
			++predSnpIt;
		TSnpIt currPredGenomeEnd = predSnpIt;
		if(currPredGenomeEnd - currPredGenomeBegin == 0)
			continue;
        predSnpIt = currPredGenomeBegin;
        refSnpIt = currRefGenomeBegin;

		// for each predicted SNP, check if matched
        for(; predSnpIt != currPredGenomeEnd; ++predSnpIt)
		{
            while(refSnpIt != currRefGenomeEnd && (*predSnpIt).originalPos > (*refSnpIt).originalPos) //refSnp iterator is behind
            {
                if(!empty(options.outputSnpFN))// SNP was not matched
                    write(fileFN,*predSnpIt,genomeIDs[i],tagAppendEmpty,options);
                ++refSnpIt;
                //FN++
            }
            if(refSnpIt != currRefGenomeEnd && (*predSnpIt).originalPos < (*refSnpIt).originalPos) //unmatched prediction
            {
       			write(file,*predSnpIt,genomeIDs[i],tagAppendEmpty,options);
                continue;
            }
            bool first = true;
            while(refSnpIt != currRefGenomeEnd && (*predSnpIt).originalPos == (*refSnpIt).originalPos) // there could be more than one refSnp (different genotypes)
            {
                bool res = compareSnpPair(*predSnpIt,*refSnpIt,options);
	            if(res) // match
                {
                    ++refTP; //count TP wrt ref set
                    if(first) //only write to file once
                    {
                        CharString tagAttach = ";";
                    	append(tagAttach,"matchID=");
                    	append(tagAttach,(*refSnpIt).idStr);
       		    	    write(file,*predSnpIt,genomeIDs[i],tagAttach,options);
                        ++predTP; //count TP prediction only once
                        first = false;
                    }
                }
                ++refSnpIt;
            }
            if(first && !empty(options.outputSnpFN))// SNP was not matched
                write(fileFN,*predSnpIt,genomeIDs[i],tagAppendEmpty,options);
                     
        }

    }

	if(!empty(options.outputSnpFN))
        fileFN.close(); 


	// optionally append non-overlapped reference indels to output file
	file << "###################################################" << ::std::endl;
	file << "# Total Stats: " << ::std::endl;
	file << "# Number of reference SNPs: " << length(refSnps) << ::std::endl;
	file << "# Number of predicted SNPs: " << length(predSnps) << ::std::endl;
	file << "# Number of TP predictions  : " << predTP << ::std::endl;
	file << "# Number of FP predictions  : " << length(predSnps)-predTP << ::std::endl;
	file << "# Number of FN predictions  : " << length(refSnps)-refTP << ::std::endl;

	file.close();

	if(options._debugLevel > 0)
	{
	    ::std::cout << "# Total Stats: " << ::std::endl;
	    ::std::cout << "# Number of reference SNPs: " << length(refSnps) << ::std::endl;
	    ::std::cout << "# Number of predicted SNPs: " << length(predSnps) << ::std::endl;
	    ::std::cout << "# Number of TP predictions  : " << predTP << ::std::endl;
	    ::std::cout << "# Number of FP predictions  : " << length(predSnps)-predTP << ::std::endl;
	    ::std::cout << "# Number of FN predictions  : " << length(refSnps)-refTP << ::std::endl;
	}
	
	return 0;
	
}






}

#endif
