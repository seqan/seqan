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

#ifndef SEQAN_HEADER_PARAMCHOOSER_H
#define SEQAN_HEADER_PARAMCHOOSER_H

#include <iostream>
#include <fstream>
#include <sstream>
//#include <sys/types.h>
#include <errno.h>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include "razers.h"
#include "readSimulator.h"

namespace seqan
{
// ls in directory dir, store filenames in files
template<typename TPath, typename TFilenameString>
int getDir(TPath path, TFilenameString &files)
{
	Directory dir(path);
	if (!dir)
	{
		::std::cout << "Error(" << errno << ") opening " << dir << ::std::endl;
		return errno;
	}

	for (; !atEnd(dir); goNext(dir))
		appendValue(files, value(dir));
	return 0;
}


struct ParamChooserOptions
{
	typedef double TFloat;
	unsigned minThreshold;					// minimum value for threshold parameter 
	unsigned maxWeight;                                           // maximum value of q
	bool chooseOneGappedOnly;      // choose onegapped (or ungapped) shape (discard all other gapped shapes)
	bool chooseUngappedOnly;      // choose ungapped shape (discard all gapped shapes)
	bool useDefaultShapes;

	// global input parameters
	unsigned totalN;				// sequence length
	unsigned totalK;					// errors
	TFloat optionLossRate;		// in
	TFloat chosenLossRate;			// out
	TFloat optionErrorRate;		// 
	bool optionHammingOnly;

	bool extrapolate;
	unsigned extrapolN;
	unsigned extrapolK;
	unsigned maxComputedEditN;
	unsigned maxComputedHammingN;

	TFloat optionProbINSERT;
	TFloat optionProbDELETE;

	CharString fparams;
	CharString fgparams;

	bool fnameCount0;
	bool fnameCount1;
	bool prefixCount;
	const char *fname[2];
	const char *fprefix[1];
	CharString paramFolderPath;
	CharString paramFolder;
	const char *shapeFile;

	int qualityCutoff;
	bool solexaQual;

	bool appendToPrevious;
	bool verbose;

	ParamChooserOptions()
	{
		minThreshold = 1;					// minimum value for threshold parameter 
		maxWeight = 14;                                           // maximum value of q

		// global input parameters
		totalN = 32;				// sequence length
		totalK = 2;					// errors
		optionLossRate = (TFloat)0.01;		// in
		chosenLossRate = (TFloat)0.0;		// out
		optionErrorRate = (TFloat)0.05;		// 
		optionHammingOnly = false;

		extrapolate = false;
		extrapolN = 65;
		extrapolK = 4;
		maxComputedEditN = 59;
		maxComputedHammingN = 75;

		chooseOneGappedOnly = false;      // choose onegapped (or ungapped) shape (discard all other gapped shapes)
		chooseUngappedOnly = false;      // choose ungapped shape (discard all gapped shapes)
		useDefaultShapes = true;
		
		optionProbINSERT = (TFloat)0.0;
		optionProbDELETE = (TFloat)0.0;

		qualityCutoff = 20;
		fnameCount0 = 0;
		fnameCount1 = 0;
		prefixCount = 0;
		fname[0] = "";
		fname[1] = "";
		fprefix[0] =  "" ;
		
		shapeFile = "";
		paramFolderPath = "";
		paramFolder = "";
	
		appendToPrevious = false;
		verbose = true;
		solexaQual = true;
	}
};




template<typename TValue>
inline TValue
_convertSolexaQual2ErrProb(TValue sq)
{
	return pow((TValue)10, sq / (TValue)-10) / ((TValue)1 + pow((TValue)10, sq / (TValue)-10));
}

template<typename TValue>
inline TValue
_convertPhredQual2ErrProb(TValue sq)
{
	return pow((TValue)10, sq / (TValue)-10);
}

template<typename TValue>
inline TValue
_convertSolexaQual2PhredQual(TValue sq)
{
	return (TValue)10 * log((TValue)1 + pow((TValue)10, sq / (TValue)10)) / log((TValue)10);
}



//compute average position dependent error distribution (assumes solexa qualtiy values in prb.txt format)
template<typename TFile, typename TDistribution>
void
qualityDistributionFromPrbFile(TFile & file, TDistribution & avg, ParamChooserOptions & pm_options)
{
//IOREV see comments in other paramChooser.h
	typedef typename Value<TDistribution>::Type TFloat;

	String<TFloat> qualitySum;
	String<int> count;
	resize(qualitySum,pm_options.totalN,(TFloat)0.0);
	resize(count,pm_options.totalN,0);

    typename DirectionIterator<TFile, Input>::Type reader(file);
    if (atEnd(reader))
        return;

    skipUntil(reader, NotFunctor<IsWhitespace>());

	int kickout = 0;
	String<int> tempReadQual;
	resize(tempReadQual,pm_options.totalN);

    seqan::CharString buffer;
    while (!atEnd(reader))
	{
		int avgReadQual = 0;
		for (unsigned pos = 0; !atEnd(reader) && (pos < pm_options.totalN); ++pos)
		{
            int quals[4] = {0, 0, 0, 0};

            for (int i = 0; i < 4; ++i)
            {
                clear(buffer);
                skipUntil(reader, NotFunctor<IsWhitespace>());
                readUntil(buffer, reader, IsWhitespace());
                quals[i] = (int) lexicalCast<double>(buffer);
            }

			int qual = std::max(std::max(quals[0], quals[1]), std::max(quals[2], quals[3]));

			avgReadQual += qual;
			tempReadQual[pos] = qual;

//			::std::cout << qual << " ";
//			f = _convertSolexaQual2ErrProb(f);

//			qualitySum[pos] += _convertSolexaQual2ErrProb((TFloat)qual);
//			++count[pos];
		}
		if((int)(avgReadQual/pm_options.totalN) < pm_options.qualityCutoff) 
		{
			++kickout;
            skipLine(reader);
			continue;
		}
		else{
			for (unsigned pos = 0; (pos < pm_options.totalN); ++pos)
			{
	//			qualitySum[pos] += _convertSolexaQual2ErrProb((TFloat)qual);
				qualitySum[pos] += tempReadQual[pos];
				++count[pos];
			}
		}
//		::std::cout << ::std::endl;
        skipLine(reader);
	}
	if (pm_options.verbose)
    {
        ::std::cout << " Readcount = " << count[0] << "\t";
        ::std::cout << " kicked out " << kickout << " low quality reads." << std::endl;
    }

	resize(avg,pm_options.totalN,(TFloat)0.0);
	for(unsigned t = 0; t < pm_options.totalN; ++t)
	{
		TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
 		f = _convertSolexaQual2ErrProb(f);
		avg[t] = f;
	}
}


template <typename TFile, typename TDistribution>
void
qualityDistributionFromFastQFile(TFile & file, TDistribution & avg, ParamChooserOptions & pm_options)
{
//IOREV _duplicate_ this should use the existing fastq implementation
    typedef typename Value<TDistribution>::Type TFloat;

    String<int> qualitySum, count;
    resize(qualitySum, pm_options.totalN, 0);
    resize(count, pm_options.totalN, 0);

    SeqFileIn reader(file);
    CharString fastaId, seq, qual;
    while (!atEnd(reader))
    {
        readRecord(fastaId, seq, qual, reader);

        for (unsigned i = 0; i != pm_options.totalN && i != length(qual); ++i)
        {
            qualitySum[i] += qual[i] - 33;
            count[i]++;
        }
    }
    if (pm_options.verbose)
        std::cout << " Readcount = " << count[0] << std::endl;

    resize(avg, pm_options.totalN, (TFloat)0.0);
    for (unsigned t = 0; t < pm_options.totalN; ++t)
    {
        TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
        if (pm_options.solexaQual)
            f = _convertSolexaQual2PhredQual(f);
        f = _convertPhredQual2ErrProb(f);
        avg[t] = f;
    }
}

// find all *_prb.txt files in directory prbPath and compute average position dependent quality distribution
// compute average over all averages and store in errorDistribution
template <typename TPath, typename TError>
void
getAvgFromPrbDirectory(TPath prbPath, TError & errorDistribution, ParamChooserOptions & pm_options)
{
    typedef typename Value<TError>::Type TFloat;

    resize(errorDistribution, pm_options.totalN, (TFloat)0.0);

    String<std::string> files;
    getDir(prbPath, files);
    unsigned countPrbs = 0;
    for (unsigned int i = 0; i < length(files); i++)
    {
        if (suffix(files[i], length(files[i]) - 6) == ".fastq")
        {
            std::cout << "Processing " << files[i] << "..." << std::endl;
            TError avg_act;
            resize(avg_act, pm_options.totalN);
            std::fstream filestrm;
            std::stringstream sstrm;
            sstrm << prbPath << files[i];
            filestrm.open(sstrm.str().c_str(), std::ios_base::in);
            qualityDistributionFromFastQFile(filestrm, avg_act, pm_options);
            filestrm.close();
            for (unsigned j = 0; j < pm_options.totalN; ++j)
            {
//				std::cout << " " << avg_act[j];
                errorDistribution[j] += avg_act[j];
            }
            ++countPrbs;
            continue;
        }
        if (suffix(files[i], length(files[i]) - 8) == "_prb.txt")
        {
            std::cout << "Processing " << files[i] << "..." << std::endl;
            TError avg_act;
            resize(avg_act, pm_options.totalN);
            std::fstream filestrm;
            std::stringstream sstrm;
            sstrm << prbPath << files[i];
            filestrm.open(sstrm.str().c_str(), std::ios_base::in);
            qualityDistributionFromPrbFile(filestrm, avg_act, pm_options);
            filestrm.close();
            for (unsigned j = 0; j < pm_options.totalN; ++j)
            {
//				std::cout << " " << avg_act[j];
                errorDistribution[j] += avg_act[j];
            }
            ++countPrbs;
            continue;
        }
    }
    for (unsigned j = 0; j < pm_options.totalN; ++j)
        errorDistribution[j] /= (TFloat)countPrbs;
    std::cout << "Writing average error probabilities to " << pm_options.fprefix[0] << "_errorProb.dat" << std::endl;
    std::fstream out;
    std::stringstream avgOut;
    avgOut << pm_options.fprefix[0] << "_errorProb.dat";
    out.open(avgOut.str().c_str(), std::ios_base::out);
    if (!out.is_open())
        std::cout << "Couldn't write to file " << avgOut.str() << std::endl;
    else
        for (unsigned j = 0; j < pm_options.totalN; ++j)
            out << errorDistribution[j] << std::endl;
    out.close();


}

//////////////////////////////////////////////////////////////////////////////
// Returns the approximated minimum coverage of a one-gapped shape with weight q, span s at threshold t
template<typename TValueQ, typename TValueS, typename TValueT>
inline TValueS getMinCov(TValueQ q, TValueS s, TValueT t)
{
	TValueS mincov;
	if(t > s - q + 1){
		mincov = q + 2 * (t - 1) - (t - (s - q + 1));
	}
	else mincov = q + 2 * (t - 1);

	return mincov;
}

template <typename TShape, typename TIter>
inline void
readShape(TShape &shape, TIter &iter)
{
    typedef AssertFunctor<OrFunctor<EqualsChar<'0'>, EqualsChar<'1'> >, ParseError> TShapeAsserter;

    clear(shape);
    readUntil(shape, iter, IsWhitespace(), TShapeAsserter());
}

template <typename TShapes, typename TFile>
int
parseShapesFromFile(TShapes & shapeStrings,
                    TFile & file,
                    ParamChooserOptions &)
{
    DirectionIterator<std::fstream, Input>::Type reader(file);
    skipUntil(reader, NotFunctor<IsWhitespace>());

    CharString shape;
    while (!atEnd(reader))
    {
        readShape(shape, reader);
        appendValue(shapeStrings, shape);
        skipUntil(reader, NotFunctor<IsWhitespace>());
    }
    return length(shapeStrings);
}

template<typename TError>
void
interpolateErrorDistr(TError & errorDistr, ParamChooserOptions & pm_options)
{
	if(!pm_options.extrapolate) return;

	unsigned totalN = pm_options.extrapolN;
	unsigned totalLargeN = pm_options.totalN;
	
	// prepare log error distribution 
	TError shorterErrorDistr;
	resize(shorterErrorDistr, totalN);
	float x =(float) (totalLargeN-2)/(totalN-2);
	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	for(unsigned j = 0; j < totalN; ++j) 
	{
		typename Value<TError>::Type newVal;
		float index = j * x;
		float add = index - (int)index;
		if(j<totalN-1)
			newVal = errorDistr[(int)index] + add * (errorDistr[(int)index+1] - errorDistr[(int)index]) ;
		else newVal = errorDistr[totalLargeN-1];
//		std::cout << newVal << std::endl;
		shorterErrorDistr[j] = newVal;
	}
	clear(errorDistr);
	errorDistr = shorterErrorDistr;

}

template<typename TError>
void
makeSelectedStatsFile(TError & errorDistr, ParamChooserOptions & pm_options)
{
//IOREV see comments in other paramChooser.h

	typedef typename Value<TError>::Type TFloat;
	
	unsigned totalN = pm_options.totalN;
	// unsigned totalK = pm_options.totalK;
	if(pm_options.extrapolate == true)
	{
		totalN = pm_options.extrapolN;
		// totalK = pm_options.extrapolK;
		if(totalN != length(errorDistr)/4)
			interpolateErrorDistr(errorDistr,pm_options);
	}

	unsigned maxErrors =  2 + (unsigned) totalN / 10;
	//unsigned maxErrors = 1 + (unsigned) totalN / 10;
	unsigned minErrors = 0;// (unsigned) totalN / 10;
	if(maxErrors<5 && totalN >= 30) maxErrors = 5;
	
	typedef typename Value<TError>::Type TErrorValue;
	String<TErrorValue> logErrorDistribution;
	
	String<CharString> shapeStrings;
	
	if(*pm_options.shapeFile != 0)
	{
		::std::fstream filestrm;
		filestrm.open(pm_options.shapeFile,::std::ios_base::in);
		int result = parseShapesFromFile(shapeStrings,filestrm,pm_options);
		if(result == 0) std::cerr << "0 shapes parsed." << std::endl;
		else if(pm_options.verbose) std::cerr << result <<" shapes parsed." << std::endl;
		filestrm.close();
	}
	
	if(pm_options.useDefaultShapes)
	{
		//q=14
		if(pm_options.optionHammingOnly)
		{
			appendValue(shapeStrings,"1111111111100000111");
			appendValue(shapeStrings,"11101110110001110110001");
			appendValue(shapeStrings,"1111011010001110011011"); 
		}
		else
		{
			appendValue(shapeStrings,"1111111111100111");
			appendValue(shapeStrings,"111111111111101");
		}
		appendValue(shapeStrings,"11111111111111");
	
		//q=13
		if(pm_options.optionHammingOnly)
		{
			appendValue(shapeStrings,"11111111110000111");
			appendValue(shapeStrings,"110101111001100010111");  
		}
		else
		{
			appendValue(shapeStrings,"111111111100111");
			appendValue(shapeStrings,"11111111111101");
		}
		appendValue(shapeStrings,"1111111111111");
	
		//q=12
		if(pm_options.optionHammingOnly)
		{
			appendValue(shapeStrings,"111111111000111");
			appendValue(shapeStrings,"1110100111010011101");
		}
		else
		{
			appendValue(shapeStrings,"11111111100111");
			appendValue(shapeStrings,"1111111111101");
		}
		appendValue(shapeStrings,"111111111111");
	
		//q=11
		if(pm_options.optionHammingOnly)
		{
			appendValue(shapeStrings,"11111110001111");
			appendValue(shapeStrings,"11111101110101");  //median shape
		}
		else
		{
			appendValue(shapeStrings,"1111111100111");
			appendValue(shapeStrings,"111111111101");
		}
		appendValue(shapeStrings,"11111111111");
		
		//q=10
		if(pm_options.optionHammingOnly)
		{
			appendValue(shapeStrings,"1111111000111");
			appendValue(shapeStrings,"111001001010011101");
		}
		else
		{
			appendValue(shapeStrings,"111111100111");
			appendValue(shapeStrings,"11111111101");
		}
		appendValue(shapeStrings,"1111111111");
		
		
		if(totalN < 50)
		{
			//q=9
			appendValue(shapeStrings,"111111111");
			if(pm_options.optionHammingOnly)
			{
				appendValue(shapeStrings,"111111100011");
				appendValue(shapeStrings,"111001001010001011");
			}
			else
			{
				appendValue(shapeStrings,"11111110011");
				appendValue(shapeStrings,"1111111101");
			}
		}
		
		if(totalN < 40)
		{
			//q=8
			appendValue(shapeStrings,"11111111");
			if(pm_options.optionHammingOnly)
			{
				appendValue(shapeStrings,"11111100011");
				appendValue(shapeStrings,"101001111000101");  //median shape
			}
			else
			{
				appendValue(shapeStrings,"1111110011");
				appendValue(shapeStrings,"111111101");
			}
		}
		
		if(totalN < 36)
		{
			//q=7
			appendValue(shapeStrings,"1111111");
			if(pm_options.optionHammingOnly)
			{
				appendValue(shapeStrings,"1111100011");
				appendValue(shapeStrings,"10110000001100101");
			}
			else
			{
				appendValue(shapeStrings,"111110011");
				appendValue(shapeStrings,"11111101");
			}
		}
		
		if(totalN < 32)
		{
			//q=6
			appendValue(shapeStrings,"111111");
			if(pm_options.optionHammingOnly)
			{
				appendValue(shapeStrings,"1111100001");
				appendValue(shapeStrings,"11000000100100101");
			}
			else
			{
				appendValue(shapeStrings,"11111001");
				appendValue(shapeStrings,"1111101");
			}
		}
	}

	unsigned minT = 0;

	String<unsigned> weights;
	resize(weights,length(shapeStrings),0);
	for(unsigned i = 0; i < length(shapeStrings) ; ++i)
		for(unsigned pos = 0; pos < length(shapeStrings[i]) ; ++pos)
			if(shapeStrings[i][pos] == '1')
				++weights[i];
		
	// prepare log error distribution 
	resize(logErrorDistribution, 4*totalN);
	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	double remainingProb = 1.0 - pm_options.optionProbINSERT - pm_options.optionProbDELETE;
	for(unsigned j = 0; j < totalN; ++j) 
	{
		logErrorDistribution[SEQAN_MISMATCH*totalN+j] = _transform(errorDistr[j]);
		logErrorDistribution[SEQAN_INSERT*totalN+j]   = _transform(pm_options.optionProbINSERT);
		logErrorDistribution[SEQAN_DELETE*totalN+j]   = _transform(pm_options.optionProbDELETE);
		logErrorDistribution[SEQAN_MATCH*totalN+j]    = _transform(remainingProb - errorDistr[j]);
	}

#ifdef RUN_RAZERS
	// generate genome and reads
	::std::cout << "Simulate reads..."<<::std::endl;
	TReadSet testReads;
	StringSet<Dna5String> testGenome;
	//StringSet<Dna5String> testReads;
	StringSet<CharString> dummyIDs;
	resize(testGenome, 1);
	std::mt19937 rng(/*seed=*/time(NULL));
	simulateGenome(rng, testGenome[0], 500000);					// generate 1Mbp genomic sequence
	simulateReads(rng,
		testReads, dummyIDs, testGenome, 
		50000, maxErrors+1, logErrorDistribution, 0, 0, 0.5, true);	// generate 50K reads

#endif


	bool first = true;
	
	for(int i = length(shapeStrings)-1; i >= 0; --i)
	{
		if(length(shapeStrings[i])>totalN) continue;
		
		unsigned maxT = totalN-length(shapeStrings[i])+2;
		String<TFloat> sensMat;

		//if(pm_options.verbose)::std::cout << "do DP\n";
//		if(pm_options.verbose)::std::cerr << "do loss rate DP" << std::endl;
		try 
		{
			Shape<Dna, GenericShape> shape;
			stringToShape(shape, shapeStrings[i]);
			if (pm_options.optionHammingOnly)
				qgramFilteringSensitivity(sensMat, shape, totalN, maxErrors - 1, maxT - 1, HammingDistance(), ThreshExact(), logErrorDistribution);
			else
				qgramFilteringSensitivity(sensMat, shape, totalN, maxErrors - 1, maxT - 1, EditDistance(), ThreshExact(), logErrorDistribution);
//			initPatterns(states, shapeStrings[i], maxErrors-1, logErrorDistribution, pm_options.optionHammingOnly);
//			computeFilteringLoss(found, states, length(shapeStrings[i]), maxT, maxErrors,  logErrorDistribution);
		}
		catch (std::bad_alloc&) 
		{
			std::cout << shapeStrings[i] << " threw bad_alloc exception, skipping this shape." << std::endl;
			continue;
		}
	
		for(unsigned e = minErrors; e < maxErrors; ++e) {
			bool highestOptimalFound = false;
			for(unsigned t = maxT-1; t > minT; --t) {
				TFloat lossrate = 1.0 - (TFloat) _transformBack(sensMat[e*maxT+t]);
				if(lossrate <= 0.0){
					if(highestOptimalFound) break;
					else highestOptimalFound = true;
				}
				if(lossrate > 0.2) continue;

//				unsigned gminCov = getMinCov(weights[i], length(shapeStrings[i]), t);

				// create the whole file name
				::std::stringstream datName;
				datName << pm_options.fgparams;
				datName << pm_options.fprefix[0]<<"_N" << totalN << "_";
				if(!pm_options.optionHammingOnly) datName << "L.dat";
				else datName <<"H.dat";
				
				// if datName-file doesnt exist, write the title on
				if(!pm_options.appendToPrevious && first){
					first = false;
					::std::ofstream fout(datName.str().c_str(), ::std::ios::out);
					fout << "errors\tshape\t\tt\t\tlossrate";
					fout << "\tPM";
					fout << ::std::endl << ::std::endl;
					fout.close();
				}
				
#ifdef RUN_RAZERS
				// count verifications
				String<ReadMatch<int> > matches;
				RazerSOptions<RazerSSpec<false, true> > razersOptions;
				razersOptions.errorRate = (double)e / (double)totalN;
				razersOptions.errorRate += 0.0000001;
				razersOptions.threshold = t;
				razersOptions._debugLevel = 2;
				razersOptions.hammingOnly = pm_options.optionHammingOnly;
				int dummy=0;
				assign(razersOptions.shape, shapeStrings[i]);
				mapReads(matches, testGenome, testReads,dummy, dummy, razersOptions);
#endif

				// write shape with its properties into file
				::std::ofstream fout(datName.str().c_str(), ::std::ios::app | ::std::ios::out);
				if(!fout.is_open())
					std::cerr << "Couldn't write to file " << datName.str() << std::endl;
				fout << e << "\t";
				fout << shapeStrings[i] << "\t\t";
				fout << t << "\t\t";
				fout << lossrate;
#ifdef RUN_RAZERS
				fout << "\t\t" << razersOptions.FP + razersOptions.TP;
#else
				fout << "\t\t0";
#endif
				fout << ::std::endl;
				fout.close();
				
			} // t-loop
		}
	}
}


template<typename TSStr>
void
getParamsFilename(TSStr & paramsfile, ParamChooserOptions & pm_options)
{
//IOREV see comments in other paramChooser.h
	int N = pm_options.totalN;
	if(pm_options.extrapolate)
		N = pm_options.extrapolN;
	paramsfile.str("");
	paramsfile << pm_options.fgparams << pm_options.fprefix[0];

	paramsfile << "_N" << N ;
	//if(prefixCount) paramsfile << fgparams<< fprefix[0]<<"_N" << totalN;
	//else paramsfile << fgparams<<"userdef_N" << totalN;
	if(pm_options.optionHammingOnly) paramsfile << "_H";
	else paramsfile << "_L";
	paramsfile << ".dat";
}



//////////////////////////////////////////////////////////////////////////////


template<typename TShape>
inline int
numGaps(TShape & currShape)
{
    int count = 0;
    unsigned j=0;
    bool ingap = false;
    while(j<length(currShape))
    {
        if (currShape[j]=='0')
        {
            if(ingap) ++j;
            else ++count;
            ingap = true;
        }
        else ingap = false;
        ++j;
    }

    return count;

}


//////////////////////////////////////////////////////////////////////////////
// Get parameters q and t optimal for given loss rate
template <typename TFile, typename TSpec>
bool
parseGappedParams(RazerSOptions<TSpec> & r_options, TFile & file, ParamChooserOptions & pm_options)
{
    typedef float TFloat;
    String<CharString> shapes;
    resize(shapes, 14); //best shape for each possible value of q
    String<unsigned> thresholds;
    resize(thresholds, 14); //corresponding t
    String<unsigned> measure;
    resize(measure, 14); //potential matches
    String<TFloat> lossrates;
    resize(lossrates, 14); //lossrates
    double extrapolFactor = 1.0; // no extrapolation
    unsigned errorsWanted = (int)(pm_options.optionErrorRate * pm_options.totalN);
    if (pm_options.extrapolate)
    {
        extrapolFactor = (double)pm_options.totalN / pm_options.extrapolN;
        errorsWanted = pm_options.extrapolK;
    }

    DirectionIterator<std::fstream, Input>::Type reader(file);

    if (atEnd(reader))
    {
        if (pm_options.verbose)
            std::cerr << "Loss rate file is empty!" << std::endl;
        return false;
    }
    if (*reader == 's' || *reader == 'e') //header line
    {
        skipLine(reader);
        skipLine(reader);
    }

    bool atLeastOneFound = false;
    CharString buffer;
    CharString currShape;
    while (!atEnd(reader))
    {
        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        unsigned numErrors = lexicalCast<unsigned>(buffer);
        skipUntil(reader, NotFunctor<IsWhitespace>());
        if (numErrors != errorsWanted)
        {
            skipLine(reader);
            continue;
        }

        readShape(currShape, reader);
        if ((pm_options.chooseUngappedOnly && numGaps(currShape) > 0) || (pm_options.chooseOneGappedOnly && numGaps(currShape) > 1))
        {
            skipLine(reader);
            continue;
        }
        skipUntil(reader, NotFunctor<IsWhitespace>());

        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        unsigned currThreshold = (unsigned)(lexicalCast<unsigned>(buffer) * extrapolFactor); //when extrapolating from shorter read lengths, threshold can be at least linearly increased
        skipUntil(reader, NotFunctor<IsWhitespace>());

        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        TFloat currLossrate = lexicalCast<double>(buffer);
        skipUntil(reader, NotFunctor<IsWhitespace>());

        clear(buffer);
        readUntil(buffer, reader, IsWhitespace());
        unsigned currMeasure = lexicalCast<unsigned>(buffer); // potential matches measured on simulated reads

        //std::cout << numErrors << "\t" << currShape << "\t" << currThreshold << "\t" << currLossrate << "\t" << currMeasure << std::endl;
        if (currThreshold >= pm_options.minThreshold && currLossrate <= pm_options.optionLossRate /*&& val > bestSoFar*/)
        {

            unsigned weight = 0;
            for (unsigned pos = 0; pos < length(currShape); ++pos)
                if (currShape[pos] == '1')
                    ++weight;
            if (length(shapes[weight - 1]) > 0)  // if this is not the first shape with weight weight
            {
                // compare currShape to the best one found so far
                if (currMeasure <= measure[weight - 1])
                {
                    if (currMeasure == measure[weight - 1])
                    {
                        bool undecided = false;
                        //next measure: threshold
                        if (thresholds[weight - 1] > currThreshold)
                        {
                            skipLine(reader);
                            continue;
                        }
                        else if (thresholds[weight - 1] == currThreshold)
                            undecided = true;

                        //if still undecided: next measure: span
                        if (undecided && length(shapes[weight - 1]) > length(currShape))
                        {
                            skipLine(reader);
                            continue;
                        }
                        else if (undecided && length(shapes[weight - 1]) < length(currShape))
                            undecided = false;

                        //if still undecided: next measure: lossrate
                        if (undecided && lossrates[weight - 1] < currLossrate)
                        {
                            skipLine(reader);
                            continue;
                        }
                    }
                    shapes[weight - 1] = currShape;
                    measure[weight - 1] = currMeasure;
                    thresholds[weight - 1] = currThreshold;
                    lossrates[weight - 1] = currLossrate;
                    atLeastOneFound = true;
                }

            }
            else
            {
                shapes[weight - 1] = currShape;
                measure[weight - 1] = currMeasure;
                thresholds[weight - 1] = currThreshold;
                lossrates[weight - 1] = currLossrate;
                atLeastOneFound = true;

            }
        }
        skipLine(reader);

    }
    if (!atLeastOneFound)
    {
        if (pm_options.verbose)
            std::cerr << std::endl << "!!! Something wrong with file? !!!" << std::endl;
        return false;
    }
    int i;
    for (i = pm_options.maxWeight - 1; i > 0; --i)
        if (length(shapes[i]) > 0) // if a shape of weight i+1 has been found
            break;
    if (i == 0)
    {
        if (pm_options.verbose)
            std::cerr << std::endl << "!!! Something wrong with file? !!!" << std::endl;
        return false;
    }
    if (thresholds[i] == 1 && length(shapes[i - 1]) > 0 && thresholds[i - 1] > 2)
        --i;
    pm_options.chosenLossRate = lossrates[i];
    assign(r_options.shape, shapes[i]);
    r_options.threshold = thresholds[i];
    // suggest a suitable combination of q and t

    return true;
}

// extrapolate if n is large
template <typename TOptions>
void
extrapolateNK(TOptions & pm_options)
{
	
	pm_options.extrapolate = true;

	double recErrorRatio = (double)pm_options.totalN/pm_options.totalK;
	int bestN, currN;
	if(pm_options.optionHammingOnly) bestN = pm_options.maxComputedHammingN;
	else bestN = pm_options.maxComputedEditN;
	currN = bestN;

	for(; currN > 49; --currN)
	{
		double bestRecErrorRatio = (double)bestN/(ceil((double)bestN * 1.0/recErrorRatio));
		double currRecErrorRatio = (double)currN/(ceil((double)currN * 1.0/recErrorRatio));
		if(currRecErrorRatio > bestRecErrorRatio)
			bestN = currN;
	}
	pm_options.extrapolN = bestN;
	pm_options.extrapolK = (unsigned)ceil(((double)bestN * 1.0/recErrorRatio)-0.00001);
	

}


template<typename TSpec>
bool
chooseParams(RazerSOptions<TSpec> & r_options, ParamChooserOptions & pm_options)
{
	typedef float TFloat;
	static const TFloat epsilon = (TFloat)0.00000000001;	
	pm_options.optionLossRate += epsilon;


#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(r_options.maqMapping)
	{
		//artSeedLength denotes the length of the 3' end of the read for which 2-error matches are found with 100% sensitvity
		if(r_options.artSeedLength==28) r_options.shape = "1101001110100111010011";//q13
		if(r_options.artSeedLength==29) r_options.shape = "11101001110100111010011";//q14
		if(r_options.artSeedLength==30) r_options.shape = "111010011101001110100101";//q14
		if(r_options.artSeedLength==31) r_options.shape = "111010011101001110100101";//q14
		if(r_options.lowMemory && r_options.artSeedLength==29) r_options.shape = "10101001110100111010011";
		r_options.threshold = 1; 
		if(r_options.lowMemory && r_options.artSeedLength==30)		
		{
			r_options.shape = "111110101000011111";
			r_options.threshold = 2;		
		}
		if(r_options.lowMemory && r_options.artSeedLength==31)		
		{
			r_options.shape = "1101000110100011010001101";
			r_options.threshold = 2;		
		}  
		if(r_options.artSeedLength==32)		//perfect setting for queueing system
		{
			r_options.shape = "11110010001111001000111";
			r_options.threshold = 2;		
		}
		
		if(r_options.artSeedLength<28 || r_options.artSeedLength> 32)
			::std::cerr << "Warning: This Version of RazerS nly supports quality-seed-lengths of 28 to 32. Using default filter parameters." << ::std::endl;
		else return true;
	}
#endif

	pm_options.fgparams = pm_options.paramFolderPath;
	if( length(pm_options.paramFolder) > 0)
		append(pm_options.fgparams, pm_options.paramFolder);
	else
		append(pm_options.fgparams, "gapped_params/");

	
	if(pm_options.optionProbINSERT <= epsilon && pm_options.optionProbDELETE <= epsilon)
		pm_options.optionHammingOnly=true;


	pm_options.totalK = (int)(pm_options.optionErrorRate * pm_options.totalN);
	
	if ((pm_options.optionHammingOnly && pm_options.totalN > pm_options.maxComputedHammingN )
		|| (!pm_options.optionHammingOnly && pm_options.totalN > pm_options.maxComputedEditN ))
		extrapolateNK(pm_options);
	
    // compute data specific loss rates
	if (pm_options.fnameCount0 || pm_options.fnameCount1)
	{
		if(!pm_options.prefixCount)
		{
			pm_options.fprefix[0] = "results";
//			pm_options.fprefix[0] = "userdef";
//			::std::cerr << "\nNo session id given, using prefix 'userdef'"<<::std::endl;
		}
		String<TFloat> errorDistribution;
		resize(errorDistribution,pm_options.totalN);
		//error distribution given --> read file containing error distr and compute loss rates
		if(pm_options.fnameCount1)
		{
			::std::fstream file;
			file.open(pm_options.fname[1],::std::ios_base::in | ::std::ios_base::binary);
			if(!file.is_open())
			{
				::std::cerr << "Couldn't open file "<<pm_options.fname[1]<<std::endl;
				return false;
			}
			unsigned count = 0;
            DirectionIterator<std::fstream, Input>::Type reader(file);
            CharString buffer;
            while (!atEnd(reader) && count < pm_options.totalN)
            {
                clear(buffer);
                skipUntil(reader, NotFunctor<IsWhitespace>());
                readUntil(buffer, reader, IsWhitespace());
                lexicalCastWithException(errorDistribution[count], buffer); // + (TFloat) 1.0/maxN;
                ++count;
            }
            file.close();
			if(count != pm_options.totalN + 1)
			{
				::std::cerr << "Error distribution file must contain at least " << pm_options.totalN << " probability values (one value per line)." << std::endl;
				return false;
			}
		}
		else // read qualtiy files and compute position dependent avg error probabilites
		{
			getAvgFromPrbDirectory(pm_options.fname[0],errorDistribution,pm_options);
		}

		::std::fstream file;
		//if(prefixCount)
		makeSelectedStatsFile(errorDistribution,pm_options);
	}
	else if(!pm_options.prefixCount) pm_options.fprefix[0] = "results";


	// get name of loss rate file
	::std::stringstream paramsfile;
	getParamsFilename(paramsfile,pm_options);
        if (pm_options.verbose)
        {
               ::std::cerr << ::std::endl;
               ::std::cerr << "Read length      = " << pm_options.totalN << "bp" << std::endl;
               ::std::cerr << "Max num errors   = " << pm_options.totalK << std::endl;
               ::std::cerr << "Recognition rate = " << 100.0*(1.0-pm_options.optionLossRate) << "%" << std::endl;
		if(pm_options.extrapolate) ::std::cerr << "Extrapolating from read length " << pm_options.extrapolN << " and " << pm_options.extrapolK << "errors." << std::endl;
        }
	
	// parse loss rate file and find appropriate filter criterium
	if(pm_options.verbose) ::std::cerr << std::endl << "--> Reading " <<  paramsfile.str()<<::std::endl;
	::std::fstream file;
	file.open(paramsfile.str().c_str(),::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
	{
		if(pm_options.verbose)::std::cerr << "Couldn't open file "<<paramsfile.str()<<::std::endl;
		return false;
	}
	else
	{
		parseGappedParams(r_options,file,pm_options);
		if(pm_options.verbose) ::std::cout << std::endl << " Choose "<< std::endl << "shape: " << r_options.shape << std::endl << " and " << std::endl << "threshold: " << r_options.threshold << std::endl <<" to achieve optimal performance for expected recognition rate >= " << (100.0-100.0*pm_options.optionLossRate) << "% (expected recognition = " << (100.0-pm_options.chosenLossRate*100.0) <<"%)" <<std::endl << std::endl;
		file.close();
	}

	return true;
}
}

#endif
