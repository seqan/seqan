#ifndef SEQAN_HEADER_INDELSIM_H
#define SEQAN_HEADER_INDELSIM_H


#include "seqan/platform.h"
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_parsing.h>
#include <iostream>
#include <fstream>
#include <cmath>

#define EPS_PREC 0.00000001
using namespace std;

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default options


//____________________________________________________________________________
// Global Parameters



struct IndelCheck{};

template<typename TSpec = IndelCheck>
struct IndelSimOptions
{
	//const char *output;	// one file for manipulated reference
	//const char *outputInfo;	// and one for indel information
	CharString output;	// one file for manipulated reference
	CharString outputIndel;	// and one for indel information
	CharString outputSnp;	// and one for indel information
	CharString inputIndel;	// and one for indel information

	unsigned numIndels;
	String<Pair<int,int> > ranges;	// indel sizes that are supposed to be simulated
					// probability among and within ranges is uniform
	String<double> rangeProbs;	// corresponding probabilities
	
	unsigned numSnps;
	int minDistance;		// minimal distance between two implanted indels
	double duplicationProb;
	int noNsInRange;		// window around randomly generated position +- noNsInRange 
					// may not contain Ns in original sequence

    bool diploid;    
	int _debugLevel;
	
	IndelSimOptions()
	{
		output = "";		
		outputIndel = "";	
		outputSnp = "";	
		inputIndel = "";	
		
		numIndels = 100;
//		appendValue(ranges,Pair<int,int>(-30,30));	//simulated uniformly in interval [-30,30]
//		appendValue(rangeProbs,1.0);	//simulated uniformly in interval [-30,30]
		
		minDistance  = 2*30;
		duplicationProb = 0.5;
		noNsInRange = 10;
		
        diploid = false;
		_debugLevel = 0;
		
	}
	
};


struct IndelInfo{
	
	unsigned genomeId;
	unsigned originalPos;
	unsigned simPos;
	int indelSize;
	bool duplication;
    bool haplotypeA;
    bool haplotypeB;
    bool heterozygote;
	CharString idStr;
	
};

struct SnpInfo{
	
	unsigned genomeId;
	unsigned originalPos;
	unsigned simPos;
	Dna base;
    bool haplotypeA;
    bool haplotypeB;
    bool heterozygote;
	CharString idStr;
	
};

//____________________________________________________________________________
// helper parsing functions

template<typename TChar>
inline bool
parse_isDigit(TChar const c)
{
//IOREV _duplicate_ use ctype.h's isdigit() instead
	//return (((int ) c >=  48) && ((int) c <=  57));
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
		    (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}


template<typename TFile, typename TChar>
inline long double
parse_readDouble(TFile & file, TChar& c)
{
//IOREV _duplicate_ _hasCRef_ see _parseReadDouble in misc_parsing.h
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!parse_isDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

template<typename TFile, typename TChar>
inline void 
parse_skipWhitespace(TFile& file, TChar& c)
{
//IOREV _hasCRef_ _duplicate_ see _parseSkipWhitespace() in misc_parsing.h
	if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) break;
	}
}


template<typename TFile, typename TChar, typename TString>
void
_parseReadWordUntilWhitespace(TFile& file, TString& str, TChar& c)
{
//IOREV _hasCRef_ _duplicate_ see equally named function in misc_parsing.h
        append(str,c);
        if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
                c = _streamGet(file);
                return;
        }
        while (!_streamEOF(file)) {
                c = _streamGet(file);
                if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
                append(str, c);
        }
        return;
}   




template <typename TOptions>
bool loadGenomes(const char* fileName, 
		StringSet<Dna5String> &genomes,
		StringSet<CharString> & genomeIDs, 
		::std::map<CharString,unsigned> &gIdStringToIdNumMap,
		TOptions &)
{
	
	MultiFasta multiFasta;
	if (!open(multiFasta.concat,fileName,OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
	resize(genomeIDs,seqCount);
	resize(genomes,seqCount);
	
	for(unsigned i = 0; i < seqCount; ++i)
	{
		CharString temp;
		assignSeq(genomes[i], multiFasta[i], Fasta());
		assignSeqId(temp, multiFasta[i], Fasta());
		for (unsigned pos = 0; pos < length(temp); ++pos)
		{
			if(temp[pos]=='\t' || temp[pos]=='\b' || temp[pos]==' ')
			{
				resize(temp,pos);
				break;
			}
		}
		genomeIDs[i] = temp;
		gIdStringToIdNumMap.insert(std::make_pair(temp, i)); 
	}
	return (seqCount > 0);
}



//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TStringSet, typename TNameSet>
bool loadFasta(TStringSet &genomes, TNameSet &fastaIDs, char const *fileName)
{
	// count sequences
	unsigned seqCount = 0;

	ifstream file;
	file.open(fileName, ios_base::in | ios_base::binary);
	if (!file.is_open()) return false;
	while (!_streamEOF(file)) {
		goNext(file, Fasta());
		++seqCount;
	}

	// import sequences
	file.clear();
	file.seekg(0, ios_base::beg);
	resize(fastaIDs, seqCount);
	resize(genomes, seqCount);
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 
	{
		readShortID(file, fastaIDs[i], Fasta());		// read Fasta id
		read(file, genomes[i], Fasta());		// read sequence
	}
	file.close();
	return (seqCount > 0);
}



template < typename TSimGenomes, typename TSimIDs, typename TOptions>
void saveFasta(
	TSimGenomes const &simGenomes,		// generated read sequences
	TSimIDs const &simIDs,			// corresponding Fasta ids
	TOptions &options,
	string const &genomeFName)
{
	ostringstream fileName;
	//if (*options.output != 0)
	//	fileName << options.output;
	if (*toCString(options.output) != 0)
		fileName << options.output;
	else
		fileName << genomeFName << ".indeled";

	ofstream file;
	
	file.open(fileName.str().c_str(), ios_base::out | ios_base::trunc);
	if (!file.is_open()) {
		cerr << "\nFailed to open output file" << endl;
		return;
	}
	else
		if(options._debugLevel > 0)cout << "\nWriting simulated genome to " << fileName.str() << "\n";

	unsigned numGenomes = length(simGenomes);
    if(options._debugLevel > 1) std::cout << "Simulated " << numGenomes << "genomes\n";
	for(unsigned i = 0; i < numGenomes; ++i)
	{
		file << '>' << simIDs[i] << std::endl;
		file << simGenomes[i] << std::endl;
	}
	file.close();
}


template <typename TOptions, typename TGenomeIDs>
void saveSnpInfos(
	std::map<int,SnpInfo> &snpMap,
	TGenomeIDs &genomeIDs,
	TOptions &options,
	string const &genomeFName)
{
	ostringstream fileName;
//	if (*options.outputInfo != 0)
//		fileName << options.outputInfo;
	if (*toCString(options.outputSnp) != 0)
		fileName << options.outputSnp;
	else
	{
		//if (*options.output != 0)
		//	fileName << options.output << ".info";
		if (*toCString(options.output) != 0)
			fileName << options.output << ".info";
		else
			fileName << genomeFName << ".snped.info";
	}
	ofstream file;
	
	file.open(fileName.str().c_str(), ios_base::out | ios_base::trunc);
	if (!file.is_open()) {
		cerr << "\nFailed to open output file" << endl;
		return;
	}
	else
		if(options._debugLevel) cout << "\nWriting SNP infos to " << fileName.str() << "\n";
	
	std::map<int,SnpInfo>::iterator it = snpMap.begin();
	int count = 0;
	while(it != snpMap.end())
	{
		file << genomeIDs[(it->second).genomeId] << '\t';
		file << "simulated\tsnv\t";
		file << (it->second).originalPos + 1 << '\t';
		file << (it->second).originalPos + 1 << '\t';
		file << ".\t+\t.\t";
		file << "ID=" << count<<";simPos=" << (it->second).simPos + 1;
		if((it->second).heterozygote) file << ";geno=het";
        else file << ";geno=homo";
		file << std::endl;
		++count;
		++it;
	}
	file.close();
}



template <typename TOptions, typename TGenomeIDs>
void saveIndelInfos(
	std::map<int,IndelInfo> &indelMap,
	TGenomeIDs &genomeIDs,
	TOptions &options,
	string const &genomeFName)
{
	ostringstream fileName;
//	if (*options.outputInfo != 0)
//		fileName << options.outputInfo;
	if (*toCString(options.outputIndel) != 0)
		fileName << options.outputIndel;
	else
	{
		//if (*options.output != 0)
		//	fileName << options.output << ".info";
		if (*toCString(options.output) != 0)
			fileName << options.output << ".info";
		else
			fileName << genomeFName << ".indeled.info";
	}
	ofstream file;
	
	file.open(fileName.str().c_str(), ios_base::out | ios_base::trunc);
	if (!file.is_open()) {
		cerr << "\nFailed to open output file" << endl;
		return;
	}
	else
		if(options._debugLevel) cout << "\nWriting indel infos to " << fileName.str() << "\n";
	
	std::map<int,IndelInfo>::iterator it = indelMap.begin();
	int count = 0;
	while(it != indelMap.end())
	{
		file << genomeIDs[(it->second).genomeId] << '\t';
		if((it->second).indelSize > 0)
			file << "simulated\tdeletion\t";
		else
			file << "simulated\tinsertion\t";
		file << (it->second).originalPos + 1 << '\t';
		if((it->second).indelSize > 0)
			file << (it->second).originalPos + (it->second).indelSize << '\t';
		else
			file << (it->second).originalPos + 1 << '\t';
		file << ".\t+\t.\t";
		file << "ID=" << count<<";size=" <<(it->second).indelSize << ";simPos=" << (it->second).simPos + 1;
		if((it->second).duplication) file << ";duplication=1";
		if((it->second).heterozygote) file << ";geno=het";
        else file << ";geno=homo";
		file << std::endl;
		++count;
		++it;
	}
	file.close();
}



/////////////////////////////////////////////////////////////
// read Gff input file containing indels
template <
	typename TCString,
	typename TIndelSet,
	typename TGenomes,
	typename TGenomeMap,
	typename TOptions
>
int readGFF(
	TCString				filename,
	TIndelSet 				&indelSet,
	TGenomes				&genomes,
	TGenomeMap				&gIdStringToIdNumMap,
	TOptions				&options)
{
//IOREV _nodoc_ what is this doing here? isn't there GFF file format support in store/store_io_gff.h ?
	typedef typename Value<TIndelSet>::Type	TIndel;
	typedef int				TId;
	typedef int				TContigPos;
	
	
	::std::ifstream file;
	file.open(filename, ::std::ios_base::in | ::std::ios_base::binary);
	if (!file.is_open()) return 1;
		
	TIndel indel = {0,0,0,0,0,0,0,0,0};
	int absMaxValue = 1;
	clear(indelSet);
	char c = _streamGet(file);
	while (!_streamEOF(file))
	{
		
		if(c == '#')
			_parseSkipLine(file,c);	
	
		// skip whitespaces just in case (actually there shouldnt be a whitespace at the beginning of a line)
		_parseSkipWhitespace(file, c);
	
		if(c == '#')
			_parseSkipLine(file,c);	
		// and read entry in column 1  --> genomeID
		CharString temp_str;
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(prefix(temp_str,3)=="chr")
			temp_str = suffix(temp_str,3);
		
		TId contigId;
		//check if the genomeID is in our map of relevant genomeIDs, otherwise skip match
		typename TGenomeMap::iterator it = gIdStringToIdNumMap.find(temp_str);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if(it != gIdStringToIdNumMap.end()) contigId = it->second;
		else
		{
			_parseSkipLine(file,c);
			continue;
		}
		
		// skip whitespaces and read entry in column 2
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		
		// skip whitespaces and read entry in column 3
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		
		// skip whitespaces and read entry in column 4  --> genomic begin position
		_parseSkipWhitespace(file, c);
		indel.originalPos = (TContigPos) _parseReadNumber(file,c) - 1;
		if(options._debugLevel > 1) 
			::std::cout << indel.originalPos << "\t";
		
		// skip whitespaces and read entry in column 5  --> genomic end position // not needed here
		_parseSkipWhitespace(file, c);
		indel.indelSize = _parseReadNumber(file,c) - indel.originalPos;
		indel.duplication = 0;
		if(indel.originalPos + indel.indelSize > length(genomes[contigId]) )
		{
			_parseSkipLine(file,c);
			continue;
		}

		// skip whitespaces and read entry in column 6  --> score (percent identity or mapping quality) or a '.'
		int readSupport = 1000; //  --> no information about read support (reference indel)
		_parseSkipWhitespace(file, c);
		if(c=='.')
			c = _streamGet(file);               // 
		else 
			readSupport = (TContigPos) _parseReadDouble(file,c); // number of supporting reads
			
		if(options._debugLevel > 1) 
			::std::cout << readSupport << "\t";
		
		// skip whitespaces and read entry in column 7  --> strand information: '+' or '-' // not needed here
		_parseSkipWhitespace(file, c);
		c = _streamGet(file);
		
		// skip whitespaces and read entry in column 8  --> always '.' here
		_parseSkipWhitespace(file, c);
		c = _streamGet(file);
		
		// skip whitespaces and read entry in column 9  --> tags, extra information. first tag is always "ID"
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadIdentifier(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\n";
		if(temp_str!="ID") ::std::cout << "first feature field should be 'ID'"<<::std::endl;
		
		// skip the "="
		c = _streamGet(file);
		
		// read the ID
		clear(temp_str);
		clear(indel.idStr);
		CharString indelID;
		_parseReadIdentifier(file,indel.idStr,c);
		if(options._debugLevel > 1) 
			::std::cout << "myID = "<< indel.idStr << "\n";
		
		// process tags in a loop
		CharString current_tag;
		_parseSkipWhitespace(file,c); 
		while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n'))) // while in same line
		{
			// different tags are separated by ';'  
			while(c != ';')
			{
				if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
					break;
				c = _streamGet(file);
			}
			if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
				break;
			
			// get the current tag
			clear(current_tag);
			c = _streamGet(file);
			_parseReadIdentifier(file,current_tag,c);
			if(options._debugLevel > 1) 
				::std::cout << current_tag << " in features\n";
			if(current_tag=="size")
			{
		//		if(c == '-') 
					c = _streamGet(file);
				indel.indelSize = _parseReadNumber(file,c); 
				if(options._debugLevel > 1) 
					::std::cout << indel.indelSize << " indel size\n";
			}
			else
			{
				if(current_tag=="duplication") indel.duplication = true;
		//		else _parseSkipLine(file,c);
			}
		}
		appendValue(indelSet,indel);
		if(abs(indel.indelSize) > absMaxValue) absMaxValue = abs(indel.indelSize);
		_parseSkipWhitespace(file, c);
		
	}

	file.close();
	options.minDistance = 2 * absMaxValue;
	if(!empty(options.ranges))
		options.minDistance = 2 * options.ranges[length(options.ranges)-1].i2;

	if(options._debugLevel > 0) std::cout << "MinDistance = " << options.minDistance << std::endl;
	
	if(options._debugLevel > 0) ::std::cout << "Parsed "<<length(indelSet)<<" indels." << ::std::endl;
	
	return 0;
}




//____________________________________________________________________________




template<typename TOptions, typename TGenomes>
bool
simulateSnps( ::std::map<int,SnpInfo>	&snpMap,
						 TGenomes &genomes,
						 String<unsigned> & limits,
						 TOptions & options)
{
	typedef ::std::map<int,SnpInfo> TSnpMap;
	typedef typename TSnpMap::iterator TSnpMapIt;

	unsigned sumLen = limits[length(limits)-1];

	SnpInfo info;
    unsigned count = 0, trials = 0;
	while(count < options.numSnps)
	{
		++trials;
		
        // get random position
		int randPos = (int) (rand() % sumLen);
//		std::cout << "randPos = " << randPos << std::endl;
	
	
		// is position already in map?
		TSnpMapIt it = snpMap.find(randPos);
        if(it != snpMap.end())
            continue;
	
        // which chromosome are we on?
		int limit_counter = 0;	
		while(randPos >= (int)limits[limit_counter])
			++limit_counter;
		int chr = limit_counter - 1;
		
		// is position in a stretch of or too close to Ns?
		unsigned relPos = randPos-limits[chr];
		unsigned pos = randPos-limits[chr]-options.noNsInRange;
		for (; pos < randPos-limits[chr]+options.noNsInRange; ++pos)
			if(genomes[chr][pos] == 'N') break;
		if(pos != randPos-limits[chr]+options.noNsInRange) continue;
		
		// now generate base 
        Dna refBase = genomes[chr][relPos];
        (void)refBase;
		Dna randBase = (Dna) (rand() % 4);
        while(randBase == genomes[chr][relPos])
            randBase = (Dna) (rand() % 4);
		info.base = randBase;

		// store simulation info
		info.genomeId = chr;
		info.originalPos = relPos;
		info.simPos = 0; // not yet determinable, we will never know anyway..
		
        // determine genotype if diploid simulation
        if(options.diploid)
        {
            double randomNum = ((double)rand()/(double)RAND_MAX);
            
            // haplotype 1
            if (randomNum <= 0.666666666)
                info.haplotypeA = true;
            else info.haplotypeA = false;

            // haplotype 2
            if (randomNum >= 0.333333333)
                info.haplotypeB = true;
            else info.haplotypeB = false;

            // remember that it is a het
            if (randomNum <= 0.666666666 && randomNum >= 0.333333333)
                info.heterozygote = false;
            else info.heterozygote = true;

        }
 
		snpMap.insert(std::make_pair(randPos, info));
		++count;
		if(options._debugLevel > 1 && count%100==0) std::cout << "." << std::flush;
		if(options._debugLevel > 1 && count%1000==0) std::cout << count<< std::flush;
		
	}
	return true;
}






template<typename TOptions, typename TGenomes>
bool
simulateIndelsFromRanges( ::std::map<int,IndelInfo>	&indelMap,
						 TGenomes &genomes,
						 String<unsigned> & limits,
						 TOptions & options)
{
	typedef ::std::map<int,IndelInfo> TIndelMap;
	typedef typename TIndelMap::iterator TIndelMapIt;

	unsigned sumLen = limits[length(limits)-1];

	// make sure range probs are given 
	if(length(options.ranges) != length(options.rangeProbs))
	{
		if(options._debugLevel > 0)std::cout << "No range probabilities given. Using same probability for each bucket.\n";
		clear(options.rangeProbs);
		resize(options.rangeProbs,length(options.ranges),0.2);
	}
	
	// make sure range probs add up to 1
	double sumProbs = 0.0;
	for(unsigned i = 0; i < length(options.rangeProbs); ++i)
		sumProbs += options.rangeProbs[i];
	if(sumProbs < EPS_PREC)
	{
		std::cout << "Something wrong with range probabilities.\nExiting...\n";
		return 1;
	}
	for(unsigned i = 0; i < length(options.rangeProbs); ++i)
		options.rangeProbs[i] /= sumProbs;
	
	// cumulative probabilities
	String<double> rangeProbLimit;
	resize(rangeProbLimit,length(options.ranges));
	rangeProbLimit[0] = options.rangeProbs[0];
	for(unsigned i = 1; i < length(options.ranges); ++i)
		rangeProbLimit[i] = options.rangeProbs[i] + rangeProbLimit[i-1];
	
	if(options._debugLevel > 1)
	{
		::std::cout << "RangeProbs:" << std::endl;
		for(unsigned i = 0; i < length(options.rangeProbs); ++i)
			::std::cout << options.rangeProbs[i] << '\t';
		::std::cout << std::endl;
		::std::cout << "RangeProbLimits:" << std::endl;
		for(unsigned i = 0; i < length(rangeProbLimit); ++i)
			::std::cout << rangeProbLimit[i] << '\t';
		::std::cout << std::endl;
	}


	// generate indels to implant
	if(options._debugLevel > 0) std::cout << "Start generating positions." << std::endl;


	unsigned count = 0;
	unsigned trials = 0;
	IndelInfo info;
	while(count < options.numIndels)
	{
		++trials;
		if(trials > options.numIndels * options.numIndels * options.numIndels) return false;
		if(options._debugLevel > 2)
		{
			std::cout << "currently in map: ";
			for (TIndelMapIt itt = indelMap.begin(); itt != indelMap.end(); ++itt)
				::std::cout << itt->first << "\t";
			std::cout << std::endl;
		}
		// get random position
		int randPos = (int) (rand() % sumLen);
//		std::cout << "randPos = " << randPos << std::endl;
	
		// which chromosome are we on?
		int limit_counter = 0;	
		while(randPos >= (int)limits[limit_counter])
			++limit_counter;
		int chr = limit_counter - 1;
		
		// is position too close to chromosome borders?
		if(randPos > (int) limits[limit_counter] - options.minDistance 
			|| randPos < (int) limits[limit_counter - 1] + options.minDistance)
			continue;
		
	
		// is position too close to previously generated indel position?
		TIndelMapIt it = indelMap.upper_bound(randPos);
		
		if(options._debugLevel > 1)std::cout << "compare with upper it->first=" << it->first<< ::std::endl;
		//if(it != indelMap.end()) std::cout<< it->first << " <- upperbound " << std::endl;
		if(it != indelMap.end() && it->first - randPos < options.minDistance)
			continue;
		
		if(it != indelMap.begin())
		{
			--it;
			if(options._debugLevel > 1)std::cout << "compare with lower it->first=" << it->first<< ::std::endl;
//			TIndelMapIt it2 = indelMap.lower_bound(randPos);
//			//if(it2 != indelMap.end() && randPos > it2->first)  --it;
//			if(it2 != indelMap.end()) std::cout<< it2->first << " <- lowerbound " << std::endl;
//			if(!(indelMap.size() == 0) && it != indelMap.end()) std::cout<< it->first << " <- lowerbound " << std::endl;
			if(!(indelMap.size() == 0) && randPos - it->first < options.minDistance)
				continue;
		}
		
		// is position in a stretch of or too close to Ns?
		unsigned pos = randPos-limits[chr]-options.noNsInRange;
		for (; pos < randPos-limits[chr]+options.noNsInRange; ++pos)
			if(genomes[chr][pos] == 'N') break;
		if(pos != randPos-limits[chr]+options.noNsInRange) continue;
		
		// now generate indel size
		int indelSize = 0;
		
		// which range?
		double prob = (double)rand()/RAND_MAX;
		int rangeBucket = 0;
		while(prob >= rangeProbLimit[rangeBucket])
			++rangeBucket;
		
		// which size?
		int rangeLen = (options.ranges[rangeBucket]).i2 - (options.ranges[rangeBucket]).i1;
		while(indelSize == 0)  	// 0 -> either not yet simulated or actually indelsize 0  was generated (which is not allowed)
		{
			int randSize = (int) (rand() % rangeLen);
			indelSize = (options.ranges[rangeBucket]).i1 + randSize;
		}
		
		// safety check
		if(indelSize >= (options.ranges[rangeBucket]).i2  || indelSize < (options.ranges[rangeBucket]).i1)
		{
			std::cerr << "Simulated indelSize outside bucket!!! H�?\n";
			continue;
		}
		
		bool duplication = false;
		if(indelSize < 0) // insertion
		{
			double prob = (double)rand()/RAND_MAX;
			if(prob < options.duplicationProb)
				duplication = true;
			
		}
		else // deletion
		{
			// is end position in a stretch of or too close to Ns?
			unsigned endpos=randPos-limits[chr]+indelSize-options.noNsInRange;
			for (; endpos < randPos-limits[chr]+indelSize+options.noNsInRange; ++endpos)
				if(genomes[chr][endpos] == 'N') break;
			if(endpos != randPos-limits[chr]+indelSize+options.noNsInRange) continue;
		}
		
		
		// store simulation info
		info.genomeId = chr;
		info.originalPos = randPos - limits[chr];
		info.simPos = 0; // not yet determinable
		info.indelSize = indelSize;
		info.duplication = duplication; 
		
        // determine genotype if diploid simulation
        if(options.diploid)
        {
            double randomNum = ((double)rand()/(double)RAND_MAX);
            // haplotype 1
            if (randomNum <= 0.666666666)
                info.haplotypeA = true;
            else info.haplotypeA = false;

            // haplotype 2
            if (randomNum >= 0.333333333)
                info.haplotypeB = true;
            else info.haplotypeB = false;

            // remember that it is a het
            if (randomNum <= 0.666666666 && randomNum >= 0.333333333)
                info.heterozygote = false;
            else info.heterozygote = true;
        }
 
		indelMap.insert(std::make_pair(randPos, info));
		++count;
		if(options._debugLevel > 1 && count%100==0) std::cout << "." << std::flush;
		if(options._debugLevel > 1 && count%1000==0) std::cout << count<< std::flush;
		
	}
	return true;
}


template<typename TOptions, typename TGenomeMap, typename TGenomes>
bool
simulateIndelsFromGff( ::std::map<int,IndelInfo> &indelMap, 
					  String<unsigned> & limits,
					  TGenomes &genomes,
					  TGenomeMap &gIdStringToIdNumMap,
					  TOptions &options)
{
	typedef ::std::map<int,IndelInfo> TIndelMap;
	typedef typename TIndelMap::iterator TIndelMapIt;

	// draw indels from given gff file
	if(options._debugLevel > 0) std::cout << "Start generating positions." << std::endl;

	StringSet<IndelInfo> indelSet;
	if (readGFF(toCString(options.inputIndel), indelSet, genomes, gIdStringToIdNumMap, options) > 0) 
	{
		std::cerr << "Template indels " << options.inputIndel << " can't be loaded." << std::endl;
		return 0;
	}
	if(options._debugLevel > 0) std::cout << "Number of indels to sample from: " << length(indelSet) << std::endl;
	
	// if also ranges are given, shrink the indel set to only the desired indel sizes 
	if(!empty(options.ranges))
		reduceIndelSetAccordingToProbs(indelSet,options);

	int sumLen = length(indelSet);
	unsigned count = 0;
	unsigned trials = 0;
	IndelInfo info;
	while(count < options.numIndels)
	{
		++trials;
		if(trials > options.numIndels * options.numIndels * options.numIndels) return false;
		if(options._debugLevel > 1) std::cout << "Trial" << trials << "\t" << std::flush;
		if(options._debugLevel > 2)
		{
			std::cout << "currently in map: ";
			for (TIndelMapIt itt = indelMap.begin(); itt != indelMap.end(); ++itt)
				::std::cout << itt->first << "\t";
			std::cout << std::endl;
		}
	
		// get random position in indelSet
		int indelsetPos = (int) (rand() % sumLen);
		//		std::cout << "randPos = " << randPos << std::endl;


		unsigned chr = indelSet[indelsetPos].genomeId;
		unsigned randPos = limits[chr] + indelSet[indelsetPos].originalPos;
				
		// is position too close to previously generated indel position?
		TIndelMapIt it = indelMap.find(randPos);
		if(it != indelMap.end()) continue;

		it = indelMap.upper_bound(randPos);
		if(options._debugLevel > 1)std::cout << "randPos=" << randPos << "\t";
		if(options._debugLevel > 1)std::cout << "compare with upper it->first=" << it->first<< ::std::endl;
		//if(it != indelMap.end()) std::cout<< it->first << " <- upperbound " << std::endl;
		if(it != indelMap.end() && it->first - randPos < (unsigned)options.minDistance)
			continue;
		
		if(it != indelMap.begin())
		{
			--it;
			if(options._debugLevel > 1)std::cout << "compare with lower it->first=" << it->first<< ::std::endl;
			if(!(indelMap.empty()) && randPos - it->first < (unsigned)options.minDistance)
				continue;
		}
		
		// is position in a stretch of or too close to Ns?
		unsigned pos=randPos-limits[chr]-options.noNsInRange;
		for (; pos < randPos-limits[chr]+options.noNsInRange; ++pos)
			if(genomes[chr][pos] == 'N') break;
		if(pos != randPos-limits[chr]+options.noNsInRange) continue;
		
		
		if(options._debugLevel > 1) std::cout << "randPos = " << randPos << std::endl;
		indelMap.insert(std::make_pair(randPos, indelSet[indelsetPos]));
		++count;

		if(options._debugLevel > 1 && count%100==0) std::cout << "." << std::flush;
		if(options._debugLevel > 1 && count%1000==0) std::cout << count<< std::flush;
		
	}
	return true;
}

    template <typename TIndelInfo>
    struct LessSize : public ::std::binary_function < TIndelInfo, TIndelInfo, bool >
    {
        inline bool operator() (TIndelInfo const &a, TIndelInfo const &b) const 
        {
            // indel size
            return (a.indelSize < b.indelSize);

        }
    };

    template <typename T>
    struct Less : public ::std::binary_function < T, T, bool >
    {
        inline bool operator() (T const &a, T const &b) const 
        {
            // indel size
            return (a < b);

        }
    };



template<typename TString>
void
generateRandomNumberString(TString &deleteStr,
			int first,
			int last,
			int numDel)
{

	//std::cout << "first=" << first << " last=" << last << " numDel=" << numDel << std::endl;
	clear(deleteStr);
	resize(deleteStr,numDel);
	std::vector<int> helpSet;
	
	for(int i = first; i < last; ++i)
		helpSet.push_back(i);
	
	for(int i = 0; i < numDel; ++i)
	{
		int randPos = (int) (rand() % helpSet.size());
		deleteStr[i] = helpSet[randPos];
		helpSet.erase(helpSet.begin()+randPos);
	}
	std::sort(begin(deleteStr,Standard()),end(deleteStr,Standard()), Less<unsigned>());
	return;
}



template<typename TOptions, typename TIndelInfo>
bool
reduceIndelSetAccordingToProbs( StringSet<TIndelInfo> & indelSet,
				TOptions &options)
{

	typedef typename Iterator<StringSet<TIndelInfo> >::Type TIndelIterator;

	// draw indels from given gff file
	if(options._debugLevel > 0) std::cout << "Reducing indel set to fit desired range probablities." << std::endl;


	// make sure range probs are given 
	if(length(options.ranges) != length(options.rangeProbs))
	{
		if(options._debugLevel > 0)std::cout << "No range probabilities given. Using same probability for each bucket.\n";
		clear(options.rangeProbs);
		resize(options.rangeProbs,length(options.ranges),0.2);
	}
	
	// make sure range probs add up to 1
	double sumProbs = 0.0;
	for(unsigned i = 0; i < length(options.rangeProbs); ++i)
		sumProbs += options.rangeProbs[i];
	if(sumProbs < EPS_PREC)
	{
		std::cout << "Something wrong with range probabilities.\nExiting...\n";
		return 1;
	}

	for(unsigned i = 0; i < length(options.rangeProbs); ++i)
		options.rangeProbs[i] /= sumProbs;
	
	if(options._debugLevel > 1)
	{
		::std::cout << "RangeProbs:" << std::endl;
		for(unsigned i = 0; i < length(options.rangeProbs); ++i)
			::std::cout << options.rangeProbs[i] << '\t';
		::std::cout << std::endl;
	}

	::std::sort(begin(indelSet,Standard()),end(indelSet,Standard()),LessSize<TIndelInfo>());

	TIndelIterator it = begin(indelSet,Standard());
	TIndelIterator endIt = end(indelSet,Standard());

	String<unsigned> rangeCounts;
	resize(rangeCounts,length(options.ranges));
	unsigned shrinkRange = 0;
	unsigned totalShrinkSize = length(indelSet);
	for (unsigned i = 0; i < length(options.rangeProbs); ++i)
	{

		while(it != endIt && (*it).indelSize < options.ranges[i].i1) //smallest indel size allowed (largest insertion or smallest del)
			++it;	//skip large insertions
		TIndelIterator keepItBegin = it;
		while(it != endIt && (*it).indelSize < options.ranges[i].i2) //smallest indel size allowed (largest insertion or smallest del)
			++it;	//skip large insertions
		rangeCounts[i] = it - keepItBegin;
		if(rangeCounts[i] / options.rangeProbs[i]  < totalShrinkSize )
		{
			totalShrinkSize = (unsigned int)(rangeCounts[i] / options.rangeProbs[i]);
			shrinkRange = i;
		}
	}
	if(totalShrinkSize < options.numIndels) 
		totalShrinkSize = _min(length(indelSet),options.numIndels);

	if(options._debugLevel > 0)
		std::cout << "Mincount range " << options.ranges[shrinkRange].i1 << " to " << options.ranges[shrinkRange].i2 << ":" << rangeCounts[shrinkRange] << "\n";
	
	if(options._debugLevel > 0)
		std::cout << "Total shrink size " << totalShrinkSize << "\n";

	it = begin(indelSet,Standard());
	TIndelIterator keepIt = begin(indelSet,Standard());
	for (unsigned i = 0; i < length(options.rangeProbs); ++i)
	{
		while(it != endIt && (*it).indelSize < options.ranges[i].i1) //smallest indel size allowed (largest insertion or smallest del)
			++it;	//skip large insertions
		TIndelIterator rangeItBegin = it;
		while(it != endIt && (*it).indelSize < options.ranges[i].i2) //smallest indel size allowed (largest insertion or smallest del)
			++it;	//skip large insertions
		TIndelIterator rangeItEnd = it;
		int numIndelsToKeepFromRange = _min(int(totalShrinkSize * options.rangeProbs[i])+1, (int)rangeCounts[i]);
			
		// generate a string of random numbers that indicates which indels are to be deleted from the indelset
		String<unsigned> deleteStr;
		generateRandomNumberString(deleteStr,(int)0,(int)rangeCounts[i],(int)rangeCounts[i]-numIndelsToKeepFromRange);
		
		//
		typename Iterator<String<unsigned> >::Type deleteIt = begin(deleteStr);
		typename Iterator<String<unsigned> >::Type deleteItEnd = end(deleteStr);
		unsigned count = 0;
		it = rangeItBegin;
//		std::cout << "Range = " << options.ranges[i].i1 << " to " << options.ranges[i].i2 << std::endl;
//		for (unsigned f = 0 ; f < length(deleteStr); ++f)
//			std::cout << deleteStr[f] << ",";
		while(it != rangeItEnd)
		{
			if(deleteIt != deleteItEnd && count == *deleteIt)
				++deleteIt;
			else{
				*keepIt = *it;
		//		std::cout << "keep: " << (*keepIt).indelSize <<std::endl;
				++keepIt;
			}
			++count;
			++it;
		}
	}
	resize(indelSet,keepIt-begin(indelSet),Exact());
	if(options._debugLevel > 0)::std::cout << "IndelSet resized to " << length(indelSet) << std::endl;
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// 1) generate a set of indels to be implanted
// 2) apply to original sequence in order to obtain manipulated sequence
template <
	typename TGenomeSet,
	typename TGenomeIDs,
	typename TOptions
>
int simulateIndels(
	TGenomeSet			&genomes,		// original reference sequence
	TGenomeIDs 			&genomeIDs,		// genome ids
	::std::map<CharString,unsigned> &gIdStringToIdNumMap,
	TGenomeSet			&simGenomes,		// manipulated sequence
	TGenomeIDs 			&simIDs,		// and corresponding ids
	::std::map<int,IndelInfo>	&indelMap,		// map of all implanted indels
	::std::map<int,SnpInfo>	    &snpMap,		// map of all implanted snps
	TOptions 			&options)		// options
{
	
	typedef typename Value<TGenomeSet>::Type 		TGenome;
	typedef typename Value<TGenome>::Type 			TAlphabet;
	//typedef typename Value<TGenomeIDs>::Type 		TGenomeID;
	typedef Dna TSimAlphabet;
	
	typedef ::std::map<int,IndelInfo> TIndelMap;
	typedef typename TIndelMap::iterator TIndelMapIt;
	
    typedef ::std::map<int,SnpInfo> TSnpMap;
    typedef typename TSnpMap::iterator TSnpMapIt;
	
	if(empty(genomes))
	{
		std::cerr << "Nothing to do, genome empty.\n";
		return 1;
	}
	if(empty(options.ranges) && *toCString(options.inputIndel) == 0)
	{
		std::cerr << "Nothing to do, neither ranges nor indel information given.\n";
		return 1;
	}
	if(options.numIndels == 0)
	{
		std::cerr << "Nothing to do, number of indels is zero.\n";
		return 1;
	}	

	// limits string to know where one sequence ends and the next starts
	String<unsigned> limits;
	resize(limits,length(genomes)+1);
	unsigned sumLen = 0;
	for(unsigned i = 0; i < length(genomes); ++i)
	{
		limits[i] = sumLen;
		sumLen += length(genomes[i]);
	}
	limits[length(genomes)] = sumLen;

	if(options._debugLevel > 0) std::cout << "Creating genome names." << std::endl;
	if(options.diploid)
    {
        resize(simIDs,length(genomeIDs)*2);
	    resize(simGenomes,length(genomes)*2);
    }
    else
    {
        resize(simIDs,length(genomeIDs));
	    resize(simGenomes,length(genomes));
    }
	
	for (unsigned i = 0; i < length(genomeIDs); ++i)
	{
		simIDs[i] = genomeIDs[i];
		append(simIDs[i], "simA");
	}
    if(options.diploid)
    {
    	for (unsigned i = 0; i < length(genomeIDs); ++i)
	    {
		    simIDs[i+length(genomes)] = genomeIDs[i];
		    append(simIDs[i+length(genomes)], "simB");
    	}
    }


	
	// get simulation instructions
	bool result = true;
	if(*toCString(options.inputIndel) != 0)
		result = simulateIndelsFromGff(indelMap,limits,genomes,gIdStringToIdNumMap,options);
	else
		result = simulateIndelsFromRanges(indelMap,genomes,limits,options);

	if(!result) 
	{
		std::cout << "WARNING! Only " << indelMap.size() << " indels could be implanted.\n";
	}

// TODO: implement sampling from input snp file    
//    if(*toCString(options.inputSnps) != 0)
//		result = simulateSnpsFromGff(snpMap,limits,genomes,gIdStringToIdNumMap,options);
//	else
		result = simulateSnps(snpMap,genomes,limits,options);
	if(!result) 
	{
		std::cout << "WARNING! Only " << snpMap.size() << " SNPs could be implanted.\n";
	}


	// first implant snps
	if(options._debugLevel > 0) std::cout << "Start implanting SNPs." << std::endl;
    TGenomeSet simGenomes1; // temporary, will have just snps, used as template in second step (indel sim)
    if(!options.diploid) resize(simGenomes1,length(genomes));
    else resize(simGenomes1,2*length(genomes));
    for (unsigned i = 0; i < length(genomes); ++i)
        simGenomes1[i] = genomes[i];
    if(options.diploid)
        for (unsigned i = 0; i < length(genomes); ++i)
            simGenomes1[length(genomes)+i] = genomes[i];
	TSnpMapIt sit = snpMap.begin();
	while(sit != snpMap.end())
	{
		unsigned chromosome = (sit->second).genomeId;
		unsigned chromosomalPos = (sit->second).originalPos; // current original position
		Dna base = (sit->second).base;
        if(options._debugLevel > 1)
        {
            std::cout << "base = " << base << "\tchromosome = " << chromosome << "\tpos = " << chromosomalPos << std::endl;
        }
        if(options.diploid)
        {
            // haplotype 1
            if ((sit->second).haplotypeA)
                simGenomes1[chromosome][chromosomalPos] = base;

            // haplotype 2
            if ((sit->second).haplotypeB)
                simGenomes1[chromosome+length(genomes)][chromosomalPos] = base;

        }
        else
            simGenomes1[chromosome][chromosomalPos] = base;

	    ++sit;	
	}


	// now implant indels
	if(options._debugLevel > 0) std::cout << "Start implanting indels." << std::endl;
	TIndelMapIt it = indelMap.begin();
	unsigned lastOriPos = 0; // lastOriPos stores the last position on the original sequence
	unsigned lastSimPos = 0; // lastSimPos stores the last position on the simulated sequence
	unsigned lastChromosome = 0;
	while(it != indelMap.end())
	{
        if(options.diploid && !(it->second).haplotypeA)// its a het, dont implant this one here
        {
            ++it; continue;
        }
		unsigned chromosome = (it->second).genomeId;
		if(chromosome != lastChromosome) // reset to local position 0 if we are on a new chromosome
		{
			// append sequence up to end
			append(simGenomes[lastChromosome],infix(simGenomes1[lastChromosome],lastOriPos,length(simGenomes1[lastChromosome])));
			lastOriPos = lastSimPos = 0;
		}
		unsigned chromosomalPos = (it->second).originalPos; // current original position
		int indelSize = (it->second).indelSize;
		bool duplication = (it->second).duplication;
		if(lastOriPos > chromosomalPos) std::cout <<"lastOriPos > chromosomalPos ! darf nicht" << std::endl;
		if(options._debugLevel > 2)
		{
			std::cout << "Position "<<chromosomalPos<<" on sequence "<<chromosome<<"." << std::endl;
			std::cout << "lastOriPos = "  << lastOriPos << " chromosomalPos = " << chromosomalPos << std::endl; 
		}
		// append sequence up to current position
		append(simGenomes[chromosome],infix(simGenomes1[chromosome],lastOriPos,chromosomalPos));
		
		
		if(options._debugLevel > 2) std::cout << "IndelSize "<<indelSize<<"." << std::endl;
		if(indelSize > 0) // deletion: remove sequence starting at chromosomalPos
		{
			lastSimPos += chromosomalPos-lastOriPos; // but not in the simulated one
			(it->second).simPos = lastSimPos;	 // remember position of indel in simulated sequence (nur so erstmal)
			lastOriPos = chromosomalPos + indelSize; // we move indelSize many positions ahead in the original seq
		}
		else	// insertion
		{
			if(duplication)	// duplicate sequence ending at current position
				append(simGenomes[chromosome],infix(simGenomes1[chromosome],chromosomalPos+indelSize,chromosomalPos));
			else		// randomly generate new insert sequence
				for(int i = 0; i < -indelSize; ++i)
					append(simGenomes[chromosome],(TAlphabet)(rand() % ValueSize<TSimAlphabet>::VALUE));
					
			lastSimPos = lastSimPos + chromosomalPos-lastOriPos - indelSize; 	// simulated one moves -indelSize many positions further
			(it->second).simPos = lastSimPos + indelSize;		// remember position of indel in simulated sequence (nur so erstmal)
			lastOriPos = chromosomalPos; 				// move to current position in original seq
		}
		lastChromosome = chromosome;
		++it;
		
	}
	// append sequence up to end
	if(lastOriPos != length(genomes[lastChromosome]))
		append(simGenomes[lastChromosome],infix(simGenomes1[lastChromosome],lastOriPos,length(genomes[lastChromosome])));
	for(unsigned i = 0; i < length(genomes); ++i)
		if(empty(simGenomes[i])) simGenomes[i] = simGenomes1[i];

    if(options.diploid)
    {
        unsigned offset = length(genomes);
        it = indelMap.begin();
    	unsigned lastOriPos = 0; // lastOriPos stores the last position on the original sequence
    	unsigned lastSimPos = 0; // lastSimPos stores the last position on the simulated sequence
	    unsigned lastChromosome = 0;
	    while(it != indelMap.end())
	    {
            if(options.diploid && !(it->second).haplotypeB)// its a het, dont implant this one here
            {
                ++it; continue;
            }
		    unsigned chromosome = (it->second).genomeId;
	    	if(chromosome != lastChromosome) // reset to local position 0 if we are on a new chromosome
	    	{
		    	// append sequence up to end
			    append(simGenomes[lastChromosome+offset],infix(simGenomes1[lastChromosome+offset],lastOriPos,length(simGenomes1[lastChromosome])));
    			lastOriPos = lastSimPos = 0;
	    	}
		    unsigned chromosomalPos = (it->second).originalPos; // current original position
		    int indelSize = (it->second).indelSize;
		    bool duplication = (it->second).duplication;
		    if(lastOriPos > chromosomalPos) std::cout <<"lastOriPos > chromosomalPos ! darf nicht" << std::endl;
		    if(options._debugLevel > 2)
		    {
			    std::cout << "Position "<<chromosomalPos<<" on sequence "<<chromosome<<"." << std::endl;
			    std::cout << "lastOriPos = "  << lastOriPos << " chromosomalPos = " << chromosomalPos << std::endl; 
		    }
	    	// append sequence up to current position
		    append(simGenomes[chromosome+offset],infix(simGenomes1[chromosome+offset],lastOriPos,chromosomalPos));
		
		
    		if(options._debugLevel > 2) std::cout << "IndelSize "<<indelSize<<"." << std::endl;
	    	if(indelSize > 0) // deletion: remove sequence starting at chromosomalPos
		    {
			    lastSimPos += chromosomalPos-lastOriPos; // but not in the simulated one
    		    (it->second).simPos = lastSimPos;	 // remember position of indel in simulated sequence (nur so erstmal)
	    		lastOriPos = chromosomalPos + indelSize; // we move indelSize many positions ahead in the original seq
		    }
    		else	// insertion
	    	{
		    	if(duplication)	// duplicate sequence ending at current position
			    	append(simGenomes[chromosome+offset],infix(simGenomes1[chromosome+offset],chromosomalPos+indelSize,chromosomalPos));
    			else		// randomly generate new insert sequence
	    			for(int i = 0; i < -indelSize; ++i)
		    			append(simGenomes[chromosome+offset],(TAlphabet)(rand() % ValueSize<TSimAlphabet>::VALUE));
			    		
    		    lastSimPos = lastSimPos + chromosomalPos-lastOriPos - indelSize; 	// simulated one moves -indelSize many positions further
	        	(it->second).simPos = lastSimPos + indelSize;		// remember position of indel in simulated sequence (nur so erstmal)
		    	lastOriPos = chromosomalPos; 				// move to current position in original seq
		    }
    		lastChromosome = chromosome;
	    	++it;
		
    	}
	    // append sequence up to end
	    if(lastOriPos != length(genomes[lastChromosome]))
		    append(simGenomes[lastChromosome+offset],infix(simGenomes1[lastChromosome+offset],lastOriPos,length(genomes[lastChromosome])));
    	for(unsigned i = 0; i < length(genomes); ++i)
	    	if(empty(simGenomes[i+offset])) simGenomes[i+offset] = simGenomes1[i+offset];

        if(options._debugLevel > 1)  std::cout << "Simulated diploid genome.\n";
    
    }

	
	return 0;
	
}




}

#endif
