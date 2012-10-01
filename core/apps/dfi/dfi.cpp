 /*==========================================================================
                     Dfi - The Deferred Frequency Index
                   http://www.seqan.de/projects/dfi.html

 ============================================================================
  This is an application of the Dfi algorithm in
  "Efficient string mining under constraints via the deferred frequency index"

 ============================================================================
  Copyright (C) 2008 by David Weese and Marcel H. Schulz

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/index.h>
#include <../../extras/include/seqan/math.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace seqan;

//#define DEBUG_ENTROPY

typedef double TFloat;
#define DFI_PLUS_EPSILON  + 0.0000001
#define DFI_MINUS_EPSILON - 0.0000001

//typedef Rational<__int64> TFloat;
//#define DFI_PLUS_EPSILON
//#define DFI_MINUS_EPSILON

//////////////////////////////////////////////////////////////////////////////
// predicates for the Frequent Pattern Mining Problem

	// minimal frequency predicate
	struct PredMinFreq 
	{	
		String<unsigned> minFreq;

		template <typename TDataSet>
		PredMinFreq(String<unsigned> _minFreq, TDataSet const &):
			minFreq(_minFreq) {}
			
		inline bool operator()(DfiEntry_ const &entry) const {
			for (unsigned i = 0; i < length(entry.freq); ++i)
				if (entry.freq[i] < minFreq[i])
					return false;
			return true;
		}
	};

	// maximal frequency predicate
	struct PredMaxFreq 
	{	
		String<unsigned> maxFreq;

		template <typename TDataSet>
		PredMaxFreq(String<unsigned> _maxFreq, TDataSet const &):
			maxFreq(_maxFreq) {}
			
		inline bool operator()(DfiEntry_ const &entry) const {
			for (unsigned i = 0; i < length(entry.freq); ++i)
				if (entry.freq[i] > maxFreq[i])
					return false;
			return true;
		}
	};


//////////////////////////////////////////////////////////////////////////////
// predicates for the Emerging Substring Mining Problem

	// minimal support predicate for D0
	struct PredMinSupp
	{	
		unsigned minFreq;

		template <typename TDataSet>
		PredMinSupp(TFloat _minSupp, TDataSet const &ds)
		{
			// emerging substring mode
			
			if ((_minSupp * (TFloat)ds[1]) < (TFloat)1) {
				cerr << "Support must be at least 1/|db_1|... exit!" << endl;
				exit(1);
			}
			// adapt parameters from support to frequency
			minFreq = (unsigned) ceil(_minSupp * (TFloat)ds[1] DFI_MINUS_EPSILON);
		}
			
		inline bool operator()(DfiEntry_ const &entry) const {
			return entry.freq[0] >= minFreq;
		}
	};

	// minimal growth predicate from D1->D0
	struct PredEmerging
	{
		TFloat growthRate;

		template <typename TDataSet>
		PredEmerging(TFloat _growthRate, TDataSet const &ds) {
			growthRate = _growthRate * ((TFloat) ds[1] / (TFloat) (ds[2] - ds[1]) DFI_MINUS_EPSILON);
		}
			
		// HINT: here growthRate is frequency-related, not support-related
		inline bool operator()(DfiEntry_ const &entry) const {
			return (TFloat)entry.freq[0] >= (TFloat)entry.freq[1] * growthRate;
		}
	};


//////////////////////////////////////////////////////////////////////////////
// predicates for the Maximum Entropy Mining Problem

	// minimal support predicate for at least one dataset
	struct PredMinAllSupp
	{	
		String<unsigned> minFreq;

		template <typename TDataSet>
		PredMinAllSupp(TFloat _minSupp, TDataSet const &ds)
		{
			resize(minFreq, length(ds) - 1, Exact());
			// adapt parameters from support to frequency
			for (unsigned i = 1; i < length(ds); ++i)
				minFreq[i - 1] = (unsigned) ceil(_minSupp * (TFloat)(ds[i] - ds[i - 1]) DFI_MINUS_EPSILON);
		}
			
		inline bool operator()(DfiEntry_ const &entry) const {
			for (unsigned i = 0; i < length(entry.freq); ++i)
				if (entry.freq[i] >= minFreq[i])
					return true;
			return false;
		}
	};

	// maximum entropy predicate
	struct PredEntropy
	{
		double maxEntropy;
		String<TFloat> dsLengths;

		template <typename TDataSet>
		PredEntropy(double _maxEntropy, TDataSet const &ds):
			maxEntropy(_maxEntropy DFI_PLUS_EPSILON)
		{
			resize(dsLengths, length(ds) - 1, Exact());
			for (unsigned i = 1; i < length(ds); ++i)
				dsLengths[i - 1] = ds[i] - ds[i - 1];
		}

		// support based entropy
		inline double
		getEntropy(DfiEntry_ const &entry) const
		{
			TFloat sum = 0;
			double H = 0;

			for (unsigned i = 0; i < length(entry.freq); ++i)
				sum += (TFloat)entry.freq[i] / (TFloat)dsLengths[i];
			
			double lSum = log((double)sum);					// sum cannot be zero
				
			for (unsigned i = 0; i < length(entry.freq); ++i)
				if (entry.freq[i])
				{
					double freq = (TFloat)entry.freq[i] / (TFloat)dsLengths[i];
					H += freq * (log(freq) - lSum);
				}
			H /= (double)-sum * log((double)length(dsLengths));		// normalize by datasets (divide by log m)
			return H;
		}
			
		inline bool operator()(DfiEntry_ const &entry) const 
		{
			return getEntropy(entry) <= maxEntropy;
		}
	};




//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta files
//
// seq ........ StringSet containing all sequences
// fileName ... array of database file names
// fileCount .. number of databases
// ds ......... seq[ds[i]],seq[ds[i]+1],...,seq[ds[i+1]-1] are the seqs. of database i
//
template <typename TSequences, typename TFileNames, typename TDatasets>
bool loadDatasets(
	TSequences		&seqs, 
	TFileNames		const &fileNames, 
	TDatasets		&ds)
{
	resize(ds, length(fileNames) + 1);
	ds[0] = 0;

	CharString seq;
	MultiFasta multiFasta;
	for(unsigned s = 0; s < length(fileNames); ++s)
	{
		if (!open(multiFasta.concat, toCString(fileNames[s]), OPEN_RDONLY)) return false;
		AutoSeqFormat format;
		guessFormat(multiFasta.concat, format);
		split(multiFasta, format);		
		unsigned seqCount = length(multiFasta);
		
		ds[s + 1] = ds[s] + seqCount;
		resize(seqs, ds[s + 1]);
		for(unsigned i = 0; i < seqCount; ++i)
			assignSeq(seqs[ds[s] + i], multiFasta[i], format);

		close(multiFasta.concat);
	}
	return (back(ds) > 0);
}


template <typename TSize>
struct SubstringEntry
{
	Pair<unsigned>	lPos;
	unsigned		len, freqSum;
	Pair<TSize>		range;
};

template <typename TSubstringEntry>
struct LessSubstringEnd : public binary_function<TSubstringEntry, TSubstringEntry, bool >
{
	inline bool operator() (TSubstringEntry const &a, TSubstringEntry const &b) const 
	{
		// sequence number
		if (a.lPos.i1 < b.lPos.i1) return true;
		if (a.lPos.i1 > b.lPos.i1) return false;

		// end position
		unsigned x = a.lPos.i2 + a.len;
		unsigned y = b.lPos.i2 + b.len;		
		if (x < y) return true;
		if (x > y) return false;
		
		x = a.range.i2 - a.range.i1;
		y = b.range.i2 - b.range.i1;
		if (x < y) return true;
		if (x > y) return false;
		
		return a.len > b.len;
	}
};

template <typename TSubstringEntry>
struct LessRange : public binary_function<TSubstringEntry, TSubstringEntry, bool >
{
	inline bool operator() (TSubstringEntry const &a, TSubstringEntry const &b) const 
	{
		// left border
		if (a.range.i1 < b.range.i1) return true;
		if (a.range.i1 > b.range.i1) return false;

		// right border
		if (a.range.i2 < b.range.i2) return true;
		if (a.range.i2 > b.range.i2) return false;

		return a.len > b.len;
	}
};

template <typename TSubstringEntry>
struct LessLex : public binary_function<TSubstringEntry, TSubstringEntry, bool >
{
	inline bool operator() (TSubstringEntry const &a, TSubstringEntry const &b) const 
	{
		if (a.range.i1 < b.range.i1) return true;
		if (a.range.i1 > b.range.i1) return false;
		return a.len < b.len;
	}
};

template <typename TMatchString>
inline void compactMatches(TMatchString &matches)
{
	typedef typename Iterator<TMatchString, Standard>::Type TIter;
	TIter src = begin(matches, Standard());
	TIter dst = src;
	TIter srcEnd = end(matches, Standard());
	
	unsigned lastSeq = ~0;
	unsigned lastPos = 0;
	unsigned lastRange = 0;
	
	if (src == srcEnd) return;
	for (; src != srcEnd; ++src)
	{
		if ((*src).len == ~0u)
		{
			*dst = *src;
			++dst;
			continue;
		}
		if (((*src).lPos.i1 == lastSeq) && ((*src).lPos.i2 + (*src).len == lastPos) && ((*src).range.i2 - (*src).range.i1 == lastRange))
			continue;
		lastSeq = (*src).lPos.i1;
		lastPos = (*src).lPos.i2 + (*src).len;
		lastRange = (*src).range.i2 - (*src).range.i1;
		*dst = *src;
		++dst;
	}
	resize(matches, dst - begin(matches, Standard()));
}

template <typename TMatchString>
inline void compactSameParentFreqMatches(TMatchString &matches)
{
	typedef typename Iterator<TMatchString, Standard>::Type TIter;
	TIter src = begin(matches, Standard());
	TIter dst = src;
	TIter srcEnd = end(matches, Standard());
	
	unsigned r1 = ~0;
	unsigned r2 = ~0;
	
	if (src == srcEnd) return;
	for (; src != srcEnd; ++src)
	{
		if (((*src).range.i1 == r1) && ((*src).range.i2 == r2))
			continue;
		r1 = (*src).range.i1;
		r2 = (*src).range.i2;
		if ((*src).len != ~0u)
		{
			*dst = *src;
			++dst;
		}
	}
	resize(matches, dst - begin(matches, Standard()));
}

template <typename TMatchString, typename TIndex, typename TDBLookup, typename TSeen, typename TEntry>
inline void compactSameSuffLinkFreqMatches(TMatchString &matches, TIndex const &index, TDBLookup const &dbLookup, TSeen &seen, TEntry &entry)
{
	typedef typename Iterator<TMatchString, Standard>::Type TIter;
	typedef typename Fibre<TIndex, FibreSA>::Type TSA;
	typedef typename Iterator<TSA, Standard>::Type TSAIter;

	TIter src = begin(matches, Standard());
	TIter dst = src;
	TIter srcEnd = end(matches, Standard());
	unsigned lastSeq = ~0;
	unsigned lastPos = 0;
	unsigned lastFreqSum = ~0;
	
	for (; src != srcEnd; ++src)
	{
		// count frequencies (debug)
		TSAIter oc = begin(indexSA(index), Standard()) + (*src).range.i1;
		TSAIter ocEnd = begin(indexSA(index), Standard()) + (*src).range.i2;
		arrayFill(begin(seen, Standard()), end(seen, Standard()), false);
		arrayFill(begin(entry.freq, Standard()), end(entry.freq, Standard()), 0);
		unsigned freqSum = 0;
		for (; oc != ocEnd; ++oc)
		{
			unsigned seqNo = getSeqNo(*oc, stringSetLimits(index));
			if (!seen[seqNo])
			{
				seen[seqNo] = true;
				++entry.freq[dbLookup[seqNo]];
				++freqSum;
			}
		}
		if (((*src).lPos.i1 == lastSeq) && ((*src).lPos.i2 + (*src).len == lastPos) && (lastFreqSum == freqSum))
			continue;
		lastSeq = (*src).lPos.i1;
		lastPos = (*src).lPos.i2 + (*src).len;
		lastFreqSum = freqSum;
		*dst = *src;
		++dst;
	}
	resize(matches, dst - begin(matches, Standard()));
}

namespace seqan {
/*
	template < typename TObject, typename TPredHull, typename TPred >
	struct Fibre< Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > >, FibreSA> 
	{
		typedef Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > > TIndex;
		typedef String< 
			typename SAValue<TIndex>::Type,
			MMap<>
		> Type;
	};
	template < typename TObject, typename TPredHull, typename TPred >
	struct Fibre< Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > > const, FibreSA>:
		public struct Fibre< Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > >, FibreSA> {};
*/

/*	template < typename TObject, typename TPredHull, typename TPred >
	struct Fibre< Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > >, FibreDir> 
	{
		typedef Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > > TIndex;
		typedef String< 
			typename Size<TIndex>::Type,
			MMap<>
		> Type;
	};
	template < typename TObject, typename TPredHull, typename TPred >
	struct Fibre< Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > > const, FibreDir>: 
		public struct Fibre< Index<TObject, IndexWotd< Dfi<TPredHull, TPred> > >, FibreDir> {};
*/

	template <typename TInt>
	inline bool
	_convertOptionValue(CommandLineOption const & opt, Rational<TInt> & dst, CharString const & src)
	{
		if (!isDoubleOption(opt)) return false;
		std::istringstream stream(toCString(src));
		return !(stream >> dst).fail();
	}

}

//////////////////////////////////////////////////////////////////////////////
// Create Dfi and output substrings within constraints band
//
// TPred .......... frequency predicate
// TPredHull ...... monotonic hull of the predicate above
// TAlphabet ...... sequence alphabet (Dna, Dna5, AminoAcid, char, ...)
//
// fileName ....... array of database file names
// paramPred ... .. parameter for the frequency predicate
// paramPredHull .. parameter for the monotonic hull
//
template <
	typename TPredHull, 
	typename TPred, 
	typename TAlphabet,
	typename TParamPredHull,
	typename TParamPred,
	typename TFileNames
>
int runDFI(
	TFileNames		const &fileNames,
	TParamPredHull	paramPredHull,
	TParamPred		paramPred,
	bool			maximal)
{
	typedef String<TAlphabet, Alloc<> >								TString;
	typedef StringSet<TString>										TStringSet;
	typedef Index<TStringSet, IndexWotd<
		Dfi<TPredHull, TPred> > >									TIndex;
	typedef Iter<TIndex, VSTree<TopDown<ParentLinks<> > > >			TIter;
	typedef SubstringEntry<typename Size<TIndex>::Type>				TSubstringEntry;

	String<unsigned>	ds;
	TStringSet			mySet;

	if (!loadDatasets(mySet, fileNames, ds)) {
		cerr << "Database read error... exit!" << endl;
		return 1;
	}

	TPred					pred(paramPred, ds);
	TPredHull				predHull(paramPredHull, ds);
	TIndex					index(mySet, predHull, pred);
	Pair<unsigned>			lPos;
	String<TSubstringEntry>	matches;

	// set index partition of sequences into datasets
	index.ds = ds;

	// database lookup table
	typedef typename Fibre<TIndex, FibreSA>::Type TSA;
	typedef typename Iterator<TSA, Standard>::Type TSAIter;
	String<unsigned>	dbLookup;
	String<bool>		seen;
	DfiEntry_			entry;
	
	resize(dbLookup, length(mySet));
	resize(seen, length(mySet));
	resize(entry.freq, length(ds) - 1);
	for (unsigned d = 0, i = 0; i < length(mySet); ++i)
	{
		while (ds[d + 1] == i) ++d;
		dbLookup[i] = d;
	}
#ifdef DEBUG_ENTROPY	
	PredEntropy			entrp(0, ds);
	unsigned			freqSumLast = ~0;
#endif

	TIter it(index);
	goBegin(it);
	
	
	if (maximal)
	{
		for (unsigned m = 0; !atEnd(it); goNext(it), ++m)
		{
			resize(matches, m + 1, Generous());
			posLocalize(matches[m].lPos, getOccurrence(it), stringSetLimits(container(it)));
			matches[m].range = range(it);
			matches[m].len = repLength(it);
			if (dirAt(value(it).node, index) & TIndex::DFI_PARENT_FREQ)
			{
				++m;
				resize(matches, m + 1, Generous());
				matches[m].lPos.i1 = 0;
				matches[m].lPos.i2 = 0;
				matches[m].range = range(index, nodeUp(it));
				matches[m].len = ~0u;
			}
		}
		
		sort(begin(matches, Standard()), end(matches, Standard()), LessSubstringEnd<TSubstringEntry>());
		compactMatches(matches);
		sort(begin(matches, Standard()), end(matches, Standard()), LessRange<TSubstringEntry>());
		compactSameParentFreqMatches(matches);
		sort(begin(matches, Standard()), end(matches, Standard()), LessSubstringEnd<TSubstringEntry>());
		compactSameSuffLinkFreqMatches(matches, index, dbLookup, seen, entry);
		sort(begin(matches, Standard()), end(matches, Standard()), LessLex<TSubstringEntry>());
		
		typedef typename Iterator<String<TSubstringEntry>, Standard>::Type TMatchIter;
		TMatchIter mit = begin(matches, Standard());
		TMatchIter mitEnd = end(matches, Standard());
		
		for (; mit != mitEnd; ++mit)
		{
#ifdef DEBUG_ENTROPY
			// count frequencies (debug)
			TSAIter oc = begin(indexSA(index), Standard()) + (*mit).range.i1;
			TSAIter ocEnd = begin(indexSA(index), Standard()) + (*mit).range.i2;
			arrayFill(begin(seen, Standard()), end(seen, Standard()), false);
			arrayFill(begin(entry.freq, Standard()), end(entry.freq, Standard()), 0);
			unsigned freqSum = 0;
			for (; oc != ocEnd; ++oc)
			{
				unsigned seqNo = getSeqNo(*oc, stringSetLimits(index));
				if (!seen[seqNo])
				{
					seen[seqNo] = true;
					++entry.freq[dbLookup[seqNo]];
					++freqSum;
				}
			}
				
			double H = entrp.getEntropy(entry);
			if (H <= 0.0) H = 0.0;
//			if (freqSum != freqSumLast)
#endif
			{
#ifdef DEBUG_ENTROPY
				cout << left << setw(14) << H << "[";
				for (unsigned i = 0; i < length(entry.freq); ++i)
					cout << right << setw(6) << entry.freq[i];
				cout << "]      \"";
#endif
				cout << infix(
					mySet[getSeqNo((*mit).lPos)], 
					getSeqOffset((*mit).lPos),
					getSeqOffset((*mit).lPos) + (*mit).len);
#ifdef DEBUG_ENTROPY
				cout << "\"";
				freqSumLast = freqSum;
#endif
				cout << endl;
			}
		}
	} 
	else 
	{
		for (; !atEnd(it); goNext(it))
		{
			posLocalize(lPos, getOccurrence(it), stringSetLimits(container(it)));
			unsigned len = repLength(it);
			for(unsigned l = parentRepLength(it) + 1; l <= len; ++l)
			{
#ifdef DEBUG_ENTROPY
				// count frequencies (debug)
				typedef typename Infix< typename Fibre<TIndex, FibreSA>::Type const >::Type TOccs;
				typedef typename Iterator<TOccs, Standard>::Type TOccIter;
	            TOccs occs = getOccurrences(it);
	            TOccIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
				arrayFill(begin(seen, Standard()), end(seen, Standard()), false);
				arrayFill(begin(entry.freq, Standard()), end(entry.freq, Standard()), 0);
				for (; oc != ocEnd; ++oc)
				{
					unsigned seqNo = getSeqNo(*oc, stringSetLimits(index));
					if (!seen[seqNo])
					{
						seen[seqNo] = true;
						++entry.freq[dbLookup[seqNo]];
					}
				}
					
				double H = entrp.getEntropy(entry);
				if (H <= 0.0) H = 0.0;
				cout << left << setw(14) << H << "[";
				for (unsigned i = 0; i < length(entry.freq); ++i)
					cout << right << setw(6) << entry.freq[i];
				cout << "]      \"";
#endif
				cout << infix(
	                mySet[getSeqNo(lPos)],
	                getSeqOffset(lPos),
	                getSeqOffset(lPos) + l);
#ifdef DEBUG_ENTROPY
			cout << "\"";
#endif
			cout << endl;
			}
		}
	}
	
	return 0;
}

int main(int argc, const char *argv[])
{
	int optionAlphabet = 0;   // 0..char, 1..protein, 2..dna
	int optionPredicate = -1; // 0..minmax, 1..growth, 2..entropy
	String<unsigned> optionMinFreq;
	String<unsigned> optionMaxFreq;
	TFloat optionMinSupp = 0;
	TFloat optionGrowthRate = 0;
	double optionEntropy = 0;
	bool optionMaximal = false;
		
	CommandLineParser parser;
	string rev = "$Revision$";
	addVersionLine(parser, "Dfi version 2.0 20100107 [" + rev.substr(11, 4) + "]");

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "**************************************************************");
	addTitleLine(parser, "***           Dfi - The Deferred Frequency Index           ***");
	addTitleLine(parser, "*** (c) Copyright 2010 by David Weese and Marcel H. Schulz ***");
	addTitleLine(parser, "**************************************************************");
	addUsageLine(parser, "[OPTION]... --minmax  <min_1> <max_1> <database 1> ... --minmax <min_m> <max_m> <database m>");
	addUsageLine(parser, "[OPTION]... --growth  <rho_s> <rho_g> <database 1> <database 2>");
	addUsageLine(parser, "[OPTION]... --entropy <rho_s> <alpha> <database 1> <database 2> ... <database m>");
	
	addOption(parser, CommandLineOption("f",  "minmax",  2, "solve Frequent Pattern Mining Problem", OptionType::Int | OptionType::Label | OptionType::List));
	addOption(parser, CommandLineOption("g",  "growth",  2, "solve Emerging Substring Mining Problem", OptionType::Double | OptionType::Label));
	addOption(parser, CommandLineOption("e",  "entropy", 2, "solve Entropy Mining Problem", OptionType::Double | OptionType::Label));
	addHelpLine(parser, "");
	addOption(parser, CommandLineOption("p",  "protein",    "use AminoAcid alphabet (for proteomes)", OptionType::Boolean));
	addOption(parser, CommandLineOption("d",  "dna",        "use DNA alphabet (for genomes)", OptionType::Boolean));
	addOption(parser, CommandLineOption("m",  "maximal",    "output only left and right maximal substrings", OptionType::Boolean));
	addHelpLine(parser, "The default is byte alphabet");

	if (argc == 1)
	{
		shortHelp(parser, cerr);	// print short help and exit
		return 0;
	}

	bool stop = !parse(parser, argc, argv, cerr);

	//////////////////////////////////////////////////////////////////////////////
	// Extract options
	getOptionValueLong(parser, "maximal", optionMaximal);
	if (isSetLong(parser, "protein")) optionAlphabet = 1;
	if (isSetLong(parser, "dna")) optionAlphabet = 2;
	if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit
	
	//////////////////////////////////////////////////////////////////////////////
	// Check options
	if (isSetLong(parser, "minmax"))
	{
		optionPredicate = 0;
		unsigned cons = length(getOptionValuesLong(parser, "minmax")) / 2;
		resize(optionMinFreq, cons);
		resize(optionMaxFreq, cons);
		for (unsigned d = 0; d < cons; ++d)
		{
			getOptionValueLong(parser, "minmax", 2*d,   optionMinFreq[d]);
			getOptionValueLong(parser, "minmax", 2*d+1, optionMaxFreq[d]);
			if ((optionMinFreq[d] < 1) && (stop = true))
				cerr << "Minimum frequency threshold must be at least 1." << endl;
			if ((optionMaxFreq[d] < 1) && (stop = true))
				cerr << "Maximum frequency threshold must be greater than 0." << endl;
		}
		if ((argumentCount(parser) < cons) && (stop = true))
			cerr << "Please specify " << cons << " databases." << endl;
		if ((argumentCount(parser) > cons) && (stop = true))
			cerr << "Please specify " << argumentCount(parser) << " min/max constraints." << endl;
	}
	if (isSetLong(parser, "growth"))
	{
		optionPredicate = 1;
		getOptionValueLong(parser, "growth", 0, optionMinSupp);
		getOptionValueLong(parser, "growth", 1, optionGrowthRate);
		if ((optionMinSupp <= (TFloat)0 || optionMinSupp > (TFloat)1) && (stop = true))
			cerr << "Support threshold must be greater than 0 and less than or equal to 1." << endl;
		if ((optionGrowthRate < (TFloat)1) && (stop = true))
			cerr << "Growth rate must not be less than 1." << endl;
		if ((argumentCount(parser) != 2) && (stop = true))
			cerr << "Please specify 2 databases." << endl;
	}
	if (isSetLong(parser, "entropy"))
	{
		optionPredicate = 2;
		getOptionValueLong(parser, "entropy", 0, optionMinSupp);
		getOptionValueLong(parser, "entropy", 1, optionEntropy);
		if ((optionMinSupp <= (TFloat)0 || optionMinSupp > (TFloat)1) && (stop = true))
			cerr << "Support threshold must be greater than 0 and less than or equal to 1." << endl;
		if ((optionEntropy <= 0.0 || optionEntropy > 1.0) && (stop = true))
			cerr << "Entropy must not be grater than 0 and less or equal to 1." << endl;
		if ((argumentCount(parser) < 2) && (stop = true))
			cerr << "Please specify at least 2 databases." << endl;
	}

	if (stop)
	{
		cerr << "Exiting ..." << endl;
		return -1;
	}

	switch (optionPredicate)
	{
		case 0:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinFreq, PredMaxFreq, unsigned char> (getArgumentValues(parser), optionMinFreq, optionMaxFreq, optionMaximal);
				case 1: return runDFI<PredMinFreq, PredMaxFreq, AminoAcid> (getArgumentValues(parser), optionMinFreq, optionMaxFreq, optionMaximal);
				case 2: return runDFI<PredMinFreq, PredMaxFreq, Dna> (getArgumentValues(parser), optionMinFreq, optionMaxFreq, optionMaximal);
			}
		case 1:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinSupp, PredEmerging, unsigned char> (getArgumentValues(parser), optionMinSupp, optionGrowthRate, optionMaximal);
				case 1: return runDFI<PredMinSupp, PredEmerging, AminoAcid> (getArgumentValues(parser), optionMinSupp, optionGrowthRate, optionMaximal);
				case 2: return runDFI<PredMinSupp, PredEmerging, Dna> (getArgumentValues(parser), optionMinSupp, optionGrowthRate, optionMaximal);
			}
		case 2:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinAllSupp, PredEntropy, unsigned char> (getArgumentValues(parser), optionMinSupp, optionEntropy, optionMaximal);
				case 1: return runDFI<PredMinAllSupp, PredEntropy, AminoAcid> (getArgumentValues(parser), optionMinSupp, optionEntropy, optionMaximal);
				case 2: return runDFI<PredMinAllSupp, PredEntropy, Dna> (getArgumentValues(parser), optionMinSupp, optionEntropy, optionMaximal);
			}
	}
	cerr << "Please choose a mining problem." << endl;
	cerr << "Exiting ..." << endl;
	return -1;
}
