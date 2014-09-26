// ==========================================================================
//                    STELLAR - SwifT Exact LocaL AligneR
//                   http://www.seqan.de/projects/stellar/
// ==========================================================================
// Copyright (C) 2010-2012 by Birte Kehr
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STELLAR_OUTPUT_H
#define SEQAN_HEADER_STELLAR_OUTPUT_H

#include <iostream>
#include <seqan/align.h>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
// Computes the length adjustment for E-value computation
// Based on the NCBI BLAST code by Tom Madden.
template<typename TSize>
TSize
_computeLengthAdjustment(TSize dbLength, TSize queryLength) {

	const double K = 0.34;
	const double logK = log(K);
	const double alphaByLambda = 1.8/1.19;
	const double beta = -3;
	const TSize maxIterations = 20;

	double n = (double)dbLength;
	double m = (double)queryLength;
	double totalLen;

	double val = 0, val_min = 0, val_max;
	bool converged = false;

	 /* Choose val_max to be the largest nonnegative value that satisfies
      *    K * (m - val) * (n - N * val) > max(m,n)
      * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */

    { // scope of mb, and c, the coefficients in the quadratic formula (the variable mb is -b, a=1 ommited)
        double mb = m + n;
        double c  = n * m - _max(m, n) / K;

        if(c < 0) {
            return 0;
        } else {
            val_max = 2 * c / (mb + sqrt(mb * mb - 4 * c));
        }
    } // end scope of mb and c

	for(TSize i = 1; i <= maxIterations; i++) {  
        totalLen = (m - val) * (n - val);
        double val_new  = alphaByLambda * (logK + log(totalLen)) + beta;  // proposed next value of val
        if(val_new >= val) { // val is no bigger than the true fixed point
            val_min = val;
            if(val_new - val_min <= 1.0) {
                converged = true;
                break;
            }
            if(val_min == val_max) { // There are no more points to check
                break;
            }
        } else { // val is greater than the true fixed point
            val_max = val;
        }
        if(val_min <= val_new && val_new <= val_max) { // ell_new is in range. Accept it.
            val = val_new;
        } else { // else val_new is not in range. Reject it.
            val = (i == 1) ? val_max : (val_min + val_max) / 2;
        }
    }

	if(converged) { // the iteration converged
        // If val_fixed is the (unknown) true fixed point, then we wish to set lengthAdjustment to floor(val_fixed).
		// We assume that floor(val_min) = floor(val_fixed)
        return (TSize) val_min;

        // But verify that ceil(val_min) != floor(val_fixed)
        val = ceil(val_min);
        if( val <= val_max ) {
          totalLen = (m - val) * (n - val);
          if(alphaByLambda * (logK + log(totalLen)) + beta >= val) {
            // ceil(val_min) == floor(val_fixed)
            return (TSize) val;
          }
        }
    } else { // the iteration did not converge
        // Use the best value seen so far.
        return (TSize) val_min;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Computes a CIGAR string and mutations from rows of StellarMatch.
template<typename TRow, typename TString>
void
_getCigarLine(TRow const & row0, TRow const & row1, TString & cigar, TString & mutations) { 
    typedef typename Size<TRow>::Type TSize;

    TSize pos = 0;

    SEQAN_ASSERT_EQ(length(row0), length(row1));
    TSize dbEndPos = length(row0);
    TSize queryEndPos = length(row1);

    bool first = true;
    TSize readBasePos = pos + clippedBeginPosition(row1);
    TSize readPos = 0;
	while (pos < dbEndPos || pos < queryEndPos) {
		int matched = 0;
		int inserted = 0;
		int deleted = 0;
		while (pos != dbEndPos && pos != queryEndPos &&
               !isGap(row0, pos) && !isGap(row1, pos)) {
            ++readPos;
			if (value(row0, pos) != value(row1, pos)) {
				if (first) first = false;
				else mutations << ",";
				mutations << readPos << value(source(row1), readBasePos);
			}
			++readBasePos;
			++pos;
			++matched;
		}
		if (matched > 0) cigar << matched << "M" ;
		while (pos < dbEndPos && isGap(row1, pos)) {
			++pos;
			++deleted;
		}
		if (deleted > 0) cigar << deleted << "D";
		while (pos < queryEndPos && isGap(row0, pos)) {
			++pos;
			++readPos;
			if (first) first = false;
			else mutations << ",";
			mutations << readPos << value(source(row1), readBasePos);
			++readBasePos;
			++inserted;
		}
		if (inserted > 0) cigar << inserted << "I";
	}
}

///////////////////////////////////////////////////////////////////////////////
// Determines the length and the number of matches of two alignment rows
template<typename TRow, typename TSize>
inline void
_analyzeAlignment(TRow const & row0, TRow const & row1, TSize & aliLen, TSize & matches) {
	TSize pos = 0;
    SEQAN_ASSERT_EQ(length(row0), length(row1));
    TSize end0 = length(row0);
    TSize end1 = length(row1);

	matches = 0;
    while ((pos < end0) && (pos < end1)) {
        if (!isGap(row0, pos) && !isGap(row1, pos)) {
            if (value(row0, pos) == value(row1, pos)) {
                ++matches;
            }
        }
        ++pos;
    }

	aliLen = _max(end0, end1);
}

///////////////////////////////////////////////////////////////////////////////
// Calculates the identity of two alignment rows (percentage of matching positions).
template<typename TRow>
double
_computeIdentity(TRow const & row0, TRow const & row1) {
    typedef typename Size<TRow>::Type TSize;
    TSize matches, aliLen;
	_analyzeAlignment(row0, row1, aliLen, matches);

    return floor(1000000.0 * matches / aliLen) / 10000.0;
}

///////////////////////////////////////////////////////////////////////////////
// Calculates the E-value from two alignment rows and a specified length adjustment
template<typename TRow, typename TSize>
double
_computeEValue(TRow & row0, TRow & row1, TSize lengthAdjustment) {
	TSize m = length(source(row0)) - lengthAdjustment;
	TSize n = length(source(row1)) - lengthAdjustment;
	double minusLambda = -1.19; // -lambda
	double K = 0.34;

	TSize matches, aliLen;
	_analyzeAlignment(row0, row1, aliLen, matches);
	// score = 1 * matches - 2 * errors (mismatches or gaps)
	//       = matches - 2 * (aliLen - matches)
	TSize score = matches - 2 * (aliLen - matches);

	return K * (double)m * (double)n * exp(minusLambda * (double)score);
}

///////////////////////////////////////////////////////////////////////////////
// Calculates the E-value for an alignment with the specified score and number
//  of matches, and the specified length of query and database sequence
template<typename TSize>
double
_computeEValue(TSize score, TSize len0, TSize len1) {
	double minusLambda = -1.19; // -lambda
	double K = 0.34;

	TSize lengthAdjustment = _computeLengthAdjustment(len0, len1);
	TSize m = len0 - lengthAdjustment;
	TSize n = len1 - lengthAdjustment;

	return K * (double)m * (double)n * exp(minusLambda * (double)score);
}

///////////////////////////////////////////////////////////////////////////////
// Writes rows of a StellarMatch in gff format to a file.
template<typename TId, typename TSize, typename TRow, typename TFile>
void
_writeMatchGff(TId const & databaseID,
              TId const & patternID,
              bool const databaseStrand,
			  TSize lengthAdjustment,
              TRow const & row0,
              TRow const & row1,
              TFile & file) {
//IOREV _recordreading_ unclear how this is related to GFF support from store_io
    typedef typename Value<typename Source<TRow>::Type>::Type TAlphabet;

    for (typename Position<TId>::Type i = 0; i < length(databaseID) && value(databaseID, i) > 32; ++i) {
        file << value(databaseID, i);
    }

    file << "\tStellar";
    file << "\teps-matches";

    if (databaseStrand) {
        file << "\t" << beginPosition(row0) + beginPosition(source(row0)) + 1;
        file << "\t" << endPosition(row0) + beginPosition(source(row0));
    } else {
        file << "\t" << length(source(row0)) - (endPosition(row0) + beginPosition(source(row0))) + 1;
        file << "\t" << length(source(row0)) - (beginPosition(row0) + beginPosition(source(row0)));
    }

    file << "\t" << _computeIdentity(row0, row1);

    file << "\t" << (databaseStrand ? '+' : '-');

    file << "\t.\t";
    for (typename Position<TId>::Type i = 0; i < length(patternID) && value(patternID, i) > 32; ++i) {
        file << value(patternID, i);
    }

	//file << ";seq2Length=" << length(source(row1));

    file << ";seq2Range=" << beginPosition(row1) + beginPosition(source(row1)) + 1;
    file << "," << endPosition(row1) + beginPosition(source(row1));

    if (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)
        file << ";eValue=" << _computeEValue(row0, row1, lengthAdjustment);

    std::stringstream cigar, mutations;
    _getCigarLine(row0, row1, cigar, mutations);
    file << ";cigar=" << cigar.str();
    file << ";mutations=" << mutations.str();
    file << "\n";
}

///////////////////////////////////////////////////////////////////////////////
// Writes rows of a StellarMatch in human readable format to file.
template<typename TId, typename TSize, typename TRow, typename TFile>
void
_writeMatch(TId const & databaseID,
            TId const & patternID,
            bool const databaseStrand,
			TSize lengthAdjustment,
            TRow const & row0,
            TRow const & row1,
            TFile & file) {
//IOREV _recordreading_ _stub_ 
    typedef typename Value<typename Source<TRow>::Type>::Type TAlphabet;

	// write database ID
	file << "Database sequence: " << databaseID;
	if (!databaseStrand) file << " (complement)" << std::endl;
	else file << std::endl;

	// write database positions
	file << "Database positions: ";
	if (databaseStrand) {
		file << beginPosition(row0) + beginPosition(source(row0));
		file << ".." << endPosition(row0) + beginPosition(source(row0));
	} else {
		file << length(source(row0)) - beginPosition(row0) + beginPosition(source(row0));
		file << ".." << length(source(row0)) - endPosition(row0) + beginPosition(source(row0));
	}
	file << std::endl;

	// write query ID
	file << "Query sequence: " << patternID << std::endl;

	// write query positions
	file << "Query positions: ";
	file << beginPosition(row1) + beginPosition(source(row1));
	file << ".." << endPosition(row1) + beginPosition(source(row1));
	file << std::endl;
    
    if (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)
    {
	    // write e-value
	    file << "E-value: " << _computeEValue(row0, row1, lengthAdjustment) << std::endl;
    }

	file << std::endl;

	// write match
	Align<typename Source<TRow>::Type> align;
	appendValue(align.data_rows, row0);
	appendValue(align.data_rows, row1);
	file << align;
	file << "----------------------------------------------------------------------\n" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Calls _writeMatchGff for each match in StringSet of String of matches.
//   = Writes matches in gff format to a file.
template<typename TInfix, typename TQueryId, typename TIds, typename TDatabases, typename TMode, typename TFile,
         typename TString>
bool
_outputMatches(StringSet<QueryMatches<StellarMatch<TInfix, TQueryId> > > & matches,
			   TIds const & ids,
			   TDatabases & databases,
			   TMode verbose,
			   TFile & fileName,
			   TString & format) {
	typedef StellarMatch<TInfix, TQueryId> TMatch;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;
	typedef typename Iterator<String<TMatch> >::Type TIterator;
    typedef typename Value<TInfix>::Type TAlphabet;

	TSize numMatches = 0;
	TSize totalLength = 0;
    TSize maxLength = 0;

	std::ofstream file;
	file.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::app);
	if (!file.is_open()) {
		std::cerr << "Could not open output file." << std::endl;
		return 1;
	}

	// output matches on positive database strand
	for (TSize i = 0; i < length(matches); i++) {
		QueryMatches<TMatch> &queryMatches = value(matches, i);

		TIterator it = begin(queryMatches.matches);
		TIterator itEnd = end(queryMatches.matches);

		if (it != itEnd && (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)) {
			queryMatches.lengthAdjustment = _computeLengthAdjustment(length(source((*it).row1)),
                                                                     length(source((*it).row2)));
		}

		while (it < itEnd) {
			TSize len = _max(length((*it).row1), length((*it).row2));
			totalLength += len;
			if(len > maxLength) maxLength = len;

			if ((*it).orientation) {
				if (format == "gff")
					_writeMatchGff((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment,
                                   (*it).row1, (*it).row2, file);
				else
					_writeMatch((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment,
                                (*it).row1, (*it).row2, file);
			}

			++it;
		}
		numMatches += length(queryMatches.matches);
	}

    if (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)
    {
	    reverseComplement(databases);

	    // output matches on negative database strand
	    for (TSize i = 0; i < length(matches); i++) {
		    QueryMatches<TMatch> &queryMatches = value(matches, i);

		    TIterator it = begin(queryMatches.matches);
		    TIterator itEnd = end(queryMatches.matches);

		    while (it < itEnd) {
			    if (!(*it).orientation) {
				    if (format == "gff")
					    _writeMatchGff((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment, (*it).row1, (*it).row2, file);
				    else 
					    _writeMatch((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment, (*it).row1, (*it).row2, file);
			    }
			    ++it;
		    }
	    }
    }

	file.close();

	std::cout << "# Eps-matches     : " << numMatches << std::endl;
	if (verbose) {
		if (numMatches > 0) {
			std::cout << "Longest eps-match : " << maxLength << std::endl;
			std::cout << "Avg match length  : " << totalLength / numMatches << std::endl;
		}
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Calls _writeMatchGff for each match in StringSet of String of matches.
//   = Writes matches in gff format to a file.
// Writes disabled query sequences to disabledFile.
template<typename TInfix, typename TQueryId, typename TQueries, typename TDatabases, typename TIds,
         typename TMode, typename TFile, typename TString>
bool 
_outputMatches(StringSet<QueryMatches<StellarMatch<TInfix, TQueryId> > > & matches, 
			   TQueries & queries,
			   TIds const & ids,
			   TDatabases & databases,
			   TMode verbose,
			   TFile & fileName,
			   TString & format,
			   TString & disabledFile) {
	typedef StellarMatch<TInfix, TQueryId> TMatch;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;
	typedef typename Iterator<String<TMatch> >::Type TIterator;
    typedef typename Value<TInfix>::Type TAlphabet;

    TSize maxLength = 0;
    TSize totalLength = 0;
    TSize numMatches = 0;
    TSize numDisabled = 0;

    std::ofstream daFile, file;
	file.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::app);
	if (!file.is_open()) {
		std::cerr << "Could not open output file." << std::endl;
		return 1;
	}

	daFile.open(toCString(disabledFile), ::std::ios_base::out | ::std::ios_base::app);
	if (!daFile.is_open()) {
		std::cerr << "Could not open file for disabled queries." << std::endl;
		return 1;
	}

	// output matches on positive database strand
    for (TSize i = 0; i < length(matches); i++) {
		QueryMatches<TMatch> &queryMatches = value(matches, i);
		if (queryMatches.disabled) {
			daFile << ">" << ids[i] << "\n";
			daFile << queries[i] << "\n\n";
			++numDisabled;
		}

		TIterator it = begin(queryMatches.matches);
		TIterator itEnd = end(queryMatches.matches);

		if (it != itEnd && (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)) {
			queryMatches.lengthAdjustment = _computeLengthAdjustment(length(source((*it).row1)),
                                                                     length(source((*it).row2)));
		}

		while (it < itEnd) {
            TSize len = _max(length((*it).row1), length((*it).row2));
            totalLength += len;
            if(len > maxLength) maxLength = len;

			if ((*it).orientation) {
				if (format == "gff")
					_writeMatchGff((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment,
                                   (*it).row1, (*it).row2, file);
				else 
					_writeMatch((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment,
                                (*it).row1, (*it).row2, file);
			}
			++it;
        }
        numMatches += length(queryMatches.matches);
	}

    if (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)
    {
	    reverseComplement(databases);
	
	    // output matches on positive database strand
	    for (TSize i = 0; i < length(matches); i++) {
		    QueryMatches<TMatch> &queryMatches = value(matches, i);
		    TIterator it = begin(queryMatches.matches);
		    TIterator itEnd = end(queryMatches.matches);

		    while (it < itEnd) {
			    if (!(*it).orientation) {
				    if (format == "gff")
					    _writeMatchGff((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment,
                                       (*it).row1, (*it).row2, file);
				    else 
					    _writeMatch((*it).id, ids[i], (*it).orientation, queryMatches.lengthAdjustment,
                                    (*it).row1, (*it).row2, file);
			    }
			    ++it;
            }
        }
    }
	daFile.close();
	file.close();

	std::cout << "# Eps-matches     : " << numMatches << std::endl;
	if (verbose) {
		if (numMatches > 0) {
			std::cout << "Longest eps-match : " << maxLength << std::endl;
			std::cout << "Avg match length  : " << totalLength / numMatches << std::endl;
		}
		std::cout << "# Disabled queries: " << numDisabled << std::endl;
	}

	return 0;
}

#endif
