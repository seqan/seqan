// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011-2018, Kathrin Trappe, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_APPS_GUSTAF_CREATE_STELLARMATCHES_FROM_FILE_H_
#define SEQAN_APPS_GUSTAF_CREATE_STELLARMATCHES_FROM_FILE_H_

#include <iostream>
#include <fstream>
#include <seqan/file.h>
#include <seqan/parse_lm.h>
#include "../stellar/stellar.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Function _getShortId()
// ----------------------------------------------------------------------------

// Creates a short Id out of a long one (i.e. it takes the prefix til the first white space)
template <typename TId>
inline void _getShortId(TId & shortId, TId const & longId)
{
    clear(shortId);
    for (typename Position<TId>::Type i = 0; i < length(longId) && isgraph(value(longId, i)); ++i)
    {
        appendValue(shortId, value(longId, i));
    }
}

// ----------------------------------------------------------------------------
// Function _createStellarMatches()
// ----------------------------------------------------------------------------

// Takes the values from the localMatchStore and creates a Stellar match out of them.
// The Stellar match is, sorted by read, appended to stQueryMatches. Reads are identified by their short query IDs
// (sQueryIds), the correct database by their short database IDs (databaseIds), which are both given in the GFF file
template <typename TSequence, typename TId>
bool _createStellarMatches(StringSet<TSequence> & queries,
                           StringSet<TId> const & sQueryIds,
                           StringSet<TSequence> & databases,
                           StringSet<TId> const & databaseIds,
                           LocalMatchStore<> & lmStore,
                           StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & stQueryMatches,
                           unsigned numThreads)
{
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;
    typedef typename Row<TAlign>::Type TRow;
    typedef typename LocalMatchStore<>::TMatchStore TMatchStore;
    typedef typename Iterator<TMatchStore, Standard>::Type TMatchStoreIterator;

    omp_set_num_threads(numThreads);
    Splitter<TMatchStoreIterator> setSplitter(begin(lmStore.matchStore, Standard()), end(lmStore.matchStore, Standard())); 

    SEQAN_OMP_PRAGMA(parallel for shared(stQueryMatches))
    for (int jobId = 0; jobId < static_cast<int>(length(setSplitter)); ++jobId)
    {
            for (TMatchStoreIterator it = setSplitter[jobId]; it != setSplitter[jobId + 1]; ++it)
	    {
		TInfix dbInf;
		TInfix queryInf;

		// Take match query ID and find the right index in queries
		unsigned iDB = std::numeric_limits<unsigned>::max();     // position of database sequence in databases
		unsigned iQuery = std::numeric_limits<unsigned>::max();  // position of query/read sequence in queries

		// Takes the short chromosome Id (from the Stellar match file) and looks up the corresponding long chromosome
		// Id entry from the reference input file
		for (unsigned j = 0; j < length(databaseIds); ++j)
		{
		    if (lmStore.sequenceNameStore[(*it).subjectId] == databaseIds[j])
		    {
			iDB = j;
			break;
		    }
		}

		// Takes the short read Id (from the Stellar match file) and looks up the corresponding long read Id entry
		// from the read input file
		for (unsigned j = 0; j < length(sQueryIds); ++j)
		{
		    if (lmStore.sequenceNameStore[(*it).queryId] == sQueryIds[j])
		    {
			iQuery = j;
			break;
		    }
		}
		// Sanity check for read and query Id:
		// skips entry if no corresponding entry in the input file could not be found, else creates StellarMatch object
		if (iDB == std::numeric_limits<unsigned>::max() || iQuery == std::numeric_limits<unsigned>::max())
		{
		    std::cerr << "Read or database does not exist for match: " << it - begin(lmStore.matchStore, Standard())
			      << " subjectId: " << lmStore.sequenceNameStore[(*it).subjectId]
			      << " queryId: " << lmStore.sequenceNameStore[(*it).queryId] << std::endl;
		    std::cerr << "Skipping entry" << std::endl;
		    continue;
		}

		// Sanity checks for database and read positions in case a wrong database or read file have been used
		// Checking for valid begin and end positions within identified database
		SEQAN_ASSERT_LEQ_MSG((*it).subjectBeginPos, length(databases[iDB]),
				     "Match begin position exceeds database length! Wrong genome?");
		SEQAN_ASSERT_LEQ_MSG((*it).subjectEndPos, length(databases[iDB]),
				     "Match end position exceeds database length! Wrong genome?");
		// Checking for valid begin and end positions within identified query
		SEQAN_ASSERT_LEQ_MSG((*it).queryBeginPos, length(queries[iQuery]),
				     "Match begin position exceeds query length! Wrong read/contig sequence?");
		SEQAN_ASSERT_LEQ_MSG((*it).queryEndPos, length(queries[iQuery]),
				     "Match end position exceeds query length! Wrong read/contig sequence?");

		if (((*it).subjectBeginPos  > length(databases[iDB]))
		   || ((*it).subjectEndPos > length(databases[iDB])))
		{
		    std::cerr << "Match begin or end position exceeds database length! Wrong genome?" << std::endl;
		    //return 1;
		}
		if (((*it).queryBeginPos > length(queries[iQuery]))
		   || ((*it).queryEndPos > length(queries[iQuery])))
		{
		    std::cerr << "Match begin or end position exceeds query length! Wrong read/contig sequence?" << std::endl;
		    std::cerr <<  (*it).queryBeginPos << " " << (*it).queryEndPos << " " <<
		    length(queries[iQuery]) << " " << sQueryIds[iQuery] << std::endl;
		    //return 1;
		}

		// Checking orientation and swapping positions for reverse matches to apply them to stellar format
		bool orientation = true;
		if ((*it).subjectBeginPos > (*it).subjectEndPos)
		{
		    orientation = false;
                    std::swap((*it).subjectBeginPos, (*it).subjectEndPos);
		}
		// Computing infices for alignment rows for stellar matches
		queryInf = infix(queries[iQuery], (*it).queryBeginPos, (*it).queryEndPos);
		dbInf = infix(databases[iDB], (*it).subjectBeginPos, (*it).subjectEndPos);

		// Creating align object for stellar format
		TAlign localAlign;
		resize(rows(localAlign), 2);
		setSource(row(localAlign, 0), host(dbInf));
		setSource(row(localAlign, 1), host(queryInf));
		TRow & row1 = row(localAlign, 0);
		TRow & row2 = row(localAlign, 1);


		// Set begin and end positions of align
		setBeginPosition(row2, (*it).queryBeginPos);
		setBeginPosition(row1, (*it).subjectBeginPos);

		// setBeginPosition(row1, 0);
		// setBeginPosition(row2, 0);
		setEndPosition(row1, (*it).subjectEndPos);
		setEndPosition(row2, (*it).queryEndPos);

		unsigned gapIndex = 0;
		// Inserting gaps into rows according to cigar line
		if (length(lmStore.cigarStore) > (*it).id)
		{
		    String<CigarElement<> > const & cigar = lmStore.cigarStore[(*it).id];
		    for (unsigned j = 0; j < length(cigar); ++j)
		    {
			// std::cout << cigar[j].count << cigar[j].operation;
			if (cigar[j].operation == 'I')
			{
			    for (unsigned gap = 0; gap < cigar[j].count; ++gap)
			    {
				insertGap(row1, gapIndex);
				++gapIndex;
			    }
			}
			else if (cigar[j].operation == 'D')
			{
			    for (unsigned gap = 0; gap < cigar[j].count; ++gap)
			    {
				insertGap(row2, gapIndex);
				++gapIndex;
			    }
			}
			else
                        {
			    gapIndex += cigar[j].count;
                        }
		    }
		}
		// Create Stellar match and append it to stQueryMatches
		StellarMatch<TSequence, TId> match(localAlign, databaseIds[iDB], orientation);
		appendValue(stQueryMatches[iQuery].matches, match);
	    }
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function _getStellarMatchesFromFile()
// ----------------------------------------------------------------------------

// Reads in a file with Stellar matches in gff format and creates StellarMatch object from the entries
template <typename TSequence, typename TId, typename TMatches>
bool _getStellarMatchesFromFile(StringSet<TSequence> & queries,
                                StringSet<TId> & sQueryIds,
                                StringSet<TSequence> & databases,
                                StringSet<TId> & databaseIDs,
                                CharString const & smFileName,
                                TMatches & stQueryMatches,
                                unsigned numThreads)
{
    // Open file with Stellar matches
    std::fstream inStreamMatches(toCString(smFileName), std::ios::in | std::ios::binary);
    if (!inStreamMatches.good())
    {
        std::cerr << "Could not open Stellar file " << smFileName << std::endl;
        return 1;
    }
    // Read local matches in GFF Stellar format.
    seqan::DirectionIterator<std::fstream, seqan::Input>::Type iter(inStreamMatches);
    LocalMatchStore<> lmStore;
    unsigned i = 0;
    while (!atEnd(iter))
    {
        try
        {
            readRecord(lmStore, iter, StellarGff());
        }
        catch (seqan::IOError const & ioErr)
        {
            std::cerr << "Invalid Stellar GFF record #" << i << "(" << ioErr.what() << ")\n";
        }
        i += 1;
    }
    // Creating Stellar Matches from input
    resize(stQueryMatches, length(queries));
    
    if (!_createStellarMatches(queries, sQueryIds, databases, databaseIDs, lmStore, stQueryMatches, numThreads))
        return 1;

    return 0;
}

#endif // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_CREATE_STELLARMATCHES_FROM_FILE_H_
