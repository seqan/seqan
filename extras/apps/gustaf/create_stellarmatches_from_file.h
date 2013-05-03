// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011, Knut Reinert, FU Berlin
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

#ifndef SEQAN_EXTRAS_APPS_GUSTAF_CREATE_STELLARMATCHES_FROM_FILE_H_
#define SEQAN_EXTRAS_APPS_GUSTAF_CREATE_STELLARMATCHES_FROM_FILE_H_

#include <iostream>
#include <fstream>
#include <seqan/file.h>
#include <seqan/parse_lm.h>
#include "../../../core/apps/stellar/stellar.h"

using namespace seqan;

// Creates a short Id out of a long one (i.e. it takes the prefix til the first white space)
template <typename TId>
void _getShortId(TId & shortId, TId const & longId)
{
    clear(shortId);
    for (typename Position<TId>::Type i = 0; i < length(longId) && isgraph(value(longId, i)); ++i)
    {
        appendValue(shortId, value(longId, i));
    }
}

// Gets a Set of Ids and creates a set of short Ids from the longer ones
template <typename TId>
void _getShortIds(StringSet<TId> & shortIds, StringSet<TId> & longIds)
{
    for (unsigned index = 0; index < length(longIds); ++index)
    {
        TId & longId = longIds[index];
        TId sId;
        _getShortId(sId, longId);
        // std::cerr << sId << std::endl;
        appendValue(shortIds, sId);
    }
}

// Takes the values from the localMatchStore and creates a Stellar match out of them.
// The Stellar match is, sorted by read, appended to stQueryMatches. Reads are identified by their short query IDs
// (sQueryIds), the coorect database by their short database IDs (sDBIds), which are both given in the GFF file
template <typename TSequence, typename TId>
bool _createStellarMatches(StringSet<TSequence> & queries,
                           StringSet<TId> const & sQueryIds,
                           StringSet<TSequence> & databases,
                           StringSet<TId> const & sDBIds,
                           StringSet<TId> const & databaseIds,
                           LocalMatchStore<> & lmStore,
                           StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & stQueryMatches)
{
    typedef typename Infix<TSequence>::Type TInfix;
    typedef Segment<TInfix, InfixSegment> TSegment;
    typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;
    typedef typename Row<TAlign>::Type TRow;
    typedef typename LocalMatchStore<>::TPosition TPosition;

    for (unsigned i = 0; i < length(lmStore.matchStore); ++i)
    {
        TInfix dbInf;
        TInfix queryInf;

        // Take match query ID and find the right index in queries
        unsigned iDB = maxValue<unsigned>();     // position of database sequence in databases
        unsigned iQuery = maxValue<unsigned>();  // position of query/read sequence in queries

        // Takes the short chromosome Id (from the Stellar match file) and looks up the corresponding long chromosome
        // Id entry from the reference input file
        for (unsigned j = 0; j < length(sDBIds); ++j)
        {
            if (lmStore.sequenceNameStore[lmStore.matchStore[i].subjectId] == sDBIds[j])
            {
                iDB = j;
                // std::cerr << "iDB" <<  sDBIds[j] << " " << lmStore.sequenceNameStore[lmStore.matchStore[i].subjectId] << std::endl;
                break;
            }
        }

        // Takes the short read Id (from the Stellar match file) and looks up the corresponding long read Id entry
        // from the read input file
        for (unsigned j = 0; j < length(sQueryIds); ++j)
        {
            if (lmStore.sequenceNameStore[lmStore.matchStore[i].queryId] == sQueryIds[j])
            {
                iQuery = j;
                // std::cerr << "iQuery" <<  sQueryIds[j] << " " << lmStore.sequenceNameStore[lmStore.matchStore[i].queryId] << std::endl;
                break;
            }
        }
        // Sanity check for read and query Id:
        // skips entry if no corresponding entry in the input file could not be found, else creates StellarMatch object
        if (iDB == maxValue<unsigned>() || iQuery == maxValue<unsigned>())
        {
            std::cerr << "Read or database does not exist for match: " << i
                      << " subjectId: " << lmStore.sequenceNameStore[lmStore.matchStore[i].subjectId]
                      << " queryId: " << lmStore.sequenceNameStore[lmStore.matchStore[i].queryId] << std::endl;
            std::cerr << "Skipping entry" << std::endl;
            continue;
        }

        // Sanity checks for database and read positions in case a wrong database or read file have been used
        // Checking for valid begin and end positions within identified database
        SEQAN_ASSERT_LEQ_MSG(lmStore.matchStore[i].subjectBeginPos, length(databases[iDB]),
                             "Match begin position exceeds database length! Wrong genome?");
        SEQAN_ASSERT_LEQ_MSG(lmStore.matchStore[i].subjectEndPos, length(databases[iDB]),
                             "Match end position exceeds database length! Wrong genome?");
        // Checking for valid begin and end positions within identified query
        SEQAN_ASSERT_LEQ_MSG(lmStore.matchStore[i].queryBeginPos, length(queries[iQuery]),
                             "Match begin position exceeds query length! Wrong read/contig sequence?");
        SEQAN_ASSERT_LEQ_MSG(lmStore.matchStore[i].queryEndPos, length(queries[iQuery]),
                             "Match end position exceeds query length! Wrong read/contig sequence?");

        if ((lmStore.matchStore[i].subjectBeginPos  > length(databases[iDB]))
           || (lmStore.matchStore[i].subjectEndPos > length(databases[iDB])))
        {
            std::cerr << "Match begin or end position exceeds database length! Wrong genome?" << std::endl;
            return 1;
        }
        if ((lmStore.matchStore[i].queryBeginPos > length(queries[iQuery]))
           || (lmStore.matchStore[i].queryEndPos > length(queries[iQuery])))
        {
            std::cerr << "Match begin or end position exceeds query length! Wrong read/contig sequence?" << std::endl;
            std::cerr <<  lmStore.matchStore[i].queryBeginPos << " " << lmStore.matchStore[i].queryEndPos << " " <<
            length(queries[iQuery]) << " " << sQueryIds[iQuery] << std::endl;
            return 1;
        }

        // Checking orientation and swapping positions for reverse matches to apply them to stellar format
        bool orientation = true;
        if (lmStore.matchStore[i].subjectBeginPos > lmStore.matchStore[i].subjectEndPos)
        {
            orientation = false;
            TPosition tmp = lmStore.matchStore[i].subjectBeginPos;
            lmStore.matchStore[i].subjectBeginPos = lmStore.matchStore[i].subjectEndPos;
            lmStore.matchStore[i].subjectEndPos = tmp;
        }
        // Computing infices for alignment rows for stellar matches
        queryInf = infix(queries[iQuery],
                         lmStore.matchStore[i].queryBeginPos,
                         lmStore.matchStore[i].queryEndPos);
        dbInf = infix(databases[iDB],
                      lmStore.matchStore[i].subjectBeginPos,
                      lmStore.matchStore[i].subjectEndPos);

        // Creating align object for stellar format
        TAlign localAlign;
        resize(rows(localAlign), 2);
        setSource(row(localAlign, 0), host(dbInf));
        setSource(row(localAlign, 1), host(queryInf));
        TRow & row1 = row(localAlign, 0);
        TRow & row2 = row(localAlign, 1);


        // Set begin and end positions of align
        setBeginPosition(row2, lmStore.matchStore[i].queryBeginPos);
        setBeginPosition(row1, lmStore.matchStore[i].subjectBeginPos);

        // setBeginPosition(row1, 0);
        // setBeginPosition(row2, 0);
        setEndPosition(row1, lmStore.matchStore[i].subjectEndPos);
        setEndPosition(row2, lmStore.matchStore[i].queryEndPos);

        unsigned gapIndex = 0;
        // Inserting gaps into rows according to cigar line
        if (length(lmStore.cigarStore) > lmStore.matchStore[i].id)
        {
            String<CigarElement<> > const & cigar = lmStore.cigarStore[lmStore.matchStore[i].id];
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
                    gapIndex += cigar[j].count;
            }
        }
        // Create Stellar match and append it to stQueryMatches
        StellarMatch<TSequence, TId> match(localAlign, databaseIds[iDB], orientation);
        appendValue(stQueryMatches[iQuery].matches, match);
    }
    return 0;
}

// Reads in a file with Stellar matches in gff format and creates StellarMatch object from the entries
template <typename TSequence, typename TId, typename TMatches>
bool _getStellarMatchesFromFile(StringSet<TSequence> & queries,
                                StringSet<TId> & queryIDs,
                                StringSet<TSequence> & databases,
                                StringSet<TId> & databaseIDs,
                                CharString const & smFileName,
                                TMatches & stQueryMatches)
{
    StringSet<TId> sQueryIds; // Allocator for short query Ids, needed bc Stellar only prints these short Ids to file
    StringSet<TId> sDBIds;    // Allocator for short db Ids
    _getShortIds(sQueryIds, queryIDs);
    _getShortIds(sDBIds, databaseIDs);

    // Open file with Stellar matches
    std::fstream inStreamMatches(toCString(smFileName), std::ios::in | std::ios::binary);
    if (!inStreamMatches.good())
    {
        std::cerr << "Could not open Stellar file " << smFileName << std::endl;
        return 1;
    }
    // Read local matches in GFF Stellar format.
    RecordReader<std::fstream, SinglePass<> > recordReader(inStreamMatches);
    LocalMatchStore<> lmStore;
    unsigned i = 0;
    while (!atEnd(recordReader))
    {
        int res = readRecord(lmStore, recordReader, StellarGff());
        if (res != 0)
            std::cerr << "Invalid Stellar GFF record #" << i << '\n';
        i += 1;
    }
    // Creating Stellar Matches from input
    resize(stQueryMatches, length(queries));

    // Infices/segments of db and query sequence, needed for stellar match
    typedef typename Infix<TSequence>::Type TInfix;
    typedef Segment<TInfix, InfixSegment> TSegment;

    if (!_createStellarMatches(queries, sQueryIds, databases, sDBIds, databaseIDs, lmStore, stQueryMatches))
        return 1;

    return 0;
}

#endif // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_CREATE_STELLARMATCHES_FROM_FILE_H_