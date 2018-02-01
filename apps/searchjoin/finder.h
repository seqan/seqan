// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the DbFinder class.
// ==========================================================================

#ifndef SANDBOX_ESIRAGUSA_APPS_SEARCHJOIN_DBFINDER_H_
#define SANDBOX_ESIRAGUSA_APPS_SEARCHJOIN_DBFINDER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/parallel.h>

#include "verifier.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Online
// ----------------------------------------------------------------------------

struct Online_;
typedef Tag<Online_>    Online;

// ----------------------------------------------------------------------------
// Tag Top
// ----------------------------------------------------------------------------

struct Top_;
typedef Tag<Top_>       Top;

// ----------------------------------------------------------------------------
// Tag Bottom
// ----------------------------------------------------------------------------

struct Bottom_;
typedef Tag<Bottom_>    Bottom;

// ----------------------------------------------------------------------------
// Class DbFinder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSpec = void>
struct DbFinder
{
    typedef Db<TText>                                               TDb;
    typedef Db<TText, TDbQuerySpec>                                 TDbQuery;
    typedef DbIndex<TIndex>                                         TDbIndex;
    typedef DbIndex<TIndex, Query>                                  TDbQueryIndex;
    typedef Verifier<TText, TDbQuerySpec>                           TVerifier;
    typedef typename Size<TText>::Type                              TTextSize;

    TDb /* const */     & db;
    TDbIndex            dbIndex;
    Holder<TDbQuery>    query;
    TDbQueryIndex       queryIndex;
    TDelegate           & delegate;
    TVerifier           verifier;
    TTextSize           minSeedLength;

    DbFinder(TDb /* const */ & db, TDelegate & delegate) :
        db(db),
        dbIndex(),
        delegate(delegate),
        verifier(db),
        minSeedLength(0)
    {}

    template <typename TFinder>
    void operator()(TFinder const & finder)
    {
        onFind(*this, finder);
    }
};

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
struct DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Parallel>
{
    typedef Db<TText>                                           TDb;
    typedef Db<TText, TDbQuerySpec>                             TDbQuery;
    typedef DbIndex<TIndex>                                     TDbIndex;
    typedef DbIndex<TIndex, Query>                              TDbQueryIndex;
    typedef Verifier<TText, TDbQuerySpec>                       TVerifier;
    typedef typename Size<TText>::Type                          TTextSize;
    typedef typename Size<TIndex>::Type                         TDepth;

    typedef Backtracking<HammingDistance, Bottom>               TBacktrackingExt;
    typedef Backtracking<EditDistance, Bottom>                  TBacktrackingApx;
    typedef Finder_<TIndex, TIndex, TBacktrackingExt>           TFinderExt;
    typedef Finder_<TIndex, TIndex, TBacktrackingApx>           TFinderApx;
    typedef String<TFinderApx>                                  TFindersApx;
    typedef String<TFinderExt>                                  TFindersExt;

    TDb /* const */     & db;
    TDbIndex            dbIndex;
    Holder<TDbQuery>    query;
    TDbQueryIndex       queryIndex;
    TDelegate           & delegate;
    TVerifier           verifier;
    TTextSize           minSeedLength;

    TDepth              parallelDepth;
    TFindersExt         findersExt;
    TFindersApx         findersApx;

    DbFinder(TDb /* const */ & db, TDelegate & delegate) :
        db(db),
        dbIndex(),
        delegate(delegate),
        verifier(db),
        minSeedLength(0),
        parallelDepth(2)
    {}

    template <typename TFinder>
    void operator()(TFinder const & finder)
    {
        onFind(*this, finder);
    }
};

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
struct DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Online>
{
    typedef Db<TText>                                               TDb;
    typedef Db<TText, TDbQuerySpec>                                 TDbQuery;
    typedef Verifier<TText, TDbQuerySpec>                           TVerifier;
    typedef typename Size<TText>::Type                              TTextSize;

    TDb /* const */     & db;
    Holder<TDbQuery>    query;
    TDelegate           & delegate;
    TVerifier           verifier;

    DbFinder(TDb /* const */ & db, TDelegate & delegate) :
        db(db),
        delegate(delegate),
        verifier(db)
    {}
};

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
struct DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Exact>
{
    typedef Db<TText>                                               TDb;
    typedef Db<TText, TDbQuerySpec>                                 TDbQuery;
    typedef DbIndex<TIndex>                                         TDbIndex;

    TDb /* const */     & db;
    TDbIndex            dbIndex;
    TDelegate           & delegate;

    DbFinder(TDb /* const */ & db, TDelegate & delegate) :
        db(db),
        delegate(delegate)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction PatternIterator_
// ----------------------------------------------------------------------------

namespace seqan {

template <typename TDistance, typename TSpec>
struct PatternIterator_<TDbDnaSaSmall, TDbDnaSaSmall, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<TDbDnaSaSmall, TopDown<Truncated<Preorder> > >::Type  Type;
};

template <typename TDistance, typename TSpec>
struct PatternIterator_<TDbGeoSaSmall, TDbGeoSaSmall, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<TDbGeoSaSmall, TopDown<Truncated<Preorder> > >::Type  Type;
};

template <typename TDistance, typename TSpec>
struct PatternIterator_<TDbDnaSaHuge, TDbDnaSaHuge, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<TDbDnaSaHuge, TopDown<Truncated<Preorder> > >::Type  Type;
};

template <typename TDistance, typename TSpec>
struct PatternIterator_<TDbGeoSaHuge, TDbGeoSaHuge, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<TDbGeoSaHuge, TopDown<Truncated<Preorder> > >::Type  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _inTerminalState()                          [Finder_<Backtracking>]
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_inTerminalState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, Top> > const & finder,
                 StageFinal_ const & /* tag */)
{
    return _inActiveState(finder, StageFinal_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_inTerminalState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, Top> > const & finder,
                 StageLower_ const & /* tag */)
{
    return _inActiveState(finder, StageLower_());
}

// ----------------------------------------------------------------------------
// Function _moveToNextStage()                          [Finder_<Backtracking>]
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, Top> > const & /* finder */,
                 StageLower_ const & /* tag */)
{
    return false;
}

// ----------------------------------------------------------------------------
// Function _stayInCurrentStage()                       [Finder_<Backtracking>]
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, Top> > const & finder,
                    StageInitial_ const & /* tag */)
{
    return !_moveToNextStage(finder, StageInitial_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, Top> > const & finder,
                    StageExact_ const & /* tag */)
{
    return !_moveToNextStage(finder, StageExact_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, Top> > const & finder,
                    StageInitial_ const & /* tag */)
{
    return !_moveToNextStage(finder, StageInitial_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, Top> > const & finder,
                    StageUpper_ const & /* tag */)
{
    return !_moveToNextStage(finder, StageUpper_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, Top> > const & finder,
                    StageDiagonal_ const & /* tag */)
{
    return !_moveToNextStage(finder, StageDiagonal_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, Top> > const & /* finder */,
                    StageLower_ const & /* tag */)
{
    return true;
}

} // namespace seqan

// ----------------------------------------------------------------------------
// Function setMinSeedLength()                                       [DbFinder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSpec, typename TSeedLength>
void setMinSeedLength(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, TSpec> & dbFinder, TSeedLength minSeedLength)
{
    dbFinder.minSeedLength = minSeedLength;
}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSeedLength>
void setMinSeedLength(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Online> & /* dbFinder */, TSeedLength /* minSeedLength */)
{}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSeedLength>
void setMinSeedLength(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Exact> & /* dbFinder */, TSeedLength /* minSeedLength */)
{}

// ----------------------------------------------------------------------------
// Function index()                                                  [DbFinder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSpec>
void index(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, TSpec> & dbFinder)
{
    // Build database index.
    build(dbFinder.dbIndex, dbFinder.db, Default());

    // Free some memory if we don't need a q-gram directory
    clear(dbFinder.dbIndex.dir);
    shrinkToFit(dbFinder.dbIndex.dir);
}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
void index(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Exact> & dbFinder)
{
    // Build database index.
    build(dbFinder.dbIndex, dbFinder.db, seqan::Exact());
}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
void index(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Online> & /* dbFinder */)
{}

// ----------------------------------------------------------------------------
// Function prepare()                                                [DbFinder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSpec>
void prepare(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, TSpec> & dbFinder, Db<TText, TDbQuerySpec> & query)
{
    // Set query.
    setValue(dbFinder.query, query);
    setValue(dbFinder.verifier.query, query);

    // Build query index.
    buildQuery(dbFinder.queryIndex, value(dbFinder.query), dbFinder.minSeedLength);
}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
void prepare(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Online> & dbFinder, Db<TText, TDbQuerySpec> & query)
{
    // Set query.
    setValue(dbFinder.query, query);
    setValue(dbFinder.verifier.query, query);
}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
void prepare(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Exact> &, Db<TText, TDbQuerySpec> &)
{}

// ----------------------------------------------------------------------------
// Function onFind()                                                 [DbFinder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSpec, typename TBacktracking>
inline void
onFind(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, TSpec> & dbFinder,
       Finder_<TIndex, TIndex, TBacktracking> const & finder)
{
    typedef typename TextIterator_<TIndex, TIndex, TBacktracking>::Type     TTextIterator;
    typedef typename PatternIterator_<TIndex, TIndex, TBacktracking>::Type  TPatternIterator;
    typedef typename Size<TIndex>::Type                                     TSize;
    //typedef typename Size<TText>::Type                                    TErrors;
    typedef typename Fibre<TIndex, FibreSA>::Type const                     TSAFibre;
    typedef typename Infix<TSAFibre>::Type                                  TOccurrences;

    TTextIterator const & textIt = back(finder.textStack);
    TPatternIterator const & patternIt = back(finder.patternStack);

    TOccurrences textOccurrences = getOccurrences(textIt);
    TOccurrences patternOccurrences = getEmptyEdges(patternIt);

    TSize textOccurrencesCount = length(textOccurrences);
    TSize patternOccurrencesCount = length(patternOccurrences);

    for (TSize i = 0; i < textOccurrencesCount; ++i)
        for (TSize j = 0; j < patternOccurrencesCount; ++j)
        {
            TSize dbId1 = getSeqNo(textOccurrences[i]);
            TSize dbId2 = getSeqNo(patternOccurrences[j]);
            TSize pos1 = getSeqOffset(textOccurrences[i]);
            TSize pos2 = getSeqOffset(patternOccurrences[j]);

            if (_validHit(dbFinder, dbId1, dbId2, pos1, pos2))
                dbFinder.verifier(dbId1, dbId2, dbFinder.delegate);
        }
}

// ----------------------------------------------------------------------------
// Function _validHit()                                              [DbFinder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSpec>
inline bool
_validHit(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, TSpec> const & dbFinder,
          typename Size<TIndex>::Type dbId1,
          typename Size<TIndex>::Type dbId2,
          typename Size<TIndex>::Type pos1,
          typename Size<TIndex>::Type pos2)
{
    if (dbId1 >= dbId2)
        return false;

    int maxErrors = getErrors(value(dbFinder.query), dbId2);
    int dbLen1 = length(dbFinder.db.text[dbId1]);
    int dbLen2 = length(dbFinder.db.text[dbId2]);

    // Check maximal length difference.
    if (dbLen1 + maxErrors < dbLen2 || dbLen2 + maxErrors < dbLen1)
        return false;

    // Check band according to Lemma 1.
    int diag = pos1 - pos2;
    int leftDiag = -((dbLen2 - dbLen1 + maxErrors) / 2);
    int rightDiag =  (dbLen1 - dbLen2 + maxErrors) / 2;

    return (leftDiag <= diag && diag <= rightDiag);
}

template <typename TText, typename TIndex, typename TDelegate, typename TSpec>
inline bool
_validHit(DbFinder<TText, TIndex, Query, TDelegate, TSpec> const & dbFinder,
          typename Size<TIndex>::Type dbId1,
          typename Size<TIndex>::Type dbId2,
          typename Size<TIndex>::Type pos1,
          typename Size<TIndex>::Type pos2)
{
    int maxErrors = getErrors(value(dbFinder.query), dbId2);
    int dbLen1 = length(dbFinder.db.text[dbId1]);
    int dbLen2 = length(value(dbFinder.query).text[dbId2]);

    // Check maximal length difference.
    if (dbLen1 + maxErrors < dbLen2 || dbLen2 + maxErrors < dbLen1)
        return false;

    // Check band according to Lemma 1.
    int diag = pos1 - pos2;
    int leftDiag = -((dbLen2 - dbLen1 + maxErrors) / 2);
    int rightDiag =  (dbLen1 - dbLen2 + maxErrors) / 2;

    if (!(leftDiag <= diag && diag <= rightDiag))
        return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function onFind()                                       [DbFinder<Parallel>]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
inline void
onFind(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Parallel> & dbFinder,
       Finder_<TIndex, TIndex, Backtracking<EditDistance, Top> > const & finder)
{
    typedef Backtracking<EditDistance, Bottom>                          TBacktracking;
    typedef Finder_<TIndex, TIndex, TBacktracking>                      TFinder;

    TFinder finderBottom;
    _setScoreThreshold(finderBottom, finder.maxScore);
    appendValue(finderBottom.textStack, back(finder.textStack));
    appendValue(finderBottom.patternStack, back(finder.patternStack));
    appendValue(finderBottom.scoreStack, back(finder.scoreStack));
    back(finderBottom.patternStack).depth = dbFinder.queryIndex.seedLength;

    appendValue(dbFinder.findersApx, finderBottom);
}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
inline void
onFind(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Parallel> & dbFinder,
       Finder_<TIndex, TIndex, Backtracking<HammingDistance, Top> > const & finder)
{
    typedef Backtracking<HammingDistance, Bottom>                       TBacktracking;
    typedef Finder_<TIndex, TIndex, TBacktracking>                      TFinder;

    TFinder finderBottom;
    _setScoreThreshold(finderBottom, finder.maxScore);
    appendValue(finderBottom.textStack, back(finder.textStack));
    appendValue(finderBottom.patternStack, back(finder.patternStack));
    appendValue(finderBottom.scoreStack, back(finder.scoreStack));
    back(finderBottom.patternStack).depth = dbFinder.queryIndex.seedLength;

    appendValue(dbFinder.findersExt, finderBottom);
}

// ----------------------------------------------------------------------------
// Function execute()                                                [DbFinder]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate, typename TSpec>
void execute(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, TSpec> & dbFinder)
{
    typedef Backtracking<EditDistance>                                          TBacktrackingApx;
    typedef Backtracking<HammingDistance>                                       TBacktrackingExt;
    typedef Finder_<TIndex, TIndex, TBacktrackingApx>                           TFinderApx;
    typedef Finder_<TIndex, TIndex, TBacktrackingExt>                           TFinderExt;
    typedef typename TextIterator_<TIndex, TIndex, TBacktrackingApx>::Type      TTextIterator;
    typedef typename PatternIterator_<TIndex, TIndex, TBacktrackingApx>::Type   TPatternIterator;

    // Instantiate a finder.
    TFinderApx finderApx;
    TFinderExt finderExt;

    unsigned seedSetsCount = length(dbFinder.queryIndex.errors);
    unsigned seedSet = 0;

    // Find exactly query index in database index.
    if (seedSetsCount > 0 && dbFinder.queryIndex.errors[seedSet] == 0)
    {
        TTextIterator textIt(dbFinder.dbIndex.index);
        TPatternIterator patternIt(dbFinder.queryIndex.index[seedSet]);
        patternIt.depth = dbFinder.queryIndex.seedLength;

        _setScoreThreshold(finderExt, 0);
        _initState(finderExt, textIt, patternIt);
        _find(finderExt, dbFinder, StageInitial_());
        clear(finderExt);

        seedSet++;
    }

    // Find approximately query index in database index.
    for (; seedSet < seedSetsCount; ++seedSet)
    {
        unsigned seedErrors = dbFinder.queryIndex.errors[seedSet];

        TTextIterator textIt(dbFinder.dbIndex.index);
        TPatternIterator patternIt(dbFinder.queryIndex.index[seedSet]);
        patternIt.depth = dbFinder.queryIndex.seedLength;

        _setScoreThreshold(finderApx, seedErrors);
        _initState(finderApx, textIt, patternIt);
        _find(finderApx, dbFinder, StageInitial_());
        clear(finderApx);
    }
}

// ----------------------------------------------------------------------------
// Function execute()                                      [DbFinder<Parallel>]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
void execute(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Parallel> & dbFinder)
{
    typedef Backtracking<EditDistance, Top>                                     TBacktrackingApx;
    typedef Backtracking<HammingDistance, Top>                                  TBacktrackingExt;
    typedef Backtracking<EditDistance, Bottom>                                  TBacktrackingApxBottom;
    typedef Backtracking<HammingDistance, Bottom>                               TBacktrackingExtBottom;
    typedef Finder_<TIndex, TIndex, TBacktrackingApx>                           TFinderApx;
    typedef Finder_<TIndex, TIndex, TBacktrackingExt>                           TFinderExt;
    typedef Finder_<TIndex, TIndex, TBacktrackingApxBottom>                     TFinderApxBottom;
    typedef Finder_<TIndex, TIndex, TBacktrackingExtBottom>                     TFinderExtBottom;
    typedef typename TextIterator_<TIndex, TIndex, TBacktrackingApx>::Type      TTextIterator;
    typedef typename PatternIterator_<TIndex, TIndex, TBacktrackingApx>::Type   TPatternIterator;

    // Instantiate a finder.
    TFinderApx finderApx;
    TFinderExt finderExt;

    unsigned seedSetsCount = length(dbFinder.queryIndex.errors);
    unsigned seedSet = 0;

    // Find exactly query index in database index.
    if (seedSetsCount > 0 && dbFinder.queryIndex.errors[seedSet] == 0)
    {
        TTextIterator textIt(dbFinder.dbIndex.index);
        TPatternIterator patternIt(dbFinder.queryIndex.index[seedSet]);
        patternIt.depth = _min(dbFinder.queryIndex.seedLength, dbFinder.parallelDepth);
        
        _setScoreThreshold(finderExt, 0);
        _initState(finderExt, textIt, patternIt);
        _find(finderExt, dbFinder, StageInitial_());
        clear(finderExt);

        seedSet++;
    }

    // Find approximately query index in database index.
    for (; seedSet < seedSetsCount; ++seedSet)
    {
        unsigned seedErrors = dbFinder.queryIndex.errors[seedSet];

        TTextIterator textIt(dbFinder.dbIndex.index);
        TPatternIterator patternIt(dbFinder.queryIndex.index[seedSet]);
        patternIt.depth = _min(dbFinder.queryIndex.seedLength, dbFinder.parallelDepth);

        _setScoreThreshold(finderApx, seedErrors);
        _initState(finderApx, textIt, patternIt);
        _find(finderApx, dbFinder, StageInitial_());
        clear(finderApx);
    }

    unsigned findersExtCount = length(dbFinder.findersExt);
    unsigned findersApxCount = length(dbFinder.findersApx);

    std::cout << "Ext jobs count:\t\t\t\t" << findersExtCount << std::endl;
    std::cout << "Apx jobs count:\t\t\t\t" << findersApxCount << std::endl;

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (int finderExt = 0; finderExt < (int)findersExtCount; ++finderExt)
    {
        TFinderExtBottom & finderBottom = dbFinder.findersExt[finderExt];
        _find(finderBottom, dbFinder, StageInitial_());
        clear(finderBottom);
    }

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (int finderApx = 0; finderApx < (int)findersApxCount; ++finderApx)
    {
        TFinderApxBottom & finderBottom = dbFinder.findersApx[finderApx];
        _find(finderBottom, dbFinder, StageInitial_());
        clear(finderBottom);
    }
}

// ----------------------------------------------------------------------------
// Function execute()                                         [DbFinder<Exact>]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
void execute(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Exact> & dbFinder)
{
    typedef typename Fibre<TIndex, FibreSA>::Type const      TSA;
    typedef typename Size<TIndex>::Type                      TSize;
    typedef typename Value<TText const>::Type                TString;
    typedef typename Iterator<TString const, Standard>::Type TStringIter;

    TText const &text = dbFinder.db.text;
    TSA const &sa = indexSA(dbFinder.dbIndex.index);
    String<TSize> &dir = dbFinder.dbIndex.dir;
    TSize sortedPrefix = dbFinder.dbIndex.sortedPrefix;

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic,1))
	for (int i = 1; i < (int)length(dir); ++i)
	{
        size_t bktBegin = dir[i - 1];
        size_t bktEnd = dir[i];

        // We need at least 2 suffixes.
		if (bktBegin + 2 > bktEnd)
            continue;

        TSize runBegin = bktBegin;
        TSize seqNo1 = getSeqNo(sa[runBegin]);
        for (size_t i = bktBegin + 1; i < bktEnd; ++i)
        {
            TSize seqNo2 = getSeqNo(sa[i]);
            SEQAN_ASSERT_EQ(getSeqOffset(sa[i]), 0u);
            if (length(text[seqNo1]) == length(text[seqNo2]))
            {
                TStringIter it1Beg = begin(text[seqNo1], Standard()) + sortedPrefix - 1;
                TStringIter it1End = end(text[seqNo1], Standard()) - 1;
                TStringIter it2End = end(text[seqNo2], Standard()) - 1;

                for (; it1End != it1Beg; --it1End, --it2End)
                {
                    if (*it1End != *it2End)
                        break;
                }

                if (it1End == it1Beg)
                    continue;
            }

            if (runBegin + 1 < i)
            {
                for (unsigned j = runBegin; j < i; ++j)
                    for (unsigned k = j + 1; k < i; ++k)
                        dbFinder.delegate(getSeqNo(sa[j]), getSeqNo(sa[k]));
            }
            runBegin = i;
            seqNo1 = getSeqNo(sa[runBegin]);
        }

        if (runBegin + 1 < bktEnd)
        {
            for (unsigned j = runBegin; j < bktEnd; ++j)
                for (unsigned k = j + 1; k < bktEnd; ++k)
                    dbFinder.delegate(getSeqNo(sa[j]), getSeqNo(sa[k]));
        }
    }
}

// NOTE(esiragusa): This is computing exact join using a dfs traversal.
//namespace seqan
//{
//    struct SAInfix;
//
//    template <typename TText>
//    struct Fibre<Index<TText, IndexSa<SAInfix> >, FibreSA>
//    {
//        typedef Segment<typename Fibre<Index<TText, IndexSa<> >, FibreSA>::Type const, InfixSegment> Type;
//    };
//
//}
//
//template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
//void execute(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Exact> & dbFinder)
//{
//    typedef Index<TText, IndexSa<SAInfix> >                                     TSAInfixIndex;
//    typedef typename Iterator<TSAInfixIndex, TopDown<ParentLinks<> > >::Type    TIter;
//
//    TText const &text = dbFinder.db.text;
//    TSA const &sa = indexSA(dbFinder.dbIndex.index);
//    String<TSize> &dir = dbFinder.dbIndex.dir;
//    TSize sortedPrefix = dbFinder.dbIndex.sortedPrefix;
//
//    TSAInfixIndex subIndex(text);
//    indexSA(subIndex) = infix(sa, bktBegin, bktEnd);
//
//    TIter iter(subIndex);
//    value(iter).range.i1 = 0;
//    _setSizeInval(value(iter).range.i2);
//    value(iter).parentRight = value(iter).range.i2;
//    value(iter).repLen = sortedPrefix;
//    value(iter).lastChar = suffix(text, sa[bktBegin])[sortedPrefix - 1];
//
//    while (!atEnd(iter))
//    {
//        if (isRightTerminal(iter))
//        {
//            typename Infix<TSA>::Type emptyEdges = getEmptyEdges(iter);
//            TSize len = length(emptyEdges);
//            if (len < 2)
//                continue;
//            std::cout << len << '\t' << representative(iter) << std::endl;
//            for (unsigned i = 0; i < len; ++i)
//                for (unsigned j = i + 1; j < len; ++j)
//                    dbFinder.delegate(getSeqNo(sa[j]), getSeqNo(sa[k]));
//        }
//        goNext(iter);
//    }
//}

// ----------------------------------------------------------------------------
// Function execute()                                        [DbFinder<Online>]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
void execute(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Online> & dbFinder)
{
    typedef Db<TText>                           TDb;
    typedef typename Size<TDb>::Type            TDbSize;
    typedef typename MakeSigned<TDbSize>::Type  TDbSSize;

    TDbSSize dbSize = (TDbSSize)length(dbFinder.db.text) - 1;
    TDbSSize dbQuerySize = length(value(dbFinder.query).text);

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (TDbSSize dbId1 = 0; dbId1 < dbSize; ++dbId1)
        for (TDbSSize dbId2 = dbId1 + 1; dbId2 < dbQuerySize; ++dbId2)
            dbFinder.verifier(dbId1, dbId2, dbFinder.delegate);
}

template <typename TText, typename TIndex, typename TDelegate>
void execute(DbFinder<TText, TIndex, Query, TDelegate, Online> & dbFinder)
{
    typedef Db<TText>                               TDb;
    typedef Db<TText, Query>                        TDbQuery;
    typedef typename Size<TDb>::Type                TDbSize;
    typedef typename MakeSigned<TDbSize>::Type      TDbSSize;
    typedef typename Size<TDbQuery>::Type           TDbQuerySize;
    typedef typename MakeSigned<TDbQuerySize>::Type TDbQuerySSize;

    TDbSSize dbSize = length(dbFinder.db.text);
    TDbQuerySSize dbQuerySize = length(value(dbFinder.query).text);

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (TDbSSize dbId = 0; dbId < dbSize; ++dbId)
        for (TDbQuerySSize queryId = 0; queryId < dbQuerySize; ++queryId)
            dbFinder.verifier(dbId, queryId, dbFinder.delegate);
}

#endif  // #ifndef SANDBOX_ESIRAGUSA_APPS_SEARCHJOIN_DBFINDER_H_

