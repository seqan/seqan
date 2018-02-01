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
// This file contains database type definitions.
// ==========================================================================

#ifndef SEQAN_APPS_SEARCHJOIN_DB_H_
#define SEQAN_APPS_SEARCHJOIN_DB_H_

#include <numeric>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

#define SEQAN_DEBUG

using namespace seqan;

// ============================================================================
// Database Type Definitions
// ============================================================================

// ----------------------------------------------------------------------------
// Database Basic Types
// ----------------------------------------------------------------------------

struct SmallDb_;
typedef Tag<SmallDb_>                           SmallDb;

struct HugeDb_;
typedef Tag<HugeDb_>                            HugeDb;

struct StringOfMaxSize256_;
typedef Tag<StringOfMaxSize256_>                StringOfMaxSize256;

typedef Owner<ConcatDirect<> >                  TDbSSetSpec;
typedef Alloc<StringOfMaxSize256>               TDbStringSpec;

typedef String<Dna5, TDbStringSpec>             TDbDnaString;
typedef StringSet<TDbDnaString, TDbSSetSpec>    TDbDna;

typedef String<char, TDbStringSpec>             TDbGeoString;
typedef StringSet<TDbGeoString, TDbSSetSpec>    TDbGeo;

// ----------------------------------------------------------------------------
// Database Basic Metafunctions
// ----------------------------------------------------------------------------

namespace seqan
{
//template <>
//struct Size<TDbDnaString>
//{
//    typedef unsigned char               Type;
//};
template <>
struct Size<TDbDna>
{
    typedef unsigned int                Type;
};
template <>
struct StringSetLimits<TDbDna>
{
    typedef String<Size<TDbDna>::Type>  Type;
};

//template <>
//struct Size<TDbGeoString>
//{
//    typedef unsigned char               Type;
//};
template <>
struct Size<TDbGeo>
{
    typedef unsigned int                Type;
};
template <>
struct StringSetLimits<TDbGeo>
{
    typedef String<Size<TDbGeo>::Type>  Type;
};
}

// ============================================================================
// Database Index Type Definitions
// ============================================================================

// ----------------------------------------------------------------------------
// Database Suffix Array Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TDbDna, IndexSa<SmallDb> >            TDbDnaSaSmall;
typedef Index<TDbGeo, IndexSa<SmallDb> >            TDbGeoSaSmall;

typedef Index<TDbDna, IndexSa<HugeDb> >             TDbDnaSaHuge;
typedef Index<TDbGeo, IndexSa<HugeDb> >             TDbGeoSaHuge;

namespace seqan
{
template <>
struct Fibre<TDbDnaSaSmall, FibreSA>
{
    typedef String<Pair<unsigned int, unsigned char, BitPacked<24, 8> >, DefaultIndexStringSpec<TDbDnaSaSmall>::Type> Type;
};

template <>
struct Fibre<TDbGeoSaSmall, FibreSA>
{
    typedef String<Pair<unsigned int, unsigned char, BitPacked<24, 8> >, DefaultIndexStringSpec<TDbGeoSaSmall>::Type> Type;
};

template <>
struct Fibre<TDbDnaSaHuge, FibreSA>
{
    typedef String<Pair<unsigned int, unsigned char, Pack>, DefaultIndexStringSpec<TDbDnaSaHuge>::Type> Type;
};

template <>
struct Fibre<TDbGeoSaHuge, FibreSA>
{
    typedef String<Pair<unsigned int, unsigned char, Pack>, DefaultIndexStringSpec<TDbDnaSaHuge>::Type> Type;
};
}

// ----------------------------------------------------------------------------
// Database Shape Length Definitions
// ----------------------------------------------------------------------------

//template <typename TText>
//struct ShapeLength
//{
//    static const unsigned VALUE = std::numeric_limits<typename Size<TText>::Type>::max();
//};
//
//template <>
//struct ShapeLength<TDbDna>
//{
//    static const unsigned VALUE = 10u;
//};
//
//template <>
//struct ShapeLength<TDbGeo>
//{
//    static const unsigned VALUE = 3u;
//};

// ----------------------------------------------------------------------------
// Database QGram Index with Bucket Refinement Type Definitions
// ----------------------------------------------------------------------------

//typedef UngappedShape<10>                               TShapeDna;
//typedef UngappedShape<3>                                TShapeGeo;
//
//typedef IndexQGram<TShapeDna, BucketRefinement>         TDbDnaQGramSpec;
//typedef IndexQGram<TShapeGeo, BucketRefinement>         TDbGeoQGramSpec;
//
//typedef Index<TDbDna, TDbDnaQGramSpec>                  TDbDnaQGram;
//typedef Index<TDbGeo, TDbGeoQGramSpec>                  TDbGeoQGram;

//namespace seqan
//{
//template <>
//struct Fibre<TDbDnaQGram, FibreDir>
//{
//    typedef String<unsigned int, DefaultIndexStringSpec<TDbDnaQGram>::Type>             Type;
//};
//
//template <>
//struct Fibre<TDbGeoQGram, FibreDir>
//{
//    typedef String<unsigned int, DefaultIndexStringSpec<TDbGeoQGram>::Type>             Type;
//};
//}

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Query
// ----------------------------------------------------------------------------

struct Query_;
typedef Tag<Query_> Query;

// ----------------------------------------------------------------------------
// Class DbParser
// ----------------------------------------------------------------------------

template <typename TDb, typename TSpec = void>
struct DbParser
{
    TDb         & db;
    CharString  _id;
    CharString  _text;

    EqualsChar<','> delim;
//    AssertFunctor<TFunctor, ParseError> asserter(functor);

    DbParser(TDb & db) :
        db(db)
    {}
};

template <typename TDb>
struct DbParser<TDb, Query>
{
    TDb         & db;
    CharString  _id;
    CharString  _text;
    CharString  _errors;

    EqualsChar<':'> delim;
    EqualsChar<','> delim2;

    DbParser(TDb & db) :
        db(db)
    {}
};

// ----------------------------------------------------------------------------
// Class Db
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec = void>
struct Db
{
    typedef Owner<ConcatDirect<> >                  TIdsSpec;
    typedef StringSet<CharString, TIdsSpec>         TIds;
    typedef typename Size<TText>::Type              TTextSize;
    typedef TTextSize                               TErrors;

    TText       text;
    TIds        ids;

    TTextSize   minLength;
    TTextSize   maxLength;
    TTextSize   avgLength;

    TErrors     errors;

    Db() :
        minLength(0),
        maxLength(0),
        avgLength(0),
        errors(0)
    {}
};

template <typename TText>
struct Db<TText, Query>
{
    typedef Owner<ConcatDirect<> >                  TIdsSpec;
    typedef StringSet<CharString, TIdsSpec>         TIds;
    typedef typename Size<TText>::Type              TTextSize;
    typedef TTextSize                               TErrors;

    TText       text;
    TIds        ids;

    TTextSize   minLength;
    TTextSize   maxLength;
    TTextSize   avgLength;

    String<TErrors> errors;

    TErrors     minErrors;
    TErrors     maxErrors;
    TErrors     avgErrors;

    Db() :
        minLength(0),
        maxLength(0),
        avgLength(0),
        minErrors(0),
        maxErrors(0),
        avgErrors(0)
    {}
};

// ----------------------------------------------------------------------------
// Class DbIndex
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = void>
struct DbIndex
{
    typedef typename Size<TIndex>::Type TSize;

    TIndex          index;
    String<TSize>   dir;
    TSize           sortedPrefix;
};

template <typename TIndex>
struct DbIndex<TIndex, Query>
{
    typedef typename Size<TIndex>::Type TSize;
    typedef String<TIndex>              TIndices;
    typedef String<TSize>               TErrors;

    TIndices    index;
    TErrors     errors;
    TSize       seedLength;

    DbIndex() :
        seedLength(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size<Db>::Type                                             [Db]
// ----------------------------------------------------------------------------

namespace seqan
{
template <typename TText, typename TSpec>
struct Size<Db<TText, TSpec> >
{
    typedef typename Size<TText>::Type         Type;
};

template <typename TIndex, typename TSpec>
struct Size<DbIndex<TIndex, TSpec> >
{
    typedef typename DbIndex<TIndex, TSpec>::TSize Type;
};
}

// ----------------------------------------------------------------------------
// Metafunction Host<Db>::Type                                             [Db]
// ----------------------------------------------------------------------------

namespace seqan
{
template <typename TText, typename TSpec>
struct Host<Db<TText, TSpec> >
{
    typedef TText                             Type;
};
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function parseLine()                                              [DbParser]
// ----------------------------------------------------------------------------

template <typename TDb, typename TSpec, typename TInputIt>
inline void
parseLine(DbParser<TDb, TSpec> & parser, TInputIt & inputIt)
{
    _clearBuffers(parser);
    _parseId(parser, inputIt);
    _parseText(parser, inputIt);
}

template <typename TDb, typename TInputIt>
inline void
parseLine(DbParser<TDb, Query> & parser, TInputIt & inputIt)
{
    _clearBuffers(parser);
    _parseId(parser, inputIt);
    _parseText(parser, inputIt);
    _parseErrors(parser, inputIt);
}

// ----------------------------------------------------------------------------
// Function _clearBuffers()                                          [DbParser]
// ----------------------------------------------------------------------------

template <typename TDb, typename TSpec>
inline void
_clearBuffers(DbParser<TDb, TSpec> & parser)
{
    clear(parser._id);
    clear(parser._text);
}

template <typename TDb>
inline void
_clearBuffers(DbParser<TDb, Query> & parser)
{
    clear(parser._id);
    clear(parser._text);
    clear(parser._errors);
}

// ----------------------------------------------------------------------------
// Function _parseId()                                               [DbParser]
// ----------------------------------------------------------------------------

template <typename TDb, typename TSpec, typename TInputIt>
inline void
_parseId(DbParser<TDb, TSpec> & parser, TInputIt & inputIt)
{
    // Read id.
    readUntil(parser._id, inputIt, parser.delim);

    // Add id to database.
    appendValue(parser.db.ids, parser._id);

    // Skip delim.
    skipOne(inputIt);
}

// ----------------------------------------------------------------------------
// Function _parseText()                                             [DbParser]
// ----------------------------------------------------------------------------

template <typename TDb, typename TSpec, typename TInputIt>
inline void
_parseText(DbParser<TDb, TSpec> & parser, TInputIt & inputIt)
{
    // Read text.
    readLine(parser._text, inputIt);

    // Add text to database.
    appendValue(parser.db.text, parser._text);
}

template <typename TDb, typename TInputIt>
inline void
_parseText(DbParser<TDb, Query> & parser, TInputIt & inputIt)
{
    // Read text.
    readUntil(parser._text, inputIt, parser.delim2);

    // Add text to database.
    appendValue(parser.db.text, parser._text);

    // Skip delim.
    skipOne(inputIt);
}

// ----------------------------------------------------------------------------
// Function _parseErrors()                                           [DbParser]
// ----------------------------------------------------------------------------

template <typename TDb, typename TInputIt>
inline void
_parseErrors(DbParser<TDb, Query> & parser, TInputIt & inputIt)
{
    // Read errors.
    readLine(parser._errors, inputIt);

    // Add errors to database.
    appendValue(parser.db.errors, std::atoi(toCString(parser._errors)));
}

// ----------------------------------------------------------------------------
// Function load()                                                         [Db]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TFileName>
bool load(Db<TText, TSpec> & db, TFileName const & fileName)
{
    typedef Db<TText, TSpec>                                        TDb;
    typedef DbParser<TDb, TSpec>                                    TDbParser;
    typedef typename Size<TText>::Type                              TTextSize;

    typedef std::ifstream                                           TInputStream;
    typedef typename DirectionIterator<TInputStream, Input>::Type   TInputIt;

    TInputStream inputFile;

    // Open file.
    if (!open(inputFile, toCString(fileName), OPEN_RDONLY))
        return false;

    // Instantiate an input iterator on the file.
    TInputIt inputIt = directionIterator(inputFile, Input());

    // Instantiate a parser.
    TDbParser parser(db);

    // Initialize min/max text length.
    db.minLength = std::numeric_limits<TTextSize>::max();
    db.maxLength = std::numeric_limits<TTextSize>::min();

    // Read the file.
    while (!atEnd(inputIt))
    {
        parseLine(parser, inputIt);

        // Update min/max text length.
        TTextSize textLength = length(back(db.text));
        db.minLength = _min(db.minLength, textLength);
        db.maxLength = _max(db.maxLength, textLength);
    }

    // Compute average text length.
    db.avgLength = length(db.text) ? length(db.text.concat) / length(db.text) : 0;

    _updateErrors(db);

    return true;
}

// ----------------------------------------------------------------------------
// Function _updateErrors()                                                [Db]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
void _updateErrors(Db<TText, TSpec> const & /* db */)
{}

template <typename TText>
void _updateErrors(Db<TText, Query> & db)
{
    if (length(db.text) == 0)
    {
        db.avgErrors = 0;
        db.minErrors = 0;
        db.maxErrors = 0;
    }
    else
    {
        db.avgErrors = std::accumulate(begin(db.errors, Standard()), end(db.errors, Standard()), 0) / length(db.text);
        db.minErrors = value(std::min_element(begin(db.errors, Standard()), end(db.errors, Standard())));
        db.maxErrors = value(std::max_element(begin(db.errors, Standard()), end(db.errors, Standard())));
    }
}

// ----------------------------------------------------------------------------
// Function split()                                                        [Db]
// ----------------------------------------------------------------------------

template <typename TText, typename TSeedLength>
void split(Db<TText, Query> & dbShort, Db<TText, Query> & dbLong, Db<TText, Query> /*const*/ & db, TSeedLength seedLength)
{
    typedef Db<TText, Query>                                  TDb;
    typedef typename Size<TDb>::Type                          TDbSize;
    typedef typename Value<TText>::Type                       TTextReference;
    typedef typename Size<TText>::Type                        TTextSize;

    TDbSize dbSize = length(db.text);

    // Initialize min/max text length.
    dbShort.minLength = std::numeric_limits<TTextSize>::max();
    dbShort.maxLength = std::numeric_limits<TTextSize>::min();
    dbLong.minLength = std::numeric_limits<TTextSize>::max();
    dbLong.maxLength = std::numeric_limits<TTextSize>::min();

    for (TDbSize dbId = 0; dbId < dbSize; ++dbId)
    {
        TTextReference text = db.text[dbId];
        TTextSize textLength = length(text);

        TDb *dbOut = (textLength < seedLength) ? &dbShort : &dbLong;

        (*dbOut).minLength = _min((*dbOut).minLength, textLength);
        (*dbOut).maxLength = _max((*dbOut).maxLength, textLength);

        appendValue((*dbOut).text, text);
        appendValue((*dbOut).ids, db.ids[dbId]);
        appendValue((*dbOut).errors, getErrors(db, dbId));
    }

    // Compute average text length.
    if (length(dbShort.text) == 0)
    {
        dbShort.minLength = 0;
        dbShort.maxLength = 0;
        dbShort.avgLength = 0;
    }
    else
    {
        dbShort.avgLength = length(dbShort.text.concat) / length(dbShort.text);
    }

    if (length(dbLong.text) == 0)
    {
        dbLong.minLength = 0;
        dbLong.maxLength = 0;
        dbLong.avgLength = 0;
    }
    else
    {
        dbLong.avgLength = length(dbLong.text.concat) / length(dbLong.text);
    }

    _updateErrors(dbShort);
    _updateErrors(dbLong);
}

// ----------------------------------------------------------------------------
// Function getErrors()                                                    [Db]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline typename Size<TText>::Type
getErrors(Db<TText, TSpec> const & db, typename Size<Db<TText, TSpec> >::Type /* dbId */)
{
    return db.errors;
}

template <typename TText>
inline typename Size<TText>::Type
getErrors(Db<TText, Query> const & db, typename Size<Db<TText, Query> >::Type dbId)
{
    return db.errors[dbId];
}

// ----------------------------------------------------------------------------
// Function getMinErrors()                                                 [Db]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline typename Size<TText>::Type
getMinErrors(Db<TText, TSpec> const & db)
{
    return db.errors;
}

template <typename TText>
inline typename Size<TText>::Type
getMinErrors(Db<TText, Query> const & db)
{
    return db.minErrors;
}

// ----------------------------------------------------------------------------
// Function getMaxErrors()                                                 [Db]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline typename Size<TText>::Type
getMaxErrors(Db<TText, TSpec> const & db)
{
    return db.errors;
}

template <typename TText>
inline typename Size<TText>::Type
getMaxErrors(Db<TText, Query> const & db)
{
    return db.maxErrors;
}

// ----------------------------------------------------------------------------
// Function getAvgErrors()                                                 [Db]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline typename Size<TText>::Type
getAvgErrors(Db<TText, TSpec> const & db)
{
    return db.errors;
}

template <typename TText>
inline typename Size<TText>::Type
getAvgErrors(Db<TText, Query> const & db)
{
    return db.avgErrors;
}

// ----------------------------------------------------------------------------
// Function countSeeds()                                                   [Db]
// ----------------------------------------------------------------------------

template <typename TCounts, typename TText, typename TSpec, typename TSeedLength>
TSeedLength countSeeds(TCounts & seedCounts,
                       Db<TText, TSpec> /* const */ & db,
                       TSeedLength minSeedLength)
{
    typedef Db<TText, TSpec> const                            TDb;
    typedef typename Size<TDb>::Type                          TDbSize;
    typedef typename Value<TText>::Type                       TTextReference;
    typedef typename Size<TText>::Type                        TTextSize;
    typedef TTextSize                                         TErrors;

    TDbSize dbSize = length(db.text);

    TTextSize exactSeedLength = db.avgLength / (getAvgErrors(db) + 1);
    TTextSize seedLength = _max(exactSeedLength, minSeedLength);

    TSeedLength maxSeedErrors = 0;
    resize(seedCounts, getMaxErrors(db) + 1, 0, Exact());

    for (TDbSize dbId = 0; dbId < dbSize; ++dbId)
    {
        TTextReference text = db.text[dbId];
        TTextSize textLength = length(text);

        TErrors errors = getErrors(db, dbId);

        TTextSize seedCount = _max(textLength / seedLength, 1u);
        TErrors seedErrors = errors / seedCount;
        TTextSize seedCountHigh = (errors % seedCount) + 1;
        TTextSize seedCountLow = seedCount - seedCountHigh;

        seedCounts[seedErrors] += seedCountHigh;
        if (seedErrors > 0) seedCounts[seedErrors - 1] += seedCountLow;

        maxSeedErrors = _max(maxSeedErrors, seedErrors);
    }

    resize(seedCounts, maxSeedErrors + 1);

    return seedLength;
}

// ----------------------------------------------------------------------------
// Function build()                                                   [DbIndex]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TDbIndexSpec, typename TDbSpec, typename TSpec>
void build(DbIndex<Index<TText, TIndexSpec>, TDbIndexSpec> & dbIndex,
           Db<TText, TDbSpec> /* const */ & db,
           TSpec)
{
    //typedef Db<TText, TDbSpec>                              TDb;
    //typedef typename Value<TText>::Type                     TTextReference;
    typedef typename Size<TText>::Type                      TTextSize;
    typedef Index<TText, TIndexSpec>                        TIndex;
    typedef typename Value<TIndex>::Type                    TIndexAlphabet;
    typedef typename Size<TIndex>::Type                     TIndexSize;
    typedef typename Fibre<TIndex, FibreSA>::Type           TIndexSAFibre;
    typedef typename Value<TIndexSAFibre>::Type             TIndexSAPos;

    typedef typename Iterator<TIndex, TopDown<> >::Type     TIterator;

    // we don't know the maximal seed length in advance as it depends on the query
    // hence we choose a sufficiently large number
    unsigned maxSeedLength = 10000;
    dbIndex.index = TIndex(db.text);

    // 1. create and sort q-gram buckets
    TIndexSAFibre &sa = indexSA(dbIndex.index);
    TText const &text = db.text;
    Shape<TIndexAlphabet, SimpleShape> shape;
    String<TIndexSize> &dir = dbIndex.dir;

    unsigned shapeLength;
    if (ValueSize<TIndexAlphabet>::VALUE <= 5)
        shapeLength = _min(maxSeedLength, 10u);
    else
        shapeLength = _min(maxSeedLength, 3u);

    dbIndex.sortedPrefix = shapeLength;
    TTextSize stepSize = 1;
    if (IsSameType<TSpec, Exact>::VALUE)
        stepSize = maxSeedLength;

    resize(shape, shapeLength);
    resize(sa, _qgramQGramCount(text, shape, stepSize), Exact());
    resize(dir, _fullDirLength(shape), Exact());
    Nothing nothing;

    createQGramIndex(sa, dir, nothing, text, shape, stepSize);

    // 2. refine q-gram buckets and sort up to their maxSeedLength prefix
    if (shapeLength < maxSeedLength)
    {
        maxSeedLength -= shapeLength; // delta between maxSeedLength and shapeLength is used for QGramLessOffset_
        typename Iterator<TIndexSAFibre, Standard>::Type saBegin = begin(sa, Standard());

        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic,1))
        for (int i = 1; i < (int)length(dir); ++i)
        {
            if (dir[i - 1] + 1 < dir[i])
            {
		auto infixSA = infix(sa, saBegin + dir[i - 1], saBegin + dir[i]);
	        sort(
		    infixSA,
		    QGramLessOffset_<TIndexSAPos, TText const>(text, maxSeedLength, shapeLength),
		    Parallel());
            }
//            typename Infix<TIndexSAFibre>::Type saInf = infix(sa, dir[i - 1], dir[i]);
//            _refineQGramIndexBucket(
//                saInf,
//                text,
//                shapeLength,
//                maxSeedLength);
        }
    }

    TIterator it(dbIndex.index);


//    TIndexSAFibre sa = indexSA(dbIndex.index);
//    resize(sa, lengthSum(db.text));
//
//    TDbSize dbSize = length(db.text);
//    for (TDbSize dbId = 0; dbId < dbSize; ++dbId)
//    {
//        TTextReference text = db.text[dbId];
//        TTextSize textLength = length(text);
//
//        TIndexSAPos saPos;
//        assignValueI1(saPos, dbId);
//        for (TTextSize pos = 0; pos < textLength; ++pos)
//        {
//            assignValueI2(saPos, pos);
//            appendValue(sa, saPos);
//        }
//    }
//
//    QGramLess_<TIndexSAPos, TText const> less(db.text, std::numeric_limits<TTextSize>::max());
//    sort(sa, less, Parallel());
}

// ----------------------------------------------------------------------------
// Function build()                                            [DbIndex<Query>]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TDbSpec, typename TSeedLength>
void buildQuery(DbIndex<Index<TText, IndexSa<TIndexSpec> >, Query> & dbIndex,
           Db<TText, TDbSpec> /* const */ & db,
           TSeedLength minSeedLength)
{
    typedef Db<TText, TDbSpec>                                 TDb;
    typedef typename Size<TDb>::Type                           TDbSize;
    typedef DbIndex<Index<TText, IndexSa<TIndexSpec> >, Query> TDbIndex;
    typedef typename Size<TDbIndex>::Type                      TDbIndexSize;
    typedef Index<TText, IndexSa<TIndexSpec> >                 TIndex;
    //typedef typename Fibre<TIndex, FibreSA>::Type             TIndexSAFibre;
    typedef typename Size<TText>::Type                         TTextSize;
    typedef TTextSize                                          TErrors;

    String<TDbSize> seedCounts;

    // Count seeds.
    dbIndex.seedLength = countSeeds(seedCounts, db, minSeedLength);

    // Get non-zero seed counts.
    TErrors maxSeedErrors = length(seedCounts);
    for (TErrors seedErrors = 0; seedErrors < maxSeedErrors; ++seedErrors)
        if (seedCounts[seedErrors] > 0)
            appendValue(dbIndex.errors, seedErrors);

    std::cout << "Seed length:\t\t\t\t" << dbIndex.seedLength << std::endl;
    std::cout << "Seed counts:\t\t\t\t";
    std::copy(begin(seedCounts, Standard()),
              end(seedCounts, Standard()),
              std::ostream_iterator<TDbSize>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << "Seed errors:\t\t\t\t";
    std::copy(begin(dbIndex.errors, Standard()),
              end(dbIndex.errors, Standard()),
              std::ostream_iterator<TDbIndexSize>(std::cout, ", "));
    std::cout << std::endl;

    // Resize indices.
    TErrors seedSetsCount = length(dbIndex.errors);
    resize(dbIndex.index, seedSetsCount, Exact());
    for (TErrors seedSet = 0; seedSet < seedSetsCount; ++seedSet)
    {
        TErrors seedErrors = dbIndex.errors[seedSet];
        dbIndex.index[seedSet] = TIndex(db.text);
        reserve(indexSA(dbIndex.index[seedSet]), seedCounts[seedErrors], Exact());
    }

    // Build SA fibres.
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (int seedSet = 0; seedSet < (int)seedSetsCount; ++seedSet)
    {
        TErrors seedErrors = dbIndex.errors[seedSet];
        _buildSA(indexSA(dbIndex.index[seedSet]), db, dbIndex.seedLength, seedErrors);
    }
}

// ----------------------------------------------------------------------------
// Function _buildSA()                                         [DbIndex<Query>]
// ----------------------------------------------------------------------------

template <typename TIndexSAFibre, typename TText, typename TDbSpec, typename TSeedErrors, typename TSeedLength>
void _buildSA(TIndexSAFibre & sa,
              Db<TText, TDbSpec> /* const */ & db,
              TSeedLength seedLength,
              TSeedErrors seedErrors)
{
    typedef Db<TText, TDbSpec>                          TDb;
    typedef typename Size<TDb>::Type                    TDbSize;
    typedef typename Value<TIndexSAFibre>::Type         TIndexSAPos;
    typedef typename Value<TText>::Type                 TTextReference;
    typedef typename Size<TText>::Type                  TTextSize;
    typedef TTextSize                                   TErrors;

    TDbSize dbSize = length(db.text);

    for (TDbSize dbId = 0; dbId < dbSize; ++dbId)
    {
        TTextReference text = db.text[dbId];
        TTextSize textLength = length(text);

        TErrors errors = getErrors(db, dbId);

        TTextSize seedCount = _max(textLength / seedLength, 1u);
        TSeedErrors seedErrors_ = errors / seedCount;
        TTextSize seedCountHigh = (errors % seedCount) + 1;

        TIndexSAPos seed;
        assignValueI1(seed, dbId);

        if (seedErrors_ == seedErrors)
        {
            TTextSize seedCounter = 0;
            TTextSize seedPos = 0;
            for (; seedCounter < seedCountHigh; ++seedCounter, seedPos += seedLength)
            {
                assignValueI2(seed, seedPos);
                appendValue(sa, seed);
            }
        }
        else if (seedErrors_ > 0 && seedErrors_ - 1 == seedErrors)
        {
            TTextSize seedCounter = seedCountHigh;
            TTextSize seedPos = seedCountHigh * seedLength;
            for (; seedCounter < seedCount; ++seedCounter, seedPos += seedLength)
            {
                assignValueI2(seed, seedPos);
                appendValue(sa, seed);
            }
        }
    }

    // Construct index using quicksort.
    QGramLess_<TIndexSAPos, TText const> less(db.text, seedLength);
    sort(sa, less, Parallel());
}

// ----------------------------------------------------------------------------
// Function build()                                                   [DbIndex]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Build a trie of Db using QGram bucket sort.

//template <typename TDb, typename TIndexText, typename TShape, typename TIndexSpec, typename TSpec>
//void build(DbIndex<TDb, Index<TIndexText, IndexQGram<TShape, TIndexSpec> >, TSpec> & dbIndex)
//{
//    typedef Index<TIndexText, IndexQGram<TShape, TIndexSpec> >  TIndex;
//    typedef typename Fibre<TIndex, QGramSA>::Type               TIndexSAFibre;
//    typedef typename Fibre<TIndex, QGramDir>::Type              TIndexDirFibre;
//    typedef typename Fibre<TIndex, QGramShape>::Type            TIndexShape;
//    typedef typename Fibre<TIndex, QGramBucketMap>::Type        TIndexBucketMap;
//    typedef typename Value<TIndexSAFibre>::Type                 TSAPos;
//
//    typedef typename Size<TDb>::Type                            TDbSize;
//    typedef typename Host<TDb>::Type const                      TText;
//    typedef typename Value<TText>::Type const                 TText;
//    typedef typename Iterator<TText, Standard>::Type            TTextIterator;
//
//    // NOTE(esiragusa): This is done in DbIndex constructor.
////    dbIndex.index = TIndex(dbIndex.db.text);
//
//    // NOTE(esiragusa): This is to index whole strings.
////    setStepSize(dbIndex.index, 256u);
//
//    TIndexSAFibre & sa = indexSA(dbIndex.index);
//    TIndexDirFibre & dir = indexDir(dbIndex.index);
//    TIndexShape & shape = indexShape(dbIndex.index);
//    TIndexBucketMap & bucketMap = indexBucketMap(dbIndex.index);
//
//    // Resize suffix array and directory.
//    TDbSize dbSize = length(dbIndex.db.text);
//    resize(sa, dbSize, Exact());
//    resize(dir, _fullDirLength(dbIndex.index), Exact());
//
//    // Clear directory.
//    _qgramClearDir(dir, bucketMap);
//
//    // Count qgrams.
//    for (TDbSize dbId = 0; dbId < dbSize; ++dbId)
//    {
//        TText & text = dbIndex.db.text[dbId];
//        TTextIterator textIt = begin(text, Standard());
//        ++dir[requestBucket(bucketMap, hash(shape, textIt))];
//    }
//
//    // Compute cumulative sum.
//    _qgramCummulativeSum(dir, False());
//
//    // Fill suffix array.
//    for (TDbSize dbId = 0; dbId < dbSize; ++dbId)
//    {
//        TText & text = dbIndex.db.text[dbId];
//        TTextIterator textIt = begin(text, Standard());
//
//        TSAPos saPos;
//        assignValueI1(saPos, dbId);
//        assignValueI2(saPos, 0);
//        sa[dir[getBucket(bucketMap, hash(shape, textIt)) + 1]++] = saPos;
//    }
//
//    // Refine buckets.
//    _refineQGramIndex(sa, dir, indexText(dbIndex.index), weight(shape), 256);
//    _setHost(dbIndex.index);
//}

#endif  // #ifndef SEQAN_APPS_SEARCHJOIN_DB_H_
