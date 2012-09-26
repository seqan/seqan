// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef INDEX_FM_COMPRESSED_SA_H_
#define INDEX_FM_COMPRESSED_SA_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

struct FibreSparseString_;
typedef Tag<FibreSparseString_> const FibreSparseString;

template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA;

// ==========================================================================
//Metafunctions
// ==========================================================================

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>
{
    typedef TSparseString Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Fibre<CompressedSA<TSparseString, TLfTable, TSpec> const, FibreSparseString>
{
    typedef TSparseString Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<CompressedSA<TSparseString, TLfTable, TSpec> >
{
    // TODO(singer): We actually need a proxy here.
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<const CompressedSA<TSparseString, TLfTable, TSpec> >
{
    // TODO(singer): We actually need a proxy here.
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type /*const*/ Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Value<CompressedSA<TSparseString, TLfTable, TSpec> >
{
    typedef typename Value<TSparseString>::Type Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Value<CompressedSA<TSparseString, TLfTable, TSpec> const>
{
    typedef typename Value<TSparseString>::Type const Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// forwards
template <typename TPos, typename TOffSet>
TPos addGapDistance_(TPos const & value, TOffSet const & offSet);

template <typename TSeqId, typename TSpec, typename TPos, typename TOffSet>
Pair<TSeqId, TPos> addGapDistance_(Pair<TSeqId, TPos, TSpec> const & value, TOffSet const & offSet);

/**
.Class.CompressedSA:
..cat:String
..summary:A suffix array storing only a few suffix array entries and computing the remaining on demand.
..signature:CompressedSA<TSparseString, TLfTable, TSpec>
..param.TSparseString:The string containing specific suffix array entries.
...type:Class.SparseString
..param.TLfTable:The lfTable containg an occurrence table and a prefix sum table.
...type:Class.LfTable
..param.TSpec:Possibility to specialise a compressed suffix array.
...default:void.
..remarks:The compressed suffix array can only be used with the FM index.
..include:seqan/index.h
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA
{
public:
    TSparseString   sparseString;
    TLfTable *      lfTable;

    CompressedSA(){};

    CompressedSA(TLfTable & lfTable) :
        lfTable(&lfTable)
    {}

    inline CompressedSA & operator=(CompressedSA const & other)
    {
        sparseString = other.sparseString;
        lfTable = other.lfTable;
        return *this;
    }

    typedef typename Value<typename Fibre<TSparseString, FibreValueString>::Type>::Type TCompressedSaValue;
    typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;

    template <typename TPos>
    inline TCompressedSaValue const operator[](TPos pos)
    {
        TIndicatorString const & indicatorString = getFibre(sparseString, FibreIndicatorString());
        TPos counter = 0;

        while (!getBit(indicatorString, pos))
        {
            pos = lfMapping(*lfTable, pos);
            ++counter;
        }
        return addGapDistance_(getValue(sparseString.valueString, getRank(indicatorString, pos) - 1), counter);
    }

    template <typename TPos>
    inline TCompressedSaValue operator[](TPos pos) const
    {
        TIndicatorString const & indicatorString = getFibre(sparseString, FibreIndicatorString());
        TPos counter = 0;
        while (!getBit(indicatorString, pos))
        {
            pos = lfMapping(*lfTable, pos);
            ++counter;
        }
        return addGapDistance_(getValue(sparseString.valueString, getRank(indicatorString, pos) - 1), counter);
    }

    inline bool operator==(const CompressedSA & other) const
    {
        return sparseString == other.sparseString &&
               *lfTable == *(other.lfTable);
    }

};

template <typename TPos, typename TOffSet>
TPos addGapDistance_(TPos const & value, TOffSet const & offSet)
{
    return value + offSet;
}

template <typename TSeqId, typename TSpec, typename TPos, typename TOffSet>
Pair<TSeqId, TPos> addGapDistance_(Pair<TSeqId, TPos, TSpec> const & value, TOffSet const & offSet)
{
    return Pair<TSeqId, TPos>(value.i1, value.i2 + offSet);
}

// ==========================================================================
// Functions
// ==========================================================================

/*
.Function.assignCompressionFactor
..param.container:The container holding the entries.
...type:Class.CompressedSA
*/
// template <typename TSparseString, typename TLfTable, typename TSpec, typename TValue>
// void assignCompressionFactor(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TValue value)
// {
//     assignCompressionFactor(getFibre(compressedSA, FibreSparseString()), value);
// }

/**
.Function.assignValue
..param.object
...type:Class.CompressedSA
..remarks:In the case of the compresses suffix array the new value only influences the sparse string and not the indicator string.
..include:seqan/index.h
*/
template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos, typename TValue>
void assignValue(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos, TValue value)
{
    assignValue(getFibre(compressedSA, FibreSparseString()), pos, value);
}

/**
.Function.clear
..param.object:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline void clear(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    clear(getFibre(compressedSA, FibreSparseString()));
}

/**
.Function.empty
..param.object:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool empty(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return empty(getFibre(compressedSA, FibreSparseString()));
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline bool entryStored(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, TPos const & pos)
{
    return entryStored(getFibre(compressedSA, FibreSparseString()), pos);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline bool entryStored(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos const & pos)
{
    return entryStored(getFibre(compressedSA, FibreSparseString()), pos);
}

// This function creates a compressed suffix array using a normal one.
template <typename TSparseString, typename TLfTable, typename TSpec, typename TSA, typename TCompression>
void compressedSaCreate(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, 
        TSA const & completeSA, 
        TCompression const compressionFactor, 
        unsigned offSet)
{
    typedef CompressedSA<TSparseString, TLfTable, TSpec> TCompressedSA;
    typedef typename GetValue<TSA>::Type                            TSAValue;
    typedef typename Size<TSA>::Type                                TSize;
    typedef typename Fibre<TCompressedSA, FibreSparseString>::Type  TSparseSA;
    typedef typename Fibre<TSparseSA, FibreIndicatorString>::Type   TIndicatorString;

    TSparseSA & sparseString = getFibre(compressedSA, FibreSparseString());
    TIndicatorString & indicatorString = getFibre(sparseString, FibreIndicatorString());

    TSize n = length(completeSA);

    resize(compressedSA, offSet + n);
    for (TSize i = 0; i < n; i++)
    {
        TSAValue sa = getValue(completeSA, i);
        (getSeqOffset(sa) % compressionFactor == 0) ? setBit(indicatorString, i + offSet, 1) : setBit(indicatorString, i + offSet, 0);
    }
    updateRanks_(indicatorString);

    resize(sparseString.valueString, getRank(indicatorString, length(indicatorString) - 1));

    TSize counter = 0;
    for (TSize i = 0; i < n; i++)
    {
        if (getBit(indicatorString, i + offSet))
        {
            assignValue(compressedSA.sparseString.valueString, counter, completeSA[i]);
            ++counter;
        }
    }
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSA, typename TCompression>
void compressedSaCreate(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TSA const & completeSA, TCompression const compressionFactor)
{
    compressedSaCreate(compressedSA, completeSA, compressionFactor, 0);
}

/**
.Function.getCompressionFactor
..param.container:The container holding the entries.
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Size<typename Fibre<TSparseString, FibreValueString>::Type>::Type
getCompressionFactor(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return getCompressionFactor(getFibre(compressedSA, FibreSparseString()));
}

/**
.Function.getFibre
..param.container:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type const &
getFibre(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, FibreSparseString)
{
    return compressedSA.sparseString;
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type &
getFibre(CompressedSA<TSparseString, TLfTable, TSpec>&compressedSA, FibreSparseString)
{
    return compressedSA.sparseString;
}

// This functions computes the position in the suffix array of text[sa[pos] - 1]
// iff the current position is not present in the compressed suffix array.
template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline bool getNextPos_(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, TPos & pos)
{
    typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;
    TIndicatorString const & indicatorString = compressedSA.sparseString.indicatorString;

    if (getBit(indicatorString, pos))
    {
        return true;
    }
    pos = lfMapping(*compressedSA.lfTable, pos);
    return false;
}

/**
.Function.length
..param.object:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Size<typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type>::Type
length(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return length(getFibre(compressedSA, FibreSparseString()));
}

/**
.Function.length
..param.object:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void resize(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
                   TSize size)
{
    resize(getFibre(compressedSA, FibreSparseString()), size);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void resize(CompressedSA<String<TSparseString>, TLfTable, TSpec> & compressedSA, TSize size)
{
    resize(getFibre(compressedSA, FibreSparseString()), size);
}

/**
.Function.open
..param.string:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool open(
    CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;
    if (!open(getFibre(compressedSA, FibreSparseString()), toCString(name), openMode))
    {
        return false;
    }
    return true;

}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool open(
    CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
    const char * fileName)
{
    return open(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
}

/**
.Function.save
..param.string:
...type:Class.CompressedSA
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool save(
    CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;
    if (!save(getFibre(compressedSA, FibreSparseString()), toCString(name), openMode))
    {
        return false;
    }
    return true;
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool save(
    CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
    const char * fileName)
{
    return save(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
}

// ==========================================================================
// Functions
// ==========================================================================
/**
.Function.setLfTable
..summary:Set the LfTable of the compressed suffix array..
..signature:setLfTable(CompressedSA<TSparseString, TLfTable, TSpec> compressedSa, TLfTable & lfTable)
..param.CompressedSA<TSparseString, TLfTable, TSpec>:The compressed suffix array.
...type:Class.CompressedSA
..param.lfTable
...type.Class:LfTable
..include:seqan/index.h
*/
template <typename TSparseString, typename TLfTable, typename TSpec>
void setLfTable(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TLfTable & lfTable)
{
    compressedSA.lfTable = &lfTable;
}

/**
.Function.value
..param.container:
...type:Class.CompressedSA
..remarks:Note that the compressed suffix array is read only. Therefore a const reference is return by
this function.
*/
template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<TSparseString>::Type
value(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
    return compressedSA[pos];
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<TSparseString>::Type const
value(const CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
    return compressedSA[pos];
}

}

//template <typename TSpecPairI1, typename TSpecPairI2, typename TSpecPairSpec, typename TStringSpec, typename TSparseStringSpec, typename TLfTable, typename TSpec>
//struct CompressedSA<SparseString<String<Pair<TSpecPairI1, TSpecPairI2, TSpecPairSpec>, TStringSpec>, TSparseStringSpec>, TLfTable, TSpec>
//{
//    typedef SparseString<String<Pair<TSpecPairI1, TSpecPairI2, TSpecPairSpec>, TStringSpec>, TSparseStringSpec>   TSparseString;
//    typedef typename Value<typename Fibre<SparseString<TSparseString, TStringSpec>, FibreSparseString>::Type>::Type TCompressedSaValue;
//    typedef typename Fibre<SparseString<TSparseString, TStringSpec>, FibreFibreIndicatorString>::Type TFibreIndicatorString;
//
//    TSparseString   sparseString;
//    TLfTable *        lfTable;
//
//    CompressedSA() :
//      sparseString(),
//      lfTable()
//    {}
//
//    template <typename TPos>
//    inline TCompressedSaValue const operator[](TPos pos)
//    {
//        TFibreIndicatorString const & FibreIndicatorString = getFibre(sparseString, FibreFibreIndicatorString());
//        TPos counter = 0;
//        while (!getBit(FibreIndicatorString, pos))
//        {
//            pos = lfMapping(*lfTable, pos);
//            ++counter;
//        }
//        TCompressedSaValue temp = getValue(sparseString, getRank(FibreIndicatorString, pos) - 1);
//        temp.i2 += counter;
//        return temp;
//    }
//
//    template <typename TPos>
//    inline TCompressedSaValue operator[](TPos pos) const
//    {
//        TFibreIndicatorString const & FibreIndicatorString = getFibre(sparseString, FibreFibreIndicatorString());
//        TPos counter = 0;
//        while (!getBit(FibreIndicatorString, pos))
//        {
//            pos = lfMapping(*lfTable, pos);
//            ++counter;
//        }
//        TCompressedSaValue temp = getValue(sparseString, getRank(FibreIndicatorString, pos) - 1);
//        temp.i2 += counter;
//        return temp;
//    }
//
//    inline bool operator==(const CompressedSA & b) const
//    {
//        return sparseString == b.sparseString &&
//               *lfTable == *(b.lfTable);
//    }
//
//};


#endif // INDEX_FM_COMPRESSED_SA_H_
