// ==========================================================================
//                 seqan - the library for sequence analysis
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

#ifndef INDEX_FM_LF_TABLE_H_
#define INDEX_FM_LF_TABLE_H_

namespace seqan {

// ==========================================================================
// LfTable is an object storing all necessary information for the LF-mapping.
// To be more precise, the occurrence-table data structure as well as the
// prefix-sum table are stored.
// ==========================================================================

template <typename TText, typename TSpec>
class WaveletTree;

template <typename TSpec>
struct FmiDollarSubstituted;

/**
.Tag.LF Table Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibreOccTable:The occurrence table of the lf table.
..tag.FMTablePrefixSumTable:The prefix sum table of the lf table.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/
struct FibreOccTable_;
struct FibrePrefixSumTable_;
struct FibreDollarPosition_;


typedef Tag<FibreOccTable_> const FibreOccTable;
typedef Tag<FibrePrefixSumTable_> const FMTablePrefixSumTable;
typedef Tag<FibreDollarPosition_> const FibreDollarPosition;

template <typename TOccTable, typename TPrefixSumTable>
struct LfTable;

// ==========================================================================
//Metafunctions
// ==========================================================================

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable>, FibreOccTable>
{
    typedef TOccTable Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable> const, FibreOccTable>
{
    typedef TOccTable const Type;
};

// ==========================================================================
template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>
{
    typedef TPrefixSumTable Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LfTable<TOccTable, TPrefixSumTable> const, FMTablePrefixSumTable>
{
    typedef TPrefixSumTable const Type;
};

// ==========================================================================
template <typename TOccTable, typename TPrefixSumTable>
struct Reference<LfTable<TOccTable, TPrefixSumTable> >
{
    typedef TPrefixSumTable & Type;
};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Class.LfTable:
..cat:Index
..summary:LfTable is an object storing all necessary information for the LF-mapping.
..signature:LfTable<TOccTable, TPrefixSumTable>
..param.TOccTable:The occurrence table data structure.
...type.class:WaveletTree
..param.TPrefixSumTable:The specialisation tag.
...default:String
..include:seqan/Index.h
*/
template <typename TOccTable, typename TPrefixSumTable>
struct LfTable
{
    TOccTable occTable;
    TPrefixSumTable prefixSumTable;

    LfTable() :
        occTable(),
        prefixSumTable()
    {}

    LfTable(TOccTable const & occTable, TPrefixSumTable const & prefixSumTable) :
        occTable(occTable),
        prefixSumTable(prefixSumTable)
    {}

    inline LfTable & operator=(LfTable const & other)
    {
        occTable = other.occTable;
        prefixSumTable = other.prefixSumTable;
        return *this;
    }

    inline bool operator==(const LfTable & b) const
    {
        return occTable == b.occTable &&
               prefixSumTable == b.prefixSumTable;
    }

};

// ==========================================================================
// Functions
// ==========================================================================

///.Function.clear.param.object.type:Class.LfTable
///.Function.clear.class:Class.LfTable
template <typename TOccTable, typename TPrefixSumTable>
inline void clear(LfTable<TOccTable, TPrefixSumTable> & lfTable)
{
    clear(lfTable.occTable);
    clear(lfTable.prefixSumTable);
}

// ==========================================================================
///.Function.empty.param.object.type:Class.LfTable
///.Function.empty.class:Class.LfTable
template <typename TOccTable, typename TPrefixSumTable>
inline bool empty(LfTable<TOccTable, TPrefixSumTable> & lfTable)
{
    return empty(lfTable.occTable)
           && empty(lfTable.prefixSumTable);
}

// ==========================================================================
/**
.Function.LfTable#createOccurrenceTable
..class:Class.LfTable
..cat:Index
..signature:createOccurrenceTable(lfTable, text, sentinalSub, sentinalPos)
..summary:This function creates the occurrence table data structure of a LF table.
..param.lfTable:The LF table to be constructed.
...type:Class.LfTable
..param.text:The text of which the occurrence table is constructed
...type:Class.String
..param.sentinalSub:The character used to substitute the sentinel with.
...type:Metafunction.Value type of String<TSpec>
..param.sentinalPos:The position of the sentinal(s) in the text.
...type:Class.String
...type:Class.RankSupportBitString
..returns:$bool$
*/

template <typename TText, typename TWaveletTreeSpec, typename TPrefixSumTable, typename TText2, typename TChar, typename TPos>
inline bool createOccurrenceTable(LfTable<WaveletTree<TText, TWaveletTreeSpec>, TPrefixSumTable> & lfTable, TText2 & text, TChar dollarSub, TPos const & dollarPos)
{
    waveletTreeCreate(lfTable, text, dollarSub, dollarPos);
    return true;
}

// ==========================================================================
/**
.Function.createLfTable
..summary:Creates the LF table
..signature:createLfTable(lfTable, text)
..param.lfTable:The LF table to be constructed.
...type:Class.LfTable.
..param.text:The underlying text
...type:Class.String
..include:seqan/index.h
*/
template <typename TPrefixSumTable, typename TText, typename TDollarSpec, typename TText2>
inline bool createLfTable(LfTable<WaveletTree<TText, FmiDollarSubstituted<TDollarSpec> >, TPrefixSumTable> & lfTable, TText2 const text)
{
    typedef typename SAValue<TText2>::Type   TSAValue;
    typedef typename Value<WaveletTree<TText, FmiDollarSubstituted<TDollarSpec> > >::Type     TAlphabet;
    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TDollarSpec> >, FibreDollarPosition>::Type TDollarPos;

    String<TSAValue> sa;
    resize(sa, length(text));
    createSuffixArray(sa, text, Skew7());

    createPrefixSumTable(lfTable.prefixSumTable, text);

    TAlphabet dollarSub;
    _determineDollarSubstitute(lfTable.prefixSumTable, dollarSub);

    String<TAlphabet> bwt;
    TDollarPos dollarPos = 0;
    resize(bwt, _computeBwtLength(text));
    _createBwTable(bwt, dollarPos, text, sa, dollarSub);

    clear(sa);

    createOccurrenceTable(lfTable, bwt, dollarSub, dollarPos);
    clear(bwt);

    _insertDollar(lfTable.prefixSumTable, countSequences(text));

    //here comes the dollar modification
    //addDollarNode(index.lfTable.occTable, dollarSub, dollarPos);
    return true;
}

// ==========================================================================
/**
.Function.LfTable#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.LfTable
..cat:Index
..param.container:The container holding the fibre.
...type:Class.LfTable
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.LF Table Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
*/
template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>::Type &
getFibre(LfTable<TOccTable, TPrefixSumTable>&lfTable, FMTablePrefixSumTable)
{
    return lfTable.prefixSumTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>::Type const &
getFibre(LfTable<TOccTable, TPrefixSumTable> const & lfTable, FMTablePrefixSumTable)
{
    return lfTable.prefixSumTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type &
getFibre(LfTable<TOccTable, TPrefixSumTable>&lfTable, FibreOccTable)
{
    return lfTable.occTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LfTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type const &
getFibre(LfTable<TOccTable, TPrefixSumTable> const & lfTable, FibreOccTable )
{
    return lfTable.occTable;
}

// ==========================================================================
/**
.Function.lfMapping:
..summary:Returns the position of the character L[c] in F.
..cat:Index
..signature:lfMapping(lfTable, pos)
..param.lfTable:The @Class.LfTable@ holding the occurrence and prefix sum table.
...type:Class.LfTable
..param.pos:The position in L
..returns:Returns the position of the character L[c] in F. The returned position is of the same type as pos.
..include:seqan/index.h
*/
template <typename TLfTable, typename TPos>
inline TPos lfMapping(TLfTable & lfTable,
                      TPos pos)
{
    typedef typename Fibre<TLfTable, FibreOccTable>::Type TOccTable;
    typedef typename Value<TOccTable>::Type TChar;
    TChar c = getCharacter(lfTable.occTable, pos);
    return countOccurrences(lfTable.occTable, c, pos) + getPrefixSum(lfTable.prefixSumTable, getCharacterPosition(lfTable.prefixSumTable, c)) - 1;
}

// ==========================================================================
/**
.Function.open
..param.object:
...type:Class.LfTable
*/
template <typename TOccTable, typename TPrefixSumTable>
inline bool open(
    LfTable<TOccTable, TPrefixSumTable> & lfTable,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName; 
    if (!open(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
    {
        return false;
    }
    name = fileName;    open(getFibre(lfTable, FMTablePrefixSumTable()), toCString(name), openMode);
    return true;

}

template <typename TOccTable, typename TPrefixSumTable>
inline bool open(
    LfTable<TOccTable, TPrefixSumTable> & lfTable,
    const char * fileName)
{
    return open(lfTable, fileName, DefaultOpenMode<LfTable<TOccTable, TPrefixSumTable> >::VALUE);
}

// ==========================================================================
/**
.Function.save
..param.object:
...type:Class.LfTable
*/
template <typename TOccTable, typename TPrefixSumTable>
inline bool save(
    LfTable<TOccTable, TPrefixSumTable> const & lfTable,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;
    if (!save(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
    {
        return false;
    }
    name = fileName;    save(getFibre(lfTable, FMTablePrefixSumTable()), toCString(name), openMode);
    return true;
}

template <typename TOccTable, typename TPrefixSumTable>
inline bool save(
    LfTable<TOccTable, TPrefixSumTable> const & lfTable,
    const char * fileName)
{
    return save(lfTable, fileName, DefaultOpenMode<LfTable<TOccTable, TPrefixSumTable> >::VALUE);
}

}
#endif // INDEX_FM_LF_TABLE_H_
