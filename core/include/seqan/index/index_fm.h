// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_H
#define INDEX_FM_H

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FMIndexConfig
// ----------------------------------------------------------------------------

/*!
 * @class FMIndexConfig
 * @headerfile seqan/index.h
 * @brief A configuration object that determines the data types of certain fibres of the @link FMIndex @endlink.
 * 
 * @signature template <[typename TSpec]>
 *            struct FMIndexConfig;
 *
 * @tparam TSpec The specializating type, defaults to <tt>void</tt>.
 *
 * @var unsigned FMIndexConfig::SAMPLING;
 * @brief The sampling rate determines how many suffix array entries are represented with one entry in the
 *        @link CompressedSA @endlink.
 *
 * @typedef FMIndexConfig::TValuesSpec
 * @signature typedef WaveletTree<TSpec> TValuesSpec;
 * @brief The <tt>TValuesSpec</tt> determines the type of the occurrence table. In the default @link FMIndexConfig
 *        @endlink object the type of <tt>TValuesSpec</tt> is a wavelet tree (@link WaveletTree @endlink).
 *
 * @typedef FMIndexConfig::TSentinelsSpec
 * @signature typedef TwoLevels<TSpec> TSentinelsSpec;
 * @brief The <tt>TSentinelsSpec</tt> determines the type of the sentinels in the @link FMIndex @endlink.  In the
 *        default @link FMIndexConfig @endlink object the type of <tt>TSentinelsSpec</tt> is a two level
 *        @link RankDictionary @endlink.
 */
template <typename TSpec = void>
struct FMIndexConfig
{
    typedef WaveletTree<TSpec>      TValuesSpec;
    typedef TwoLevels<TSpec>        TSentinelsSpec;

    static const unsigned SAMPLING = 10;
};

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec = void, typename TConfig = FMIndexConfig<TSpec> >
class FMIndex;

// ============================================================================
// Tags
// ============================================================================

// FM index fibres
/**
.Tag.FM Index Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibreSA:The compressed suffix array of the text.
..tag.FibreText:The original text of the index.
..tag.FibreLF:The lf table.
..tag.FibreSALF:The lf table as well as the compressed suffix array.
...remarks:This tag can only be used with the functions @Function.indexRequire@ or @Function.indexSupplied@.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/

struct FibreTempSA_;
struct FibreLF_;
struct FibreSALF_;

typedef Tag<FibreTempSA_> const         FibreTempSA;
typedef Tag<FibreLF_> const             FibreLF;
typedef Tag<FibreSALF_> const           FibreSALF;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

/*
.Metafunction.Fibre:
..summary:Type of a specific FMIndex member (fibre).
..signature:Fibre<Index<TText, FMIndex<TSentinelRankDictionary, TSpec> >, TFibreSpec>::Type
..class:Spec.FMIndex
..cat:Index
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TSentinelRankDictionary:The type of the sentinel rank dictionary.
...type:Class.SentinelRankDictionary
..param.TSpec:Tag to specify a certain variant of the FM index.
...default;$void$
..param.TFibreSpec:Tag to specify the fibre.
...type:Tag.FM Index Fibres
..returns:Fibre type.
..remarks:Some containers, such as @Spec.FMIndex@, can be seen as a bundle consisting of various fibres. Because not 
every table is a fibre we did not call them tables, however, in many cases one can think of fibres as tables. The 
fibre interface was designed to unify the access to the members of the different fibres.
To get a reference or the type of a specific fibre use @Function.getFibre@ or @Metafunction.Fibre@.		
..include:seqan/index.h
*/

/*!
 * @defgroup FMIndexFibres FM Index Fibres
 * @brief Tag to select a specific fibre of a @link FMIndex @endlink.
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a FM index.
 * 
 * @see Fibre
 * @see Index#getFibre
 * 
 * @tag FMIndexFibres#FibreText
 * @brief The original text of the index.
 * 
 * @tagFMIndexFibres#FibreSA
 * @brief The compressed suffix array of the text.
 * 
 * @tag FMIndexFibres#FibreLF
 * @brief The lf table.
 */


template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>
{
    typedef LF<TText, TSpec, TConfig>  Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreSA>
{
    typedef CompressedSA<TText, TSpec, TConfig> Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreTempSA>
{
    typedef Index<TText, FMIndex<TSpec, TConfig> >          TIndex_;
    typedef typename SAValue<TIndex_>::Type                 TSAValue_;

    // NOTE(esiragusa): External causes problems on device code.
#ifndef PLATFORM_CUDA
    typedef String<TSAValue_, External<ExternalConfigLarge<> > >                Type;
#else
    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type>     Type;
#endif
};

// ----------------------------------------------------------------------------
// Metafunction DefaultFinder
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct DefaultFinder<Index<TText, FMIndex<TSpec, TConfig> > >
{
    typedef FinderSTree Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class FMIndex 
// ----------------------------------------------------------------------------

/**
.Spec.FMIndex:
..summary:An index based on the Burrows-Wheeler transform.
..cat:Index
..general:Class.Index
..signature:Index<TText, FMIndex<TSpec, TConfig> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TOccSpec:Occurrence table specialisation. 
...type:Tag.WT
...type:Tag.SBM
...remarks:The tags are really shortcuts for the different @Class.SentinelRankDictionary@s
...default:Tag.WT
..param.TSpec:FM index specialisation.
...default:void
..include:seqan/index.h
*/

/*!
 * @class FMIndex
 * @extends Index
 * @headerfile seqan/index.h
 * @brief An index based on the Burrows-Wheeler transform.
 * 
 * @signature template <typename TText[, typename TSpec[, typename TConfig]]>
 *            class Index<TText, FMIndex<TSpec, TConfig> >;
 * 
 * @tparam TText   The text type. Types: @link String @endlink, @link StringSet @endlink
 * @tparam TSpec   FM index specialisation, defaults to <tt>void</tt>.
 * @tparam TConfig A config object which determines the data types of the different fibres, defaults to
 *                 <tt>FMIndexConfig&lt;TSpec&gt;</tt>.
 */

template <typename TText, typename TSpec, typename TConfig>
class Index<TText, FMIndex<TSpec, TConfig> >
{
public:
    typename Member<Index, FibreText>::Type         text;
    typename Fibre<Index, FibreLF>::Type            lf;
    typename Fibre<Index, FibreSA>::Type            sa;

    Index() {};

    Index(TText & text) :
        text(text)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

// Already documented
template <typename TText, typename TSpec, typename TConfig>
inline void clear(Index<TText, FMIndex<TSpec, TConfig> > & index)
{
    clear(getFibre(index, FibreText()));
    clear(getFibre(index, FibreLF()));
    clear(getFibre(index, FibreSA()));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

// This function checks whether the index is empty. Its already documented.
template <typename TText, typename TSpec, typename TConfig>
inline bool empty(Index<TText, FMIndex<TSpec, TConfig> > const & index)
{
    return empty(getFibre(index, FibreText())) &&
           empty(getFibre(index, FibreLF())) &&
           empty(getFibre(index, FibreSA()));
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#getFibre:
..summary:Returns a specific fibre of a fm index.
..signature:getFibre(index, fibreTag)
..class:Spec.FMIndex
..cat:Index
..param.index:The index holding the fibre.
...type:Spec.FMIndex
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.FM Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type &
getFibre(Index<TText, FMIndex<TSpec, TConfig> > & index, FibreLF /*tag*/)
{
    return index.lf;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type const &
getFibre(Index<TText, FMIndex<TSpec, TConfig> > const & index, FibreLF /*tag*/)
{
    return index.lf;
}

// ----------------------------------------------------------------------------
// Function indexLF()
// ----------------------------------------------------------------------------
/*!
 * @fn FMIndex#indexLF
 * @headerfile seqan/index.h
 * @brief A shortcut for <tt>getFibre(index, FibreLF())</tt>.
 *
 * @signature TFibre indexLF(index);
 * 
 * @param[in] index The FM index.
 *
 * @return TFibre A reference to the @link FMIndexFibres#FibreLF @endlink.
 */

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type &
indexLF(Index<TText, FMIndex<TSpec, TConfig> > & index)
{
    return getFibre(index, FibreLF());
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type const &
indexLF(Index<TText, FMIndex<TSpec, TConfig> > const & index)
{
    return getFibre(index, FibreLF());
}

// ----------------------------------------------------------------------------
// Function toSuffixPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn FMIndex#toSuffixPosition
 * @headerfile seqan/index.h
 * @brief This function computes the position of a specified position in the compressed suffix array (additionally 
 *        containing entries for the sentinels). The returned position corresponds to the suffix array of the original
 *        text without sentinels.
 * 
 * @signature TSAValue toSuffixPosition(index, pos, offset);
 * 
 * @param[in] index  The FM index.
 * @param[in] pos    The position in the suffix array of the FM index (with sentinels).
 *                   Types: @link UnsignedIntegerConcept @endlink
 * @param[in] offset The number of sequences in the original text.  Types: @link UnsignedIntegerConcept @endlink
 *
 * @return TSAValue The function function computes the position of a specified position in the compressed suffix array
 *                  (additionally containing entries for the sentinels).  The returned position corresponds to the
 *                  suffix array of the original text without sentinels.  The return type is @link SAValue
 *                  @endlink&lt;@link Index @endlink&lt;TText, FMIndex&lt;TSpec, TConfig&gt; &gt; &gt;::Type
 */

/**
.Function.FMIndex#toSuffixPosition
..class:Spec.FMIndex
..summary:This function computes the position of a specified position in the suffix array (additionally containing 
entries for the sentinels. The returned position correspond to the suffix array of the original text without sentinels.
..signature:toSuffixPosition(fmIndex, pos, offset)
..param.fmIndex:The FM index.
...type:Spec.FMIndex
..param.pos:The position in the suffix array of the fm index (with sentinels).
...type:Concept.UnsignedIntegerConcept
..param.offset:The number of sequences in the original text.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/

template <typename TText, typename TSpec, typename TConfig, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TSpec, TConfig> > >::Type
toSuffixPosition(Index<TText, FMIndex<TSpec, TConfig> > & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

template <typename TText, typename TSpec, typename TConfig, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TSpec, TConfig> > const>::Type
toSuffixPosition(Index<TText, FMIndex<TSpec, TConfig> > const & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#indexCreate
..summary:Creates a specific @Metafunction.Fibre@.
..signature:indexCreate(index, fibreTag)
..param.index:The index to be created.
...type:Spec.FMIndex
..param.fibreTag:The fibre of the index to be computed.
...type:Tag.FM Index Fibres.tag.FibreSALF
*/

template <typename TText, typename TSpec, typename TConfig>
inline bool indexCreate(Index<TText, FMIndex<TSpec, TConfig> > & index, FibreSALF)
{
    typedef Index<TText, FMIndex<TSpec, TConfig> >      TIndex;
    typedef typename Fibre<TIndex, FibreTempSA>::Type   TTempSA;
    typedef typename Size<TIndex>::Type                 TSize;

    TText const & text = indexText(index);

    if (empty(text))
        return false;

    TTempSA tempSA;

    // Create the full SA.
    resize(tempSA, lengthSum(text), Exact());
    createSuffixArray(tempSA, text, Skew7());

    // Create the LF table.
    createLF(indexLF(index), text, tempSA);

    // Set the FMIndex LF as the CompressedSA LF.
    setFibre(indexSA(index), indexLF(index), FibreLF());

    // Create the compressed SA.
    TSize numSentinel = countSequences(text);
    createCompressedSa(indexSA(index), tempSA, numSentinel);

    return true;
}

template <typename TText, typename TSpec, typename TConfig>
inline bool indexCreate(Index<TText, FMIndex<TSpec, TConfig> > & index, FibreSA)
{
    return indexCreate(index, FibreSALF());
}

template <typename TText, typename TSpec, typename TConfig>
inline bool indexCreate(Index<TText, FMIndex<TSpec, TConfig> > & index)
{
    return indexCreate(index, FibreSALF());
}

// ----------------------------------------------------------------------------
// Function indexSupplied()
// ----------------------------------------------------------------------------

/**
.Function.FMIndex#indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..param.fibreTag:
...type:Tag.FM Index Fibres
*/
template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline bool indexSupplied(Index<TText, FMIndex<TSpec, TConfig> > & index, FibreSALF const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLF())));
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline bool indexSupplied(Index<TText, FMIndex<TSpec, TConfig> > const & index, FibreSALF const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLF())));
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.Index#open
..class:Class.Index
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

// This function can be used to open a previously saved index.
template <typename TText, typename TSpec, typename TConfig>
inline bool open(Index<TText, FMIndex<TSpec, TConfig> > & index, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;    append(name, ".txt");
    if (!open(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!open(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    setFibre(getFibre(index, FibreSA()), getFibre(index, FibreLF()), FibreLF());

    return true;
}

// This function can be used to open a previously saved index.
template <typename TText, typename TSpec, typename TConfig>
inline bool open(Index<TText, FMIndex<TSpec, TConfig> > & index, const char * fileName)
{
    return open(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/**
.Function.Index#save
..class:Class.Index
..summary:This functions saves an index to disk.
..signature:save(index, fileName [, openMode])
..param.index:The index.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TText, typename TSpec, typename TConfig>
inline bool save(Index<TText, FMIndex<TSpec, TConfig> > const & index, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;    append(name, ".txt");
    if (!save(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!save(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    return true;
}

// This function can be used to save an index on disk.
template <typename TText, typename TSpec, typename TConfig>
inline bool save(Index<TText, FMIndex<TSpec, TConfig> > const & index, const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
}

}
#endif // INDEX_FM_H
