// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Christopher Pockrandt <christopher.pockrandt@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_BIFM_H_
#define INDEX_BIFM_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TSpec = void, typename TConfig = FMIndexConfig<TSpec> >
class BidirectionalFMIndex;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction BiFMReversedText
// ----------------------------------------------------------------------------

/*!
 * @class BiFMReversedText
 * @headerfile <seqan/index.h>
 * @brief A helper object that stores the type of the modifier of a given type @endlink.
 *
 * @signature template <TText>
 *            struct BiFMReversedText;
 */

template <typename TText>
struct BiFMReversedText
{
    typedef ModifiedString<TText, ModReverse> Type;
};

template <typename TText>
struct BiFMReversedText<ModifiedString<TText, ModReverse> >
{
    typedef TText Type;
};

template <typename TText, typename TTextConfig>
struct BiFMReversedText<StringSet<TText, TTextConfig> >
{
    typedef StringSet<typename BiFMReversedText<TText>::Type, TTextConfig> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

/*!
 * @defgroup BidirectionalFMIndexFibres Bidirectional FM Index Fibres
 * @brief Tag to select a specific fibre of a @link BidirectionalFMIndex @endlink.
 *
 * These tags can be used to get @link Fibre Fibres @endlink of a Bidirectional FM index.
 *
 * @see Fibre
 * @see Index#getFibre
 *
 * @tag BidirectionalFMIndexFibres#FibreText
 * @brief The original text of the index.
 *
 * @tag BidirectionalFMIndexFibres#FibreSA
 * @brief The compressed suffix array of the text.
 *
 * @tag BidirectionalFMIndexFibres#FibreLF
 * @brief The lf table.
 */

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, BidirectionalFMIndex<TSpec, TConfig> >, FibreLF>
{
    typedef LF<TText, TSpec, TConfig>  Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, BidirectionalFMIndex<TSpec, TConfig> >, FibreSA>
{
    typedef CompressedSA<TText, TSpec, TConfig> Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<Index<TText, BidirectionalFMIndex<TSpec, TConfig> >, FibreTempSA>
{
    typedef Index<TText, BidirectionalFMIndex<TSpec, TConfig> >          TIndex_;
    typedef typename SAValue<TIndex_>::Type                 TSAValue_;

    // NOTE(esiragusa): External causes problems on device code.
#ifndef PLATFORM_CUDA
    typedef String<TSAValue_, External<ExternalConfigLarge<> > >                Type;
#else
    typedef String<TSAValue_, typename DefaultIndexStringSpec<TText>::Type>     Type;
#endif
};



// ----------------------------------------------------------------------------
// Class BidirectionalFMIndex
// ----------------------------------------------------------------------------

/*!
 * @class BidirectionalFMIndex
 * @extends Index
 * @headerfile <seqan/index.h>
 * @brief A bidirectional index based on the Burrows-Wheeler transform.
 *
 * @signature template <typename TText[, typename TSpec[, typename TConfig]]>
 *            class Index<TText, BidirectionalFMIndex<TSpec, TConfig> >;
 *
 * @tparam TText   The text type. Types: @link String @endlink, @link StringSet @endlink
 * @tparam TSpec   Bidirectional FM index specialisation, defaults to <tt>void</tt>.
 * @tparam TConfig A config object which determines the data types of the different fibres, defaults to
 *                 <tt>FMIndexConfig&lt;TSpec&gt;</tt>. TBidirectional is only used internally to
 *                 distinguish between a stand-alone FM index and an FM index belonging to a pair of
 *                 FM indices forming a bidirectional FM index.
 *
 * @section Structure
 *
 * The Bidirectional FM index consists of two FM indices that index the original text and the reverse text.
 * Performing a regular backward search on one on the FM indices will automatically update the search
 * interval on the other index.
 */

template <typename TText, typename TSpec, typename TSpec2, typename TLengthSum, typename TBidirectional>
class Index<TText, BidirectionalFMIndex<TSpec, FMIndexConfig<TSpec2, TLengthSum, TBidirectional> > >
{
    typedef typename BiFMReversedText<TText>::Type                                                    TRevText;
    typedef Index<TRevText, FMIndex<TSpec, FMIndexConfig<TSpec2, TLengthSum, FMBidirectional> > >     TRevIndex;
    typedef Index<TText, FMIndex<TSpec, FMIndexConfig<TSpec2, TLengthSum, FMBidirectional> > >        TFwdIndex;

    public:

    TRevText    revText;
    TRevIndex        rev;
    TFwdIndex        fwd;

    Index()    {}

    Index(TText & text) :
        revText(text),
        rev(revText),
        fwd(text)
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
inline void clear(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index)
{
    clear(index.fwd);
    clear(index.rev);
    clear(index.revText);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

// This function checks whether the index is empty. Its already documented.
template <typename TText, typename TSpec, typename TConfig>
inline bool empty(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > const & index)
{
    return empty(index.fwd) && empty(index.bwd);
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type &
getFibre(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index, FibreLF /*tag*/)
{
    return index.fwd.lf;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type const &
getFibre(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > const & index, FibreLF /*tag*/)
{
    return index.fwd.lf;
}

// ----------------------------------------------------------------------------
// Function indexLF()
// ----------------------------------------------------------------------------
/*!
 * @fn BidirectionalFMIndex#indexLF
 * @headerfile <seqan/index.h>
 * @brief A shortcut for <tt>getFibre(index, FibreLF())</tt>.
 *
 * @signature TFibre indexLF(index);
 *
 * @param[in] index The FM index.
 *
 * @return TFibre A reference to the @link BidirectionalFMIndexFibres#FibreLF @endlink.
 */

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type &
indexLF(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index)
{
    return getFibre(index.fwd, FibreLF());
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreLF>::Type const &
indexLF(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > const & index)
{
    return getFibre(index.fwd, FibreLF());
}

// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
inline bool indexCreate(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index, FibreSALF)
{
    return indexCreate(index.fwd) && indexCreate(index.bwd);
}

template <typename TText, typename TSpec, typename TConfig>
inline bool indexCreate(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index, FibreSA)
{
    return indexCreate(index.fwd, FibreSALF()) && indexCreate(index.rev, FibreSALF());
}

template <typename TText, typename TSpec, typename TConfig>
inline bool indexCreate(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index)
{
    return indexCreate(index.fwd, FibreSALF()) && indexCreate(index.rev, FibreSALF());
}

// ----------------------------------------------------------------------------
// Function indexSupplied()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline bool indexSupplied(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index, FibreSALF const)
{
    return !(empty(getFibre(index.fwd, FibreSA())) || empty(getFibre(index.fwd, FibreLF()))
            ||
            empty(getFibre(index.rev, FibreSA())) || empty(getFibre(index.rev, FibreLF())));
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline bool indexSupplied(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > const & index, FibreSALF const)
{
    return !(empty(getFibre(index.fwd, FibreSA())) || empty(getFibre(index.fwd, FibreLF()))
            ||
            empty(getFibre(index.rev, FibreSA())) || empty(getFibre(index.rev, FibreLF())));
}

// This function can be used to open a previously saved index.
template <typename TText, typename TSpec, typename TConfig>
inline bool open(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index, const char * fileName)
{
    String<char> name;

    name = fileName;    append(name, ".fwd");
    bool fwdIndex = open(index.fwd, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
    if (fwdIndex)
    {
        name = fileName;    append(name, ".bwd");
        return open(index.rev, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
    }
    return false;
}

// This function can be used to save an index on disk.
template <typename TText, typename TSpec, typename TConfig>
inline bool save(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > const & index, const char * fileName)
{
    String<char> name;

    name = fileName;    append(name, ".fwd");
    bool fwdIndex = save(index.fwd, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
    if (fwdIndex)
    {
        name = fileName;    append(name, ".bwd");
        return save(index.rev, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
    }
    return false;
}

}
#endif /* INDEX_BIFM_H_ */
