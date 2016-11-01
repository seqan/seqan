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

#ifndef INDEX_BIDIRECTIONAL_H_
#define INDEX_BIDIRECTIONAL_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TIndexSpec>
class BidirectionalIndex;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction RevTextFibre
// ----------------------------------------------------------------------------

/*!
 * @class RevTextFibre
 * @headerfile <seqan/index.h>
 * @brief A helper object that stores the type of the modifier of a given type.
 *
 * @signature template <TText>
 *            struct RevTextFibre;
 */

template <typename TText>
struct RevTextFibre
{
    typedef TText Type;
};

// ----------------------------------------------------------------------------
// Class BidirectionalIndex
// ----------------------------------------------------------------------------

/*!
 * @class BidirectionalIndex
 * @headerfile <seqan/index.h>
 * @brief A bidirectional index class.
 *
 * @signature template <typename TText, typename TIndexSpec>
 *            class Index<TText, Bidirectional<TIndexSpec> >;
 *
 * @tparam TText           The text type. Types: @link String @endlink, @link StringSet @endlink
 * @tparam TIndexSpec   Unidirectional index used for constructing a bidirectional one.
 *
 * @section Structure
 *
 * The BidirectionalIndex consists of two unidirectional indices that index the original text and the reverse text.
 */

template <typename TText, typename TIndexSpec>
class Index<TText, BidirectionalIndex<TIndexSpec> >
{
    typedef typename RevTextFibre<TText>::Type  TRevText;
    typedef Index<TRevText, TIndexSpec>         TRevIndex;
    typedef Index<TText, TIndexSpec>            TFwdIndex;

    public:

    TRevText    revText;
    TRevIndex   rev;
    TFwdIndex   fwd;

    Index()    {}

    Index(TText & text) :
        revText(text),
        rev(revText),
        fwd(text)
    {
        reverse(revText);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TFibre>
inline bool indexCreate(Index<TText, BidirectionalIndex<TIndexSpec> > & index, Tag<TFibre>)
{
    return indexCreate(index.fwd, Tag<TFibre>()) && indexCreate(index.rev, Tag<TFibre>());
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

// Already documented
template <typename TText, typename TIndexSpec>
inline void clear(Index<TText, BidirectionalIndex<TIndexSpec> > & index)
{
    clear(index.fwd);
    clear(index.rev);
    clear(index.revText);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

// This function checks whether the index is empty. Its already documented.
template <typename TText, typename TIndexSpec>
inline bool empty(Index<TText, BidirectionalIndex<TIndexSpec> > const & index)
{
    return empty(index.fwd) && empty(index.rev);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec>
inline bool open(Index<TText, BidirectionalIndex<TIndexSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    if (open(index.fwd, toCString(name), openMode))
    {
        name = fileName;    append(name, ".rev");
        return open(index.rev, toCString(name), openMode);
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec>
inline bool save(Index<TText, BidirectionalIndex<TIndexSpec> > const & index, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    if (save(index.fwd, toCString(name), openMode))
    {
        name = fileName;    append(name, ".rev");
        return save(index.rev, toCString(name), openMode);
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------
// only used for testing open/save

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec>, FibreText>::Type &
getFibre(Index<TText, BidirectionalIndex<TSpec> > &index, FibreText) {
    return value(index.fwd.text);
}

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, TSpec>, const FibreText>::Type &
getFibre(Index<TText, BidirectionalIndex<TSpec> > const &index, FibreText) {
    return value(index.fwd.text);
}

}
#endif /* INDEX_BIDIRECTIONAL_H_ */
