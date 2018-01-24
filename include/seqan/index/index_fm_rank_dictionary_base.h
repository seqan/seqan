// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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

#ifndef INDEX_FM_RANK_DICTIONARY_BASE_H_
#define INDEX_FM_RANK_DICTIONARY_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TText, typename TSpec, typename TConfig>
struct LF;

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag FibreRanks
// ----------------------------------------------------------------------------

/*!
 * @defgroup RankDictionaryFibres RankDictionary Fibres
 * @brief Tag to select a specific fibre of a SequenceBitMask.
 *
 * @see Fibre
 * @see RankDictionary#getFibre
 *
 * @tag RankDictionaryFibres#FibreRanks
 * @brief The string set containing a bit string for each character.
 *
 * @tag RankDictionaryFibres#FibreSuperBlocks
 * @brief The string set containing a bit string for each character.
 *
 * @tag RankDictionaryFibres#FibreUltraBlocks
 * @brief The string set containing a bit string for each character.
 */

struct FibreRanks_;
typedef Tag<FibreRanks_> const FibreRanks;

struct FibreSuperBlocks_;
typedef Tag<FibreSuperBlocks_> const FibreSuperBlocks;

struct FibreUltraBlocks_;
typedef Tag<FibreUltraBlocks_> const FibreUltraBlocks;

// ----------------------------------------------------------------------------
// Tag RDConfig
// ----------------------------------------------------------------------------

template <
    typename TSize = size_t,
    typename TFibre = Alloc<>,
    unsigned LEVELS_ = 1,
    unsigned WORDS_PER_BLOCK_ = 1
>
struct RDConfig
{
    typedef TSize   Size;
    typedef TFibre  Fibre;

    static const unsigned WORDS_PER_BLOCK = WORDS_PER_BLOCK_;
    static const unsigned LEVELS          = LEVELS_;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class RankDictionary
// ----------------------------------------------------------------------------

/*!
 * @class RankDictionary
 * @headerfile <seqan/index.h>
 * @brief A rank dictionary is a data structure storing the rank of an element in a sequence at every position of
 *        the sequence. Only supports sequences containing @link SimpleType @endlink or bool.
 *
 * @signature template <typename TValue, typename TSpec>
 *            struct RankDictionary;
 *
 * @tparam TSpec The rank dictionary specialisation. Default: @link WaveletTree @endlink, @link Levels @endlink.
 */
template <typename TValue, typename TSpec>
struct RankDictionary;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultIndexStringSpec
// ----------------------------------------------------------------------------

template <typename TValue, template <typename, typename> class TRankDictionary, typename TSpec, typename TConfig>
struct DefaultIndexStringSpec<RankDictionary<TValue, TRankDictionary<TSpec, TConfig> > >
{
    typedef typename TConfig::Fibre Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue, template <typename, typename> class TRankDictionary, typename TSpec, typename TConfig>
struct Size<RankDictionary<TValue, TRankDictionary<TSpec, TConfig> > >
{
    typedef typename TConfig::Size  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<RankDictionary<TValue, TSpec> >
{
    typedef TValue  Type;
};

template <typename TValue, typename TSpec>
struct Value<RankDictionary<TValue, TSpec> const> :
    Value<RankDictionary<TValue, TSpec> > {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/*!
 * @fn RankDictionary#getFibre
 * @headerfile <seqan/index.h>
 * @brief Returns a specific fibre of a dictionary.
 *
 * @signature TFibre getFibre(dictionary, fibreTag);
 *
 * @param[in] dictionary The dictionary holding the fibre.
 * @param[in] fibreTag   A tag that identifies the @link Fibre @endlink. Types: @link RankDictionaryFibres @endlink.
 *
 * @return TFibre A reference to the @link Fibre @endlink object.
 */


template <typename TValue, typename TSpec>
inline typename Fibre<RankDictionary<TValue, TSpec>, FibreRanks>::Type &
getFibre(RankDictionary<TValue, TSpec> & dict, FibreRanks)
{
    return dict.ranks;
}

template <typename TValue, typename TSpec>
inline typename Fibre<RankDictionary<TValue, TSpec>, FibreRanks>::Type const &
getFibre(RankDictionary<TValue, TSpec> const & dict, FibreRanks)
{
    return dict.ranks;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn RankDictionary#clear
 * @headerfile <seqan/index.h>
 * @brief Resets the rank dictionary.
 *
 * @signature void clear(dictionary);
 *
 * @param[in,out] dictionary The rank dictionary to be cleared.
 */

template <typename TValue, typename TSpec>
inline void clear(RankDictionary<TValue, TSpec> & dict)
{
    clear(getFibre(dict, FibreRanks()));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn RankDictionary#empty
 * @headerfile <seqan/index.h>
 * @brief Returns whether or not the rank dictionary is empty.
 *
 * @signature bool empty(dictionary);
 *
 * @param[in] dictionary The rank dictionary to be checked.
 *
 * @return bool <tt>true</tt> if the dictionary is empty, <tt>false</tt> otherwise.
 */

template <typename TValue, typename TSpec>
inline bool empty(RankDictionary<TValue, TSpec> const & dict)
{
    return empty(getFibre(dict, FibreRanks()));
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()
// ----------------------------------------------------------------------------
/*!
 * @fn RankDictionary#createRankDictionary
 * @headerfile <seqan/index.h>
 * @brief This functions creates the dictionary.
 *
 * @signature void createRankDictionary(dictionary, text);
 *
 * @param[in]  text       A text to be transfered into a rank dictionary. Types: @link ContainerConcept @endlink
 * @param[out] dictionary The dictionary.
 */

template <typename TValue, typename TSpec, typename TText>
inline void
createRankDictionary(RankDictionary<TValue, TSpec> & dict, TText const & text)
{
    typedef typename Iterator<TText const, Standard>::Type      TTextIterator;

    // Resize the RankDictionary.
    resize(dict, length(text), Exact());

    // Assign the text value by value.
    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());
    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
        setValue(dict, textIt - textBegin, *textIt);

    // Update all ranks.
    updateRanks(dict);
}

template <typename TText, typename TSpec, typename TConfig, typename TOtherText, typename TSA>
inline void
createRankDictionary(LF<TText, TSpec, TConfig> & lf, TOtherText const & text, TSA const & sa)
{
    typedef typename GetValue<TSA>::Type                    TSAValue;
    typedef typename Size<TSA>::Type                        TSize;

    // Resize the RankDictionary.
    resize(lf.bwt, lengthSum(text) + 1, Exact());

    // Assign the text value by value.
    setValue(lf.bwt, 0, back(text));

    for (TSize i = 0; i < length(sa); ++i)
    {
        TSAValue pos = sa[i];

        if (pos != 0)
        {
            setValue(lf.bwt, i + 1, getValue(text, pos - 1));
        }
        else
        {
            setValue(lf.bwt, i + 1, lf.sentinelSubstitute);
            lf.sentinels = i + 1;
        }
    }

   // Update all ranks.
   updateRanks(lf.bwt);
}

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig, typename TOtherText, typename TSA>
inline void
createRankDictionary(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> & lf, TOtherText const & text, TSA const & sa)
{
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;

    // Resize the RankDictionary.
    TSize seqNum = countSequences(text);
    TSize totalLen = lengthSum(text);
    resize(lf.sentinels, seqNum + totalLen, Exact());
    resize(lf.bwt, seqNum + totalLen, Exact());

    // Fill the sentinel positions (they are all at the beginning of the bwt).
    for (TSize i = 0; i < seqNum; ++i)
    {
        if (length(text[seqNum - (i + 1)]) > 0)
        {
            setValue(lf.bwt, i, back(text[seqNum - (i + 1)]));
            setValue(lf.sentinels, i, false);
        }
    }

    // Compute the rest of the bwt.
    for (TSize i = 0; i < length(sa); ++i)
    {
        TSAValue pos;    // = SA[i];
        posLocalize(pos, sa[i], stringSetLimits(text));

        if (getSeqOffset(pos) != 0)
        {
            setValue(lf.bwt, i + seqNum, getValue(getValue(text, getSeqNo(pos)), getSeqOffset(pos) - 1));
            setValue(lf.sentinels, i + seqNum, false);
        }
        else
        {
            setValue(lf.bwt, i + seqNum, lf.sentinelSubstitute);
            setValue(lf.sentinels, i + seqNum, true);
        }
    }

   // Update all ranks.
   updateRanks(lf.bwt);
   // Update the auxiliary RankDictionary of sentinel positions.
   updateRanks(lf.sentinels);
}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------
/*!
 * @fn RankDictionary#getRank
 * @headerfile <seqan/index.h>
 * @brief Returns the rank of a specified character up to a specified position.
 *
 * @signature TSize countOccurrences(dictionary, pos[, character]);
 *
 * @param[in] dictionary The dictionary.
 * @param[in] pos The position (which is also included in the rank computation).
 * @param[in] character The character of interest. Default: <tt>true</tt>
 *
 * @return TSize The rank of a specified character up to a specified position.
 *               The result of the metafunction Size&lt;RankDictionary&gt;::Type
 */


// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------
/*!
 * @fn RankDictionary#getValue
 * @headerfile <seqan/index.h>
 * @brief Returns the character of a specified position.
 *
 * @signature TValue getValue(dictionary, pos);
 *
 * @param[in] dictionary The dictionary.
 * @param[in] pos        The position. Types: @link UnsignedIntegerConcept @endlink.
 *
 * @return TValue Returns the character of a specified position. The result is of type
 *                Value&lt;RankDictionary&gt;::Type;
 */

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

/*!
 * @fn RankDictionary#setValue
 * @headerfile <seqan/index.h>
 * @brief Sets the character at a specified position.
 *
 * @signature void setValue(dictionary, pos, character);
 *
 * @param[in] dictionary The dictionary.
 * @param[in] pos        The position. Types: @link UnsignedIntegerConcept @endlink.
 * @param[in] character  The character to be set.
 */

// ----------------------------------------------------------------------------
// Function updateRanks()
// ----------------------------------------------------------------------------

/*!
 * @fn RankDictionary#updateRanks
 *
 * @brief Updates the rank information.
 *
 * @signature void updateRanks(dict)
 *
 * @param dict The @link RankDictionary @endlink.
 */

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn RankDictionary#length
 *
 * @brief Retruns the number of element in the rank dictionary.
 *
 * @signature TSize length(dict)
 *
 * @param dict The @link RankDictionary @endlink.
 *
 * @return TSize The number of element in the rank dictionary. The result of the metafunction
 *               Size&lt;RankDictionary&gt;::Type
 */


// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------
/*!
 * @fn RankDictionary#reserve
 *
 * @brief Reserves memory for a specified number of elements in the rank dictionary
 *
 * @signature TSize reserve(dict, newCapacity, tag)
 *
 * @param dict The @link RankDictionary @endlink.
 * @param newCapacity The number of elements one wants to reserve memory for.
 * @param tag The tag that specifies the memory allocation strategy. Types: @link OverflowStrategyTags @endlink.
 *
 * @section Remarks
 *
 * This operation does not changes the content of object.
 *
 * This operation may invalidate iterators of object.
 *
 * @return TSize Retruns the number of element that could be reserved. The result of the metafunction
 *               Size&lt;RankDictionary&gt;::Type
 */

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------
/*!
 * @fn RankDictionary#resize
 *
 * @brief Increases the length of a rank dictionary to a specified value.
 *
 * @signature TSize resize(dict, newLength, tag)
 *
 * @param dict The @link RankDictionary @endlink.
 * @param newLength The number of elements the rank dictionaru contains after resizing.
 * @param tag The tag that specifies the memory allocation strategy. Types: @link OverflowStrategyTags @endlink.
 *
 * @section Remarks
 *
 * This operation may invalidate iterators of object.
 *
 * @return TSize Retruns the number of element that could be reserved. The result of the metafunction
 *               Size&lt;RankDictionary&gt;::Type
 */

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn RankDictionary#open
 * @headerfile <seqan/index.h>
 * @brief This functions loads a rank dictionary from disk.
 *
 * @signature bool open(dictionary, fileName[, openMode]);
 *
 * @param[in,out] dictionary The dictionary.
 * @param[in]     fileName   C-style character string containing the file name.
 * @param[in]     openMode   The combination of flags defining how the file should be opened.To open a file read-only,
 *                           write-only or to read and write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                           <tt>OPEN_RDWR</tt>.To create or overwrite a file add <tt>OPEN_CREATE</tt>.To append a file
 *                           if existing add <tt>OPEN_APPEND</tt>.To circumvent problems, files are always opened in
 *                           binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>
 *
 * @return bool <tt>true</tt> on success.
 */

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<TValue, TSpec> & dict, const char * fileName, int openMode)
{
    return open(getFibre(dict, FibreRanks()), fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<TValue, TSpec> & dict, const char * fileName)
{
    return open(dict, fileName, DefaultOpenMode<RankDictionary<TValue, TSpec> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------
/*!
 * @fn RankDictionary#save
 * @headerfile <seqan/index.h>
 * @brief This functions saves a dictionary to disk.
 *
 * @signature bool save(dictionary, fileName[, openMode]);
 *
 * @param[in] dictionary The dictionary.
 * @param[in] fileName   C-style character string containing the file name.
 * @param[in] openMode   The combination of flags defining how the file should be opened.To open a file read-only,
 *                       write-only or to read and write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                       <tt>OPEN_RDWR</tt>.To create or overwrite a file add <tt>OPEN_CREATE</tt>.To append a file
 *                       if existing add <tt>OPEN_APPEND</tt>.To circumvent problems, files are always opened in
 *                       binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>
 *
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */
template <typename TValue, typename TSpec>
inline bool save(RankDictionary<TValue, TSpec> const & dict, const char * fileName, int openMode)
{
    return save(getFibre(dict, FibreRanks()), fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<TValue, TSpec> const & dict, const char * fileName)
{
    return save(dict, fileName, DefaultOpenMode<RankDictionary<TValue, TSpec> >::VALUE);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_BASE_H_
