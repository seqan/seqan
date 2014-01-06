// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

/*!
 * @defgroup FMIndexRankDictionarySpec FMIndex RankDictionary Specialisations
 * @brief Tags that can be chosen to specify a certain @link RankDictionary @endlink.
 *
 * @tag FMIndexRankDictionarySpec#WT
 * @brief Tag that specifies the @link FMIndex @endlink to use a wavelet tree as the occurrence table.
 *
 * @tag FMIndexRankDictionarySpec#SBM
 * @brief Tag that specifies the @link FMIndex @endlink to use a StringSet of rank support bis strings as the occurrence table.
 *
 */

/*!
 * @defgroup FMIndexCompressionSpec FMIndex Compression Specialisations
 * @brief Tags that can be chosen to specify if the index stored the text or not.
 *
 * @tag FMIndexCompressionSpec#CompressText
 * @brief Tag to select a FM index variant that can be used such that it is not
 *        necessary to store the text after index construction. This index is
 *        very space efficient.
 *
 * @tag FMIndexCompressionSpec#void
 * @brief The text is kept and not cleared. This FM Index version is faster but more memory is required.
 */


/*!
 * @fn FMIndex#begin
 * @headerfile seqan/index.h
 * @brief Returns an iterator pointing to the root node of the virtual prefix trie of the reversed text of the index.
 * 
 * @signature TIterator begin(index, tag);
 * 
 * @param[in] index The index to be traversed.
 * @param[in] tag   The specialisation of the iterator to be returned by the function. Types: VSTree Iterator
 * 
 * @return TIterator Returns an iterator pointing to the root node of the virtual prefix trie of the reversed text
 *                   of the the index.  Types: <tt>The result of Iterator&lt;Index&lt;TText, TIndexSpec&gt;,
 *                   TSpec&gt;::Type</tt>
 */

/*!
 * @fn LF#lfMapping
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position of an character at a specified position of L in
 *        F. L corresponds to the last column of the sorted cyclic rotations of
 *        the original text, while F correspond to the first column.
 * 
 * @signature lfMapping(lfTable, pos)
 * 
 * @param lfTable The @link LF @endlink holding the occurrence and prefix
 *                sum table.
 *
 * @param pos The position in L. Types: @link UnsignedIntegerConcept @endlink
 * 
 * @return TReturn Returns the position of the character L[c] in F. The returned
 *                 position is of the same type as pos. Types: The type of the position.
 */


/*!
 * @defgroup PrefixSumTableFibres PrefixSumTable Fibres
 * 
 * @brief Tag to select a specific fibre of a @link PrefixSumTableFibres @endlink.
 * 
 * @see Fibre
 * @see Index#getFibre
 * 
 * @tag PrefixSumTableFibres#FibreEntries
 * @brief The entries in the prefix sum table.
 */

/*!
 * @class PrefixSumTable
 * @headerfile seqan/Index.h
 * @brief The prefix-sum table is a data structure which stores for each
 *        character the number of smaller lexicographic smaller characters in a
 *        given text.
 * 
 * @signature template <typename TChar, typename TSpec>
 *            class PrefixSumTable;
 * 
 * @tparam TSpec A specialisation tag, defaults to <tt>void</tt>.
 * @tparam TChar The character type
 */

/*!
 * @fn PrefixSumTable#clear
 * @headerfile seqan/index.h
 * @brief Resets the prefix sum table.
 * 
 * @signature void clear(prefixSumTable);
 *
 * @param[in,out] prefixSumTable The prefix sum table to be cleared.
 */

/*!
 * @fn PrefixSumTable#createPrefixSumTable
 * @headerfile seqan/index.h
 * @brief Creates the prefix sum table
 * 
 * @signature void createPrefixSumTable(prefixSumTable, text);
 * 
 * @param[in,out] prefixSumTable The prefix sum table to be constructed.
 * @param[in,out] text           The underlying text. Types: @link String @endlink, @link StringSet @endlink
 */

/*!
 * @fn PrefixSumTable#getAlphabetSize
 * @headerfile seqan/index.h
 * @brief Returns the number of different characters in the prefix sum table.
 * 
 * @signature TSize getAlphabetSize(prefixSumTable);
 * 
 * @param[in] prefixSumTable A prefix sum table.
 *
 * @return TSize The alphabet size.
 */

/*!
 * @fn PrefixSumTable#getCharacterPosition
 * @headerfile seqan/index.h
 * @brief Returns the position of a given character within the prefix sum table.
 * 
 * @signature TPos getCharacterPosition(prefixSumTable, character);
 * 
 * @param[in] character      A character.
 * @param[in] prefixSumTable A prefix sum table.
 *
 * @return TPos The position to return.
 */

/*!
 * @fn PrefixSumTable#getCharacter
 * @headerfile seqan/index.h
 * @brief Returns the character of a given position within the prefix sum table.
 * 
 * @signature TChar getCharacter(prefixSumTable, pos);
 * 
 * @param[in] prefixSumTable A prefix sum table.
 * @param[in] pos            A position. Types: @link UnsignedIntegerConcept @endlink
 */

/*!
 * @fn PrefixSumTable#getPrefixSum
 * @headerfile seqan/index.h
 * @brief Returns the prefix sum of a given position.
 * 
 * @signature TLength getPrefixSum(prefixSumTable, pos);
 * 
 * @param[in] prefixSumTable A prefix sum table.
 * @param[in] pos            A position. Types: @link UnsignedIntegerConcept @endlink.
 *
 * @return TLength The prefix sum of a given position.
 */

/*!
 * @fn PrefixSumTable#getValue
 * @headerfile seqan/index.h
 * @brief Returns the prefix sum of a given position.
 * 
 * @signature TValue getValue(prefixSumTable, pos);
 * 
 * @param[in] prefixSumTable A prefix sum table.
 * @param[in] pos            A position. Types: @link UnsignedIntegerConcept @endlink
 *
 * @return TValue The retrieved value.
 */

/*!
 * @fn PrefixSumTable#getFibre
 * @headerfile seqan/index.h
 * @brief Returns a specific fibre of a prefix-sum table.
 * 
 * @signature TFibre getFibre(prefixSumTable, fibreTag);
 * 
 * @param[in] fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                     @link PrefixSumTableFibres PrefixSumTanble Fibres @endlink.
 * @param[in] prefixSumTable The prefix sum table.
 * 
 * @return TFibre A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn PrefixSumTable#length
 * @headerfile seqan/index.h
 * @brief Returns the number of different characters in the prefix-sum table.
 * 
 * @signature TSize length(lfTable);
 * 
 * @param[in] lfTable The prefix-sum table.
 * 
 * @return TSSize Returns the number of different characters in the prefix-sum
 *                 table.  If the type of the characters of the prefix-sum table
 *                 consists of more than 8 bit only the characters actually
 *                 occurring in the original text are accounted for when calling
 *                 length.  Types: @link Size @endlink of the prefix-sum table.
 */

/*!
 * @fn PrefixSumTable#prefixSum
 * @headerfile seqan/index.h
 * 
 * @brief Returns a reference to the entry of the prefix sum table of a given
 *        position.
 * 
 * @signature TSize prefixSum(prefixSumTable, pos);
 * 
 * @param[in] pos            A position. Types: @link UnsignedIntegerConcept @endlink.
 * @param[in] prefixSumTable A prefix sum table.
 */

/*!
 * @fn PrefixSumTable#resize
 * @headerfile seqan/index.h
 * 
 * @brief Resize the prefix-sum table to be able to store more or less
 *        characters.
 * 
 * @signature void resize(prefixSumTable, size[, value, resizeTag]);
 * 
 * @param[in,out] prefixSumTable A prefix sum table. Types: PrefixSumTable
 * @param[in]     size           The new size. Types: @link UnsignedIntegerConcept @endlink
 * @param[in]     value          The value to be used to initialize the new storage.
 * @param[in] resizeTag         Specifies the strategy that is applied if the capacity of
 *                              <tt>object</tt> is less than <tt>newLength</tt>. (optional)
 *                              Types: @link OverflowStrategyTags @endlink Default: Specified by @link
 *                              DefaultOverflowExplicit @endlink.
 */

// TODO(holtgrew): The following seems wrong.
/*!
 * @fn PrefixSumTable#setPrefixSum
 * @headerfile seqan/index.h
 * 
 * @brief Returns a reference to the entry of the prefix-sum table of a given position.
 * 
 * @signature void setPrefixSum(prefixSumTable, value, pos);
 * 
 * @param[in,out] prefixSumTable A prefix sum table.
 * @param[in]     value          A specified value to be inserted.
 * @param[in]     pos            A position. Types: @link UnsignedIntegerConcept @endlink
 */

/*!
 * @fn PrefixSumTable#open
 * @headerfile seqan/index.h
 * @brief This functions loads a prefix-sum table from disk.
 * 
 * @signature bool open(prefixSumTable, fileName[, openMode]);
 * 
 * @param[in,out] prefixSumTable The prefix-sum table.
 * @param[in]     fileName       C-style character string containing the file name.
 * @param[in]     openMode       The combination of flags defining how the file should be
 *                               opened.  To open a file read-only, write-only or to read and
 *                               write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                               <tt>OPEN_RDWR</tt>.  To create or overwrite a file add
 *                               <tt>OPEN_CREATE</tt>.  To append a file if existing add
 *                               <tt>OPEN_APPEND</tt>.  To circumvent problems, files are always
 *                               opened in binary mode.  Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                               OPEN_APPEND</tt>.
 * 
 * @return bool <tt>true</tt> on success.
 */

/*!
 * @fn PrefixSumTable#save
 * @headerfile seqan/index.h
 * @brief This functions saves a prefix-sum table to disk.
 * 
 * @signature bool save(prefixSumTable, fileName [, openMode]);
 * 
 * @param[in] prefixSumTable The prefix-sum table.
 * @param[in] fileName       C-style character string containing the file name.
 * @param[in] openMode       The combination of flags defining how the file should be
 *                           opened.  To open a file read-only, write-only or to read and
 *                           write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                           <tt>OPEN_RDWR</tt>.  To create or overwrite a file add
 *                           <tt>OPEN_CREATE</tt>.  To append a file if existing add
 *                           <tt>OPEN_APPEND</tt>.  To circumvent problems, files are always
 *                           opened in binary mode.  Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                           OPEN_APPEND</tt>.
 * 
 * @return bool <tt>true</tt> on success.
 */









/*!
 * @class SequenceBitMask
 *
 * @extends RankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The string set bit string dictionary is a string set of rank support
 *        bit strings for constant time acces of the rank of a specified
 *        character at a specified position.
 * 
 * @signature template <typename TValue>
 *            class RankDictionary<SequenceBitMask<TValue> >;
 * 
 * @tparam TValue The value type of the text.
 * 
 * This data structure is optimized for very small alphabets, such as @link Dna @endlink
 * or  @link Dna5 @endlink.  Consider using a @link WaveletTree @endlink if
 * your alphabet size is larger.
 */

/*!
 * @defgroup SequenceBitMaskFibres SequenceBitMask Fibres
 * @brief Tag to select a specific fibre of a SequenceBitMask.
 * 
 * @see Fibre
 * @see Index#getFibre
 *
 * @tag SequenceBitMaskFibres#FibreBitStrings
 * @brief The string set containing a bit string for each character.
 */


// TODO(holtgrew): Missing return type.
/*!
 * @fn RankDictionary#getValue
 * @headerfile seqan/index.h
 * @brief Returns the character of a specified position.
 * 
 * @signature TChar getValue(dictionary, pos);
 * 
 * @param[in] dictionary The dictionary.
 * @param[in] pos        The position. Types: @link UnsignedIntegerConcept @endlink.
 */


/*!
 * @fn RankDictionary#countOccurrences
 * @headerfile seqan/index.h
 * @brief Returns the rank (number of occurrences) of a specified character up to a specified position.
 * 
 * @signature unsigned countOccurrences(dictionary, character, pos);
 * 
 * @param[in] dictionary The dictionary. 
 * @param[in] character The character of interest.
 * @param[in] pos The position (which is also included in the rank computation).
 *
 * @return unsigned The rank (number of occurrences) of a specified character up to a specified position.
 */

/*!
 * @fn RankDictionary#createRankDictionary
 * @headerfile seqan/index.h
 * @brief This functions creates the dictionary.
 * 
 * @signature void createRankDictionary(dictionary, text);
 * 
 * @param[in]  text       A text to be transfered into a wavelet tree. Types: @link String @endlink
 * @param[out] dictionary The dictionary.
 */





/*!
 * @defgroup RankSupportBitStringFibres RankSupportBitString Fibres
 * 
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link
 *        RankSupportBitString @endlink.
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a rank support
 * bit string.
 * 
 * @see Fibre
 * @see Index#getFibre
 * 
 * @tag RankSupportBitStringFibres#FibreSuperBlocks
 * @brief The super block string.
 * 
 * @tag RankSupportBitStringFibres#FibreBits
 * @brief The bit string.
 * 
 * @tag RankSupportBitStringFibres#FibreBlocks
 * @brief The block string.
 */

/*!
 * @class RankSupportBitString
 * @headerfile seqan/index.h
 * @brief A bit string supporting rank queries in constant time.
 * 
 * @signature template <typename TSpec>
 *            class RankSupportBitString;
 * 
 * @tparam TSpec Specialisation tag.  Default: void
 * 
 * The constant rank query time is achieved by evaluating precomputed
 * subsolutions.  In order to do so, the bit string is divided into blocks of
 * length l.  A super block string stores for each block of l blocks the number
 * of bits set from the beginning. In addition a block string stores the number
 * of bits set in each block from the start of the last super block block.
 * Therefore it is possible to compute the result of a rank query in constant
 * time by adding information from the bit, block and super block string.
 */

/*!
 * @fn RankSupportBitString#appendValue
 * 
 * @headerfile seqan/sequence.h
 * 
 * @brief Appends a bit to a @link RankSupportBitString @endlink.
 * 
 * @signature appendValue(bitString, bit)
 * 
 * @param target A container. Types: RankSupportBitString
 * @param bit Value that is appended to <tt>target</tt>.If the value is
 *              different from 0 it is interpreted as 1. Types: @link UnsignedIntegerConcept @endlink, bool
 */

/*!
 * @fn RankSupportBitString#clear
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets an rank-support-bit string.
 * 
 * @signature clear(bitString)
 * 
 * @param bitString The bit string to be cleared.
 */

/*!
 * @fn RankSupportBitString#getRank
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the rank (the number of bits set from the start of the bit
 *        string) of a specified position.
 * 
 * @signature getRank(bitString, pos)
 * 
 * @param bitString The bit string. Types: @link RankSupportBitString @endlink
 * @param pos Position of a bit.
 * 
 * @return TReturn @link Value @endlink of @link Fibre @endlink of the rank-support-bit string.
 */

/*!
 * @fn RankSupportBitString#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Checks whether or not a rank-support-bit string contains any elements.
 * 
 * @signature bool empty(bitString)
 * 
 * @param bitString The rank-support-bit string to be checked.
 *
 * @return bool Returns true if the rank-support-bit string is empty and false otherwise.
 */

/*!
 * @fn RankSupportBitString#isBitSet
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether the bit with the given index is set to 1.
 * 
 * @signature isSetBit(bitString, pos)
 * 
 * @param bitString The bit string.
 * @param pos Position of the bit. Types: @link UnsignedIntegerConcept @endlink
 * 
 * @return TReturn Returns whether a specified bit is set or not.
 */

/*!
 * @fn RankSupportBitString#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a rank-support-bit string.
 * 
 * @signature getFibre(bitString, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link RankSupportBitStringFibres RankSupportBitString Fibres @endlink.
 * @param bitString The rank-support-bit string holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn RankSupportBitString#length
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of bits in the rank-support-bit string.
 * 
 * @signature length(bitString)
 * 
 * @param bitString The rank-support-bit string.
 * 
 * @return TReturn Types: @link Value @endlink of @link Fibre @endlink of the rank-support-bit string.
 */

/*!
 * @fn RankSupportBitString#resize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Resets the number of bits in the rank-support-bit string.
 * 
 * @signature TSize resize(bitString, newLenght)
 * 
 * @param bitString The rank-support-bit string.
 * @param newLength The number of elements which should be stored in the compressed suffix array. Types: @link
 * UnsignedIntegerConcept @endlink.
 *
 * @return TSize The number of elements in the rank-support-bit string. Types: The result of @link Size @endlink of the
 * rank-support-bit string.
 *
 * @section Note If the new length is smaller than the actual one then the last <tt>x<tt> items of the rank-support-bit string
 * are deleted with x = oldLength - newLength.
 */

/*!
 * @fn RankSupportBitString#setBitTo
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Set the bit with the given position to the given value.
 * 
 * @signature setBitTo(bitString, pos, value)
 * 
 * @param pos Position of the bit. Types: @link UnsignedIntegerConcept @endlink
 * @param bitString The bit string. 
 * @param bit The value of the bit. Note that values different from 0 are
 *            interpreted as 1.
 * 
 * @return TReturn <tt>void</tt>
 * 
 * @section Examples
 * 
 * @see RankSupportBitString#isBitSet
 */

/*!
 * @fn RankSupportBitString#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions opens a @link RankSupportBitString @endlink from disk.
 * 
 * @signature open(bitString, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param bitString The bit string to be opened.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn RankSupportBitString#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a @link RankSupportBitString @endlink to disk.
 * 
 * @signature save(bitString, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param bitString The bit string to be saved.
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */















































/*!
 * @defgroup SentinelRankDictionaryFibres SentinelRankDictionary Fibres
 * 
 * @brief Tag to select a specific fibre of a @link
 *        SentinelRankDictionary @endlink.
 *  
 * @see Fibre
 * @see Index#getFibre
 * 
 * @tag SentinelRankDictionaryFibres#FibreSentinelPosition
 * @brief The bit string encoding the position of the sentinel sign.
 * 
 * @tag SentinelRankDictionaryFibres#FibreRankDictionary
 * @brief The rank dictionary.
 */

/*!
 * @class SentinelRankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A rank dictionary, additional storing sentinel character which are not
 *        accounted for in a rank query.
 * 
 * @signature template <typename TRankDictionary, typename TSpec>
 *            class SentinelRankDictionary;
 * 
 * @tparam TSpec Specialisation
 * @tparam TRankDictionary The rank dictionary of a text.
 */

/*!
 * @fn SentinelRankDictionary#clear
 * @headerfile seqan/index.h
 * 
 * @brief Clears the dictionary.
 * 
 * @signature clear(dictionary)
 * 
 * @param dictionary The rank dictionary to be cleared.
 */

/*!
 * @fn SentinelRankDictionary#sentinelPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether a specified position is a sentinel position.
 * 
 * @signature sentinelPosition(dictionary, pos)
 * 
 * @param pos The position. Types: @link UnsignedIntegerConcept @endlink
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#empty
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns whether or not the dictionary is empty.
 * 
 * @signature empty(dictionary)
 * 
 * @param dictionary The rank dictionary to be checked. 
 */

/*!
 * @fn SentinelRankDictionary#getValue
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the character of a specified position.
 * 
 * @signature getCharacter(dictionary, pos)
 * 
 * @param pos The position
 * @param dictionary The rank dictionary.
 */

/*!
 * @fn SentinelRankDictionary#getFibre
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns a specific fibre of a dictionary.
 * 
 * @signature getFibre(dictionary, fibreTag)
 * 
 * @param fibreTag A tag that identifies the @link Fibre @endlink. Types:
 *                 @link SentinelRankDictionaryFibres SentinelRankDictionary Fibres @endlink.
 * @param dictionary The dictionary holding the fibre.
 * 
 * @return TReturn A reference to the @link Fibre @endlink object.
 */

/*!
 * @fn SentinelRankDictionary#countOccurrences
 * @headerfile seqan/index.h
 * @brief Returns the number of occurrences of a specified character from the start to a specified position.
 * 
 * @signature TSize countOccurrences(dictionary, character, pos);
 * 
 * @param[in] dictionary The dictionary.
 * @param[in] character  The character.
 * @param[in] pos        The position (which is included in the counting).
 *
 * @return TSize The number of occurences  (Metafunction: @link Index#Size @endlink).
 */

/*!
 * @fn SentinelRankDictionary#getSentinelSubstitute
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the character used to substitute the sentinel sign.
 * 
 * @signature getSentinelSubstitute(dictionary)
 * 
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#setSentinelSubstitute
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the character used to substitute the sentinel sign.
 * 
 * @signature setSentinelSubstitute(dictionary, character)
 * 
 * @param character The sentinel substitute.
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#setSentinelPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the sentinel position..
 * 
 * @signature setSentinelPosition(dictionary, pos)
 * 
 * @param pos The sentinel position. Types: @link UnsignedIntegerConcept @endlink
 * @param dictionary The dictionary.
 */

/*!
 * @fn SentinelRankDictionary#createSentinelRankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions creates the dictionary structure.
 * 
 * @signature void createSentinelRankDictionary(dictionary, text)
 * 
 * @param text A text to be transfered into a dictionary. Types: @link String @endlink
 * @param dictionary The dictionary. 
 */

/*!
 * @fn SentinelRankDictionary#save
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions saves a dictionary to disk.
 * 
 * @signature save(dictionary, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param dictionary The dictionary. Types: SentinelRankDictionary
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */

/*!
 * @fn SentinelRankDictionary#open
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This functions loads a dictionary from disk.
 * 
 * @signature open(dictionary, fileName [, openMode])
 * 
 * @param openMode The combination of flags defining how the file should be
 *                 opened.To open a file read-only, write-only or to read and
 *                 write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                 <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                 <tt>OPEN_CREATE</tt>.To append a file if existing add
 *                 <tt>OPEN_APPEND</tt>.To circumvent problems, files are always
 *                 opened in binary mode. Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                 OPEN_APPEND</tt>
 * @param dictionary The dictionary. Types: SentinelRankDictionary
 * @param fileName C-style character string containing the file name.
 * 
 * @return TReturn A <tt>bool</tt> which is <tt>true</tt> on success.
 */
