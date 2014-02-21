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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Code for read/write access to BAM tag dicts.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

class BamTagsDict;

inline bool hasIndex(BamTagsDict const & bamTags);
inline void buildIndex(BamTagsDict const & bamTags);

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <>
struct Host<BamTagsDict>
{
    typedef CharString Type;
};

template <>
struct Host<BamTagsDict const>
{
    typedef CharString const Type;
};

template <>
struct Size<BamTagsDict>
{
    typedef unsigned Type;
};

template <>
struct Position<BamTagsDict>
{
    typedef unsigned Type;
};


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class BamTagsDict
 * @headerfile <seqan/bam_io.h>
 * @brief Indexes start positions of BAM tags in a @link CharString @endlink and provides a dict-like API.
 *
 * @signature class BamTagsDict;
 *
 * @section Example
 *
 * @include demos/bam_io/bam_tags_dict.cpp
 *
 * Output is:
 *
 * @include demos/bam_io/bam_tags_dict.cpp.stdout
 *
 * @see getBamTypeSize
 * @see getBamTypeChar
 */

/*!
 * @fn BamTagsDict::BamTagsDict
 * @brief Constructor
 *
 * @signature BamTagsDict::BamTagsDict();
 */

/**
.Class.BamTagsDict
..cat:BAM I/O
..cat:Fragment Store
..signature:BamTagsDict
..summary:Indexes start positions of BAM tags in a @Shortcut.CharString@ and provides a dict-like API.
..example.code:
CharString samStr = "AA:Z:value1\tAB:Z:value2\tAC:i:30";
CharString bamStr;
assignSamToBam(bamStr, samStr);
BamTagsDict tags(bamStr);
std::cerr << length(tags) << std::endl;  // #=> "3"
for (unsigned i = 0; i < length(tags); ++i)
{
    std::cerr << getTagKey(tags, i) << " -> " << getTagValue(tags, i) << std::endl;
    if (getTagValue(tags, i)[0] == 'i')  // is 32 bit integer
    {
        __int32 x = 0;
        bool res = extractTagValue(x, tags, i);
        SEQAN_ASSERT_MSG(res, "Not a valid integer at pos %u!", i);
        std::cerr << "     " << x << std::endl;
    }
}
// #=> "AA -> Zvalue1"
// #=> "AB -> Zvalue2"
// #-> "AC -> i<binary representation of 30>"
#  #-> "      30"
..include:seqan/bam_io.h

.Memfunc.BamTagsDict#BamTagsDict
..class:Class.BamTagsDict
..signature:BamTagsDict()
..summary:Constructor.
..remarks:Only the default constructor is provided.
*/

class BamTagsDict
{
    typedef Host<BamTagsDict>::Type TBamTagsSequence;
    typedef Position<TBamTagsSequence>::Type TPos;

public:
    Holder<TBamTagsSequence> _host;
    mutable String<TPos> _positions;

    BamTagsDict() {}
    BamTagsDict(TBamTagsSequence & tags) : _host(tags) {}

    template <typename TPos>
    inline Infix<Host<BamTagsDict const>::Type>::Type
    operator[] (TPos pos) const
    {
        if (!hasIndex(*this))
            buildIndex(*this);
        return infix(host(*this), _positions[pos], _positions[pos + 1]);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

inline Host<BamTagsDict>::Type &
host(BamTagsDict & bamTags)
{
    return value(bamTags._host);
}

inline Host<BamTagsDict const>::Type &
host(BamTagsDict const & bamTags)
{
    return value(bamTags._host);
}

// ----------------------------------------------------------------------------
// Function hasIndex()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#hasIndex
 * @brief Returns whether the BamTagsDict has an index.
 *
 * @signature bool hasIndex(dict);
 *
 * @param[in] dict The @link BamTagsDict @endlink to query.
 *
 * @return bool true if <tt>dict</tt> has an index and false otherwise.
 */

/**
.Function.BamTagsDict#hasIndex
..class:Class.BamTagsDict
..cat:Fragment Store
..summary:Return $true$ if @Class.BamTagsDict@ has an index.
..signature:hasIndex(bamTags)
..param.bamTags:SAM Tags to query
...type:Class.BamTagsDict
..returns:$bool$
..include:<seqan/store_ex.h>
*/

inline bool
hasIndex(BamTagsDict const & bamTags)
{
    return !empty(bamTags._positions) || empty(host(bamTags));
}

// ----------------------------------------------------------------------------
// Function getBamTypeSize()
// ----------------------------------------------------------------------------

/*!
 * @fn getBamTypeSize
 * @headerfile <seqan/bam_io.h>
 * @brief Return size of the type identified by a type char.
 *
 * @signature int getBamTypeSize(c);
 *
 * @param[in] c A <tt>char</tt> that identifies a type.
 *
 * @return int The size of the type in bytes, -1 vor variable-sized types, -2 for invalid paramters.
 *
 * @see BamTagsDict
 * @see getBamTypeChar
 */

// Return sizeof() of the type identified with the given char.  Returns -2 if not
// valid, -1 if of variable length.

/**
.Function.getBamTypeSize
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:getBamTypeSize(c)
..summary:Return size of the type identified by $c$.
..param.c:The BAM type identifier
..returns:$int$ with the $sizeof()$ of the type, -1 for variable sized types, -2 for invalid parameters.
..include:seqan/bam_io.h
*/

inline int
getBamTypeSize(char c)
{
    switch (c)
    {
        case 'A':
            return 1;
        case 'f':
            return 4;
        case 'Z':
        case 'H':
        case 'B':
            return -1;
        case 'c':
        case 'C':
            return 1;
        case 's':
        case 'S':
            return 2;
        case 'i':
        case 'I':
            return 4;
    }
    return -2;
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#buildIndex
 * @brief Build index for a BamTagsDict  object.
 *
 * @signature void buildIndex(bamTags);
 *
 * @param[in,out] bamTags The BamTagsDict object to build the index for.
 */

/**
.Function.BamTagsDict#buildIndex
..class:Class.BamTagsDict
..cat:Fragment Store
..summary:Build index for a @Class.BamTagsDict@ object.
..signature:buildIndex(bamTags)
..param.bamTags:SAM Tags to build index for.
...type:Class.BamTagsDict
..returns:$void$
..include:<seqan/bam_io.h>
*/

inline void
buildIndex(BamTagsDict const & bamTags)
{
    typedef Host<BamTagsDict>::Type TTagString;
    typedef Iterator<TTagString const, Standard>::Type TIter;

    clear(bamTags._positions);
    if (empty(value(bamTags._host)))
        return;  // Done.

    appendValue(bamTags._positions, 0);
    TIter itBegin = begin(host(bamTags), Standard());
    TIter itEnd = end(host(bamTags), Standard());
    for (TIter it = itBegin; it != itEnd; )
    {
        // skip tag name (e.g. "NM")
        it += 2;

        // get tag type (e.g. 'I')
        register char c = *(it++);
        if (c == 'H' || c == 'Z')
        {
            // skip string and its end-of-string marker
            while (*it != '\0')
            {
                ++it;
                SEQAN_ASSERT(it != itEnd);
            }
            ++it;
        }
        else if (c == 'B')
        {
            // skip array of PODs
            c = *(it++);
            union {
                char raw[4];
                __uint32 len;
            } tmp;
            arrayCopyForward(it, it + 4, tmp.raw);
            it += 4 + tmp.len * getBamTypeSize(c);
        }
        else
        {
            // skip POD type (e.g. byte, int)
            it += getBamTypeSize(c);
        }
        appendValue(bamTags._positions, it - itBegin);
    }
    // if (!empty(value(bamTags._host)))
    //     appendValue(bamTags._positions, length(host(bamTags)) + 1);  // +1 since there is not tab at the end
}

// ----------------------------------------------------------------------------
// Function _dataHost()
// ----------------------------------------------------------------------------

inline Holder<Host<BamTagsDict>::Type> &
_dataHost(BamTagsDict & bamTags)
{
    return bamTags._host;
}

inline Holder<Host<BamTagsDict>::Type> const &
_dataHost(BamTagsDict const & bamTags)
{
    return bamTags._host;
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename THost>
inline void
setHost(BamTagsDict & me, THost & host_)
{
    SEQAN_CHECKPOINT;
	setValue(_dataHost(me), host_);
    clear(me._positions);
}

template <typename THost>
inline void
setHost(BamTagsDict & me, THost const & host_)
{
    SEQAN_CHECKPOINT;
	setValue(_dataHost(me), host_);
    clear(me._positions);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#length
 * @brief Returns the number of entries in a BamTagsDict.
 *
 * @signature unsigned length(tagsDict);
 *
 * @param[in] tagsDict The BamTagsDict object to query for its length.
 *
 * @return unsigned The number of entries in the BamTagsDict.
 */

///.Function.length.param.object.type:Class.BamTagsDict

inline Size<BamTagsDict>::Type
length(BamTagsDict const & tags)
{
    if (empty(value(tags._host)))
        return 0;
    if (!hasIndex(tags))
        buildIndex(tags);
    return length(tags._positions) - 1;
}

// ----------------------------------------------------------------------------
// Function getTagType()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#getTagType
 * @brief Returns the tag type char for an entry of a BamTagsDict.
 *
 * @signature char getTagType(tags, id);
 *
 * @param[in] tags The BamTagsDict to query.
 * @param[in] id   The position for which to retrieve the type <tt>__int32</tt>.
 */

/**
.Function.BamTagsDict#getTagType
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:getTagType(tagsDict, id)
..summary:Get key of a tag by index.
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.id:Index of the tag whose key to retrieve.
..returns:$char$, the SAM/BAM identifier of the type.
..include:seqan/bam_io.h
*/

template <typename TId>
inline char
getTagType(BamTagsDict const & tags, TId id)
{
    return tags[id][2];
}

// ----------------------------------------------------------------------------
// Function getTagKey()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#getTagKey
 * @brief Return key of a tag by index.
 *
 * @signature TKey getTagKey(tagsDict, id);
 *
 * @param[in] tagsDict The BamTagsDict to query.
 * @param[in] id       The index of the dict entry (<tt>__int32</tt>).
 *
 * @return TKey An infix of a @link CharString @endlink.  Will be a two-character char sequence.
 */

/**
.Function.BamTagsDict#getTagKey
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:getTagKey(tagsDict, id)
..summary:Return key of a tag by index.
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
...type:Class.BamTagsDict
..param.id:Index of the tag whose key to retrieve.
..returns:Infix of the underlying string.
..remarks:See @Class.BamTagsDict@ for an example.
..include:seqan/bam_io.h
*/

template <typename TId>
inline Infix<Host<BamTagsDict const>::Type>::Type
getTagKey(BamTagsDict const & tags, TId id)
{
    return prefix(tags[id], 2);
}

// ----------------------------------------------------------------------------
// Function findTagKey()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#findTagKey
 * @brief Find a tag by its key for a @link BamTagsDict @endlink object.
 *
 * @signature bool findTagKey(id, tagsDict, key);
 *
 * @param[out] id       The index of the tag is stored here (<tt>unsigned</tt>).
 * @param[in]  tagsDict The BamTagsDict to query.
 * @param[in]  key      The key to query for: @link CharString @endlink.
 *
 * @return bool true if the key could be found and false otherwise.
 */

/**
.Function.BamTagsDict#findTagKey
..summary:Find a tag by its key for a @Class.BamTagsDict@ object.
..class:Class.BamTagsDict
..signature:findTagKey(id, tagsDict, key)
..param.id:Index of the tag with the given key.
...type:nolink:$unsigned$
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.key:Name of the key to find.
...type:Shortcut.CharString
..returns:$bool$, indicating whether such a key could be found.
..include:seqan/bam_io.h
*/

template <typename TId, typename TKey>
inline bool
findTagKey(TId & id, BamTagsDict const & tags, TKey const & key)
{
    for (id = 0; id < (TId)length(tags); ++id)
        if (getTagKey(tags, id) == key)
            return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function extractTagValue()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#extractTagValuej
 * @brief Extract and cast "atomic" value from tags string with index <tt>id</tt>.
 *
 * @signature bool extractTagValue(dest, tags, id)
 *
 * @param[out] dest The variable to write the value to.The value is first copied in a variable of the type indicated in
 *                  the BAM file. Then it is cast into the type of <tt>dest</tt>.
 *
 * @param[in] tags The BamTagsDict object to query.
 * @param[in] id   The integer index in the dict to use (<tt>__int32</tt>).
 *
 * @return bool true if the value could be extracted.
 *
 * @section Remarks
 *
 * The function only works for atomic types such as <tt>int</tt>, not for <tt>char*</tt> or arrays.
 *
 * See @link BamTagsDict @endlink for an example.
 */

/**
.Function.BamTagsDict#extractTagValue
..class:Class.BamTagsDict
..cat:BAM I/O
..signature:extractTagValue(dest, tags, id)
..summary:Extract and cast "atomic" value from tags string with index $id$.
..param.dest:The variable to write the value to.
...remarks:The value is first copied in a variable of the type indicated in the BAM file. Then it is cast into the type of $dest$.
..param.tags:@Class.BamTagsDict@ object.
...type:Class.BamTagsDict
..params.id:Index of the tag in the tag list.
..returns:$bool$, indicating the success.
..remarks:The function only works for atomic types such as $int$, not for $char*$ or arrays.
..remarks:See @Class.BamTagsDict@ for an example.
..see:Function.BamTagsDict#getTagValue
..include:seqan/bam_io.h
*/

template <typename TResultValue, typename TId>
SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TResultValue> >, bool)
extractTagValue(TResultValue & val, BamTagsDict const & tags, TId id)
{
    typedef Infix<Host<BamTagsDict const>::Type>::Type TInfix;
    typedef Iterator<TInfix, Standard>::Type TIter;

    TInfix inf = tags[id];
    if (length(inf) < 4)
        return false;

    TIter it = begin(inf, Standard()) + 2;
    char typeC = getValue(it++);

    if (typeC == 'A' || typeC == 'c' || typeC == 'C')
    {
        val = getValue(it);
    }
    else if (typeC == 's' || typeC == 'S')
    {
        if (length(inf) < 5)
            return false;

        union {
            char raw[2];
            short i;
        } tmp;

        arrayCopyForward(it, it + 2, tmp.raw);
        val = tmp.i;
    }
    else if (typeC == 'i' || typeC == 'I' || typeC == 'f')
    {
        if (length(inf) < 7)
            return false;

        union {
            char raw[4];
            int i;
            float f;
        } tmp;

        arrayCopyForward(it, it + 4, tmp.raw);
        if (typeC == 'f')
            val = tmp.f;
        else
            val = tmp.i;
    }
    else // variable sized type or invald
    {
        return false;
    }
    return true;
}

template <typename TResultValue, typename TId>
SEQAN_FUNC_ENABLE_IF(Is<IsSequence<TResultValue> >, bool)
extractTagValue(TResultValue & val, BamTagsDict const & tags, TId id)
{
    typedef Infix<Host<BamTagsDict const>::Type>::Type TInfix;

    TInfix inf = tags[id];
    if (length(inf) < 4 || inf[2] != 'Z')
        return false;

    val = infix(inf, 3, length(inf) - 1);
    return true;
}

// ----------------------------------------------------------------------------
// Function getBamTypeChar()
// ----------------------------------------------------------------------------

/*!
 * @fn getBamTypeChar
 * @headerfile <seqan/bam_io.h>
 * @brief Return char identifying the type of the argument type.
 *
 * @signature char getBamTypeChar<T>();
 *
 * @tparam T The type to query for its type char.
 *
 * @section Remarks
 *
 * Note that this function is defined for the <tt>__int16</tt>, <tt>__uint16</tt> etc. but not for the types
 * <tt>short</tt>, <tt>int</tt> etc. An exception are 8-bit characters/char, where it is defined for <tt>__int8</tt>,
 * <tt>__uint8</tt>, and <tt>char</tt> unless <tt>char</tt> is equal to one of the other two types. This is important
 * when used in @link BamTagsDict#setTagValue @endlink etc. since BAM gives type chars for printable characters, signed
 * 8-bit numbers and unsigned 8-bit numbers.
 *
 * If <tt>__int8</tt> and <tt>__uint8</tt> are not identical to <tt>char</tt>, we can make this decision from the type,
 * otherwise we cannot and we will give the integer types a higher precedence.
 *
 * In your programs, this should not make any difference, only the written SAM/BAM will differ.
 *
 * @see BamTagsDict
 * @see getBamTypeSize
 */

/**
.Function.getBamTypeChar
..class:Class.BamTagsDict
..cat:BAM I/O
..summary:Return char identifying the type of the atomic argument.
..signature:getBamTypeChar<T>()
..param.T:The type to get the BAM char for.
..returns:$char$ describing the BAM type. One of $ACcSsIifZ$.
..remarks:Note that this function is defined for the $__int16$, $__uint16$ etc. but not for the types $short$, $int$ etc. An exception are 8-bit characters/char, where it is defined for $__int8$, $__uint8$, and $char$ unless $char$ is equal to one of the other two types. This is important when used in @Function.BamTagsDict#setTagValue@ etc. since BAM gives type chars for printable characters, signed 8-bit numbers and unsigned 8-bit numbers.
..remarks:If $__int8$ and $__uint8$ are not identical to $char$, we can make this decision from the type, otherwise we cannot and we will give the integer types a higher precedence.
..remarks:In your programs, this should not make any difference, only the written SAM/BAM will differ.
..include:seqan/bam_io.h
*/

template <typename TValue>
struct BamTypeChar
{
    static const char VALUE;
};

template <>
const char BamTypeChar<char>::VALUE = 'A';
template <>
const char BamTypeChar<__int8>::VALUE = 'C';
template <>
const char BamTypeChar<__uint8>::VALUE = 'c';
template <>
const char BamTypeChar<__int16>::VALUE = 's';
template <>
const char BamTypeChar<__uint16>::VALUE = 'S';
template <>
const char BamTypeChar<__int32>::VALUE = 'i';
template <>
const char BamTypeChar<__uint32>::VALUE = 'I';
template <>
const char BamTypeChar<float>::VALUE = 'f';
template <>
const char BamTypeChar<double>::VALUE = 'f';
template <>
const char BamTypeChar<char *>::VALUE = 'Z';

template <typename T>
inline char getBamTypeChar(T const &)
{
    return BamTypeChar<T>::VALUE;
}

// ----------------------------------------------------------------------------
// Function setTagValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Test me!

/*!
 * @fn BamTagsDict#setTagValue
 *
 * @headerfile seqan/bam_io.h
 *
 * @brief Set the value of a tag through a @link BamTagsDict @endlink.
 *
 * @signature bool setTagValue(tags, key, val[, typeC]);
 *
 * @param[in,out] tags  The BamTagsDict to modify.
 * @param[in]     key   The key of the tag.Must be a string of length 2. Types: CharString
 * @param[in]     val   The value to set the the tag to.

 * @param[in]     typeC BAM type char to use.For portability (so the generated files are the same on all platforms), use
 *                      a signed/unsigned qualified type for <tt>val</tt> or give <tt>typeC</tt>. Also see the remarks
 *                      for @link getBamTypeChar @endlink. Types: getBamTypeChar@.

 *
 * @return bool true on success, false on failure.  This function can fail if the key is not a valid tag id (e.g. does
 *              not have length 2) or if the type of <tt>val</tt> is not an atomic value or a string (anything but
 *              <tt>char *</tt>, <tt>char const *</tt>, a character, integer or float type is invalid).
 *
 * @section Remarks
 *
 * Note that <tt>setTagValue</tt> does not cast the type, so <tt>typeC</tt> only influences the type character written
 * out but <tt>val</tt> is written out in binary without modification.
 *
 * @section Examples
 *
 * An example setting some atomic tag values.
 *
 * @code{.cpp}
 * CharString rawTagsText;
 * BamTagsDict tags(rawTagsText);
 * setTagValue(tags, "XA", 9);    // int
 * setTagValue(tags, "XB", 9u);   // unsigned int
 * setTagValue(tags, "XC", 'X');  // char
 * @endcode
 *
 * If <tt>char</tt> is equal to <tt>__int8</tt> or <tt>__uint8</tt> then the last line produces an entry with type 'c'
 * or 'C'. To make sure that the type char 'A' (for "printable character") is written to the file, give it explicitely:
 *
 * @code{.cpp}
 * setTagValue(tags, "XC", 'X', 'A');  // Overrwrite XC, enforce type 'printable character'.
 * @endcode
 *
 * Note that on most systems <tt>int</tt>s have a width of 32 bytes, but the C++ standard leaves this open. For all
 * types but characters, you should not give an explicit type char but use one of the types with explicit width and
 * signed/unsigned qualifier such as <tt>__int32</tt>, <tt>__uint32</tt> etc.
 *
 * @code{.cpp}
 * // The following is not recommended since the type of <tt>x</tt> is not "unsigned 32 bit int."
 * __int32 x = -1;
 * setTagValue(tags, "XB", x, 'I');
 * // Instead, explicitely use an unsigned type if you need one.  Note that your compiler
 * // might warn you about assigning -1 to an unsigned variable so you know that you are
 * // probably doing something unintended.
 * __uint32 y = -1;
 * setTagValue(tags, "XB", y);
 *
 * // Do not do this!
 * setTagValue(tags, "XA", 9, 'f');    // BOGUS since 9 is not a floating point number.
 * @endcode
 *
 * @see getBamTypeChar
 */

/**
.Function.BamTagsDict#setTagValue
..class:Class.BamTagsDict
..cat:BAM I/O
..summary:Set the value of a tag through a @Class.BamTagsDict@.
..signature:setTagValue(tags, key, val[, typeC])
..param.tags:The dict to modify.
...type:Class.BamTagsDict
..param.key:The key of the tag.
...type:Shortcut.CharString
...remarks:Must be a string of length 2.
..param.val:The value to set the the tag to.
..param.typeC:BAM type char to use.
...type:nolink:By default, the type is inflected using @Function.getBamTypeChar@.
...remarks:For portability (so the generated files are the same on all platforms), use a signed/unsigned qualified type for $val$ or give $typeC$. Also see the remarks for @Function.getBamTypeChar@.
..returns:$bool$ indicating the success. This function can fail if the key is not a valid tag id (e.g. does not have length 2) or if the type of $val$ is not an atomic value or a string (anything but $char *$, $char const *$, a character, integer or float type is invalid).
..see:Function.getBamTypeChar
..remarks:Note that $setTagValue$ does not cast the type, so $typeC$ only influences the type character written out but $val$ is written out in binary without modification.
..include:seqan/bam_io.h
..example.text:An example setting some atomic tag values.
..example.code:
CharString rawTagsText;
BamTagsDict tags(rawTagsText);
setTagValue(tags, "XA", 9);    // int
setTagValue(tags, "XB", 9u);   // unsigned int
setTagValue(tags, "XC", 'X');  // char
..example.text:If $char$ is equal to $__int8$ or $__uint8$ then the last line produces an entry with type 'c' or 'C'. To make sure that the type char 'A' (for "printable character") is written to the file, give it explicitely:
..example.code:
setTagValue(tags, "XC", 'X', 'A');  // Overrwrite XC, enforce type 'printable character'.
..example.text:Note that on most systems $int$s have a width of 32 bytes, but the C++ standard leaves this open. For all types but characters, you should not give an explicit type char but use one of the types with explicit width and signed/unsigned qualifier such as $__int32$, $__uint32$ etc.
..example.code:
// The following is not recommended since the type of $x$ is not "unsigned 32 bit int."
__int32 x = -1;
setTagValue(tags, "XB", x, 'I');
// Instead, explicitely use an unsigned type if you need one.  Note that your compiler
// might warn you about assigning -1 to an unsigned variable so you know that you are
// probably doing something unintended.
__uint32 y = -1;
setTagValue(tags, "XB", y);

// Do not do this!
setTagValue(tags, "XA", 9, 'f');    // BOGUS since 9 is not a floating point number.
*/

// Convert "atomic" value to BAM tag.  Return whether val was atomic.
template <typename TBamValueSequence, typename TValue>
SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TValue> >, bool)
_toBamTagValue(TBamValueSequence & result, TValue const & val, char typeC)
{
    appendValue(result, typeC);

    if (typeC == 'A' || typeC == 'c' || typeC == 'C')
    {
        appendValue(result, (char)val);
    }
    else if (typeC == 's' || typeC == 'S')
    {
        union {
            char raw[2];
            short i;
        } tmp;

        tmp.i = val;

        append(result, toRange(&tmp.raw[0], &tmp.raw[2]));
    }
    else if (typeC == 'i' || typeC == 'I' || typeC == 'f')
    {
        union {
            char raw[4];
            int i;
            float f;
        } tmp;

        if (typeC == 'f')
            tmp.f = val;
        else
            tmp.i = val;
            
        append(result, toRange(&tmp.raw[0], &tmp.raw[4]));
    }
    else // non-string and variable sized type or invald
    {
        resize(result, length(result) - 1);
        return false;
    }
    return true;
}

template <typename TBamValueSequence, typename TValue>
SEQAN_FUNC_ENABLE_IF(IsSequence<TValue>, bool)
_toBamTagValue(TBamValueSequence & result, TValue const & val, char typeC)
{
    if (typeC != 'Z')
        return false;
    
    appendValue(result, typeC);
    append(result, val);
    appendValue(result, '\0');
    return true;
}


// Sets an atomic value in a BamTagsDict.
// Returns true successful, can fail if val not atomic or key is not a valid tag id (2 chars).

template <typename TKey, typename TValue>
inline bool
setTagValue(BamTagsDict & tags, TKey const & key, TValue const & val, char typeC)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    // Build value to insert/append.
    if (length(key) != 2u)
        return false;

    Position<BamTagsDict>::Type id = 0;
    if (findTagKey(id, tags, key))
    {
        CharString bamTagVal;
        if (!_toBamTagValue(bamTagVal, val, typeC))
            return false;

        replace(host(tags), tags._positions[id] + 2, tags._positions[id + 1], bamTagVal);
        clear(tags._positions);
    }
    else
    {
        append(host(tags), key);
        if (!_toBamTagValue(host(tags), val, typeC))
        {
            resize(host(tags), length(host(tags)) - length(key));
            return false;
        }
        appendValue(tags._positions, length(host(tags)));
    }

    return true;
}

template <typename TKey, typename TValue>
inline bool
setTagValue(BamTagsDict & tags, TKey const & key, TValue const & val)
{
    return setTagValue(tags, key, val, BamTypeChar<TValue>::VALUE);
}



template <typename TSequence, typename TKey, typename TValue>
inline bool
appendTagValue(TSequence & tags, TKey const & key, TValue const & val, char typeC)
{
    // Build value to insert/append.
    if (length(key) != 2u)
        return false;

    append(tags, key);
    return _toBamTagValue(tags, val, typeC);
}

template <typename TKey, typename TValue>
inline bool
appendTagValue(BamTagsDict & tags, TKey const & key, TValue const & val, char typeC)
{
    if (appendTagValue(host(tags), key, val, typeC))
    {
        appendValue(tags._positions, length(host(tags)));
        return true;
    }
    return false;
}


template <typename TDictOrString, typename TKey, typename TValue>
inline bool
appendTagValue(TDictOrString & tags, TKey const & key, TValue const & val)
{
    return appendTagValue(tags, key, val, BamTypeChar<TValue>::VALUE);
}


// ----------------------------------------------------------------------------
// Function eraseTag()
// ----------------------------------------------------------------------------

/*!
 * @fn BamTagsDict#eraseTag
 * @brief Erase a tag from BamTagsDict.
 *
 * @signature bool eraseTag(tagsDict, key);
 *
 * @param[in,out] tagsDict The BamTagsDict to erase the tag from.
 * @param[in]     key      The key of the tag to erase, of type @link CharString @endlink.
 *
 * @return bool true if the tag was present for erasing, false if not.
 */

/**
.Function.BamTagsDict#eraseTag
..class:Class.BamTagsDict
..summary:Erase tag from @Class.BamTagsDict@.
..cat:BAM I/O
..signature:eraseTag(tagsDict, key)
..param.tags:The dict to erase from.
...type:Class.BamTagsDict
..param.key:The key of the entry to remove.
...type:Shortcut.CharString
..returns:$bool$, indicating whether the key was present.
..include:seqan/bam_io.h
 */

template <typename TKey>
inline bool
eraseTag(BamTagsDict & tags, TKey const & key)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    Position<BamTagsDict>::Type id = 0;
    if (!findTagKey(id, tags, key))
        return false;

    erase(host(tags), tags._positions[id], tags._positions[id + 1]);
    clear(tags._positions);
    return true;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
