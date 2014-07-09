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
// Code to convert between SAM and BAM format tags (textual <-> binary).
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_SAM_TAGS_TO_BAM_TAGS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assignTagsSamToBam()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TIter, typename TType1, typename TType2>
void _assignTagsToBamOneTagHelper(TTarget & target,
                                  CharString & buffer2,
                                  TIter & it,
                                  __int32 nEntries,
                                  TType1 /*tag*/,
                                  TType2 /*tag*/)
{
    for (int i = 0; i < nEntries; ++i)
    {
        clear(buffer2);
        for (; !atEnd(it) && *it != ',' && *it != '\t'; skipOne(it))
            appendValue(buffer2, *it);
        TType1 x = 0;  // short to avoid textual interpretation in lexicalCast<> below.
        x = lexicalCast<TType1>(buffer2);
        if(IsSameType<TType2, __int8>::VALUE)
        {
            appendValue(target, static_cast<__int8>(x));
            if (!atEnd(it) && *it == ',')
                skipOne(it);  // Skip ','.
            else
                break;  // End of field or end of string.
        }
        else
        {
            char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
            for (int i = 0; i < BytesPerValue<TType1>::VALUE; ++i, ++ptr)
                appendValue(target, *ptr);
            if (!atEnd(it) && *it == ',')
                goNext(it);  // Skip ','.
            else
                break;  // End of field or end of string.
        }
    }
}

template <typename TTarget, typename TForwardIter>
void _assignTagsSamToBamOneTag(TTarget & target, TForwardIter & iter, CharString & buffer)
{
    clear(buffer);
    readUntil(buffer, iter, CountDownFunctor<>(2));
    append(target, buffer);
    
    clear(buffer);
    readUntil(buffer, iter, CountDownFunctor<>(3));
    SEQAN_ASSERT_EQ(buffer[0], ':');
    SEQAN_ASSERT_EQ(buffer[2], ':');
    char t = buffer[1];
    appendValue(target, t);
    
    switch (t)
    {
    case 'A':
        {
            clear(buffer);
            readOne(buffer, iter);
            append(target, buffer);
        }
        break;
    case 'i':
        {
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
            __int32 x;
            x = lexicalCast<__int32>(buffer);
            char const * ptr = reinterpret_cast<char const *>(&x);
            for (int i = 0; i < 4; ++i, ++ptr)
                appendValue(target, *ptr);
        }
        break;
    case 'f':
        {
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
            float x;
            x = lexicalCast<float>(buffer);
            char const * ptr = reinterpret_cast<char const *>(&x);
            // TODO(singer): Why not just append buffer?
            for (int i = 0; i < 4; ++i, ++ptr)
                appendValue(target, *ptr);
        }
        break;
    case 'H':
    case 'Z':
        {
            // TODO(holtgrew): Could test on even length in case of 'H'.
            readUntil(target, iter, OrFunctor<IsTab, IsNewline>());
            appendValue(target, '\0');
        }
        break;
    case 'B':
        {
            CharString buffer2; // TODO(holtgrew): Also give from outside.
            char const c = 'a';

            // Read type.
            char t2;
            readOne(t2, iter);
            appendValue(target, t2);

            // Read array contents.
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
            typename Iterator<CharString, Rooted>::Type it, it2;
            // Search first non-comma position.
            it = begin(buffer, Rooted());
            for (;!atEnd(it) && *it == ','; ++it)
                continue;
            // Count number of entries.
            __int32 nEntries = !atEnd(it);  // At least one if array not empty.
            for (it2 = it; !atEnd(it2); ++it2)
                nEntries += (*it2 == ',');
            // Write out array length to result.
            char const * ptr = reinterpret_cast<char const *>(&nEntries);
            for (int i = 0; i < 4; ++i, ++ptr)
                appendValue(target, *ptr);

            // Now, write out the arrays, depending on the entry type.
            // TODO(holtgrew): Whee, this could be a bit more compact...
            switch (t2)
            {
            case 'c':
                _assignTagsToBamOneTagHelper(target, buffer2, it, nEntries, __int16(), __int8());
                /*for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __int16 x = 0;  // short to avoid textual interpretation in lexicalCast<> below.
                    bool b = lexicalCast2<__int16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    appendValue(target, static_cast<__int8>(x));
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }*/
                break;
            case 'C':
                _assignTagsToBamOneTagHelper(target, buffer2, it, nEntries, __uint16(), __int8());
                /*for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __uint16 x = 0;  // short to avoid textual interpretation in lexicalCast<> below.
                    bool b = lexicalCast2<__uint16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    appendValue(target, static_cast<__int8>(x));
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }*/
                break;
            case 's':
                _assignTagsToBamOneTagHelper(target, buffer2, it, nEntries, __int16(), &c);
                /*for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __int16 x = 0;
                    bool b = lexicalCast2<__int16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 2; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }*/
                break;
            case 'S':
                _assignTagsToBamOneTagHelper(target, buffer2, it, nEntries, __uint16(), &c);
                /*for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __uint16 x = 0;
                    bool b = lexicalCast2<__uint16>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 2; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }*/
                break;
            case 'i':
                _assignTagsToBamOneTagHelper(target, buffer2, it, nEntries, __int32(), &c);
                /*for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __int32 x = 0;
                    bool b = lexicalCast2<__int32>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 4; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }*/
                break;
            case 'I':
                _assignTagsToBamOneTagHelper(target, buffer2, it, nEntries, __uint32(), &c);
                /*for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    __uint32 x = 0;
                    bool b = lexicalCast2<__uint32>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 4; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }*/
                break;
            case 'f':
                _assignTagsToBamOneTagHelper(target, buffer2, it, nEntries, float(), &c);
                /*for (int i = 0; i < nEntries; ++i)
                {
                    SEQAN_ASSERT_NOT(atEnd(it));
                    clear(buffer2);
                    for (; !atEnd(it) && *it != ',' && *it != '\t'; goNext(it))
                        appendValue(buffer2, *it);
                    float x = 0;
                    bool b = lexicalCast2<float>(x, buffer2);
                    (void)b;
                    SEQAN_ASSERT(b);
                    char const * ptr = reinterpret_cast<char const *>(&x);  // write out byte-wise
                    for (int i = 0; i < 4; ++i, ++ptr)
                        appendValue(target, *ptr);
                    if (!atEnd(it) && *it == ',')
                        goNext(it);  // Skip ','.
                    else
                        break;  // End of field or end of string.
                }*/
                break;
            default:
                SEQAN_FAIL("Invalid array type: %c!", t2);
            }
        }
        break;
    default:
        SEQAN_ASSERT_FAIL("Invalid tag type: %c!", t);
    }
}

/*!
 * @fn assignTagsSamToBam
 * @headerfile <seqan/bam_io.h>
 * @brief Assign tags in SAM format to tags in BAM format.
 *
 * @signature void assignTagsBamToSam(bamTags, samTags);
 *
 * @param bamTags[out] A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the target BAM tags.
 * @param samTags[in] A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the source SAM tags.
 *
 * @see assignTagsBamToSam
 */

/**
.Function.assignTagsSamToBam
..cat:BAM I/O
..summary:Assign tags in SAM format to tags in BAM format.
..signature:assignTagsSamToBam(bamTags, samTags)
..param.bamTags:Destination BAM tags.
...type:Shortcut.CharString
..param.samTags:Source SAM tags.
...type:Shortcut.CharString
..returns:$void$
..include:seqan/bam_io.h
*/

template <typename TTarget, typename TSource>
void assignTagsSamToBam(TTarget & target, TSource & source)
{
    // Handle case of empty source sequence.
    if (empty(source))
        clear(target);

    typedef typename Iterator<TSource, Rooted>::Type TSourceIter;
    TSourceIter it = begin(source, Rooted());

    CharString buffer;

    while (!atEnd(it))
    {
        if (value(it) == '\t')
            skipOne(it);

        _assignTagsSamToBamOneTag(target, it, buffer);
    }
}

// ----------------------------------------------------------------------------
// Function assignTagsBamToSam()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSourceIter>
void _assignTagsBamToSamOneTag(TTarget & target, TSourceIter & it)
{
    // Copy tag name.
    SEQAN_ASSERT_NOT(atEnd(it));
    appendValue(target, *it++);
    SEQAN_ASSERT_NOT(atEnd(it));
    appendValue(target, *it++);
    unsigned char t = *it;

    // Add ':'.
    appendValue(target, ':');

    // Add type.
    SEQAN_ASSERT_NOT(atEnd(it));
    if (*it == 'c' || *it == 'C' || *it == 's' || *it == 'S' || *it == 'i' || *it == 'I')
        appendValue(target, 'i');
    else
        appendValue(target, *it);
    ++it;

    // Add ':'.
    appendValue(target, ':');

    // Convert the payload, depending on the field's type.

    switch (t)
    {
    case 'A':
        appendValue(target, *it++);
        break;
    case 'c':
        {
            SEQAN_ASSERT_NOT(atEnd(it));
            __int8 x = *it++;
            char buffer[4];
            snprintf(buffer, 4, "%d", x);
            append(target, buffer);
        }
        break;
    case 'C':
        {
            SEQAN_ASSERT_NOT(atEnd(it));
            char buffer[4];
            __uint8 x = *it++;
            snprintf(buffer, 4, "%u", x);
            append(target, buffer);
        }
        break;
    case 's':
        {
            __int16 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 2; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            char buffer[32];
            snprintf(buffer, 32, "%d", x);
            append(target, buffer);
        }
        break;
    case 'S':
        {
            __uint16 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 2; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            char buffer[32];
            snprintf(buffer, 32, "%u", x);
            append(target, buffer);
        }
        break;
    case 'i':
        {
            __int32 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            char buffer[32];
            snprintf(buffer, 32, "%d", x);
            append(target, buffer);
        }
        break;
    case 'I':
        {
            __uint32 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            char buffer[32];
            snprintf(buffer, 32, "%u", x);
            append(target, buffer);
        }
        break;
    case 'f':
        {
            float x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            char buffer[32];
            snprintf(buffer, 32, "%g", x);
            append(target, buffer);
        }
        break;
    case 'Z':
        {
            while (*it != '\0')
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                appendValue(target, *it++);
            }
            SEQAN_ASSERT_NOT(atEnd(it));
            it++;
        }
        break;
    case 'H':
        {
            while (*it != '\0')
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                appendValue(target, *it++);
            }
            SEQAN_ASSERT_NOT(atEnd(it));
            it++;
        }
        break;
    case 'B':
        {
            // Read type.
            char t2 = *it++;
            appendValue(target, t2);
            // Read array length.
            __int32 x = 0;
            char * ptr = reinterpret_cast<char *>(&x);
            for (int i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                *ptr++ = *it++;
            }
            // Depending on t2, read array.
            // TODO(holtgrew): Whee, this could be a bit more compact...
            switch (t2)
            {
            case 'c':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    __int8 y = *it++;
                    char buffer[32];
                    snprintf(buffer, 32, "%d", y);
                    append(target, buffer);
                }
                break;
            case 'C':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    __uint8 y = *it++;
                    char buffer[32];
                    snprintf(buffer, 32, "%u", y);
                    append(target, buffer);
                }
                break;
            case 's':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    __int16 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 2; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    char buffer[32];
                    snprintf(buffer, 32, "%d", y);
                    append(target, buffer);
                }
                break;
            case 'S':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    __uint16 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 2; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    char buffer[32];
                    snprintf(buffer, 32, "%d", y);
                    append(target, buffer);
                }
                break;
            case 'i':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    __int32 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 4; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    char buffer[32];
                    snprintf(buffer, 32, "%d", y);
                    append(target, buffer);
                }
                break;
            case 'I':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    __uint32 y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 4; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    char buffer[32];
                    snprintf(buffer, 32, "%u", y);
                    append(target, buffer);
                }
                break;
            case 'f':
                for (__int32 i = 0; i < x; ++i)
                {
                    appendValue(target, ',');
                    float y = 0;
                    char * ptr = reinterpret_cast<char *>(&y);
                    for (int i = 0; i < 4; ++i)
                    {
                        SEQAN_ASSERT_NOT(atEnd(it));
                        *ptr++ = *it++;
                    }
                    char buffer[32];
                    snprintf(buffer, 32, "%g", y);
                    append(target, buffer);
                }
                break;
            default:
                SEQAN_FAIL("Invalid array type: %c!", t2);
            }
        }
        break;
    default:
        SEQAN_ASSERT_FAIL("Invalid tag type: %c!", t);
    }
}

/*!
 * @fn assignTagsBamToSam
 * @headerfile <seqan/bam_io.h>
 * @brief Assign tags in BAM format to tags in SAM format.
 *
 * @signature void assignTagsBamToSam(samTags, bamTags);
 *
 * @param samTags[out] A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the target SAM tags.
 * @param bamTags[in] A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the source BAM tags.
 *
 * @see assignTagsSamToBam
 */

/**
.Function.assignTagsBamToSam
..cat:BAM I/O
..summary:Assign tags in BAM format to tags in SAM format.
..signature:assignTagsSamToBam(bamTags, samTags)
..param.samTags:Destination SAM tags.
...type:Shortcut.CharString
..param.bamTags:Source BAM tags.
...type:Shortcut.CharString
..returns:$void$
..include:seqan/bam_io.h
..see:Function.assignTagsSamToBam
*/

template <typename TTarget, typename TSource>
void assignTagsBamToSam(TTarget & target, TSource const & source)
{
    // Handle case of empty source sequence.
    if (empty(source))
        clear(target);

    clear(target);

    typedef typename Iterator<TSource const, Rooted>::Type TSourceIter;
    TSourceIter it = begin(source, Rooted());

    bool first = true;
    while (!atEnd(it))
    {
        if (!first)
            appendValue(target, '\t');
        first = false;
        _assignTagsBamToSamOneTag(target, it);
    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_SAM_TAGS_TO_BAM_TAGS_H_
