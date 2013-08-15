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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Functions for tokenizing streams
// ==========================================================================

// TODO(holtgrew): Some skip*, skipUntil* functions are missing.

#ifndef SEQAN_STREAM_TOKENIZE_H
#define SEQAN_STREAM_TOKENIZE_H

#include <cctype>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

/*!
 * @enum TokenizeResult
 * @headerfile <seqan/stream.h>
 * @brief Enum with return values for Tokenization operations.
 *
 * @signature enum TokenizeResult;
 *
 * @var TokenizeResult SUCCESS = 0;
 * @brief Reading the specified data succeeded.
 *
 * @var TokenizeResult EOF_BEFORE_SUCCESS = 1024;
 * @brief End of file was reached before the read/skip operation succeeded.
 *
 * @var TokenizeResult NO_SUCCESS = 1025;
 * @brief The pattern was not found but we are not at EOF (may be returned when tokenizing on limited scope.
 */

/**
.Enum.TokenizeResult
..cat:Input/Output
..summary:Enum with return values for Tokenizing operations.
..value.SUCCESS:Reading the specified data succeeded.
..value.1-1023:File Error passed through
..value.EOF_BEFORE_SUCCESS:End of file was reached before the pattern was found.
..value.NO_SUCCESS:The pattern was not found, but we are not at EOF (may be returned when tokenizing on limited scope)
..include:seqan/stream.h
 */

enum TokenizeResult {
    SUCCESS = 0,
    EOF_BEFORE_SUCCESS = 1024,
    NO_SUCCESS = 1025
};

// ----------------------- Helper structs / private tags -----------------------
struct Whitespace__;
struct Blank__;
struct Char__;
struct Digit__;
struct Alpha__;
struct AlphaNum__;
struct UnixEOL__;
struct BackslashR__;
struct Graph__;
struct TabOrLineBreak__;
struct Identifier__;

typedef Tag<Whitespace__> Whitespace_;
typedef Tag<Blank__> Blank_;
typedef Tag<Char__> Char_;
typedef Tag<Digit__> Digit_;
typedef Tag<Alpha__> Alpha_;
typedef Tag<AlphaNum__> AlphaNum_;
typedef Tag<UnixEOL__> UnixEOL_;
typedef Tag<BackslashR__> BackslashR_;
typedef Tag<Graph__> Graph_;
typedef Tag<TabOrLineBreak__> TabOrLineBreak_;
typedef Tag<Identifier__> Identifier_;

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

/*!
 * @defgroup FileFormatTokenization File Format Tokenization
 * @brief Helper code for tokenization when reading from @link RecordReader @endlink.
 */

// ----------------------------------------------------------------------------
// helper functions for tags
// ----------------------------------------------------------------------------

inline int
_charCompare(int const c, Whitespace_ const & /* tag*/)
{
    return isspace(c);
}

inline int
_charCompare(int const c, Blank_ const & /* tag*/)
{
    return isblank(c);
}

inline int
_charCompare(int const c, Alpha_ const & /* tag*/)
{
    return isalpha(c);
}

inline int
_charCompare(int const c, AlphaNum_ const & /* tag*/)
{
    return isalnum(c);
}

inline int
_charCompare(int const c, Identifier_ const & /* tag*/)
{
    return isalnum(c) || c == '-' || c == '_';
}

inline int
_charCompare(int const c, Digit_ const & /* tag*/)
{
    return isdigit(c);
}

inline int
_charCompare(int const c, Graph_ const & /* tag*/)
{
    return isgraph(c);
}

inline int
_charCompare(int const c, UnixEOL_ const & /* tag*/)
{
    return (c == '\n');
}

inline int
_charCompare(int const c, BackslashR_ const & /* tag*/)
{
    return (c == '\r');
}

inline int
_charCompare(int const c, Tag<AminoAcid_> const & /* tag*/)
{
    switch (c)
    {
        case '*':
        case 'A':
        case 'B':
        case 'C':
        case 'D':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'I':
        case 'K':
        case 'L':
        case 'M':
        case 'N':
        case 'P':
        case 'Q':
        case 'R':
        case 'S':
        case 'T':
        case 'V':
        case 'W':
        case 'X':
        case 'Y':
        case 'Z':
        case 'a':
        case 'b':
        case 'c':
        case 'd':
        case 'e':
        case 'f':
        case 'g':
        case 'h':
        case 'i':
        case 'k':
        case 'l':
        case 'm':
        case 'n':
        case 'p':
        case 'q':
        case 'r':
        case 's':
        case 't':
        case 'v':
        case 'w':
        case 'x':
        case 'y':
        case 'z':
            return true;
    }
    return false;
}

inline int
_charCompare(int const c, Tag<Dna_> const & /* tag*/)
{
    switch (c)
    {
        case 'a':
        case 'c':
        case 'g':
        case 't':
        case 'A':
        case 'C':
        case 'G':
        case 'T': return true;
    }
    return false;
}

inline int
_charCompare(int const c, Tag<Rna_> const & /* tag*/)
{
    switch (c)
    {
        case 'a':
        case 'c':
        case 'g':
        case 'u':
        case 'A':
        case 'C':
        case 'G':
        case 'U': return true;
    }
    return false;
}

inline int
_charCompare(int const c, Tag<DnaQ_> const & /* tag*/)
{
    return _charCompare(c, Tag<Dna_>());
}

inline int
_charCompare(int const c, Tag<Dna5_> const & /* tag*/)
{
    switch (c)
    {
        case 'a':
        case 'c':
        case 'g':
        case 't':
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'n':
        case 'N':
            return true;
    }
    return false;
}

inline int
_charCompare(int const c, Tag<Rna5_> const & /* tag*/)
{
    switch (c)
    {
        case 'a':
        case 'c':
        case 'g':
        case 'u':
        case 'A':
        case 'C':
        case 'G':
        case 'U':
        case 'n':
        case 'N':
            return true;
    }
    return false;
}

inline int
_charCompare(int const c, Tag<Dna5Q_> const & /* tag*/)
{
    return _charCompare(c, Tag<Dna5_>());
}

template <typename TSpec>
inline int
_charCompare(int const c,
             Tag<SimpleType<unsigned char, TSpec> > const & /* tag */ )
{
    return _charCompare(c, Tag<TSpec>());
}

inline int
_charCompare(int const c, TabOrLineBreak_ const & /* tag*/)
{
    return c == '\t' || c == '\r' || c == '\n';
}

//TODO(h4nn3s): add for AminoAcid and Rna-tags


// ----------------------------------------------------------------------------
// Function _readHelper() [other read functions use this]
// ----------------------------------------------------------------------------

// read chars from record reader depending on condition
template <typename TTagSpec, // specialization of character comparison
          typename TRecordReader, // record reader
          typename TBuffer> // usually charstring, but maybe DnaString or so
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TTagSpec> const & tag,
            bool const desiredOutcomeOfComparison) 
/*   desired behaviour of loop -> "readUntil()" or "readwhile()"  */
{
    typedef char TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (bool( _charCompare(c, tag)) == desiredOutcomeOfComparison)
            return 0;
        appendValue(buffer, c, Generous());
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

template <typename TSpec, // specialization of character comparison
          typename TRecordReader, // record reader
          typename TBuffer> // usually charstring, but maybe DnaString or so
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TSpec> const & tag)
/*   default behaviour of loop is "readUntil()" */
{
    return _readHelper(buffer, reader, tag, true);
}

// same as above but allows for a second character type to be ignored
template <typename TTagSpec,
          typename TTagSpec2, // specialization of character class to be ignored
          typename TRecordReader,
          typename TBuffer>
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TTagSpec> const & compTag,
            Tag<TTagSpec2> const & skipTag,
            bool const desiredOutcomeOfComparison)
{
    typedef char TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (!_charCompare(c, skipTag))
        {
            if (bool (_charCompare(c, compTag)) == desiredOutcomeOfComparison)
                return 0;
            appendValue(buffer, c, Generous());
        }
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}


template <typename TTagSpec,
          typename TTagSpec2, // specialization of character class to be ignored
          typename TRecordReader,
          typename TBuffer>
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TTagSpec> const & compTag,
            Tag<TTagSpec2> const & skipTag)

{
    return _readHelper(buffer, reader, compTag, skipTag, true);
}


// ----------------------------------------------------------------------------
// Function _skipHelper() [other skip functions use this]
// ----------------------------------------------------------------------------


// same as above, just don't save characters read
template <typename TTagSpec, typename TRecordReader>
inline int
_skipHelper(TRecordReader & reader,
            Tag<TTagSpec> const & tag,
            bool const desiredOutcomeOfComparison)
{
    typedef char TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (bool(_charCompare(c, tag)) == desiredOutcomeOfComparison)
            return 0;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}


template <typename TTagSpec, typename TRecordReader>
inline int
_skipHelper(TRecordReader & reader,
            Tag<TTagSpec> const & tag)
{
    return _skipHelper(reader, tag, true);
}

// ----------------------------------------------------------------------------
// Function _countHelper()
// ----------------------------------------------------------------------------

// same as skip, but count the characters read (for doublepass-io)
template <typename TTagSpec, typename TRecordReader>
inline int
_countHelper(unsigned & count,
            TRecordReader & reader,
            Tag<TTagSpec> const & tag,
            bool const desiredOutcomeOfComparison)
{
    typedef char TChar;
    count = 0;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        count++;
        if (bool(_charCompare(c, tag)) == desiredOutcomeOfComparison)
            return 0;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

template <typename TTagSpec, typename TRecordReader>
inline int
_countHelper(unsigned & count,
            TRecordReader & reader,
            Tag<TTagSpec> const & tag)
{
    return _countHelper(count, reader, tag);
}

// same as above but allows for a second character type to be ignored
template <typename TTagSpec,
          typename TTagSpec2, // specialization of character class to be ignored
          typename TRecordReader>
inline int
_countHelper(unsigned & count,
            TRecordReader & reader,
            Tag<TTagSpec> const & compTag,
            Tag<TTagSpec2> const & skipTag,
            bool const desiredOutcomeOfComparison)
{
    count = 0;
    typedef char TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (!_charCompare(c, skipTag))
        {
            if (bool (_charCompare(c, compTag)) == desiredOutcomeOfComparison)
                return 0;
            ++count;
        }
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}


template <typename TTagSpec,
          typename TTagSpec2, // specialization of character class to be ignored
          typename TRecordReader>
inline int
_countHelper(unsigned & count,
            TRecordReader & reader,
            Tag<TTagSpec> const & compTag,
            Tag<TTagSpec2> const & skipTag)

{
    return _countHelper(count, reader, compTag, skipTag, true);
}




// ----------------------------------------------------------------------------
// Function _readAndCompareWithStr()
// ----------------------------------------------------------------------------



template <typename TRecordReader, typename TString>
inline int
_readAndCompareWithStr(TRecordReader & reader,
                       TString const & str)
{
    bool win = false;
    for (unsigned i = 0; i < length(str); ++i)
    {
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;
        if (value(reader) != str[i]) //mismatch
            break;
        if (i == length(str)-1) // win
            win = true;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return win ? 0 : NO_SUCCESS;
}



// ----------------------- the real deal ------------------------------------

/*!
 * @fn FileFormatTokenization#readUntilOneOf
 * @headerfile <seqan/stream.h>
 * @brief Read characters from RecordReader into buffer until one of the given character is encountered.
 *
 * @signature int readUntilOneOf(buffer, reader, c1[, c2[, c3[, c4[, c5]]]]);
 *
 * @param[in,out] buffer The buffer to write to.  Type: @link SequenceConcept @endlink.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     c1     One of the characters to look for.  Type: <tt>char</tt>.
 * @param[in]     c2     One of the characters to look for.  Type: <tt>char</tt>.
 * @param[in]     c3     One of the characters to look for.  Type: <tt>char</tt>.
 * @param[in]     c4     One of the characters to look for.  Type: <tt>char</tt>.
 * @param[in]     c5     One of the characters to look for.  Type: <tt>char</tt>.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>on</b> the found character and the character itself is not written to the buffer.  Note that
 * even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF) occured.
 */

/**
.Function.readUntilOneOf
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream into buffer until one of the given characters is encountered
..signature:readUntilOneOf(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader, c1[, c2[, c3[, c4[, c5]]])
..param.buffer:The buffer to write to.
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *on* the found character. The character is not written to buffer.
..include:seqan/stream.h
..see:Enum.TokenizeResult
 */

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilOneOf(TBuffer & buffer, RecordReader<TStream, TPass> & reader, char c1)
{
    while (!atEnd(reader))
    {
        char c = value(reader);
        if (c == c1)
            return 0;
        appendValue(buffer, c, Generous());
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilOneOf(TBuffer & buffer, RecordReader<TStream, TPass> & reader, char c1, char c2)
{
    while (!atEnd(reader))
    {
        char c = value(reader);
        if (c == c1 || c == c2)
            return 0;
        appendValue(buffer, c, Generous());
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilOneOf(TBuffer & buffer, RecordReader<TStream, TPass> & reader, char c1, char c2, char c3)
{
    while (!atEnd(reader))
    {
        char c = value(reader);
        if (c == c1 || c == c2 || c == c3)
            return 0;
        appendValue(buffer, c, Generous());
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilOneOf(TBuffer & buffer, RecordReader<TStream, TPass> & reader, char c1, char c2, char c3, char c4)
{
    while (!atEnd(reader))
    {
        char c = value(reader);
        if (c == c1 || c == c2 || c == c3 || c == c4)
            return 0;
        appendValue(buffer, c, Generous());
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilOneOf(TBuffer & buffer, RecordReader<TStream, TPass> & reader, char c1, char c2, char c3, char c4, char c5)
{
    while (!atEnd(reader))
    {
        char c = value(reader);
        if (c == c1 || c == c2 || c == c3 || c == c4 || c == c5)
            return 0;
        appendValue(buffer, c, Generous());
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

/*!
 * @fn FileFormatTokenization#readUntilWhitespace
 * @headerfile <seqan/stream.h>
 * @brief Read characters from stream into buffer until a @link isspace whitespace @endlink is encountered.
 *
 * @signature int readUntilBlank(buffer, reader);
 *
 * @param[in,out] buffer The buffer to write to.  Type: @link SequenceConcept @endlink.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * Whitespace is more than <tt>' '</tt> and <tt>'\t'</tt>, see @link isspace @endlink.  Consider using @link
 * FileFormatTokenization#readUntilBlank @endlink for reading until space or tab.
 *
 * This function stops <b>on</b> the found character and the character itself is not written to the buffer.  Note that
 * even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF) occured.
 */

/**
.Function.readUntilWhitespace
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream into buffer until Whitespace is encountered
..signature:readUntilWhitespace(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:Whitespace is more than '' and '\t', see @Function.isspace@
..remarks:This function stops *on* the whitespace character. The whitespace is not written to buffer.
..include:seqan/stream.h
..see:Function.isspace
..see:Function.skipUntilWhitespace
..see:Enum.TokenizeResult
 */

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilWhitespace(TBuffer & buffer,
                    RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _readHelper(buffer,
                       reader,
                       Whitespace_());
}

/*!
 * @fn FileFormatTokenization#readUntilBlank
 * @headerfile <seqan/stream.h>
 * @brief Read characters from stream into buffer until a @link isblank blank @endlink is encountered.
 *
 * @signature int readUntilBlank(buffer, reader);
 *
 * @param[in,out] buffer The buffer to write to.  Type: @link SequenceConcept @endlink.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>on</b> the found character and the character itself is not written to the buffer.  Note that
 * even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF) occured.
 */

/**
.Function.readUntilBlank
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream into buffer until Blank is encountered
..signature:readUntilBlank(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:Blank is '' and '\t', see @Function.isblank@
..remarks:This function stops *on* the blank character. The blank is not written to buffer.
..include:seqan/stream.h
..see:Function.isblank
..see:Function.skipUntilBlank
..see:Enum.TokenizeResult
 */
template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilBlank(TBuffer & buffer,
               RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _readHelper(buffer,
                       reader,
                       Blank_());
}

/*!
 * @fn FileFormatTokenization#readUntilChar
 * @headerfile <seqan/stream.h>
 * @brief Read characters from stream into buffer until a given character is read.
 *
 * @signature int readUntilChar(buffer, reader, c);
 *
 * @param[in,out] buffer The buffer to write to.  Type: @link SequenceConcept @endlink.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     c      The <tt>char</tt> to look for.
 *
 * @return int 0 if there was no error and non-0 if there were errors.
 *             A special value is @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>on</b> the found character and the character itself is not written to the buffer.  Note that
 * even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF) occured.
 */

/**
.Function.readUntilChar
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream into buffer until Char is encountered
..signature:readUntilChar(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader, TChar const & x)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.x:The character to stop on
...type:nolink:$char$ or similar
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *on* the character x. It is not written to buffer.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.skipUntilChar
 */
template <typename TBuffer, typename TStream, typename TPass, typename TCharX>
inline int
readUntilChar(TBuffer & buffer,
              RecordReader<TStream, TPass> & reader,
              TCharX const & x)
{
    SEQAN_CHECKPOINT
    typedef char TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (c == x)
            return 0;
        appendValue(buffer, c, Generous());
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

/*!
 * @fn FileFormatTokenization#readNChars
 * @headerfile <seqan/stream.h>
 * @brief Read fixed number of characters from stream into buffer until a given character is read.
 *
 * @signature int readUntilChar(buffer, reader, n);
 *
 * @param[in,out] buffer The buffer to write to.  Type: @link SequenceConcept @endlink.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     n      The number of characters to read.  Type: <tt>unsigned</tt>.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops after <tt>n</tt> read characters.  Note that even when an error occurs, <tt>buffer</tt> will
 * contain the characters read until the error (or EOF) occured.
 */

/**
.Function.readNChars
..class:Class.RecordReader
..cat:Input/Output
..summary:Read exactly n characters from stream into buffer
..signature:readNChars(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader, unsigned const n)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.n:The number of characters to read
...type:nolink:$unsigned$
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..include:seqan/stream.h
..see:Enum.TokenizeResult
 */
template <typename TBuffer, typename TStream, typename TPass>
inline int
readNChars(TBuffer & buffer,
           RecordReader<TStream, TPass> & reader,
           unsigned const n)
{
    SEQAN_CHECKPOINT


    reserve(buffer, n, Exact());

    for (unsigned i = 0; i < n; ++i)
    {
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;
        appendValue(buffer, value(reader));
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return 0;
}

/*!
 * @fn FileFormatTokenization#readNCharsIgnoringWhitespace
 * @headerfile <seqan/stream.h>
 * @brief Read fixed number of non-whitespace characters.
 *
 * @signature int readUntilChar(buffer, reader, n);
 *
 * @param[in,out] buffer The buffer to write to.  Type: @link SequenceConcept @endlink.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     n      The number of non-whitespace characters to read.  Type: <tt>unsigned</tt>.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * Whitespace is more than <tt>' '</tt> and <tt>'\t'</tt>, see @link isspace @endlink.
 *
 * This function stops after <tt>n</tt> read non-whitespace characters.  Note that even when an error occurs,
 * <tt>buffer</tt> will contain the characters read until the error (or EOF) occured.
 */

/**
.Function.readNCharsIgnoringWhitespace
..class:Class.RecordReader
..cat:Input/Output
..summary:Read n characters from stream into buffer, but skip certain Chars
..signature:readNCharsIgnoringWhitespace(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader, unsigned const n)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.n:The number of characters to read
...type:nolink:$unsigned$
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..include:seqan/stream.h
..remarks:Whitespace characters are not counted
..see:Enum.TokenizeResult
 */


template <typename TBuffer, typename TStream, typename TPass,
          typename TIgnoredType>
inline int
_readNCharsIgnoringType(TBuffer & buffer,
                         RecordReader<TStream, TPass> & reader,
                         unsigned const n,
                         TIgnoredType /* tag */)
{
    SEQAN_CHECKPOINT


    reserve(buffer, n, Exact());

    for (unsigned i = 0; i < n; ++i)
    {
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;

        if (_charCompare(value(reader), TIgnoredType()))
            --i;
        else
            appendValue(buffer, value(reader));

        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return 0;
}

template <typename TBuffer, typename TStream, typename TPass>
inline int
readNCharsIgnoringWhitespace(TBuffer & buffer,
                             RecordReader<TStream, TPass> & reader,
                             unsigned const n)
{
    return _readNCharsIgnoringType(buffer, reader, n, Whitespace_());
}

/*!
 * @fn FileFormatTokenization#skipNChars
 * @headerfile <seqan/stream.h>
 * @brief Skip a fixed number of characters.
 *
 * @signature int skipNChars(reader, n);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     n      The number of characters to read.  Type: <tt>unsinged</tt>.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 */

/**
.Function.skipNChars
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip exactly n characters from stream
..signature:skipNChars(RecordReader<TStream, TPass> & recordReader, unsigned const n)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.n:The number of characters to skip
...type:nolink:$unsigned$
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..include:seqan/stream.h
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass>
inline int
skipNChars(RecordReader<TStream, TPass> & reader,
           unsigned const n)
{
    SEQAN_CHECKPOINT

    for (unsigned i = 0; i < n; ++i)
    {
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return 0;
}

/*!
 * @fn FileFormatTokenization#skipNCharsIgnoringWhitespace
 * @headerfile <seqan/stream.h>
 * @brief Skip a fixed number of non-whitespace characters.
 *
 * @signature int skipNChars(reader, n);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     n      The number of characters to read.  Type: <tt>unsinged</tt>.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * Whitespace is more than <tt>' '</tt> and <tt>'\t'</tt>, see @link isspace @endlink.
 */


/**
.Function.skipNCharsIgnoringWhitespace
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip n characters from stream, not counting whitespaces
..signature:skipNCharsIgnoringWhitespace(RecordReader<TStream, TPass> & recordReader, unsigned const n)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.n:The number of characters to skip
...type:nolink:$unsigned$
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..include:seqan/stream.h
..remarks:Whitespace characters are not counted towards n
..see:Enum.TokenizeResult
..see:Function.skipNChars
..see:Function.readNCharsIgnoringWhitespace
 */
template <typename TStream, typename TPass,
          typename TIgnoredType>
inline int
_skipNCharsIgnoringType(RecordReader<TStream, TPass> & reader,
                        unsigned const n,
                        TIgnoredType const & /* tag */)
{
    for (unsigned i = 0; i < n; ++i)
    {
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;
        if (_charCompare(value(reader), TIgnoredType()))
            --i;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return 0;
}

template <typename TStream, typename TPass>
inline int
skipNCharsIgnoringWhitespace(RecordReader<TStream, TPass> & reader,
                             unsigned const n)
{
    return _skipNCharsIgnoringType(reader, n, Whitespace_());
}

/*!
 * @fn FileFormatTokenization#skipUntilWhitespace
 * @headerfile <seqan/stream.h>
 * @brief Skip until the first whitespace character is encountered.
 *
 * @signature int skipUntilWhitespace(reader);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * Whitespace is more than <tt>' '</tt> and <tt>'\t'</tt>, see @link isspace @endlink.  Consider using @link
 * FileFormatTokenization#skipUntilBlank @endlink for skipping until space or tab.
 *
 * The reader will stop <b>on</b> the first whitespace character.
 */

/**
.Function.skipUntilWhitespace
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until Whitespace is encountered
..signature:skipUntilWhitespace(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:Whitespace is more than '' and '\t', see @Function.isspace@
..remarks:This function stops *on* the whitespace character. The whitespace is not skipped.
..include:seqan/stream.h
..see:Function.isspace
..see:Function.skipWhitespaces
..see:Function.readUntilWhitespace
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass>
inline int
skipUntilWhitespace(RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _skipHelper(reader, Whitespace_());
}

/*!
 * @fn FileFormatTokenization#skipUntilBlank
 * @headerfile <seqan/stream.h>
 * @brief Skip until the first blank character is encountered.
 *
 * @signature int skipUntilBlank(reader);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The reader will stop <b>on</b> the first blank character.
 */

/**
.Function.skipUntilBlank
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until Blank is encountered
..signature:skipUntilBlank(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *on* the blank character. The blank is not skipped.
..include:seqan/stream.h
..see:Function.isblank
..see:Enum.TokenizeResult
..see:Function.readUntilBlank
 */

template <typename TStream, typename TPass>
inline int
skipUntilBlank(RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _skipHelper(reader, Blank_());
}

/*!
 * @fn FileFormatTokenization#skipUntilGraph
 * @headerfile <seqan/stream.h>
 * @brief Skip until the first printable character is encountered.
 *
 * @signature int skipUntilGraph(reader);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The reader will stop <b>on</b> the first printable character.
 *
 * @see isgraph
 */

/**
.Function.skipUntilGraph
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until printable, non-' ' character is encountered
..signature:skipUntilGraph(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *on* the graph character. The graph is not skipped.
..remarks:See @Function.isgraph@ for details on the "graph"-group of characters
..include:seqan/stream.h
..see:Function.isgraph
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass>
inline int
skipUntilGraph(RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _skipHelper(reader, Graph_());
}

/*!
 * @fn FileFormatTokenization#skipUntilChar
 * @headerfile <seqan/stream.h>
 * @brief Skip until the a given character is found.
 *
 * @signature int skipUntilChar(reader, c);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     c      The <tt>char</tt> to search for.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The reader will stop <b>on</b> the first found location of <tt>c</tt>.
 */

/**
.Function.skipUntilChar
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until Char is encountered
..signature:skipUntilChar(RecordReader<TStream, TPass> & recordReader, TChar const & x)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.x:The character to stop on
...type:nolink:$char$ or similar
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *on* the character x. x is not skipped.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.readUntilChar
 */
template <typename TStream, typename TPass, typename TCharX>
inline int
skipUntilChar(RecordReader<TStream, TPass> & reader,
              TCharX const & x)
{
    SEQAN_CHECKPOINT
    typedef char TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (c == x)
            return 0;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return EOF_BEFORE_SUCCESS;
}

/*!
 * @fn FileFormatTokenization#skipUntilString
 * @headerfile <seqan/stream.h>
 * @brief Skip until the a given string is found.
 *
 * @signature int skipUntilString(reader, str);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     str    The @link SequenceConcept @endlink to find.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The reader will stop <b>behind</b> the location of <tt>str</tt>.
 */

/**
.Function.skipUntilString
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until String is encountered
..signature:skipUntilString(RecordReader<TStream, TPass> & recordReader, TString const & str)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.str:The string to stop on
...type:nolink:f$char$ or similar
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the character string
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.skipUntilChar
 */
// I think this function is stupid, but it is implemented in misc_parsing.h
// so I offer replacement; with proper buffering shiftOr would be faster
template <typename TStream, typename TPass, typename TString>
inline int
skipUntilString(RecordReader<TStream, TPass> & reader,
                             TString const & str)
{
    SEQAN_CHECKPOINT
    if (length(str) < 1)
        return -1; //TODO some better error code
    while(skipUntilChar(reader, str[0])==0)
    {
        switch(int r = _readAndCompareWithStr(reader, str))
        {
            case 0: return 0;
            case NO_SUCCESS: break;
            default: return r;
        }
    }
    return EOF_BEFORE_SUCCESS;
    /* NOTE I am not 100% sure this will get every pattern with repitions in it.
     * it should be much better than the original, since it escapes on mismatch
     * directly and doesnt read (and possibly discard) N charecters in a row,
     * which definitely leads to misses. */
}

/*!
 * @fn FileFormatTokenization#skipTabOrLineBreak
 * @headerfile <seqan/stream.h>
 * @brief Skip over tab or line break characters (<tt>'\r' or '\n'</tt>).
 *
 * @signature int skipTabOrLineBreak(reader);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The reader will stop <b>behind</b> the (possibly empty) sequence of tab or line break characters (regex:
 * <tt>[\r\n\t]</tt>).
 */

/**
.Function.readUntilTabOrLineBreak
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream until a tab or line-break occurs.
..signature:readUntilTabOrLineBreak(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to.
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *on* the tab or line-break character. This character is not written to buffer.
..include:seqan/stream.h
 */

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilTabOrLineBreak(TBuffer & buffer,
                        RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _readHelper(buffer,
                       reader,
                       TabOrLineBreak_());
}

// ---

/*!
 * @fn FileFormatTokenization#readLetters
 * @headerfile <seqan/tokenize.h>
 * @brief Read from a RecordReader as long as the characters are letters (in @link isalpha alpha @endlink class).
 *
 * @signature int readLetters(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to read the data into.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The function will stop <b>after</b> the last read character.
 * 
 * Note that even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF)
 * occured.
 */

/**
.Function.readLetters
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream as long as characters are letters
..signature:readLetters(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:Shortcut.DnaString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isalpha
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readLetters(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _readHelper(buffer,
                       reader,
                       Alpha_(),
                       false);
}

/*!
 * @fn FileFormatTokenization#readDigits
 * @headerfile <seqan/tokenize.h>
 * @brief Read from a RecordReader as long as the characters are digits (in @link isdigit digit @endlink class).
 *
 * @signature int readDigits(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to read the data into.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The function will stop <b>after</b> the last read character.
 * 
 * Note that even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF)
 * occured.
 */

/**
.Function.readDigits
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream as long as characters are digits
..signature:readDigits(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isdigit
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readDigits(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _readHelper(buffer,
                       reader,
                       Digit_(),
                       false);
}

/*!
 * @fn FileFormatTokenization#readGraphs
 * @headerfile <seqan/tokenize.h>
 * @brief Read from a RecordReader as long as the characters are printable (in @link isgraph graph @endlink class).
 *
 * @signature int readGraphs(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to read the data into.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The function will stop <b>after</b> the last read character.
 * 
 * Note that even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF)
 * occured.
 */

/**
.Function.readGraphs
..cat:Input/Output
..summary:Read characters from stream as long as characters are graph characters.
..signature:readGraphs(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isgraph
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readGraphs(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _readHelper(buffer, reader, Graph_(), false);
}

/*!
 * @fn FileFormatTokenization#readFloat
 * @headerfile <seqan/tokenize.h>
 * @brief Read characters from @link RecordReader @endlink as long as the string is a valid floating point numbers.
 *
 * @signature int readFloat(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to read the data into.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The function will stop <b>after</b> the last read character.
 * 
 * Note that even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF)
 * occured.
 */

/**
.Function.readFloat
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream as long as the number is a valid floating point numbers.
..description:Supports normal floating point numbers and scientific notation.
..signature:readFloat(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isdigit
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readFloat(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;

    // Read sign if any.
    if (value(reader) == '+' || value(reader) == '-')
        if (readNChars(buffer, reader, 1) != 0)
            return 1;

    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;

    int res = 0;

    // Either read [0-9]+(\.[0-9]*)? or \.[0-9]+.
    if (value(reader) == '.')
    {
        res = readNChars(buffer, reader, 1);
        if (res != 0)
            return res;
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;
        if (!isdigit(value(reader)))
            return 1;
        if ((res = readDigits(buffer, reader)) != 0)
            return res;
    }
    else
    {
        // A row of digits must follow.
        if (!isdigit(value(reader)))
            return 1;  // Invalid.
        res = readDigits(buffer, reader);
        // The input may end after these digits.
        if (res == EOF_BEFORE_SUCCESS)
            return 0;
        if (res != 0)
            return 1;
        // Read \.[0-9]* block.
        if (value(reader) == '.')
        {
            res = readNChars(buffer, reader, 1);
            // The input may end after the dot.
            if (res == EOF_BEFORE_SUCCESS)
                return 0;
            if (res != 0)
                return 1;
            if (isdigit(value(reader)))
                if ((res = readDigits(buffer, reader)) != 0)
                    return res;
        }
    }

    if (value(reader) != 'e' && value(reader) != 'E')
        return 0;
    if ((res = readNChars(buffer, reader, 1)) != 0)
        return res;
    if (value(reader) == '+' || value(reader) == '-')
        if ((res = readNChars(buffer, reader, 1)) != 0)
            return res;
    if (!isdigit(value(reader)))
        return 1;
    if ((res = readDigits(buffer, reader)) != 0)
        return res;

    return 0;
}

/*!
 * @fn FileFormatTokenization#readAlphaNums
 * @headerfile <seqan/tokenize.h>
 * @brief Read from a RecordReader as long as the characters are alphanumeric (in @link isalnum alnum @endlink class).
 *
 * @signature int readAlphaNums(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to read the data into.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The function will stop <b>after</b> the last read character.
 * 
 * Note that even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF)
 * occured.
 */

/**
.Function.readAlphaNums
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream as long as characters are letters
..signature:readAlphaNums(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isalnum
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readAlphaNums(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _readHelper(buffer,
                       reader,
                       AlphaNum_(),
                       false);
}

/*!
 * @fn FileFormatTokenization#readIdentifier
 * @headerfile <seqan/tokenize.h>
 * @brief Read from a RecordReader as long as the characters form an identifier (alnum, <tt>'-'</tt>, <tt>'_'</tt>).
 *
 * @signature int readIdentifier(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to read the data into.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * The function will stop <b>after</b> the last read character.
 * 
 * Note that even when an error occurs, <tt>buffer</tt> will contain the characters read until the error (or EOF)
 * occured.
 */

/**
.Function.readIdentifier
..cat:Input/Output
..summary:Read characters from stream as long as characters are identifiers (alphanumeric, $'-'$, and $'_'$).
..signature:readIdentifier(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isalnum
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readIdentifier(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    return _readHelper(buffer, reader, Identifier_(), false);
}

/*!
 * @fn FileFormatTokenization#skipWhitespaces
 * @headerfile <seqan/stream.h>
 * @brief Advance RecordReader until a whitespace occurs.
 *
 * @signature int skipWhitespaces(reader);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * Whitespace is more than <tt>' '</tt> and <tt>'\t'</tt>, see @link isspace @endlink.  Consider using @link
 * FileFormatTokenization#skipBlanks @endlink for skipping blanks.
 *
 * This function stops <b>after</b> the last skipped.
 */

/**
.Function.skipWhitespaces
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until non-Whitespace is encountered
..signature:skipWhitespaces(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:Whitespace is more than '' and '\t', see @Function.isspace@
..remarks:This function stops *behind* the last whitespace character.
..include:seqan/stream.h
..see:Function.isspace
..see:Function.readUntilWhitespace
..see:Function.skipUntilWhitespace
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass>
inline int
skipWhitespaces(RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _skipHelper(reader, Whitespace_(), false);
}

/*!
 * @fn FileFormatTokenization#skipChar
 * @headerfile <seqan/stream.h>
 * @brief Skip one character that must be equal to a given one for this function to succeed.
 *
 * @signature int skipChar(reader, c);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     c      The <tt>char</tt> to skip.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 */

/**
.Function.skipChar
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip one character that must be equal to a given one for this function to succeed.
..signature:skipChar(RecordReader<TStream, TPass> & recordReader, c)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.c:The character to skip.
...type:nolink:$char$
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..include:seqan/stream.h
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass>
inline int
skipChar(RecordReader<TStream, TPass> & reader, char const c)
{
    SEQAN_CHECKPOINT;
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (value(reader) != c)
        return 1;
    goNext(reader);
    return 0;
}

/*!
 * @fn FileFormatTokenization#skipBlanks
 * @headerfile <seqan/stream.h>
 * @brief Advance RecordReader until a blank occurs.
 *
 * @signature int skipBlanks(reader);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>after</b> the last skipped.
 */

/**
.Function.skipBlanks
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until non-Blank is encountered
..signature:skipBlanks(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops *behind* the last whitespace character.
..include:seqan/stream.h
..see:Function.isblank
..see:Function.readUntilBlank
..see:Function.skipUntilBlank
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass>
inline int
skipBlanks(RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    return _skipHelper(reader, Blank_(), false);
}

/*!
 * @fn FileFormatTokenization#readLine
 * @headerfile <seqan/stream.h>
 * @brief Read a line from a RecordReader into a buffer.
 *
 * @signature int readLine(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to append the read data to.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops on the beginning of the next line (if there is any).  End of line characters are not written to
 * buffer.  Works on ANSI, Mac, and Unix EOL.
 */

/**
.Function.readLine
..class:Class.RecordReader
..cat:Input/Output
..summary:Read a line from stream and save it to buffer
..signature:readLine(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops on the beginning of the next line, if there is a next line
..remarks:End-line characters are not written to buffer.
..remarks:Works on ANSI EOL, Mac EOL and on Unix EOL.
..include:seqan/stream.h
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readLine(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    int r = 0;

    // Loop over characters from stream, taking early exits on line breaks and errors.
    while (!atEnd(reader))
    {
        char c = value(reader);

        // Unix EOL is the simplest case.
        if (c == '\n')
        {
            if (!atEnd(reader))
                goNext(reader);
            return resultCode(reader);
        }

        // If the current character is '\r' then this can be an ANSI or a Mac line ending.
        if (c == '\r')
        {
            goNext(reader);
            if ((r = resultCode(reader)) != 0)
                return r;
            if (atEnd(reader))
                return 0;  // Assume Mac EOL at end of file.
            if (value(reader) == '\n')
            {
                // Assume Windows EOL.
                if (!atEnd(reader))
                    goNext(reader);
                return resultCode(reader);
            }
            return 0;  // Was Mac EOL, not at end of file.
        }

        appendValue(buffer, c);
        goNext(reader);
    }

    return EOF_BEFORE_SUCCESS;
}

/*!
 * @fn FileFormatTokenization#readLineStripTrailingBlanks
 * @headerfile <seqan/stream.h>
 * @brief Read a line from a RecordReader into a buffer and remove trailing blanks.
 *
 * @signature int readLineStripTrailingBlanks(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to append the read data to.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops on the beginning of the next line (if there is any).  End of line characters and all trailing
 * whitespaces are not written to buffer.  Works on ANSI, Mac, and Unix EOL.
 */

/**
.Function.readLineStripTrailingBlanks
..class:Class.RecordReader
..cat:Input/Output
..summary:Read a line from stream and save it to buffer, remove trailing blanks
..signature:readLineStripTrailingBlanks(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops on the beginning of the next line, if there is a next line
..remarks:End-line characters and all trailing blanks are not written to buffer.
..remarks:Works on ANSI EOL and on Unix EOL.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isblank
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readLineStripTrailingBlanks(TBuffer & buffer,
                            RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    int r = readLine(buffer,reader);

    if (r != 0)
        return r;   // 1234567890

    typename Size<TBuffer>::Type pos = length(buffer) -1;

    if (pos < 0)
        return 0;

    while ((isblank(value(buffer, pos)) == 0) && (pos >= 0)) --pos;

    if (pos + 1 != length(buffer))
        resize(buffer, pos+1, Exact());
    return 0;
}

/*!
 * @fn RecordReader#skipLine
 * @headerfile <seqan/stream.h>
 * @brief Skip a line in stream and go to beginning of next
 * 
 * @signature int skipLine(reader);
 * 
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * 
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 * 
 * @section Remarks
 * 
 * This function stops on the beginning of the next line, if there is a next line
 * 
 * Works on ANSI EOL and on Unix EOL.
 */

/**
.Function.skipLine
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip a line in stream and go to beginning of next
..signature:skipLine(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops on the beginning of the next line, if there is a next line
..remarks:Works on ANSI EOL and on Unix EOL.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.readLine
 */
template <typename TStream, typename TPass>
inline int
skipLine(RecordReader<TStream, TPass> & reader)
{
    int r = _skipHelper(reader, UnixEOL_());

    if (r != 0)
        return r;

    if (!atEnd(reader))
        goNext(reader); // go to beginning of next line

    return resultCode(reader);
}

/*!
 * @fn RecordReader#countLine
 * @headerfile <seqan/stream.h>
 * @brief Count characters in a line excluding <tt>'\r'</tt> and <tt>'\n'</tt>.
 * 
 * @signature int countLine(num, reader);
 *
 * @param[in,out] num    The <tt>unsigned</tt> to increment.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * 
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 * 
 * @section Remarks
 * 
 * This function stops on the beginning of the next line, if there is a next line (even though newline characters are
 * not counted).
 * 
 * Works on ANSI EOL and on Unix EOL.
 */

/**
.Function.countLine
..class:Class.RecordReader
..cat:Input/Output
..summary:count characters in a line not including \r and \n
..signature:countLine(unsigned & count, RecordReader<TStream, TPass> & recordReader)
..param.count:The variable to increment
...type:nolink:unsigned int
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function stops on the beginning of the next line, if there is a next line (even though newline characters are not counted)
..remarks:Works on ANSI EOL and on Unix EOL.
..include:seqan/stream.h
..see:Enum.TokenizeResult
 */
template <typename TStream, typename TPass>
inline int
countLine(unsigned & count, RecordReader<TStream, TPass> & reader)
{
    SEQAN_CHECKPOINT
    int r = _countHelper(count,
                        reader,
                        UnixEOL_(), // abort on Newline
                        BackslashR_()); // just ignore carriage return
    if (r != 0)
        return r;

    if (!atEnd(reader))
        goNext(reader); // go to beginning of next line

    return resultCode(reader);
}

/*!
 * @fn FileFormatTokenization#readDna5IgnoringWhitespaces
 * @headerfile <seqan/stream.h>
 * @brief Read characters from stream as long as they are DNA5 characters, skipping whitespaces.
 *
 * @signature int readDna5IgnoringWhitespaces(buffer, reader);
 *
 * @param[in,out] buffer The @link String @endlink to append the read data to.
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>on</b> the first non-matching character, whitespaces such as newlines are discarded.
 */

/**
.Function.readDna5IgnoringWhitespaces
..class:Class.RecordReader
..cat:Input/Output
..summary:Read characters from stream, as long as they are DNA5 characters. Skip over whitespaces.
..signature:readDna5IgnoringWhitespaces(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.Dna5String
...type:Shortcut.CharString
...type:nolink:or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function read Dna5-Characters and simply discards whitespaces, including newlines.
..remarks:It stops on the first character that is not Dna5 and not whitespace.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.isspace
..see:Spec.Dna5
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readDna5IgnoringWhitespaces(TBuffer & buffer,
                            RecordReader<TStream, TPass> & reader)
{
    return _readHelper(buffer, reader, Tag<Dna5_>(), Whitespace_(), false);
} 
// this would read a fasta or fastq sequence, since meta and qualities begin 
// with special chars

/*!
 * @fn FileFormatTokenization#sipUntilLineBeginsWithChar
 * @headerfile <seqan/stream.h>
 * @brief Skip input until the first graphical (see @link isgraph @endlink) character of a line is equal to c.
 *
 * @signature int skipUntilLineBeginsWithChar(reader, c);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     c      The <tt>char</tt> to look for.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>on</b> the first occurence of <tt>c</tt> at the beginning of a line (ignoring leading
 * whitespace).
 */

/**
.Function.skipUntilLineBeginsWithChar
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip input until the first graphical(see @Function.isgraph@) character of a line is equal to c
..signature:skipUntilLineBeginsWithChar(RecordReader<TStream, TPass> & recordReader, TChar const & c)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.c:The character to look for, must be in @Function.isgraph@
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function skips lines, and then non-graphical chars at the line beginnings,
..remarks:until it finds a line where the first graph-character is c. It stops on c, not behind.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.skipUntilLineBeginsWithStr
 */
// also skips non-graph characters at linestart, c must not be non-graph
// NOTE that position is on c, like in old function
// NOTE also, that we DO NOT reset in case of search failure as that would
// require seek
template <typename TStream, typename TPass, typename TChar>
inline int
skipUntilLineBeginsWithChar(RecordReader<TStream, TPass> & reader,
                            TChar const & c )
{
    int r = 0;
    while (!atEnd(reader) && (r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        if (value(reader) == c)
            return 0;
    }
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    return r;
}

/*!
 * @fn FileFormatTokenization#sipUntilLineBeginsWithStr
 * @headerfile <seqan/stream.h>
 * @brief Skip input until a line begins with str (str itself must begin with a printable character).
 *
 * @signature int skipUntilLineBeginsWithChar(reader, str);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     c      The @link SequenceConcept @endlink to look for.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>behind</b> the first occurence of <tt>str</tt> at the beginning of a line.  The function skips
 * over leading whitespace.
 */

/**
.Function.skipUntilLineBeginsWithStr
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip input until a line begins with str.
..signature:skipUntilLineBeginsWithStr(RecordReader<TStream, TPass> & recordReader, TString const & str)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.str:The string to look for, must begin with a char in @Function.isgraph@
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function skips lines, and then non-graphical chars (see @Function.isgraph@) at the line beginnings,
..remarks:until it finds a line where the first string beginning with graph-characters is equal to str.
..remarks:It stops behaind the string.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.readLine
 */
// also skips non-graph characters at linestart, str must not begin with non-graph
// NOTE that position is behind str, like in old function, if a char behind
// string exists. NOTE also, that we DO NOT reset in case of search failure
// as that would require seek
template <typename TStream, typename TPass, typename TString>
inline int
skipUntilLineBeginsWithStr(RecordReader<TStream, TPass> & reader,
                           TString const & str )
{
    int r = 0;
    while (!atEnd(reader) && (r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        switch(r = _readAndCompareWithStr(reader, str))
        {
            case 0: return 0;
            case NO_SUCCESS: break;
            default: return r;
        }
    }
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    return r;
}

/*!
 * @fn FileFormatTokenization#sipUntilLineBeginsWithOneCharOfStr
 * @headerfile <seqan/stream.h>
 * @brief Skip input until a line begins with one of the character in str.
 *
 * @signature int skipUntilLineBeginsWithOneCharOfStr(reader, str);
 *
 * @param[in,out] reader The @link RecordReader @endlink to read from.
 * @param[in]     str    The @link SequenceConcept @endlink with the characters to look for.
 *
 * @return int 0 if there was no error on reading or non-0 if there were errors.  A special value is
 *             @link TokenizeResult EOF_BEFORE_SUCCESS @endlink.
 *
 * @section Remarks
 *
 * This function stops <b>on</b> the first found occurence.  All characters in <tt>str</tt> must be printable (see @link
 * isgraph @endlink).
 */

/**
.Function.skipUntilLineBeginsWithOneCharOfStr
..class:Class.RecordReader
..cat:Input/Output
..summary:Skip input until a line begins with a one of the characters in str
..signature:skipUntilLineBeginsWithOneCharOfStr(RecordReader<TStream, TPass> & recordReader, TString const & str)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.str:The set of characters to look for, must be in @Function.isgraph@
..returns:0 if there was no error reading
..returns:non-zero value on errors, especially EOF_BEFORE_SUCCESS
...type:nolink:$int$
...type:Enum.TokenizeResult
..remarks:This function skips lines, and then non-graphical chars (see @Function.isgraph@) at the line beginnings,
..remarks:until it finds a line where the first graph-char is one of the characters in str.
..remarks:It stops on the char.
..include:seqan/stream.h
..see:Enum.TokenizeResult
..see:Function.readLine
 */
// also skips non-graph characters at linestart, c must not be non-graph
// NOTE that position is on c, like in old function
// NOTE also, that we DO NOT reset in case of search failure as that would
// require seek
template <typename TStream, typename TPass, typename TString>
inline int
skipUntilLineBeginsWithOneCharOfStr(RecordReader<TStream, TPass> & reader,
                                    TString const & str )
{

    int r = 0;
    while (!atEnd(reader) && (r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        for (unsigned i = 0; i < length(str); ++i)
        {
            if (value(reader) == value(str,i))
                return 0;
        }

    }
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    return r;
}

// TODO missing functions from misc_parsing.h:
// _parseUntilBeginLine with max num lines
//_parseUntilBeginLineOneOf(
// all functions that require lexical_cast<>() which isnt there yet

// TODO stuff from other files (shouldnt be that much and much is duplicate)

}  // namespace seqan

#endif // #ifdef SEQAN_STREAM_RECORD_READER_SINGLE_H_
