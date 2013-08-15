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
// Collection of char checking functions
// ==========================================================================

// ==========================================================================
// We add documentation for cctype here for completeness.
// ==========================================================================

/*!
 * @fn isalnum
 * @headerfile <cctype>
 * @brief Check if character is alphanumeric.
 * 
 * @signature bool isalnum(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is alphanumeric.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, plain c and listed for completeness.
 * 
 * @see http://www.cplusplus.com/reference/clibrary/cctype/isalnum/
 * @see http://www.cppreference.com/wiki/string/character_classes
 * @see FileFormatTokenization#readAlphaNums
 * @see FileFormatTokenization#readIdentifier
 */

/**
.Function.isalnum:
..summary:Check if character is alpha-numeric
..cat:Input/Output
..signature:isalnum(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks: This is non-seqan, plain c and listed for completeness
..see:http://www.cplusplus.com/reference/clibrary/cctype/isalnum/
..see:http://www.cppreference.com/wiki/string/character_classes
..include:<cctype>
*/

/*!
 * @fn isalpha
 * @headerfile <cctype>
 * @brief Check if character is a upper or lowercase letter
 * 
 * @signature bool isalpha(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is a character.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, plain c and listed for completeness.
 * 
 * @see http://www.cplusplus.com/reference/clibrary/cctype/isalpha/
 * @see http://www.cppreference.com/wiki/string/character_classes
 * @see FileFormatTokenization#readLetters
 */

/**
.Function.isalpha:
..summary:Check if character is a upper or lowercase letter
..cat:Input/Output
..signature:isalpha(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks: This is non-seqan, plain c and listed for completeness
..see:http://www.cplusplus.com/reference/clibrary/cctype/isalpha/
..see:http://www.cppreference.com/wiki/string/character_classes
..include:<cctype>
*/

/*!
 * @fn isdigit
 * @headerfile <cctype>
 * @brief Check if character is a digit
 * 
 * @signature bool isdigit(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is a digit.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, plain c and listed for completeness.
 * 
 * @see http://www.cplusplus.com/reference/clibrary/cctype/isdigit/
 * @see http://www.cppreference.com/wiki/string/character_classes
 * @see FileFormatTokenization#readDigits
 * @see FileFormatTokenization#readFloat
 */

/**
.Function.isdigit:
..summary:Check if character is a digit
..cat:Input/Output
..signature:isdigit(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks: This is non-seqan, plain c and listed for completeness
..see:http://www.cplusplus.com/reference/clibrary/cctype/isdigit/
..see:http://www.cppreference.com/wiki/string/character_classes
..include:<cctype>
*/

/*!
 * @fn iscntrl
 * @headerfile <cctype>
 * @brief Check if character is a control character
 * 
 * @signature bool iscntrl(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is a control character.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, plain c and listed for completeness.
 * 
 * @see http://www.cplusplus.com/reference/clibrary/cctype/iscntrl/
 * @see http://www.cppreference.com/wiki/string/character_classes
 */

/**
.Function.iscntrl:
..summary:Check if character is a control character
..cat:Input/Output
..signature:iscntrl(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks: This is non-seqan, plain c and listed for completeness
..see:http://www.cplusplus.com/reference/clibrary/cctype/iscntrl/
..see:http://www.cppreference.com/wiki/string/character_classes
..include:<cctype>
*/

/*!
 * @fn isprint
 * @headerfile <cctype>
 * @brief Check if character is printable, i.e. not a control character
 * 
 * @signature bool isprint(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is a printable character.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, plain c and listed for completeness.
 * 
 * @see http://www.cplusplus.com/reference/clibrary/cctype/isprint/
 * @see http://www.cppreference.com/wiki/string/character_classes
 * @see isgraph
 */

/**
.Function.isprint:
..summary:Check if character is printable, i.e. not a control character
..cat:Input/Output
..signature:isprint(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks: This is non-seqan, plain c and listed for completeness
..see:http://www.cplusplus.com/reference/clibrary/cctype/isprint/
..see:http://www.cppreference.com/wiki/string/character_classes
..see:Function.isgraph
..include:<cctype>
*/

/*!
 * @fn isgraph
 * @headerfile <cctype>
 * @brief Check if character is printable and not white space
 * 
 * @signature bool isgraph(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is a printable character.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, plain c and listed for completeness.
 * 
 * @see http://www.cplusplus.com/reference/clibrary/cctype/isgraph/
 * @see http://www.cppreference.com/wiki/string/character_classes
 * @see isprint
 * @see isspace
 * @see FileFormatTokenization#readGraphs
 * @see FileFormatTokenization#skipUntilGraph
 */

/**
.Function.isgraph:
..summary:Check if character is printable and not white space
..cat:Input/Output
..signature:isgraph(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks: This is non-seqan, plain c and listed for completeness
..see:http://www.cplusplus.com/reference/clibrary/cctype/isgraph/
..see:http://www.cppreference.com/wiki/string/character_classes
..see:Function.isprint
..see:Function.isspace
..include:<cctype>
*/

/*!
 * @fn isspace
 * @headerfile <cctype>
 * @brief Check if character is a white-space character
 * 
 * @signature bool isspace(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is a printable character.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, plain c and listed for completeness.
 * 
 * NOTE: White-Space contains more than space and tab characters, if you want to check for that, use Function.isblank
 * instead!
 * 
 * @see http://www.cplusplus.com/reference/clibrary/cctype/isspace/
 * @see http://www.cppreference.com/wiki/string/character_classes
 * @see isblank
 * @see FileFormatTokenization#readUntilWhitespace
 * @see FileFormatTokenization#skipUntilWhitespace
 * @see isgraph
 * @see FileFormatTokenization#readDna5IgnoringWhitespaces
 * @see FileFormatTokenization#skipWhitespaces
 */

/**
.Function.isspace:
..summary:Check if character is a white-space character
..cat:Input/Output
..signature:isspace(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks: This is non-seqan, plain c and listed for completeness
..remarks: NOTE: White-Space contains more than space and tab characters, if you want to check for that, use Function.isblank instead!
..see:http://www.cplusplus.com/reference/clibrary/cctype/isspace/
..see:http://www.cppreference.com/wiki/string/character_classes
..see:Function.isblank
..include:<cctype>
*/

/*!
 * @fn isblank
 * @headerfile seqan/stream.h
 * @brief Check if character is either ' ' or '\t'.
 * 
 * @signature bool isblank(c);
 * 
 * @param c The character to be checked. Type: <tt>int</tt>
 * 
 * @return bool <tt>true</tt> (non-zero) if <tt>c</tt> is a blank character.
 * 
 * @section Remarks
 * 
 * This is non-SeqAn, specified in POSIX and listed for completeness.
 * 
 * For visual studio, we define it ourselves.
 * 
 * @see http://www.opengroup.org/onlinepubs/009695399/functions/isblank.html
 * @see http://www.cppreference.com/wiki/string/character_classes
 * @see isspace
 * @see FileFormatTokenization#readLineStripTrailingBlanks
 * @see FileFormatTokenization#skipUntilBlank
 * @see FileFormatTokenization#skipBlanks
 * @see FileFormatTokenization#readUntilBlank
 */

/**
.Function.isblank:
..summary:Check if character is either ' ' or '\t'
..cat:Input/Output
..signature:isblank(int c)
..param.c:the character to be checked
...type:nolink:$int$
..returns:$true$(non-zero) if the check is successful and $false$ (zero) otherwise
..remarks:This is non-seqan, specified in POSIX and listed for completeness.
..remarks:For visual studio, we define it ourselves.
..see:http://pubs.opengroup.org/onlinepubs/009695399/functions/isblank.html
..see:http://www.cppreference.com/wiki/string/character_classes
..see:Function.isspace
..include:seqan/stream.h
*/

#ifdef PLATFORM_WINDOWS_VS
// TODO(holtgrew): Adding basic/POSIX functions for OS that do not support them theirselves should go into platform.
inline
int isblank(int c)
{
    return (c == ' ') || (c == '\t');
}
#endif  // #ifdef PLATFORM_WINDOWS_VS
