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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// This file contains routines to read BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_BLAST_BLAST_TABULAR_LOWLEVEL_H_
#define SEQAN_BLAST_BLAST_TABULAR_LOWLEVEL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BlastTabularLL
// ----------------------------------------------------------------------------

/*!
 * @class BlastTabularLL
 * @signature typedef Tag<BlastTabularLL_> BlastTabularLL;
 * @headerfile <seqan/blast.h>
 * @brief Low-Level support for Blast Tabular file formats
 *
 * There are three blast format related tags in SeqAn:
 *
 * <li> @link BlastReport @endlink with the FormattedFile output specialization @link BlastReportFileOut @endlink</li>
 * <li> @link BlastTabular @endlink with the FormattedFile output and input specializations
 * @link BlastTabularFileOut @endlink and @link BlastTabularFileIn @endlink</li>
 * <li> @link BlastTabularLL @endlink which provides light-weight, but very basic tabular IO </li>
 *
 * This is the third tag, it offers <b>low-level</b> support for reading and writing NCBI Blast compatible
 * <b>tabular</b> files, <b>without comment lines</b> -- although files with comment lines can be read if the comment
 * lines are skipped. These are the formats that are available in legacy Blast
 * (<tt>blastall</tt> executable) with the parameters <tt>-m 8</tt> and <tt>-m 9</tt> (with comment lines)
 * and in BLAST+ (<tt>blastx</tt>, <tt>blastn</tt>...) with
 * the parameters <tt>-outfmt 6</tt> and <tt>-outfmt 7</tt> respectively.
 *
 * For most situations @link BlastTabular @endlink is more adequate. Use this tag's interface only for quick parsing
 * of matches in a file, e.g counting and filtering purposes. This interface does not offer a FormattedFile
 * abstraction and no convenience data structures, it does no transformations on the data.
 *
 * The reference Blast implementation used for developing the SeqAn support is NCBI Blast+ 2.2.26 and
 * NCBI Blast 2.2.26 for the legacy support.
 *
 * @section Input example
 *
 * The following example program extracts the list of matching query-subject-pairs from a blast tabular file and prints
 * it to std::out:
 *
 * @include demos/blast/blast_in_lowlevel.cpp
 *
 * The output looks like this:
 *
 * @include demos/blast/blast_in_lowlevel.cpp.stdout_
 *
 */

struct BlastTabularLL_;
typedef Tag<BlastTabularLL_> BlastTabularLL;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function onMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularLL#onMatch
 * @brief Returns whether the iterator is on the beginning of a match line.
 * @signature bool onMatch(stream, blastTabularLL)
 * @headerfile seqan/blast.h
 *
 * @param[in] iter              An input iterator over a stream or any fwd-iterator over a string.
 * @param[in] blastTabularLL    The @link BlastTabularLL @endlink tag.
 *
 * @throw IOError On low-level I/O errors.
 *
 * @return bool true or false
 */

template <typename TFwdIterator>
inline bool
onMatch(TFwdIterator & iter,
        BlastTabularLL const &)
{
    return (value(iter) != '#');
}

// ----------------------------------------------------------------------------
// Function skipUntilMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularLL#skipUntilMatch
 * @brief Skip arbitrary number of comment lines until the beginning of a match is reached.
 * @signature void skipUntilMatch(stream, blastTabularLL);
 * @headerfile seqan/blast.h
 *
 * @param[in,out] stream         An input iterator over a stream or any fwd-iterator over a string.
 * @param[in]     blastTabularLL The @link BlastTabularLL @endlink tag.
 *
 * @section Remarks
 *
 * This is also part of the low-level IO and not required if you use readRecord.
 * Call this function whenever you are not @link BlastTabularLL#onMatch @endlink, but want to be, e.g. to
 * @link BlastTabularLL#readMatch @endlink.
 *
 * Since it is legal for files to end with comment lines, this function does not throw if end-of-file is reached.
 * You need to check that after calling.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator>
inline void
skipUntilMatch(TFwdIterator & iter,
               BlastTabularLL const & /*tag*/)
{
    while ((!atEnd(iter)) && value(iter) == '#') // skip comments
        skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function readMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularLL#readMatch
 * @brief Low-level BlastTabular file reading.
 * @signature void readMatch(stream, blastTabularLL, args ...);
 * @headerfile seqan/blast.h
 *
 * @param[in,out] stream         An input iterator over a stream or any fwd-iterator over a string.
 * @param[in]     blastTabularLL The @link BlastTabularLL @endlink tag.
 * @param[out]    args           Arbitrary typed variables able to hold the fields.
 *
 * @section Remarks
 *
 * Use this signature only if you do not or cannot use @link BlastMatch
 * @endlinkes. You can specify any number of arguments that are expected
 * to be able to hold the values in the columns read, i.e. if you pass a
 * double as argument and the value in the column cannot be successfully cast
 * to double, an exception will be thrown. If you want to be on the safe side,
 * you can pass CharStrings and evaluate them in another way.
 *
 * You may specify less columns than are available in the file, all but the first
 * n will be discarded.
 *
 * No transformations are made on the data, e.g. the positions are still
 * one-indexed and flipped for reverse strand matches.
 *
 * See @link BlastTabularLL @endlink for an example of low-level IO.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

// arbitrary columns
template <typename TTarget>
inline SEQAN_FUNC_ENABLE_IF(IsSequence<TTarget>)
_assignOrCast(TTarget & target, std::string const & source)
{
    assign(target, source);
}

template <typename TTarget>
inline SEQAN_FUNC_ENABLE_IF(Is<NumberConcept<TTarget>>)
_assignOrCast(TTarget & target, std::string const & source)
{
    target = lexicalCast<TTarget>(source);
}

template <typename TFwdIterator,
          typename TArg>
inline void
_readMatchImplBlastTab(TFwdIterator & iter,
                       TArg & arg)
{
    static std::string buffer;
    clear(buffer);

    readUntil(buffer, iter, OrFunctor<IsTab,IsNewline>());
    _assignOrCast(arg, buffer);

    // as this is the last requested field, go to beginning of next line
    skipLine(iter);
}

template <typename TFwdIterator,
          typename TArg,
          typename... TArgs>
inline void
_readMatchImplBlastTab(TFwdIterator & iter,
                       TArg & arg,
                       TArgs & ... args)
{
    static std::string buffer;
    clear(buffer);

    readUntil(buffer, iter, IsTab());
    skipOne(iter, IsTab());
    _assignOrCast(arg, buffer);

    // recurse to next argument
    _readMatchImplBlastTab(iter, args...);
}

// custom arguments
template <typename TFwdIterator,
          typename... TArgs>
inline void
readMatch(TFwdIterator & iter,
           BlastTabularLL const &,
           TArgs & ... args)
{
    // comment lines should have been read or skipped
    if (SEQAN_UNLIKELY(!onMatch(iter, BlastTabularLL())))
        SEQAN_THROW(ParseError("ERROR: Not on beginning of Match (you should have skipped comments)."));

    _readMatchImplBlastTab(iter, args...);
}

// ----------------------------------------------------------------------------
// Function writeMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularLL#writeMatch
 * @headerfile seqan/blast.h
 * @brief Low-level file-writing for blast tabular formats
 * @signature void writeMatch(stream, blastTabularLL, columns...)
 *
 * @section Remarks
 *
 * This is a very leight-weight alternative to @link BlastTabularFileOut#writeRecord @endlink. It doesn't require
 * @link BlastMatch @endlinkes, @link BlastRecord @endlinks or the use of @link FormattedFile @endlink.
 * It supports an arbitrary amount of and arbitrary typed columns to be printed.
 *
 * Use this only if you do not require comment lines and you are prepared to do all transformations on the data
 * yourself, i.e. this function does none of the match adjustments mentioned in
 * @link BlastTabularFileOut#writeRecord @endlink.
 *
 * @param[in,out] stream         The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...).
 * @param[in]     blastTabularLL The @link BlastTabularLL @endlink tag.
 * @param[in]     columns...     Any number of printable parameters.
 *
 * @throw IOError On low-level I/O errors.
 */

template <typename TFwdIterator>
inline void
_writeFields(TFwdIterator & /**/,
             BlastTabularLL const & /*tag*/)
{
}

template <typename TFwdIterator, typename TField, typename... TFields>
inline void
_writeFields(TFwdIterator & stream,
             BlastTabularLL const & /*tag*/,
             TField const & field1,
             TFields const & ... fields)
{
    write(stream, '\t');
    write(stream, field1);
    _writeFields(stream, BlastTabularLL(), fields... );
}

// Function for arbitrary number and typed fields
template <typename TFwdIterator, typename TField, typename... TFields>
inline void
writeMatch(TFwdIterator & stream,
           BlastTabularLL const & /*tag*/,
           TField const & field1,
           TFields const & ... fields)
{
    write(stream, field1);

    _writeFields(stream, BlastTabularLL(), fields...);
    write(stream, '\n');
}

}

#endif // SEQAN_BLAST_BLAST_TABULAR_LOWLEVEL_H_