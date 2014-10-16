// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Smart file for reading/writing files in UCSC format.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_UCSC_IO_UCSC_FILE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_UCSC_IO_UCSC_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

typedef Tag<Ucsc_<> > Ucsc;


// ----------------------------------------------------------------------------
// Typedef UcscFileIn
// ----------------------------------------------------------------------------

/*!
 * @class UcscFileIn
 * @extends SmartFile
 * @headerfile <seqan/ucsc_io.h>
 * @brief @link SmartFile @endlink for reading UCSC <tt>knownGenes.txt</tt> and <tt>knownIsoforms.txt</tt> files.
 *
 * @signature typedef SmartFile<Ucsc, Input> UcscFileIn;
 */
typedef SmartFile<Ucsc, Input>   UcscFileIn;

// ----------------------------------------------------------------------------
// Typedef UcscFileOut
// ----------------------------------------------------------------------------

/*!
 * @class UcscFileInOut
 * @extends SmartFile
 * @headerfile <seqan/ucsc_io.h>
 * @brief @link SmartFile @endlink for reading UCSC <tt>knownGenes.txt</tt> and <tt>knownIsoforms.txt</tt> files.
 *
 * @signature typedef SmartFile<Ucsc, Output> UcscFileOut;
 */
typedef SmartFile<Ucsc, Output> UcscFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct SmartFileContext<SmartFile<Ucsc, TDirection, TSpec>, TStorageSpec>
{
    typedef UcscIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<SmartFile<Ucsc, TDirection, TSpec> >
{
    typedef TagSelector<
                TagList<UcscKnownGene,
                TagList<UcscKnownIsoforms
                > >
            > Type;
};

// ----------------------------------------------------------------------------
// Function readRecord(); UcscRecord
// ----------------------------------------------------------------------------

/*!
 * @fn UcscFileIn#readRecord
 * @brief Reading records from a UcscFileIn.
 *
 * @signature void readRecord(record, file);
 *
 * @param[out] record The resulting @link UcscRecord @endlink.
 * @param[out] file   The UcscFileIn to read from.
 *
 * @throw IOError in case of problems.
 */

// support for dynamically chosen file formats
template <typename TForwardIter>
inline void
readRecord(UcscRecord & /* record */,
           UcscIOContext & /* context */,
           TForwardIter & /* iter */,
           TagSelector<> const & /* format */)
{
    SEQAN_FAIL("UcscFileIn: File format not specified.");
}

template <typename TForwardIter, typename TTagList>
inline void
readRecord(UcscRecord & record,
           UcscIOContext & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(record, context, iter, TFormat());
    else
        readRecord(record, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
void readRecord(UcscRecord & record, SmartFile<Ucsc, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); UcscRecord
// ----------------------------------------------------------------------------

/*!
 * @fn UcscFileIn#writeRecord
 * @brief Writing records to a UcscFileIn.
 *
 * @signature void readRecord(record, file);
 *
 * @param[out] file   The UcscFileIn to write to.
 * @param[out] record The @link UcscRecord @endlink to write out.
 *
 * @throw IOError in case of problems.
 */

// support for dynamically chosen file formats
template <typename TTarget>
inline void
writeRecord(TTarget & /* target */,
            UcscRecord const & /* record */,
            TagSelector<> const & /* format */)
{
    SEQAN_FAIL("UcscFileOut: File format not specified.");
}

template <typename TTarget, typename TTagList>
inline void
writeRecord(TTarget & target,
            UcscRecord const & record,
            TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        writeRecord(target, record, TFormat());
    else
        writeRecord(target, record, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
writeRecord(SmartFile<Ucsc, Output, TSpec> & file, UcscRecord const & record)
{
    writeRecord(file.iter, record, file.format);
}

}  // namespace seqan

#endif // SEQAN_CORE_INCLUDE_SEQAN_UCSC_IO_UCSC_FILE_H_
