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
// Authors: Joerg Winkler <j.winkler@fu-berlin.de>
//          Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================
// Class for reading/writing files in Ebpseq format.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_EBPSEQ_FILE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_EBPSEQ_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Ebpseq_;
typedef Tag<Ebpseq_> Ebpseq;

// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef EbpseqFileIn
// ----------------------------------------------------------------------------

/*!
 * @class EbpseqFileIn
 * @signature typedef FormattedFile<Ebpseq, Input> EbpseqFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/bpseq_io.h>
 * @brief Class for reading BPSEQ files.
 *
 * @see EbpseqHeader
 * @see EbpseqRecord
 */

typedef FormattedFile<Ebpseq, Input>   EbpseqFileIn;

// ----------------------------------------------------------------------------
// Typedef EbpseqFileOut
// ----------------------------------------------------------------------------

/*!
 * @class EbpseqFileOut
 * @signature typedef FormattedFile<Ebpseq, Output> EbpseqFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/vcf_io.h>
 * @brief Class for writing BPSEQ files.
 *
 * @see EbpseqHeader
 * @see EbpseqRecord
 */

typedef FormattedFile<Ebpseq, Output>  EbpseqFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Ebpseq, T>
{
    static unsigned char const VALUE[19];
};

template <typename T>
unsigned char const MagicHeader<Ebpseq, T>::VALUE[19] =
{
    '#', '#', 'f', 'i', 'l', 'e', 'f', 'o', 'r', 'm', 'a', 't', '=', 'E', 'B', 'P', 'S', 'E', 'Q'
};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Ebpseq, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Ebpseq, T>::VALUE[1] =
{
    ".ebpseq"     // default output extension
};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Ebpseq, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    //typedef EbpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
    typedef RnaIOContext                                             Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------
template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Ebpseq, TDirection, TSpec> >
{
    typedef Ebpseq Type;
};

// ----------------------------------------------------------------------------
// Function readHeader(); EbpseqHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(RnaHeader & header, FormattedFile<Ebpseq, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(); EbpseqRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RnaRecord & record, FormattedFile<Ebpseq, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeHeader(); EbpseqHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<Ebpseq, Output, TSpec> & file, RnaHeader & header)
{
    writeHeader(file.iter, header, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); EbpseqRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Ebpseq, Output, TSpec> & file, RnaRecord & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

}  // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_EBPSEQ_FILE_H_
