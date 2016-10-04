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
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================
// Class for reading/writing files in Bpseq format.
// ==========================================================================
// TODO(weese:) add Bcf I/O and integrate it

#ifndef SEQAN_BPSEQ_IO_BPSEQ_FILE_H_
#define SEQAN_BPSEQ_IO_BPSEQ_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Bpseq_;
typedef Tag<Bpseq_> Bpseq;

//struct Bcf_;
//typedef Tag<Bcf_> Bcf;

// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef BpseqFileIn
// ----------------------------------------------------------------------------

/*!
 * @class BpseqFileIn
 * @signature typedef FormattedFile<Bpseq, Input> BpseqFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/bpseq_io.h>
 * @brief Class for reading BPSEQ files.
 *
 * @see BpseqHeader
 * @see BpseqRecord
 */

typedef FormattedFile<Bpseq, Input>   BpseqFileIn;

// ----------------------------------------------------------------------------
// Typedef BpseqFileOut
// ----------------------------------------------------------------------------

/*!
 * @class BpseqFileOut
 * @signature typedef FormattedFile<Bpseq, Output> BpseqFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/vcf_io.h>
 * @brief Class for writing BPSEQ files.
 *
 * @see BpseqHeader
 * @see BpseqRecord
 */

typedef FormattedFile<Bpseq, Output>  BpseqFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Bpseq, T>
{
    static unsigned char const VALUE[18];
};

template <typename T>
unsigned char const MagicHeader<Bpseq, T>::VALUE[18] =
{
    '#', '#', 'f', 'i', 'l', 'e', 'f', 'o', 'r', 'm', 'a', 't', '=', 'B', 'P', 'S', 'E', 'Q'  // BPSEQ's magic header
};

//template <typename T>
//struct MagicHeader<Bcf, T>
//{
//    static unsigned char const VALUE[5];
//};
//
//template <typename T>
//unsigned char const MagicHeader<Bcf, T>::VALUE[5] = { 'B', 'C', 'F', '\2', '\1' };  // BCF2's magic header

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Bpseq, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Bpseq, T>::VALUE[1] =
{
    ".bpseq"     // default output extension
};

//template <typename T>
//struct FileExtensions<Bcf, T>
//{
//    static char const * VALUE[1];    // default is one extension
//};
//
//template <typename T>
//char const * FileExtensions<Bcf, T>::VALUE[1] =
//{
//    ".bcf"     // default output extension
//};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Bpseq, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    //typedef BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
    typedef RnaIOContext                                             Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------
// Riferirsi al file top chiamato sequence_file.h nel modulo seq_io per generare il sistema chiamato rna_io
template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Bpseq, TDirection, TSpec> >
{
// TODO(weese:) Enable this, as soon as someone implements BCF

//#if SEQAN_HAS_ZLIB
//    typedef TagSelector<
//                TagList<Bcf, /// TODO aggiungere il tag selector con gli altri file format
//                TagList<Bpseq
//                > >
//            > Type;
//#else
    typedef Bpseq Type;
//#endif
};

// --------------------------------------------------------------------------
// Function _mapBamFormatToCompressionFormat()
// --------------------------------------------------------------------------

//inline BgzfFile
//_mapFileFormatToCompressionFormat(Bcf)
//{
//    return BgzfFile();
//}

// ----------------------------------------------------------------------------
// Function readRecord(); BpseqRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RnaRecord & record, FormattedFile<Bpseq, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); BpseqRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Bpseq, Output, TSpec> & file, RnaRecord & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

}  // namespace seqan

#endif // SEQAN_BPSEQ_IO_BPSEQ_FILE_H_
