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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Smart file for reading/writing files in Fasta or Fastq format.
// ==========================================================================

#ifndef SEQAN_SEQ_IO_SEQUENCE_FILE_H_
#define SEQAN_SEQ_IO_SEQUENCE_FILE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

typedef SmartFile<Fastq, Input>     SeqFileIn;
typedef SmartFile<Fastq, Output>    SeqFileOut;

// --------------------------------------------------------------------------
// Tag AutoSeqFormat
// --------------------------------------------------------------------------
// if TagSelector is set to -1, the file format is auto-detected

/*!
 * @class AutoSeqFormat
 * @extends TagSelector
 * @headerfile <seqan/file.h>
 * @brief Auto-detects and stores a file format.
 *
 * @signature typedef TagList<Fastq, TagList<Fasta, TagList<Raw> > > SeqFormats;
 * @signature typedef TagSelector<SeqFormat> AutoSeqFormat;
 */

/**
.Class.AutoSeqFormat
..summary:Auto-detects and stores a file format.
..cat:Input/Output
..general:Class.TagSelector
..signature:AutoSeqFormat
..remarks:Currently, it is defined as $TagSelector<SeqFormats>$, with:
...code:
	typedef
		TagList<Fastq,
		TagList<Fasta,
		TagList<QSeq,
		TagList<Raw> > > > 						SeqFormats;
..include:seqan/file.h
*/

typedef
    TagList<Fastq,
    TagList<Fasta,
//    TagList<QSeq,   // doesn't work as it uses STL strings and parsers
    TagList<Raw
    > > > //>
    SeqFormats;

typedef TagSelector<SeqFormats> AutoSeqFormat;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormats<SmartFile<Fastq, TDirection, TSpec> >
{
    typedef AutoSeqFormat Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function read(); Without qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename SmartFile<Fastq, TDirection, TSpec>::TStream> >, void)
readRecord(TIdString & meta, TSeqString & seq, SmartFile<Fastq, TDirection, TSpec> & file)
{
    readRecord(meta, seq, file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function read(); With qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename SmartFile<Fastq, TDirection, TSpec>::TStream> >, void)
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, SmartFile<Fastq, TDirection, TSpec> & file)
{
    readRecord(meta, seq, qual, file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); Without qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename SmartFile<Fastq, TDirection, TSpec>::TStream> >, void)
writeRecord(SmartFile<Fastq, TDirection, TSpec> & file, TIdString const & meta, TSeqString const & seq)
{
    writeRecord(file.iter, meta, seq, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); With qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename SmartFile<Fastq, TDirection, TSpec>::TStream> >, void)
writeRecord(SmartFile<Fastq, TDirection, TSpec> & file, TIdString const & meta, TSeqString const & seq, TQualString const & qual)
{
    writeRecord(file.iter, meta, seq, qual, file.format);
}

}  // namespace seqan

#endif // SEQAN_SEQ_IO_SEQUENCE_FILE_H_