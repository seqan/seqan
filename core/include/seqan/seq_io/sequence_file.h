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
// ==========================================================================
// TODO(esiragusa): move the content of file/file_format_mmap.h (e.g. AutoSeqFormat) inside this file

#ifndef SEQAN_SEQ_IO_SEQUENCE_FILE_H_
#define SEQAN_SEQ_IO_SEQUENCE_FILE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SequenceFile
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec = void>
struct SequenceFile
{
//    typedef VirtualStream<char, TDirection>             TStream;
//    typedef typename Iterator<TStream, Standard>::Type  TIter;

    typedef std::ofstream                               TStream;
    typedef Iter<TStream, StreamIterator<Output> >      TIter;

    AutoSeqFormat   format;
    TStream         stream;
    TIter           iter;

    SequenceFile() {}

    SequenceFile(const char *fileName, int openMode = DefaultOpenMode<SequenceFile>::VALUE)
    {
        if (!open(*this, fileName, openMode))
            throw std::runtime_error(std::string("Could not open file ") + fileName);
    }

    ~SequenceFile()
    {
        close(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultOpenMode
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct DefaultOpenMode<SequenceFile<TDirection, TSpec> > :
    DefaultOpenMode<typename SequenceFile<TDirection, TSpec>::TStream> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setFormat()
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TFormat>
inline void setFormat(SequenceFile<TDirection, TSpec> & file, TFormat)
{
    assign(file.format, TFormat());
}

// ----------------------------------------------------------------------------
// Function open(fileName)
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
inline bool open(SequenceFile<TDirection, TSpec> & file,
                 const char *fileName,
                 int openMode = DefaultOpenMode<SequenceFile<TDirection, TSpec> >::VALUE)
{
    typedef typename SequenceFile<TDirection, TSpec>::TIter TIter;

    if (!guessFormatFromFilename(fileName, file.format))
        return false;

    if (!open(file.stream, fileName, openMode))
        return false;

    file.iter = TIter(file.stream);

    return true;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
inline bool close(SequenceFile<TDirection, TSpec> & file)
{
    return close(file.stream);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename SequenceFile<TDirection, TSpec>::TStream> >, bool)
atEnd(SequenceFile<TDirection, TSpec> const & file)
{
    return atEnd(file.iter);
}

// ----------------------------------------------------------------------------
// Function read(); Without qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename SequenceFile<TDirection, TSpec>::TStream> >, void)
read(SequenceFile<TDirection, TSpec> & file, TIdString & meta, TSeqString & seq)
{
    switch (value(file.format))
    {
        case Find<AutoSeqFormat, Fasta>::VALUE:
            readRecord(meta, seq, file.iter, Fasta());
            break;

        case Find<AutoSeqFormat, Fastq>::VALUE:
            readRecord(meta, seq, file.iter, Fastq());
            break;
    }
}

// ----------------------------------------------------------------------------
// Function read(); With qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename SequenceFile<TDirection, TSpec>::TStream> >, void)
read(SequenceFile<TDirection, TSpec> & file, TIdString & meta, TSeqString & seq, TQualString & qual)
{
    readRecord(meta, seq, qual, file.iter, Fastq());
}

// ----------------------------------------------------------------------------
// Function write(); Without qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename SequenceFile<TDirection, TSpec>::TStream> >, void)
write(SequenceFile<TDirection, TSpec> & file, TIdString const & meta, TSeqString const & seq)
{
    switch (value(file.format))
    {
        case Find<AutoSeqFormat, Fasta>::VALUE:
            writeRecord(file.iter, meta, seq, Fasta());
            break;

        case Find<AutoSeqFormat, Fastq>::VALUE:
            writeRecord(file.iter, meta, seq, Fastq());
            break;
    }
}

// ----------------------------------------------------------------------------
// Function write(); With qualities
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline SEQAN_FUNC_ENABLE_IF(Is<OutputStreamConcept<typename SequenceFile<TDirection, TSpec>::TStream> >, void)
write(SequenceFile<TDirection, TSpec> & file, TIdString const & meta, TSeqString const & seq, TQualString const & qual)
{
    writeRecord(meta, seq, qual, file.iter, Fastq());
}

}  // namespace seqan

#endif // SEQAN_SEQ_IO_SEQUENCE_FILE_H_