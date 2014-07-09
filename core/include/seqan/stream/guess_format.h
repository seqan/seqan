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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// TODO(esiragusa): tests

#ifndef SEQAN_STREAM_GUESS_FORMAT_
#define SEQAN_STREAM_GUESS_FORMAT_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================
// TODO(esiragusa): move all tags into a common header.

// --------------------------------------------------------------------------
// Tag Fasta
// --------------------------------------------------------------------------

/**
.Tag.File Format.tag.Fasta:
    FASTA file format for sequences.
..include:seqan/file.h
*/
struct TagFasta_;
typedef Tag<TagFasta_> Fasta;

// --------------------------------------------------------------------------
// Tag Fastq
// --------------------------------------------------------------------------

/**
.Tag.File Format.tag.Fastq:
    FASTQ file format for sequences.
..include:seqan/file.h
*/
struct TagFastq_;
typedef Tag<TagFastq_> Fastq;

// --------------------------------------------------------------------------
// Tag QSeq
// --------------------------------------------------------------------------

/**
.Tag.File Format.tag.QSeq:
	QSeq format, used for most of the Illumina read files.
..include:seqan/file.h
*/
struct QSeq_;
typedef Tag<QSeq_> QSeq;

// --------------------------------------------------------------------------
// Tag Raw
// --------------------------------------------------------------------------

struct Raw_;
typedef Tag<Raw_> Raw;

// --------------------------------------------------------------------------
// Tag FileFormatExtensions
// --------------------------------------------------------------------------

template <typename TFormat, typename T = void>
struct FileFormatExtensions;


template <typename T>
struct FileFormatExtensions<Fasta, T>
{
    static char const * VALUE[6];
};

template <typename T>
char const * FileFormatExtensions<Fasta, T>::VALUE[6] =
{
    ".fa",      // default output extension
    ".fasta",
    ".faa",     // FASTA Amino Acid file
    ".ffn",     // FASTA nucleotide coding regions file
    ".fna",     // FASTA Nucleic Acid file
    ".frn"
};

template <typename T>
struct FileFormatExtensions<Fastq, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileFormatExtensions<Fastq, T>::VALUE[2] =
{
    ".fq",      // default output extension
    ".fastq"
};

template <typename T>
struct FileFormatExtensions<QSeq, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileFormatExtensions<QSeq, T>::VALUE[2] =
{
    ".txt",     // default output extension
    ".seq"
};

template <typename T>
struct FileFormatExtensions<Raw, T>
{
    static char const * VALUE[1];	// default is one extension
};

template <typename T>
char const * FileFormatExtensions<Raw, T>::VALUE[1] =
{
    ".qseq"     // default output extension
};

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
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function guessFormat()
// --------------------------------------------------------------------------
// Base case: we get here if the file format could not be determined.

template <typename TFileSeq>
inline bool guessFormat(TFileSeq &, TagSelector<> &)
{
    return false;
}

// --------------------------------------------------------------------------
// Function guessFormat(TagSelector)
// --------------------------------------------------------------------------

template <typename TFileSeq, typename TTagList>
inline bool guessFormat(TFileSeq &seq, TagSelector<TTagList> &format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == -1 || value(format) == LENGTH<TTagList>::VALUE - 1)
    {
        // if tagId is set to -1 (auto-detect) or the current format (TFormatTag) then test for TFormatTag format
        if (guessFormat(seq, TFormatTag()))
        {
            value(format) = LENGTH<TTagList>::VALUE - 1;
            return true;
        }
    }
    return guessFormat(seq, static_cast<typename TagSelector<TTagList>::Base &>(format));
}

// --------------------------------------------------------------------------
// Function guessFormatFromFilename()
// --------------------------------------------------------------------------
// Base case: we get here if the file format could not be determined.

template <typename TFilename>
inline bool guessFormatFromFilename(TFilename const &, TagSelector<>)
{
    return false;
}

// --------------------------------------------------------------------------
// Function guessFormatFromFilename(TagSelector)
// --------------------------------------------------------------------------

template <typename TFilename, typename TTagList>
inline bool guessFormatFromFilename(TFilename const &fname, TagSelector<TTagList> &format)
{
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == -1 || value(format) == LENGTH<TTagList>::VALUE - 1)
    {
        // if tagId is set to -1 (auto-detect) or the current format (TFormatTag) then test for TFormatTag format
        if (guessFormatFromFilename(fname, TFormatTag()))
        {
            value(format) = LENGTH<TTagList>::VALUE - 1;
            return true;
        }
    }
    return guessFormatFromFilename(fname, static_cast<typename TagSelector<TTagList>::Base &>(format));
}

// --------------------------------------------------------------------------
// Function guessFormatFromFilename(Tag)
// --------------------------------------------------------------------------

/*!
 * @fn guessFormatFromFilename
 * @headerfile <seqan/file.h>
 * @brief Guesses a file format from a sequence file name.
 *
 * @signature bool guessFormatFromFilename(fileName, formatTag);
 *
 * @param[in] fileName  A filename of a sequence file.
 * @param[in] formatTag A file format tag.
 *
 * @return bool <tt>true</tt> if the format represented by <tt>formatTag</tt> was recognized in <tt>fileName</tt>.
 */

/**
.Function.guessFormatFromFilename:
..summary:Guesses a file format from a sequence file name.
..cat:Input/Output
..signature:guessFormatFromFilename(fileName, formatTag)
..param.fileName:A filename of a sequence file.
...see:Class.String
..param.formatTag:A file format tag.
...type:Tag.File Format
...type:Class.AutoSeqFormat
..returns:$true$ if the format represented by $formatTag$ was recognized in $fileName$.
..see:Function.guessFormatFromFilename
..include:seqan/file.h
*/

template <typename TFilename, typename TFormat_>
inline bool guessFormatFromFilename(TFilename const & fileName, Tag<TFormat_> /*formatTag*/)
{
    typedef typename Value<TFilename>::Type                                     TValue;
    typedef ModifiedString<TFilename const, ModView<FunctorLowcase<TValue> > >	TLowcase;
    typedef Tag<TFormat_>                                                       TFormat;
    
    TLowcase lowcaseFileName(fileName);
    for (unsigned i = 0; i < sizeof(FileFormatExtensions<TFormat>::VALUE) / sizeof(char*); ++i)
        if (endsWith(lowcaseFileName, FileFormatExtensions<TFormat>::VALUE[i]))
            return true;

    return false;
}

// --------------------------------------------------------------------------
// Function _getFileFormatExtensions()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TFormat_>
inline void _getFileFormatExtensions(TStringSet &stringSet, Tag<TFormat_> /*formatTag*/)
{
    typedef Tag<TFormat_> TFormat;
    for (unsigned i = 0; i < sizeof(FileFormatExtensions<TFormat>::VALUE) / sizeof(char*); ++i)
        appendValue(stringSet, FileFormatExtensions<TFormat>::VALUE[i]);
}

template <typename TStringSet, typename TTag>
inline void _getFileFormatExtensions(TStringSet &stringSet, TagList<TTag, void> const /*formatTag*/)
{
    _getFileFormatExtensions(stringSet, TTag());
}

template <typename TStringSet, typename TTag, typename TSubList>
inline void _getFileFormatExtensions(TStringSet &stringSet, TagList<TTag, TSubList> const /*formatTag*/)
{
    _getFileFormatExtensions(stringSet, TTag());
    _getFileFormatExtensions(stringSet, TSubList());
}

template <typename TStringSet, typename TTagList>
inline void _getFileFormatExtensions(TStringSet &stringSet, TagSelector<TTagList> const /*formatTag*/)
{
    _getFileFormatExtensions(stringSet, TTagList());
}

} // namespace seqan

#endif  // #ifndef SEQAN_STREAM_GUESS_FORMAT_
