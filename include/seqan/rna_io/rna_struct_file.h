// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// Class for reading/writing RNA structure files.
// ==========================================================================

#ifndef SEQAN_RNA_IO_RNA_STRUCT_FILE_H_
#define SEQAN_RNA_IO_RNA_STRUCT_FILE_H_

namespace seqan {

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag RnaStruct
// --------------------------------------------------------------------------

/*!
 * @tag FileFormats#RnaStruct
 * @headerfile <seqan/rna_io.h>
 * @brief General RNA structure file type. The concrete type is defined by the file name suffix.
 * @signature typedef Tag<RnaStruct_> RnaStruct;
 * @see FileFormats#Connect
 * @see FileFormats#Stockholm
 * @see FileFormats#DotBracket
 * @see FileFormats#Vienna
 * @see FileFormats#Bpseq
 * @see FileFormats#Ebpseq
 *
 * This tag is used as TSpec template parameter in @link FormattedFile @endlink.
 */
struct RnaStruct_;
typedef Tag<RnaStruct_> RnaStruct;

typedef TagList<Connect,
        TagList<Stockholm,
        TagList<DotBracket,
        TagList<Vienna,
        TagList<Ebpseq,
        TagList<Bpseq
        > > > > > > RnaStructFormats;

// --------------------------------------------------------------------------
// Class RnaStructContents
// --------------------------------------------------------------------------

/*!
 * @class RnaStructContents
 * @headerfile <seqan/rna_io.h>
 * @brief Contains all the contents of an RNA structure file.
 *
 * @signature struct RnaStructContents;
 * @see RnaRecord
 * @see RnaHeader
 *
 * Container for the header and all records of an RNA structure file.
 * If no header was present in the file, a pseudo header with the required information is stored instead.
 */
struct RnaStructContents
{
    /*!
     * @var std::vector<RnaRecord> RnaStructContents::records
     * @brief All records of an RNA structure file.
     * @see RnaRecord
     */
    std::vector<RnaRecord> records;

    /*!
     * @var RnaHeader RnaStructContents::header
     * @brief The (pseudo) header of an RNA structure file.
     * @see RnaHeader
     */
    RnaHeader header;
};

// ============================================================================
// Typedefs
// ============================================================================

/*!
 * @class RnaStructFileIn
 * @signature typedef FormattedFile<RnaStruct, Input> RnaStructFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/rna_io.h>
 * @brief Class for reading DBN, DBV, CT, STH, BPSEQ, and EBPSEQ files containing RNA structures.
 */
typedef FormattedFile<RnaStruct, Input> RnaStructFileIn;

/*!
 * @class RnaStructFileOut
 * @signature typedef FormattedFile<RnaStruct, Output> RnaStructFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/rna_io.h>
 * @brief Class for writing DBN, DBV, CT, STH, BPSEQ, and EBPSEQ files containing RNA structures.
 */
typedef FormattedFile<RnaStruct, Output> RnaStructFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<RnaStruct, TDirection, TSpec>, TStorageSpec>
{
    typedef RnaIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormat
// ----------------------------------------------------------------------------

template <typename TSpec>
struct FileFormat<FormattedFile<RnaStruct, Input, TSpec> >
{
    typedef TagSelector<RnaStructFormats> Type;
};

template <typename TSpec>
struct FileFormat<FormattedFile<RnaStruct, Output, TSpec> >
{
    typedef TagSelector<RnaStructFormats> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector)
// ----------------------------------------------------------------------------

/*!
 * @fn RnaStructFileIn#readRecord
 * @brief Read one @link FormattedFileRecordConcept @endlink from an @link RnaStructFileIn @endlink object.
 * @signature void readRecord(record, fileIn);
 * @param[in,out] record   The @link RnaRecord @endlink object where to write the information into.
 * @param[in] fileIn       The @link RnaStructFileIn @endlink object to read from.
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator>
inline void
readRecord(RnaRecord &, RnaIOContext &, TFwdIterator &, TagSelector<> const &)
{
    SEQAN_FAIL("RnaStructFileIn: File format not specified.");
}

template <typename TFwdIterator, typename TTagList>
inline void
readRecord(RnaRecord & record, RnaIOContext & context, TFwdIterator & iter, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(record, context, iter, TFormat());
    else
        readRecord(record, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
readRecord(RnaRecord & record, FormattedFile<RnaStruct, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(TagSelector)
// ----------------------------------------------------------------------------

/*!
 * @fn RnaStructFileOut#writeRecord
 * @brief Write one @link FormattedFileRecordConcept @endlink into an @link RnaStructFileOut @endlink object.
 * @signature void writeRecord(fileOut, record);
 * @param[in,out] fileOut   The @link RnaStructFileOut @endlink object to write into.
 * @param[in] record        The @link RnaRecord @endlink object where to read from.
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator>
inline void
writeRecord(TFwdIterator &, RnaRecord const &, RnaIOContext &, TagSelector<> const &)
{
    SEQAN_FAIL("RnaStructFileOut: File format not specified.");
}

template <typename TFwdIterator, typename TTagList>
inline void
writeRecord(TFwdIterator & iter, RnaRecord const & record, RnaIOContext & context, TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        writeRecord(iter, record, context, TFormat());
    else
        writeRecord(iter, record, context, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

template <typename TSpec>
inline void
writeRecord(FormattedFile<RnaStruct, Output, TSpec> & file, RnaRecord const & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function readHeader(TagSelector)
// ----------------------------------------------------------------------------

/*!
 * @fn RnaStructFileIn#readHeader
 * @brief Read one @link FormattedFileHeaderConcept @endlink from an @link RnaStructFileIn @endlink object.
 * @signature void readHeader(header, fileIn);
 * @param[in,out] header  The @link RnaHeader @endlink object where to write the information into.
 * @param[in] fileIn      The @link RnaStructFileIn @endlink object to read from.
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */
template <typename TSpec>
inline void
readHeader(RnaHeader & header, FormattedFile<RnaStruct, Input, TSpec> & file)
{
    if (isEqual(file.format, Ebpseq()))  // all other files contain no header
        readHeader(header, context(file), file.iter, Ebpseq());
}

// ----------------------------------------------------------------------------
// Function writeHeader(TagSelector)
// ----------------------------------------------------------------------------

/*!
 * @fn RnaStructFileOut#writeHeader
 * @brief Write one @link FormattedFileHeaderConcept @endlink into an @link RnaStructFileOut @endlink object.
 * @signature void writeHeader(fileOut, header);
 * @param[in,out] fileOut   The @link RnaStructFileOut @endlink object to write into.
 * @param[in] record        The @link RnaHeader @endlink object where to read from.
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */
template <typename TSpec>
inline void
writeHeader(FormattedFile<RnaStruct, Output, TSpec> & file, RnaHeader const & header)
{
    if (isEqual(file.format, Ebpseq()))  // all other files contain no header
        writeHeader(file.iter, header, context(file), Ebpseq());
}

// ----------------------------------------------------------------------------
// Function readRecords
// ----------------------------------------------------------------------------

/*!
 * @fn RnaStructFileIn#readRecords
 * @brief Read @link RnaStructContents @endlink from a @link RnaStructFileIn @endlink object.
 * @signature void readRecords(contents, fileIn, maxRecords);
 * @param[in,out] contents   The @link RnaStructContents @endlink object where to write the information into.
 * @param[in] fileIn         The @link RnaStructFileIn @endlink object to read from.
 * @param[in] maxRecords     The maximum number of records to read.
 * @see RnaStructFileIn#readRecord
 * @see RnaStructFileIn#readHeader
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */
template <typename TSpec, typename TSize>
inline void readRecords(RnaStructContents & contents, FormattedFile<RnaStruct, Input, TSpec> & file, TSize maxRecords)
{
    readHeader(contents.header, file);
    RnaRecord record;
    TSize numRecords = 0;
    while (!atEnd(file.iter) && ++numRecords <= maxRecords)
    {
        if (numRecords >= UINT32_MAX)
            SEQAN_THROW(ParseError("ERROR: Capacity exceeded. File contains too many sequences: Set the maxRecords "
                                           "parameter to at most " + std::to_string(UINT32_MAX - 1u) + "."));

        readRecord(record, file);
        append(contents.records, record);
    }
    // for files without header: create pseudo header
    if (empty(contents.header.seqLabels))
        createPseudoHeader(contents.header, contents.records); // in ebpseq_read_write.h
}

// ----------------------------------------------------------------------------
// Function writeRecords
// ----------------------------------------------------------------------------

/*!
 * @fn RnaStructFileOut#writeRecords
 * @brief Write @link RnaStructContents @endlink into a @link RnaStructFileOut @endlink object.
 * @signature void writeRecords(fileOut, contents);
 * @param[in,out] fileOut   The @link RnaStructFileOut @endlink object to write into.
 * @param[in] contents      The @link RnaStructContents @endlink object where to read from.
 * @see RnaStructFileOut#writeRecord
 * @see RnaStructFileOut#writeHeader
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */
template <typename TSpec>
inline void writeRecords(FormattedFile<RnaStruct, Output, TSpec> & file, RnaStructContents const & contents)
{
    writeHeader(file, contents.header);
    for (RnaRecord const & record : contents.records)
        writeRecord(file, record);
}

}  // namespace seqan

#endif // SEQAN_RNA_IO_RNA_STRUCT_FILE_H_
