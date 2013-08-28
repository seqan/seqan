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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Simple-to-use I/O for SAM and BAM files.
//
// The actual implementation is done using the class hierarchies rooted in
// XamReader and XamWriter.  We use virtual functions and native C++
// inheritance here because the implementations requires some kind of dynamic
// lookup, and using the built-in inheritance model is probably fast.  Also,
// this API is on a very high level/layer and thus some loss in performance is
// OK.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_STREAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_STREAM_H_

#include <fstream>
#include <memory>

#include <seqan/store.h>

// TODO(holtgrew): Replace std::fstream by MMap String?
// TODO(holtgrew): Implement BAI support?

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class BamStream
 * @headerfile <seqan/bam_io.h>
 * @brief Class that provides an easy to use interface for reading and writing SAM and BAM files.
 *
 * @signature class BamStream;
 *
 * @section Example
 *
 * Read and write SAM or BAM files.
 *
 * @include demos/bam_io/bam_stream.cpp
 *
 * The output is as follows:
 *
 * @include demos/bam_io/bam_stream.cpp.stdout
 */

/*!
 * @fn BamStream::BamStream
 * @brief Constructor
 *
 * @signature BamStream::BamStream([fileName[, mode[, format]]]);
 *
 * @param[in] fileName The path to the SAM or BAM file to load, <tt>char const *</tt>.
 * @param[in] mode     The open mode, of type @link BamStream::OperationMode @endlink, defaults to <tt>READ</tt>.
 * @param[in] format   The format, of type @link BamStream::Format @endlink, defaults to <tt>AUTO</tt>.
 */

/*!
 * @var THeader BamStream::header
 * @brief The @link BamHeader @endlink of the @link BamStream @endlink object.
 *
 * SAM and BAM files have a header.  When writing SAM or BAM files, you have to fill this member before writing @link
 * BamAlignmentRecord @endlinks.  Upon writing the first record, the header will be written out.
 *
 * When reading BAM files, the header will be read upon opening the file. When reading SAM files, any header will be
 * read upon opening the file.
 *
 * Note that there is a special case when reading SAM records: If there is no header, or records refer to reference
 * sequences that are previously unknown when reading SAM then a new entry is added to @link BamHeader::sequenceInfos
 * @endlink.
 */

/*!
 * @var TBamIOContext BamStream::bamIOContext
 *
 * @brief The @link BamIOContext @endlink object to use for reading and writing @link BamAlignmentRecord @endlinks.
 *
 * When reading, the <tt>bamIOContext</tt> will be updated automatically.  When reading SAM, new reference sequences can
 * be introduced "on the fly" when a new sequence appears.  When writing, the <tt>bamIOContext</tt> is automatically
 * filled/reset when the first record is written.
 */

/**
.Class.BamStream
..cat:BAM I/O
..summary:Class that provides an easy to use interface for reading and writing SAM and BAM files.
..signature:BamStream
..example:Read and write SAM or BAM files.
..example.file:demos/bam_io/bam_stream.cpp
..example.text:The output is as follows:
..example.output:
@HD VN:1.3  SO:coordinate
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
r004    0   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
r001    83  ref 37  30  9M  =   7   -39 CAGCGCCAT   *
..include:seqan/bam_io.h

.Memfunc.BamStream#BamStream:
..class:Class.BamStream
..description:See documentation of @Class.BamStream@ for more information.
..summary:Constructor
..signature:BamStream()
..signature:BamStream(fileName[, mode[, format]])
..param.fileName:Path to the file to open.
...type:nolink:$char const *$
..param.mode:The mode to use for opening the file (read/write). Optional
...default:@Enum.BamStream\colon\colonOperationMode.value.READ@
...type:Enum.BamStream\colon\colonOperationMode
..param.format:Use this to enforce opening the file in the given format. Autodetected from file name or content if not specified. Optional.
...type:Enum.BamStream\colon\colonFormat
...default:@Enum.BamStream\colon\colonFormat.value.AUTO@

.Memvar.BamStream#header:
..class:Class.BamStream
..type:Class.BamHeader
..summary:The @Class.BamHeader@ of the @Class.BamStream@ object.
..description:
SAM and BAM files have a header.
When writing SAM or BAM files, you have to fill this member before writing @Class.BamAlignmentRecord@s.
Upon writing the first record, the header will be written out.
..description:
When reading BAM files, the header will be read upon opening the file.
When reading SAM files, any header will be read upon opening the file.
..description:
Note that there is a special case when reading SAM records:
If there is no header, or records refer to reference sequences that are previously unknown when reading SAM then a new entry is added to @Memvar.BamHeader#sequenceInfos@.\

.Memvar.BamStream#bamIOContext:
..class:Class.BamStream
..summary:The @Class.BamIOContext@ object to use for reading and writing @Class.BamAlignmentRecord@s.
..description:
When reading, the $bamIOContext$ will be updated automatically.
When reading SAM, new reference sequences can be introduced "on the fly" when a new sequence appears.
When writing, the $bamIOContext$ is automatically filled/reset when the first record is written.

.Enum.BamStream\colon\colonOperationMode:
..cat:BAM I/O
..summary:Select the operation mode of a @Class.BamStream@.
..value.READ:Open stream for reading.
..value.WRITE:Open stream for writing.
..include:seqan/bam_io.h

.Enum.BamStream\colon\colonFormat:
..cat:BAM I/O
..summary:Select the format to use for reading/writing.
..value.AUTO:Auto-detect format from file content on reading and from the file name on writing. If Auto-detection fails, SAM is used.
..value.SAM:Force reading/writing of SAM.
..value.BAM:Force reading/writing of BAM.
..include:seqan/bam_io.h
*/

/*!
 * @enum BamStream::OperationMode
 * @brief Selects teh operation mode of a @link BamStream @endlink.
 * @see BamStream
 *
 * @signature enum BamStream::OperationMode;
 *
 * @var BamStream::OperationMode BamStream::READ;
 * @brief Enum value for reading.
 *
 * @var BamStream::OperationMode BamStream::WRITE;
 * @brief Enum value for writing.
 */

/*!
 * @enum BamStream::Format
 * @brief Select the format to use for reading/writing.
 *
 * @signature enum BamStream::Format;
 *
 * @var BamStream::Format BamStream::AUTO;
 * @brief Auto-detect the format from file content on reading and from the file name on writing.  If auto-detection
 *        fails, SAM is used.
 *
 * @var BamStream::Format BamStream::SAM;
 * @brief Force reading/writing of SAM.
 *
 * @var BamStream::Format BamStream::BAM;
 * @brief Force reading/writing of BAM.
 */

class BamStream
{
public:

    // Enum for selecting read/write mode.
    enum OperationMode
    {
        READ,
        WRITE
    };

    // Enum for selecting format.  AUTO is only used as the default, after opening, only SAM and BAM are used.
    enum Format
    {
        AUTO,
        SAM,
        BAM
    };

    // Name of the BAM file.
    CharString _filename;
    // The open mode.
    OperationMode _mode;
    // The format.
    Format _format;
    // Whether or not the header was written out.
    bool _headerWritten;

    // Indicates whether stream is at end when reading.
    bool _atEnd;
    // Indicates whether there was an error when reading or writing.
    bool _isGood;

    // The BAM Header record.
    BamHeader header;
    // The BAM I/O Context and its elements.
    StringSet<CharString>                  _nameStore;
    NameStoreCache<StringSet<CharString> > _nameStoreCache;
    BamIOContext<StringSet<CharString> >   bamIOContext;

    // The actual implementation of writing SAM or BAM.
    std::SEQAN_AUTO_PTR_NAME<XamWriter_> _writer;
    // The actual implementation of reading SAM or BAM.
    std::SEQAN_AUTO_PTR_NAME<XamReader_> _reader;

    // Constructors.

    BamStream() :
        _mode(READ), _format(AUTO), _headerWritten(false), _atEnd(false), _isGood(true),
        _nameStoreCache(_nameStore), bamIOContext(_nameStore, _nameStoreCache)
    {}

    BamStream(char const * filename, OperationMode mode = READ, Format format = AUTO);

    // Write header if necessary.
    inline int _writeHeader()
    {
        if (this->_headerWritten)
            return 0;

        // Rewrite name store and cache.
        clear(_nameStore);
        for (unsigned i = 0; i < length(header.sequenceInfos); ++i)
            appendValue(_nameStore, header.sequenceInfos[i].i1);
        refresh(_nameStoreCache);

        // Write out header.
        this->_headerWritten = true;
        return this->_writer->writeHeader(header, bamIOContext);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Member Function BamStream::BamStream()
// ----------------------------------------------------------------------------

// Forward declaration is here since it refers to enum member type.
inline int open(BamStream & bamIO,
                char const * filename,
                BamStream::OperationMode mode,
                BamStream::Format format);

inline BamStream::BamStream(char const * filename, OperationMode mode, Format format) :
    _filename(filename), _mode(mode), _format(format), _headerWritten(false), _atEnd(false), _isGood(true),
    _nameStoreCache(_nameStore), bamIOContext(_nameStore, _nameStoreCache)
{
    open(*this, filename, _mode, _format);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#open
 * @brief Open a @link BamStream @endlink object for reading/writing.
 *
 * @signature int open(bamIO, fileName[, mode[, format]]);
 *
 * @param[in,out] bamIO    The @link BamStream @endlink object to open. Types: BamStream
 * @param[in]     fileName The path to the file to open, <tt>char const *</tt>.
 * @param[in]     mode     The mode to open the file in, optional, of type @link BamStream::OperationMode @endlink,
 *                         defaults to <tt>BamStream::READ</tt>.
 * @param[in]     format   The format to use, inferred from file contents (reading) or file name (writing) by default.
 *                         the path to the file to open, of type @link BamStream::Format @endlink, defaults to
 *                         <tt>AUTO</tt>.
 *
 * @return int A status code, 0 on success, a value <tt>!= 0</tt> on errors.
 */

/**
.Function.BamStream#open
..class:Class.BamStream
..summary:Open a @Class.BamStream@ object for reading/writing.
..signature:open(bamIO, fileName[, mode[, format]])
..param.bamIO:The @Class.BamStream@ object to open.
...type:Class.BamStream
..param.fileName:The path to the file to open.
...type:Shortcut.CharString
..param.mode:The mode to open the file in. Optional.
...default:$BamStream::READ$
...type:nolink:$BamStream::OperationMode$.
..param.format:The format to use, inferred from file contents (reading) or file name (writing) by default.
...default:$BamStream::AUTO$
...type:nolink:$BamStream::Format$.
..param.format:The path to the file to open.
...type:Shortcut.CharString
..returns:An $int$ status code: $0$ on success, $1$ on errors.
..include:seqan/bam_io.h
*/

inline int open(BamStream & bamIO,
                char const * fileName,
                BamStream::OperationMode mode = BamStream::READ,
                BamStream::Format format = BamStream::AUTO)
{
    bamIO._filename = fileName;
    bamIO._isGood = true;

    // Guess format if necessary.
    if (format == BamStream::AUTO)
    {
        if (mode == BamStream::READ)
        {
            format = BamStream::SAM;  // SAM is default.

            // Look whether the file is in BAM format.
            std::fstream inStream(toCString(fileName), std::ios_base::binary | std::ios_base::in);
            if (!inStream.good())
            {
                bamIO._isGood = false;
                return 1;  // Error opening the file.
            }
            char buffer[3];
            inStream.read(&buffer[0], 3);
            if (buffer[0] == '\x1F' && buffer[1] == '\x8B' && buffer[2] == '\x08')
                format = BamStream::BAM;
        }
        else  // mode == WRITE
        {
            format = BamStream::SAM;  // SAM is default.
            if (endsWith(fileName, ".bam"))
                format = BamStream::BAM;
        }
    }

#if !SEQAN_HAS_ZLIB
    // Guard against opening BAM files without zlib.
    if (format == BamStream::BAM)
    {
        std::cerr << "ERROR: Trying to open BAM file and ZLIB is not available!\n";
        bamIO._isGood = false;
        return 1;
    }
#endif  // #if !SEQAN_HAS_ZLIB

    if (mode == BamStream::READ)
    {
        if (format == BamStream::SAM)
            bamIO._reader.reset(new SamReader_());
#if SEQAN_HAS_ZLIB
        // The branch above is always taken if zlib is not available, there already is a check above.
        else
            bamIO._reader.reset(new BamReader_());
#endif  // #if !SEQAN_HAS_ZLIB
        if (bamIO._reader->open(fileName) != 0)
        {
            bamIO._isGood = false;
            return 1;
        }
    }
    else  // (format == BamStream::WRITE)
    {
        if (format == BamStream::SAM)
            bamIO._writer.reset(new SamWriter_());
#if SEQAN_HAS_ZLIB
        // The branch above is always taken if zlib is not available, there already is a check above.
        else
            bamIO._writer.reset(new BamWriter_());
#endif  // #if !SEQAN_HAS_ZLIB
        if (bamIO._writer->open(fileName) != 0)
        {
            bamIO._isGood = false;
            return 1;
        }
    }

    bamIO._mode = mode;
    bamIO._format = format;

    // Read header.
    if (bamIO._isGood && bamIO._mode == BamStream::READ)
    {
        clear(bamIO.header);
        if (bamIO._reader->readHeader(bamIO.header, bamIO.bamIOContext) != 0)
        {
            bamIO._isGood = false;
            return 1;
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function reset()
// ----------------------------------------------------------------------------

/*
 * @fn BamStream#reset
 * @brief Reset @link BamStream @endlink object to status after construction.
 *
 * @signature void reset(stream);
 *
 * @param stream The @link BamStream @endlink object to reset.
 *
 * @return int A status code, 0 on success, <tt>!= 0</tt> on error.
 */

/**
.Function.BamStream#reset
..class:Class.BamStream
..summary:Reset @Class.BamStream@ object to status after construction.
..signature:reset(bamIO)
..param.bamIO:The @Class.BamStream@ object to reset.
...type:Class.BamStream
..returns:$int$, a status code ($0$ for success, non-$0$ for error).
..include:seqan/bam_io.h
*/

inline int reset(BamStream & bamIO)
{
    return open(bamIO, toCString(bamIO._filename), bamIO._mode, bamIO._format);
}

// ----------------------------------------------------------------------------
// Function flush()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#flush
 * @brief Flush output when writing.
 *
 * @signature int flush(stream);
 *
 * @param stream The @link BamStream @endlink object to flush.
 *
 * @return int A status code, 0 on success, <tt>!= 0</tt> on errors.
 *
 * @section Remarks
 *
 * This will write out the header if no record has been written out yet.
 */

/**
.Function.BamStream#flush
..class:Class.BamStream
..summary:Flush output when writing.
..description:This will write out the header if no record has been written out yet.
..signature:flush(bamIO)
..param.bamIO:The @Class.BamStream@ object to flush.
...type:Class.BamStream
..returns:$int$ with an error code.
..include:seqan/bam_io.h
*/

inline int flush(BamStream & bamIO)
{
    if (bamIO._mode == BamStream::WRITE)
    {
        bamIO._writeHeader();
        return bamIO._writer->flush();
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#close
 * @brief Close BamStream object's underlying file.
 *
 * @signature int close(stream);
 *
 * @param stream[in,out] The @link BamStream @endlink object to close.
 *
 * @return int A status code, 0 on success, <tt>!= 0</tt> on error.
 */

/**
.Function.BamStream#close
..class:Class.BamStream
..summary:Close BamStream object's underlying file.
..signature:close(bamIO)
..param.bamIO:The @Class.BamStream@ object to close
...type:Class.BamStream
..returns:$int$ with an error code.
..include:seqan/bam_io.h
*/

inline int close(BamStream & bamIO)
{
    if (bamIO._mode == BamStream::WRITE)
    {
        bamIO._writeHeader();
        return bamIO._writer->close();
    }
    else
    {
        return bamIO._reader->close();
    }
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#atEnd
 * @brief Check whether a @link BamStream @endlink object is at end when reading.
 *
 * @signature bool atEnd(stream);
 *
 * @param[in] stream The @link BamStream @endlink object to query.
 *
 * @return bool true in case of the stream being at the end, false otherwise.
 *
 * @section Remarks
 *
 * The stream will only be guaranteed at the end after trying to read <b>after</b> the last character.
 */

/**
.Function.BamStream#atEnd
..class:Class.BamStream
..summary:Check whether a @Class.BamStream@ object is at end when reading.
..signature:atEnd(bamIO)
..param.bamIO:The @Class.BamStream@ object to query.
...type:Class.BamStream
..returns:$bool$, indicating whether the object is at the end of the file.
..include:seqan/bam_io.h
*/

inline bool atEnd(BamStream const & bamIO)
{
    SEQAN_ASSERT_EQ_MSG(bamIO._mode, BamStream::READ, "You can only call atEnd() when opened the file for reading.");
    return bamIO._reader->atEnd();
}

inline bool atEnd(BamStream & bamIO)
{
    SEQAN_ASSERT_EQ_MSG(bamIO._mode, BamStream::READ, "You can only call atEnd() when opened the file for reading.");
    return bamIO._reader->atEnd();
}

// ----------------------------------------------------------------------------
// Function isGood()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#isGood
 * @brief Check whether the @link BamStream @endlink object has is in the failure state.
 *
 * @signature bool isGood(stream);
 *
 * @param stream The @link BamStream @endlink object to query.
 *
 * @return bool true if the stream is not in an error state and false otherwise.
 */

/**
.Function.BamStream#isGood
..class:Class.BamStream
..summary:Check whether the @Class.BamStream@ object has is in the failure state.
..signature:isGood(bamIO)
..param.bamIO:The @Class.BamStream@ object to query.
...type:Class.BamStream
..returns:$bool$, indicating whether there was no error or not.
..include:seqan/bam_io.h
*/

inline bool isGood(BamStream const & bamIO)
{
    return bamIO._isGood;
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#readRecord
 * @brief Read one @link BamAlignmentRecord @endlink from a @link BamStream @endlink.
 *
 * @signature int readRecord(record, stream);
 *
 * @param[out]   record    The @link BamAlignmentRecord @endlink to read the next alignment record into.  Of type
 *                         @link BamAlignmentRecord @endlink.
 * @param[in,out] stream   The @link BamStream @endlink object to read from.
 *
 * @return int A status code, 0 on success.
 */

/**
.Function.BamStream#readRecord
..class:Class.BamStream
..summary:Read one @Class.BamAlignmentRecord@ from a @Class.BamStream@.
..signature:readRecord(record, bamIO)
..param.record:The @Class.BamAlignmentRecord@ to read the next alignment record into.
...class:Class.BamAlignmentRecord
..param.bamIO:The @Class.BamStream@ object to read from.
...type:Class.BamStream
..returns:An $int$ status code: $0$ on success, non-$0$ on failure.
..include:seqan/bam_io.h
*/

inline int readRecord(BamAlignmentRecord & record, BamStream & bamIO)
{
    int res = bamIO._reader->readRecord(record, bamIO.bamIOContext);
    bamIO._isGood = bamIO._isGood && (res == 0);
    return res;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#writeRecord
 * @brief Write one @link BamAlignmentRecord @endlink to a @link BamStream @endlink.
 *
 * @signature int writeRecord(stream, record);
 *
 * @param[in,out] bamIO  The @link BamStream @endlink object to write to.
 * @param[in]     record The @link BamAlignmentRecord @endlink to write out.
 *
 * @return int A status code, 0 on success.
 */

/**
.Function.BamStream#writeRecord
..class:Class.BamStream
..summary:Write one @Class.BamAlignmentRecord@ to a @Class.BamStream@.
..signature:writeRecord(bamIO, record)
..param.record:The @Class.BamAlignmentRecord@ to write out.
...class:Class.BamAlignmentRecord
..param.bamIO:The @Class.BamStream@ object to write to.
...type:Class.BamStream
..returns:An $int$ status code: $0$ on success, non-$0$ on failure.
..include:seqan/bam_io.h
*/

inline int writeRecord(BamStream & bamIO, BamAlignmentRecord const & record)
{
    bamIO._writeHeader();  // Does nothing if head already written out.

    int res = bamIO._writer->writeRecord(record, bamIO.bamIOContext);
    bamIO._isGood = bamIO._isGood && (res == 0);
    return res;
}

// ----------------------------------------------------------------------------
// Function fileSize()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#fileSize
 * @brief Returns the size of the file in bytes as stored on the disk.
 *
 * @signature __int64 fileSize(stream);
 *
 * @param[in] stream The @link BamStream @endlink to query.
 *
 * @return __int64 The size of the file on the disk.
 *
 * @section Remarks
 *
 * This only works when reading.
 */

// Returns size of file in bytes as stored on the disk.

inline __int64 fileSize(BamStream const & bamIO)
{
    if (bamIO._mode == BamStream::WRITE)
        return 0;
    return bamIO._reader->fileSize();
}

// ----------------------------------------------------------------------------
// Function positionInFile()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#positionInFile
 * @brief Approximate byte position in file, to be used for progress display, not for seeking.
 *
 * @signature __int64 positionInFile(stream);
 *
 * @param[in] stream The @link BamStream @endlink to query for its position in the file.
 *
 * @return __int64 The position in the file.
 *
 * @section Remarks
 *
 * This function returns the "approximate" position in the file and only works when the file is opened in BAM format.
 * It is meant for progress display in connection with @link BamStream#fileSize @endlink and not for jumping within the
 * file.  The position is approximate in the sense that it points between the block boundaries of the BGZ file.
 */

// TODO(holtgrew): Review this functionality, extend, fix.

// Returns "approximate" byte position in file.  To be used for progress display, not for seeking.  For this, we have to
// implement streamTell() and streamSeek() for BamStream.  This works for BAM only at the moment.

inline __int64 positionInFile(BamStream const & bamIO)
{
    if (bamIO._mode == BamStream::WRITE)
        return 0;
    return bamIO._reader->positionInFile();
}

// ----------------------------------------------------------------------------
// Function jumpToRegion()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#jumpToRegion
 * @brief Seek in BamStream using an index.
 *
 * You provide a region <tt>[pos, posEnd)</tt> on the reference <tt>refID</tt> that you want to jump to and the function
 * jumps to the first alignment in this region, if any.
 *
 * @signature bool jumpToRegion(stream, hasAlignments, bamIOContext, refID, pos, posEnd, index);
 *
 * @param[in,out] stream        The @link BamStream @endlink to jump with.
 * @param[out]    hasAlignments A <tt>bool</tt> that is set true if the region <tt>[pos, posEnd)</tt> has any
 *                              alignments.
 * @param[in]     refID         The reference id to jump to (<tt>__int32</tt>).
 * @param[in]     pos           The begin of the region to jump to.
 * @param[in]     posEnd        The end of the region to jump to.
 * @param[in]     index         The @link BamIndex @endlink to use for the jumping.
 *
 * @return bool true if seeking was successful, false if not.
 *
 * @section Remarks
 *
 * This function fails if <tt>refID</tt>/<tt>pos</tt> are invalid.
 *
 * @see BamIndex#jumpToRegion
 */

#if SEQAN_HAS_ZLIB
inline bool jumpToRegion(BamStream & bamIO, bool & hasAlignments, __int32 refId, __int32 pos, __int32 posEnd, BamIndex<Bai> const & index)
{
    if (bamIO._format != BamStream::BAM)
        return false;  // Can only jump in BAM files.
    if (bamIO._mode != BamStream::READ)
        return false;  // Can only jump when reading.

    BamReader_ * s = static_cast<BamReader_ *>(bamIO._reader.get());
    return s->jumpToRegion(hasAlignments, refId, pos, posEnd, index, bamIO.bamIOContext);
}
#endif  // #if SEQAN_HAS_ZLIB

// ----------------------------------------------------------------------------
// Function jumpToOrphans()
// ----------------------------------------------------------------------------

/*!
 * @fn BamStream#jumpToOrphans
 * @brief Seek to orphans block in BamStream using an index.
 *
 * @signature bool jumpToOrphans(stream, hasAlignments, index);
 *
 * @param[in,out] stream         The @link BgzfStream @endlink object to jump with.
 * @param[out]    hasAlignments  A <tt>bool</tt> that is set to true if there are any orphans.
 * @param[in]     index          The index to use for jumping.
 *
 * @see BamIndex#jumpToOrphans
 */

#if SEQAN_HAS_ZLIB
inline bool jumpToOrphans(BamStream & bamIO, BamIndex<Bai> const & index)
{
    if (bamIO._format != BamStream::BAM)
        return false;  // Can only jump in BAM files.
    if (bamIO._mode != BamStream::READ)
        return false;  // Can only jump when reading.

    BamReader_ * s = static_cast<BamReader_ *>(bamIO._reader.get());
    return s->jumpToOrphans(index, bamIO.bamIOContext);
}
#endif  // #if SEQAN_HAS_ZLIB

}  // namespace seqan;

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_STREAM_H_
