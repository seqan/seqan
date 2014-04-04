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

// TODO(holtgrew): Currently the parsers for GFF and GTF are the same. Thus, we need no file format tag for reading.  We always write GFF.

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_GFF_IO_GFF_STREAM_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_GFF_IO_GFF_STREAM_H_

#include <memory>
#include <fstream>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class GffStream
 * @headerfile <seqan/gff_io.h>
 * @brief High-level GFF/GTF I/O class.
 * 
 * @signature class GffStream;
 * 
 * The GffStream class allows to read and write GTF and GFF files.
 * 
 * @section Remarks
 * 
 * Note that the class <tt>GffStream</tt> allows to read both GFF 2, 3, and GTF.  For writing, only GFF 3 and GTF are
 * supported.
 * 
 * @section Examples
 * 
 * The following example demonstrates reading a GFF file and printing the annotation locations.
 * 
 * @include demos/gff_io/gff_stream_read.cpp
 * 
 * @see GffStream::FileFormat
 * 
 *
 * @var TCharStringSet GffStream::sequenceNames
 * @brief The names of the sequences (StringSet of CharString), updated when new sequences are seen in GFF file.
 * 
 * @var GffStream::FileFormat GffStream::fileFormat
 * @brief File format to use for writing.
 * 
 *
 * @fn GffStream::GffStream
 * @brief Constructor.
 * 
 * @signature GffStream::GffStream();
 * @signature GffStream::GffStream(fileName[, mode=READ[, fileFormat=GFF]]);
 * 
 * @param[in] fileName   The path to the file to open, <tt>char const *</tt>.
 * @param[in] mode       The open mode, GffStream::Mode, default is READ.
 * @param[in] fileFormat The open mode, GffStream::FileFormat, default is auto-detect.
 * 
 * @see GffStream::Mode
 *
 * @see GffStream#open
 */

/**
.Class.GffStream
..cat:GFF I/O
..summary:High-level GFF/GTF I/O class.
..description:The GffStream class allows to read and write GTF and GFF files.
..signature:class GffStream
..remarks:
Note that the class $GffStream$ allows to read both GFF 2, 3, and GTF.
For writing, only GFF 3 and GTF are supported.
..example:The following example demonstrates reading a GFF file and printing the annotation locations.
..example.file:demos/gff_io/gff_stream_read.cpp
..include:seqan/gff_io.h
..see:Enum.GffStream\colon\colonFileFormat

.Memfunc.GffStream#GffStream
..class:Class.GffStream
..summary:Constructor.
..signature:GffStream::GffStream()
..signature:GffStream::GffStream(fileName[, mode=READ[, fileFormat=GFF]])
..param.fileName:The path to the file to open.
...type:nolink:$char const *$
..param.mode:The open mode.
...type:Enum.GffStream\colon\colonMode
...default:$GffStream\colon\colonMode::READ$
..see:Enum.GffStream\colon\colonMode

.Memvar.GffStream#fileFormat
..class:Class.GffStream
..summary:File format to use for writing.

.Memvar.GffStream#sequenceNames
..class:Class.GffStream
..summary:The names of the sequences (@Class.StringSet@ of @Shortcut.CharString@), updated when new sequences are seen in GFF file.

.Enum.GffStream\colon\colonMode
..cat:GFF I/O
..summary:Open mode for the @Class.GffStream@ class.
..value.INVALID:Invalid open mode.
..value.READ:Open in read mode.
..value.WRITE:Open in write mode.
..include:seqan/gff_io.h

.Enum.GffStream\colon\colonFileFormat
..cat:GFF I/O
..summary:File format for writing in the @Class.GffStream@ class.
..value.GFF:GFF file format.
..value.GTF:GTF file format.
..include:seqan/gff_io.h
*/

/*!
 * @enum GffStream::FileFormat
 * @headerfile <seqan/gff_io.h>
 *
 * @signature enum GffStream::FileFormat;
 * 
 * @brief File format for writing in the @link GffStream @endlink class.
 * 
 * @see GffStream
 * 
 * @val GffStream::FileFormat GffStream::GFF
 * @brief GFF file format.
 * 
 * @val GffStream::FileFormat GffStream::GTF
 * @brief GTF file format.
 */

/*!
 * @enum GffStream::Mode
 * @headerfile <seqan/gff_io.h>
 * @brief Open mode for the GffStream class.
 * 
 * @signature enum GffStream::Mode;
 *
 * @see GffStream#open
 * @see GffStream::GffStream
 * 
 * @val GffStream::Mode GffStream::READ
 * @brief Open in read mode.
 * 
 * @val GffStream::Mode GffStream::WRITE
 * @brief Open in write mode.
 * 
 * @val GffStream::Mode GffStream::INVALID
 * @brief Invalid open mode.
 */

class GffStream
{
public:
    typedef RecordReader<std::istream, SinglePass<> > TReader_;
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > TNameStoreCache;
    typedef GffIOContext<TNameStore, TNameStoreCache> TGffIOContext;

    // The open mode.
    enum Mode
    {
        INVALID,
        READ,
        WRITE
    };

    // The file format.
    enum FileFormat
    {
        GFF,
        GTF
    };

    std::SEQAN_AUTO_PTR_NAME<std::fstream> _stream;
    std::ostream * _outStream;
    std::istream * _inStream;
    CharString _filename;
    std::SEQAN_AUTO_PTR_NAME<TReader_> _reader;
    Mode _mode;
    FileFormat fileFormat;
    int _error;
    bool _isGood;

    TNameStore sequenceNames;
    TNameStoreCache _sequenceNamesCache;
    TGffIOContext _context;

    GffStream() : _outStream(), _inStream(), _mode(INVALID), fileFormat(GFF), _error(0), _isGood(true),
                  _sequenceNamesCache(sequenceNames), _context(sequenceNames, _sequenceNamesCache)
    {}

    GffStream(char const * filename, Mode mode = READ, FileFormat fileFormat = GFF) :
            _outStream(), _inStream(), _filename(filename), _mode(mode), fileFormat(fileFormat), _error(0),
            _isGood(true), _sequenceNamesCache(sequenceNames), _context(sequenceNames, _sequenceNamesCache)
    {
        _open(filename, mode);
    }

    bool _open(char const * filename, Mode mode)
    {
        // Reset.
        _filename = filename;
        _mode = mode;
        _error = 0;
        _isGood = true;

        if (mode == READ)
        {
            if (_filename == "-")
            {
                _stream.reset();
                _inStream = &std::cin;
            }
            else
            {
                _stream.reset(new std::fstream);
                _stream->open(toCString(_filename), std::ios::binary | std::ios::in);
                if (!_stream->good())
                {
                    _isGood = false;
                    return false;
                }
                _inStream = _stream.get();
            }
            _reader.reset(new TReader_(*_inStream));
            _outStream = 0;
        }
        else if (mode == WRITE)
        {
            if (_filename == "-")
            {
                _stream.reset();
                _outStream = &std::cout;
            }
            else
            {
                _stream.reset(new std::fstream);
                _stream->open(toCString(_filename), std::ios::binary | std::ios::out);
                if (!_stream->good())
                {
                    _isGood = false;
                    return false;
                }
                _reader.reset();
                _outStream = _stream.get();
            }
            _inStream = 0;
        }
        return true;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#open
 * @brief Open a GffStream .
 * 
 * @signature bool open(gffStream, fileName[, mode[, fileFormat]]);
 * 
 * @param[in,out] gffStream  The GffStream to open.
 * @param[in]     fileName   The path to the file to open, <tt>char const *</tt>.
 * @param[in]     mode       The open mode, GffStream::Mode, defaults to <tt>GffStream::Mode::READ</tt>.
 * @param[in]     fileFormat File format to use for writing, GffStream::FileFormat.
 * 
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 * 
 * @see GffStream#isGood
 * @see GffStream::Mode
 */

/**
.Function.GffStream#open
..class:Class.GffStream
..cat:GFF I/O
..summary:Open a @Class.GffStream@.
..signature:bool open(gffStream, fileName[, mode[, fileFormat]])
..param.gffStream:The @Class.GffStream@ to open.
...type:Class.GffStream
..param.fileName:The path to the file to open.
...type:nolink:$char const *$
..param.mode:The open mode.
...type:Enum.GffStream\colon\colonMode
...default:$GffStream\colon\colonMode::READ$
..param.fileFormat:File format to use for writing.
...type:Enum.GffStream\colon\colonFileFormat
..returns:$true$ on success, $false$ on failure.
..see:Function.GffStream#isGood
..see:Enum.GffStream\colon\colonMode
..include:seqan/gff_io.h
*/

inline bool open(GffStream & stream, char const * filename, GffStream::Mode mode = GffStream::READ,
                 GffStream::FileFormat fileFormat = GffStream::GFF)
{
    stream.fileFormat = fileFormat;
    return stream._open(filename, mode);
}

// ----------------------------------------------------------------------------
// Function addSequenceName()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#addSequenceName
 * @brief Add the name of a sequence to a GffStream .
 * 
 * @signature void addSequenceName(gffStream, seqName);
 * 
 * @param[in,out] gffStream The GffStream to add the name to
 * @param[in]     seqName   The name of the sequence to append, CharString.
 */

/**
.Function.GffStream#addSequenceName
..class:Class.GffStream
..cat:GFF I/O
..summary:Add the name of a sequence to a @Class.GffStream@.
..signature:void addSequenceName(gffStream, seqName);
..param.gffStream:The @Class.GffStream@ to add the name to.
...type:Class.GffStream
..param.seqName:The name of the sequence to append.
...type:Shortcut.CharString
..include:seqan/gff_io.h
*/

inline void addSequenceName(GffStream & stream, CharString const & name)
{
    appendName(stream.sequenceNames, name, stream._sequenceNamesCache);
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#readRecord
 * @brief Read a record from a GffStream
 * 
 * @signature int readRecord(record, gffStream);
 * 
 * @param[out]    record    The GffRecord to read into, GffRecord.
 * @param[in,out] gffStream The GffStream to read from, GffStream.
 * 
 * @return int <tt>0</tt> on success, non-<tt>0</tt> on failure.
 */

/**
.Function.GffStream#readRecord
..class:Class.GffStream
..cat:GFF I/O
..summary:Read a record from a @Class.GffStream@
..signature:int readRecord(record, gffStream)
..param.record:The @Class.GffRecord@ to read into.
...type:Class.GffRecord
..param.gffStream:The @Class.GffStream@ to read from.
...type:Class.GffStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/gff_io.h
*/

inline int readRecord(GffRecord & record,
                      GffStream & stream)
{
    int res = readRecord(record, *stream._reader, stream._context, Gff());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#writeRecord
 * @brief Write a record to a GffStream.
 *
 * @signature int writeRecord(gffStream, record);
 * 
 * @param[in,out] gffStream The GffStream to write to.
 * @param[in]     record    The GffRecord to write.
 * 
 * @return int <tt>0</tt> on success, non-<tt>0</tt> on failure.
 */

/**
.Function.GffStream#writeRecord
..class:Class.GffStream
..cat:GFF I/O
..summary:Write a record to a @Class.GffStream@
..signature:int writeRecord(gffStream, record)
..param.gffStream:The @Class.GffStream@ to write to.
...type:Class.GffStream
..param.record:The @Class.GffRecord@ to write.
...type:Class.GffRecord
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/gff_io.h
*/

inline int writeRecord(GffStream & stream,
                       GffRecord const & record)
{
    int res = 0;
    if (stream.fileFormat == GffStream::GFF)
        res = writeRecord(*stream._outStream, record, stream._context, Gff());
    else
        res = writeRecord(*stream._outStream, record, stream._context, Gtf());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function flush()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#flush
 * @brief Flush to a GffStream.
 * 
 * @signature int flush(gffStream);
 * 
 * @param[in,out] gffStream The GffStream to flush. Types: GffStream
 * 
 * @return int <tt>0</tt> on success, non-<tt>0</tt> on failure.
 */

/**
.Function.GffStream#flush
..class:Class.GffStream
..cat:GFF I/O
..summary:Flush to a @Class.GffStream@
..signature:int flush(gffStream)
..param.gffStream:The @Class.GffStream@ to flush.
...type:Class.GffStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/gff_io.h
*/

inline int flush(GffStream & stream)
{
    if (stream._stream.get())
        stream._stream->flush();
    return 0;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#close
 * @brief Closes a GffStream
 * 
 * @signature int close(gffStream);
 * 
 * @param[in,out] gffStream The GffStream to close. Types: GffStream
 * 
 * @return int <tt>0</tt> on success, non-<tt>0</tt> on failure.
 */

/**
.Function.GffStream#close
..class:Class.GffStream
..cat:GFF I/O
..summary:Closes a @Class.GffStream@
..signature:int close(gffStream)
..param.gffStream:The @Class.GffStream@ to close.
...type:Class.GffStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/gff_io.h
*/

inline int close(GffStream & stream)
{
    if (stream._stream.get())
        stream._stream->close();
    return 0;
}

// ----------------------------------------------------------------------------
// Function isGood()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#isGood
 * @brief Query a GffStream for errors.
 * 
 * @signature bool isGood(gffStream);
 * 
 * @param[in] gffStream The GffStream to query.
 * 
 * @return bool <tt>true</tt> if stream is good, <tt>false</tt> otherwise.
 */

/**
.Function.GffStream#isGood
..class:Class.GffStream
..cat:GFF I/O
..summary:Query a @Class.GffStream@ for errors.
..signature:bool isGood(gffStream)
..param.gffStream:The @Class.GffStream@ to query.
...type:Class.GffStream
..returns:$true$ if stream is good, $false$ otherwise.
..include:seqan/gff_io.h
*/

inline bool isGood(GffStream const & stream)
{
    return stream._isGood;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

/*!
 * @fn GffStream#atEnd
 * @brief Query a GffStream for being at the end of the file.
 * 
 * @signature bool atEnd(gffStream);
 * 
 * @param[in] gffStream The GffStream to query. Types: GffStream
 * 
 * @return bool <tt>true</tt> if stream is at the end, <tt>false</tt> otherwise.
 */

/**
.Function.GffStream#atEnd
..class:Class.GffStream
..cat:GFF I/O
..summary:Query a @Class.GffStream@ for being at the end of the file.
..signature:bool atEnd(gffStream)
..param.gffStream:The @Class.GffStream@ to query.
...type:Class.GffStream
..returns:$true$ if stream is at the end, $false$ otherwise.
..include:seqan/gff_io.h
*/

inline bool atEnd(GffStream const & stream)
{
    return atEnd(*stream._reader);
}

inline bool atEnd(GffStream & stream)
{
    return atEnd(*stream._reader);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_GFF_IO_GFF_STREAM_H_
