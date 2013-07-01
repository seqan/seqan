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

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_BED_IO_BED_STREAM_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_BED_IO_BED_STREAM_H_

#include <memory>
#include <fstream>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.BedStream
..cat:BED I/O
..summary:High-level BED I/O class.
..signature:class BedStream
..example:The following example demonstrates reading a BED file and printing the variant locations.
..example.file:demos/bed_io/bed_stream_read.cpp
..see:Enum.BedStream\colon\colonMode
..include:seqan/bed_io.h

.Memfunc.BedStream#BedStream
..class:Class.BedStream
..summary:Constructor.
..signature:BedStream::BedStream()
..signature:BedStream::BedStream(fileName, mode=READ)
..param.fileName:The path to the file to open.
...type:nolink:$char const *$
..param.mode:The open mode.
...type:Enum.BedStream\colon\colonMode
...default:$BedStream\colon\colonMode::READ$
..see:Enum.BedStream\colon\colonMode

.Memvar.BedStream#sequenceNames
..class:Class.BedStream
..summary:The names of the sequences (@Class.StringSet@ of @Shortcut.CharString@), updated when new sequences are seen in BED file.

.Enum.BedStream\colon\colonMode
..cat:BED I/O
..summary:Open mode for the @Class.BedStream@ class.
..value.INVALID:Invalid open mode.
..value.READ:Open in read mode.
..value.WRITE:Open in write mode.
..include:seqan/bed_io.h
*/

class BedStream
{
public:
    typedef RecordReader<std::istream, SinglePass<> > TReader_;
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > TNameStoreCache;
    typedef BedIOContext<TNameStore, TNameStoreCache> TBedIOContext;

    enum Mode
    {
        INVALID,
        READ,
        WRITE
    };

    std::auto_ptr<std::fstream> _stream;
    std::ostream * _outStream;
    std::istream * _inStream;
    CharString _filename;
    std::auto_ptr<TReader_> _reader;
    Mode _mode;
    int _error;
    bool _isGood;

    TNameStore sequenceNames;
    TNameStoreCache _sequenceNamesCache;
    TBedIOContext _context;

    BedStream() : _outStream(), _inStream(), _mode(INVALID), _error(0), _isGood(true),
                  _sequenceNamesCache(sequenceNames), _context(sequenceNames, _sequenceNamesCache)
    {}

    BedStream(char const * filename, Mode mode = READ) :
            _outStream(), _inStream(), _filename(filename), _mode(mode), _error(0), _isGood(true),
            _sequenceNamesCache(sequenceNames), _context(sequenceNames, _sequenceNamesCache)
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

/**
.Function.BedStream#open
..class:Class.BedStream
..cat:BED I/O
..summary:Open a @Class.BedStream@.
..signature:bool open(bedStream, fileName, mode)
..param.bedStream:The @Class.BedStream@ to open.
...type:Class.BedStream
..param.fileName:The path to the file to open.
...type:nolink:$char const *$
..param.mode:The open mode.
...type:Enum.BedStream\colon\colonMode
...default:$BedStream\colon\colonMode::READ$
..returns:$true$ on success, $false$ on failure.
..see:Function.BedStream#isGood
..see:Enum.BedStream\colon\colonMode
..include:seqan/bed_io.h
*/

inline bool open(BedStream & stream, char const * filename, BedStream::Mode mode = BedStream::READ)
{
    return stream._open(filename, mode);
}

// ----------------------------------------------------------------------------
// Function addSequenceName()
// ----------------------------------------------------------------------------

/**
.Function.BedStream#addSequenceName
..class:Class.BedStream
..cat:BED I/O
..summary:Add the name of a sequence to a @Class.BedStream@.
..signature:void addSequenceName(bedStream, seqName);
..param.bedStream:The @Class.BedStream@ to add the name to.
...type:Class.BedStream
..param.seqName:The name of the sequence to append.
...type:Shortcut.CharString
..include:seqan/bed_io.h
*/

inline void addSequenceName(BedStream & stream, CharString const & name)
{
    appendName(stream.sequenceNames, name, stream._sequenceNamesCache);
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/**
.Function.BedStream#readRecord
..class:Class.BedStream
..cat:BED I/O
..summary:Read a record from a @Class.BedStream@
..signature:int readRecord(record, bedStream)
..param.record:The @Class.BedRecord@ to read into.
...type:Class.BedRecord
..param.bedStream:The @Class.BedStream@ to read from.
...type:Class.BedStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/bed_io.h
*/

template <typename TSpec>
inline int readRecord(BedRecord<TSpec> & record,
                      BedStream & stream)
{
    int res = readRecord(record, *stream._reader, stream._context, Bed());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/**
.Function.BedStream#writeRecord
..class:Class.BedStream
..cat:BED I/O
..summary:Write a record to a @Class.BedStream@
..signature:int writeRecord(bedStream, record)
..param.bedStream:The @Class.BedStream@ to write to.
...type:Class.BedStream
..param.record:The @Class.BedRecord@ to write.
...type:Class.BedRecord
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/bed_io.h
*/

template <typename TSpec>
inline int writeRecord(BedStream & stream,
                       BedRecord<TSpec> const & record)
{
    int res = writeRecord(*stream._outStream, record, stream._context, Bed());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function flush()
// ----------------------------------------------------------------------------

/**
.Function.BedStream#flush
..class:Class.BedStream
..cat:BED I/O
..summary:Flush to a @Class.BedStream@
..signature:int flush(bedStream)
..param.bedStream:The @Class.BedStream@ to flush.
...type:Class.BedStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/bed_io.h
*/

inline int flush(BedStream & stream)
{
    if (stream._stream.get())
        stream._stream->flush();
    return 0;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/**
.Function.BedStream#close
..class:Class.BedStream
..cat:BED I/O
..summary:Closes a @Class.BedStream@
..signature:int close(bedStream)
..param.bedStream:The @Class.BedStream@ to close.
...type:Class.BedStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/bed_io.h
*/

inline int close(BedStream & stream)
{
    if (stream._stream.get())
        stream._stream->close();
    return 0;
}

// ----------------------------------------------------------------------------
// Function isGood()
// ----------------------------------------------------------------------------

/**
.Function.BedStream#isGood
..class:Class.BedStream
..cat:BED I/O
..summary:Query a @Class.BedStream@ for errors.
..signature:bool isGood(bedStream)
..param.bedStream:The @Class.BedStream@ to query.
...type:Class.BedStream
..returns:$true$ if stream is good, $false$ otherwise.
..include:seqan/bed_io.h
*/

inline bool isGood(BedStream const & stream)
{
    return stream._isGood;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

/**
.Function.BedStream#atEnd
..class:Class.BedStream
..cat:BED I/O
..summary:Query a @Class.BedStream@ for being at the end of the file.
..signature:bool atEnd(bedStream)
..param.bedStream:The @Class.BedStream@ to query.
...type:Class.BedStream
..returns:$true$ if stream is at the end, $false$ otherwise.
..include:seqan/bed_io.h
*/

inline bool atEnd(BedStream const & stream)
{
    return atEnd(*stream._reader);
}

inline bool atEnd(BedStream & stream)
{
    return atEnd(*stream._reader);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_BED_IO_BED_STREAM_H_
