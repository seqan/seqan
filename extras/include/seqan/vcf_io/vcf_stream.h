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

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_STREAM_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_STREAM_H_

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
.Class.VcfStream
..cat:VCF I/O
..summary:High-level VCF I/O class.
..signature:class VcfStream
..example:The following example demonstrates reading a VCF file and printing the variant locations.
..example.file:demos/vcf_io/vcf_stream_read.cpp
..include:seqan/vcf_io.h

.Memfunc.VcfStream#VcfStream
..class:Class.VcfStream
..summary:Constructor.
..signature:VcfStream::VcfStream()
..signature:VcfStream::VcfStream(fileName, mode=READ)
..param.fileName:The path to the file to open.
...type:nolink:$char const *$
..param.mode:The open mode.
...type:Enum.VcfStream\colon\colonMode
...default:$VcfStream\colon\colonMode::READ$
..see:Enum.VcfStream\colon\colonMode

.Enum.VcfStream\colon\colonMode
..cat:VCF I/O
..summary:Open mode for the @Class.VcfStream@ class.
..value.INVALID:Invalid open mode.
..value.READ:Open in read mode.
..value.WRITE:Open in write mode.
..include:seqan/vcf_io.h
*/

class VcfStream
{
public:
    typedef RecordReader<std::istream, SinglePass<> > TReader_;

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
    bool _headerWritten;

    VcfHeader header;
    VcfIOContext _context;

    VcfStream() : _outStream(), _inStream(), _mode(INVALID), _error(0), _isGood(true), _headerWritten(false),
                  _context(header.sequenceNames, header.sampleNames)
    {}

    VcfStream(char const * filename, Mode mode = READ) :
            _outStream(), _inStream(), _filename(filename), _mode(mode), _error(0), _isGood(true),
            _headerWritten(false), _context(header.sequenceNames, header.sampleNames)
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
        _headerWritten = false;

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

            int res = read(header, *_reader, _context, Vcf());
            if (res != 0)
            {
                _error = res;
                _isGood = false;
            }
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
.Function.VcfStream#open
..class:Class.VcfStream
..cat:VCF I/O
..summary:Open a @Class.VcfStream@.
..signature:bool open(vcfStream, fileName, mode)
..param.vcfStream:The @Class.VcfStream@ to open.
...type:Class.VcfStream
..param.fileName:The path to the file to open.
...type:nolink:$char const *$
..param.mode:The open mode.
...type:Enum.VcfStream\colon\colonMode
...default:$VcfStream\colon\colonMode::READ$
..returns:$true$ on success, $false$ on failure.
..see:Function.VcfStream#isGood
..see:Enum.VcfStream\colon\colonMode
..include:seqan/vcf_io.h
*/

inline bool open(VcfStream & stream, char const * filename, VcfStream::Mode mode = VcfStream::READ)
{
    return stream._open(filename, mode);
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/**
.Function.VcfStream#readRecord
..class:Class.VcfStream
..cat:VCF I/O
..summary:Read a record from a @Class.VcfStream@
..signature:int readRecord(record, vcfStream)
..param.record:The @Class.VcfRecord@ to read into.
...type:Class.VcfRecord
..param.vcfStream:The @Class.VcfStream@ to read from.
...type:Class.VcfStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/vcf_io.h
*/

inline int readRecord(VcfRecord & record,
                      VcfStream & stream)
{
    int res = readRecord(record, *stream._reader, stream._context, Vcf());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/**
.Function.VcfStream#writeRecord
..class:Class.VcfStream
..cat:VCF I/O
..summary:Write a record to a @Class.VcfStream@
..signature:int writeRecord(vcfStream, record)
..param.vcfStream:The @Class.VcfStream@ to write to.
...type:Class.VcfStream
..param.record:The @Class.VcfRecord@ to write.
...type:Class.VcfRecord
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/vcf_io.h
*/

inline int writeRecord(VcfStream & stream,
                       VcfRecord const & record)
{
    if (!stream._headerWritten)
    {
        int res = write(*stream._outStream, stream.header, stream._context, Vcf());
        if (res != 0)
            stream._isGood = false;
        stream._headerWritten = true;
    }

    int res = writeRecord(*stream._outStream, record, stream._context, Vcf());
    if (res != 0)
        stream._isGood = false;
    return res;
}

// ----------------------------------------------------------------------------
// Function flush()
// ----------------------------------------------------------------------------

/**
.Function.VcfStream#flush
..class:Class.VcfStream
..cat:VCF I/O
..summary:Flush to a @Class.VcfStream@
..signature:int flush(vcfStream)
..param.vcfStream:The @Class.VcfStream@ to flush.
...type:Class.VcfStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/vcf_io.h
*/

inline int flush(VcfStream & stream)
{
    if (stream._stream.get())
        stream._stream->flush();
    return 0;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/**
.Function.VcfStream#close
..class:Class.VcfStream
..cat:VCF I/O
..summary:Closes a @Class.VcfStream@
..signature:int close(vcfStream)
..param.vcfStream:The @Class.VcfStream@ to close.
...type:Class.VcfStream
..returns:$0$ on success, non-$0$ on failure.
..include:seqan/vcf_io.h
*/

inline int close(VcfStream & stream)
{
    // Close only when not stdout/stdin.
    if (stream._stream.get())
        stream._stream->close();
    return 0;
}

// ----------------------------------------------------------------------------
// Function isGood()
// ----------------------------------------------------------------------------

/**
.Function.VcfStream#isGood
..class:Class.VcfStream
..cat:VCF I/O
..summary:Query a @Class.VcfStream@ for errors.
..signature:bool isGood(vcfStream)
..param.vcfStream:The @Class.VcfStream@ to query.
...type:Class.VcfStream
..returns:$true$ if stream is good, $false$ otherwise.
..include:seqan/vcf_io.h
*/

inline bool isGood(VcfStream const & stream)
{
    return stream._isGood;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

/**
.Function.VcfStream#atEnd
..class:Class.VcfStream
..cat:VCF I/O
..summary:Query a @Class.VcfStream@ for being at the end of the file.
..signature:bool atEnd(vcfStream)
..param.vcfStream:The @Class.VcfStream@ to query.
...type:Class.VcfStream
..returns:$true$ if stream is at the end, $false$ otherwise.
..include:seqan/vcf_io.h
*/

inline bool atEnd(VcfStream const & stream)
{
    return atEnd(*stream._reader);
}

inline bool atEnd(VcfStream & stream)
{
    return atEnd(*stream._reader);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_STREAM_H_
