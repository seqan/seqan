// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_READER_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_READER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Implementation of reading BAM files.

class BamReader_ :
    public XamReader_
{
public:
    // The BGZF Stream to read from.
    Stream<Bgzf> _stream;
    // Flag indicating whether there was an error or not.
    bool _isGood;

    BamReader_() :
        XamReader_(), _isGood(true)
    {}

    BamReader_(CharString const & filename);

    // XamReader_ interface.

    virtual int open(CharString const & filename);
    virtual bool isGood();
    virtual bool atEnd();
    virtual int readHeader(BamHeader & header, BamIOContext<StringSet<CharString> > & context);
    virtual int readRecord(BamAlignmentRecord & record, BamIOContext<StringSet<CharString> > & context);
    virtual int close();
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Member Function BamReader_::BamReader_()
// ----------------------------------------------------------------------------

BamReader_::BamReader_(CharString const & filename) :
    XamReader_(filename), _isGood(true)
{
    this->open(_filename);
}

// ----------------------------------------------------------------------------
// Member Function BamReader_::open()
// ----------------------------------------------------------------------------

int BamReader_::open(CharString const & filename)
{
    if (!seqan::open(this->_stream, toCString(filename), "r"))
    {
        _isGood = false;
        return 1;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Member Function BamReader_::isGood()
// ----------------------------------------------------------------------------

bool BamReader_::isGood()
{
    return this->_isGood;
}

// ----------------------------------------------------------------------------
// Member Function BamReader_::atEnd()
// ----------------------------------------------------------------------------

bool BamReader_::atEnd()
{
    return seqan::atEnd(this->_stream);
}

// ----------------------------------------------------------------------------
// Member Function BamReader_::readHeader()
// ----------------------------------------------------------------------------

int BamReader_::readHeader(BamHeader & header, BamIOContext<StringSet<CharString> > & context)
{
    return seqan::readRecord(header, context, this->_stream, Bam());
}

// ----------------------------------------------------------------------------
// Member Function BamReader_::readRecord()
// ----------------------------------------------------------------------------

int BamReader_::readRecord(BamAlignmentRecord & record, BamIOContext<StringSet<CharString> > & context)
{
    return seqan::readRecord(record, context, this->_stream, Bam());
}

// ----------------------------------------------------------------------------
// Member Function BamReader_::close()
// ----------------------------------------------------------------------------

int BamReader_::close()
{
    seqan::close(this->_stream);
    return 0;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_READER_H_
