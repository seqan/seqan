// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the header for the journal sequence file formats.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_HEADER_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_HEADER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class GdfRefInfo
// ----------------------------------------------------------------------------

struct GdfRefInfo
{
    unsigned    _refHash;
    CharString  _refId;
    CharString  _refFile;
    bool        _storeRef;

    GdfRefInfo() : _refHash(0), _storeRef(false)
    {}
};

// ----------------------------------------------------------------------------
// GdfFileInfo
// ----------------------------------------------------------------------------

template <typename TConfig = void>
struct GdfFileInfo
{
    unsigned    _minorFileId;
    unsigned    _majorFileId;
    unsigned    _blockSize;
    bool        _byteOrder; // true = little endian; false = big endian.
    bool        _snpCompression;  // true if 2 bit alphabet used.

    GdfFileInfo() : _minorFileId(GdfIO::FILE_VERSION_LITTLE),
                    _majorFileId(GdfIO::FILE_VERSION_BIG),
                    _blockSize(GdfFileBlockSize_<TConfig>::VALUE),
                    _byteOrder(SystemByteOrder::IS_LITTLE_ENDIAN()),
                    _snpCompression(false)
    {}
};

// ----------------------------------------------------------------------------
// Class GdfHeader
// ----------------------------------------------------------------------------

template <typename TConfig = void>
class GdfHeader
{
public:
    String<CharString>*   _nameStorePtr;  // The names of the individuals.
    GdfFileInfo<TConfig>  _fileInfos;
    GdfRefInfo            _refInfos;

    GdfHeader() : _nameStorePtr(0), _fileInfos(), _refInfos()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_HEADER_H_
