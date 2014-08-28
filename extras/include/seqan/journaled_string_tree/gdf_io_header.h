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

template <unsigned BitsPerValue>
struct GdfCompressionMode;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class GdfFileConfiguration
// ----------------------------------------------------------------------------

template <typename TSnpValue>
struct GdfFileConfiguration
{
    bool isLittleEndian;
    unsigned int refHash;
    unsigned int blockSize;
    GdfIO::CoverageCompression coverageCompression;
    GdfIO::CompressionMode compressionMode;

    template <typename TCoverageSize>
    GdfFileConfiguration(TCoverageSize coverageSize) :
        isLittleEndian(true),
        refHash(0),
        blockSize(100000),
        compressionMode(GdfCompressionMode<BitsPerValue<TSnpValue>::VALUE>::VALUE)
    {
        if (coverageSize <= MaxValue<__uint8>::VALUE)
            coverageCompression = GdfIO::COVERAGE_COMPRESSION_1_BYTE_PER_VALUE;
        else if (coverageSize <= MaxValue<__uint16>::VALUE)
            coverageCompression = GdfIO::COVERAGE_COMPRESSION_2_BYTE_PER_VALUE;
        else if (coverageSize <= MaxValue<__uint32>::VALUE)
            coverageCompression = GdfIO::COVERAGE_COMPRESSION_4_BYTE_PER_VALUE;
        else
            coverageCompression = GdfIO::COVERAGE_COMPRESSION_8_BYTE_PER_VALUE;
    }

};

// ----------------------------------------------------------------------------
// Class GdfHeader
// ----------------------------------------------------------------------------

class GdfHeader
{
public:

    GdfIO::SaveReferenceMode   referenceMode;
    CharString             referenceFilename;
    CharString             referenceId;
    String<CharString>     nameStore;  // The names of the individuals.

    GdfHeader() : referenceMode(GdfIO::SAVE_REFERENCE_MODE_ENABLED)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GdfCompressionMode
// ----------------------------------------------------------------------------

template <unsigned BitsPerValue>
struct GdfCompressionMode
{
    static const typename GdfIO::CompressionMode VALUE = GdfIO::COMPRESSION_MODE_NO_SNP_COMPRESSION;
};

template <>
struct GdfCompressionMode<2>
{
    static const typename GdfIO::CompressionMode VALUE = GdfIO::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION;
};

// ============================================================================
// Functions
// ============================================================================

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_HEADER_H_
