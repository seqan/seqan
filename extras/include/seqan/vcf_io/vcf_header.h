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

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_VCF_HEADER_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_VCF_HEADER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfHeader
// ----------------------------------------------------------------------------

/**
.Class.VcfHeader
..cat:VCF I/O
..summary:Store VCF Header information.
..signature:class VcfHeader
..include:seqan/vcf_io.h

.Memfunc.VcfHeader#VcfHeader
..class:Class.VcfHeader
..signature:VcfHeader::VcfHeader()
..summary:Only default constructor.

.Memvar.VcfHeader#sequenceNames
..class:Class.VcfHeader
..summary:Names of the sequences (@Class.StringSet@<@Shortcut.CharString@>).

.Memvar.VcfHeader#sampleNames
..class:Class.VcfHeader
..summary:Names of the samples (@Class.StringSet@<@Shortcut.CharString@>).

.Memvar.VcfHeader#headerRecords
..class:Class.VcfHeader
..summary:The meta information records (@Class.String@ of @Class.VcfHeaderRecord@).
*/

class VcfHeader
{
public:
    typedef StringSet<CharString> TNameStore_;

    // The names of the sequences.
    TNameStore_ sequenceNames;
    // The names of the samples.
    TNameStore_ sampleNames;
    // Records for the meta information lines.
    String<VcfHeaderRecord> headerRecords;

    VcfHeader()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/**
.Function.VcfHeader#clear
..class:Class.VcfHeader
..summary:Clear a @Class.VcfHeader@.
..signature:void clear(header)
..param.header:@Class.VcfHeader@ to clear.
...type:Class.VcfHeader
..include:seqan/vcf_io.h
*/

inline void clear(VcfHeader & header)
{
    clear(header.sequenceNames);
    clear(header.sampleNames);
    clear(header.headerRecords);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_VCF_HEADER_H_
