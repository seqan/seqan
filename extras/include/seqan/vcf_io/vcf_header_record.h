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

// TODO(holtgrew): Parse more than just the key/value pair.

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_HEADER_RECORD_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_HEADER_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfHeaderRecord
// ----------------------------------------------------------------------------

/**
.Class.VcfHeaderRecord
..cat:VCF I/O
..signature:class VcfHeaderRecord
..summary:Store key/value pair for VCF header records.
..include:seqan/vcf_io.h

.Memfunc.VcfHeaderRecord#VcfHeaderRecord
..class:Class.VcfHeaderRecord
..summary:Constructor
..description:The default constructor and construction from key/value pair are provided.
..signature:VcfHeaderRecord::VcfHeaderRecord()
..signature:VcfHeaderRecord::VcfHeaderRecord(key, value)
..param.key:Key of the header record.
...type:Shortcut.CharString
..param.key:Key of the header record.
...type:Shortcut.CharString

.Memvar.VcfHeaderRecord#key
..class:Class.VcfHeaderRecord
..summary:Key of the header record (@Shortcut.CharString@).

.Memvar.VcfHeaderRecord#value
..class:Class.VcfHeaderRecord
..summary:Value of the header record (@Shortcut.CharString@).
*/

class VcfHeaderRecord
{
public:
    // Record's key.
    CharString key;
    // Record's value.
    CharString value;

    // Default constructor.
    VcfHeaderRecord()
    {}

    // Construct directly with key/value.
    VcfHeaderRecord(CharString const & key, CharString const & value) :
            key(key), value(value)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/**
.Function.VcfHeaderRecord#clear
..class:Class.VcfHeaderRecord
..summary:Clear a @Class.VcfHeaderRecord@.
..signature:void clear(record)
..param.record:@Class.VcfHeaderRecord@ to clear.
...type:Class.VcfHeaderRecord
..include:seqan/vcf_io.h
*/

inline void clear(VcfHeaderRecord & record)
{
    clear(record.key);
    clear(record.value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_HEADER_RECORD_H_
