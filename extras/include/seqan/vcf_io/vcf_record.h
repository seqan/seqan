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

// TODO(holtgrew): Get out more than just strings...

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_RECORD_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfRecord
// ----------------------------------------------------------------------------

/**
.Class.VcfRecord
..cat:VCF I/O
..summary:Information for one VCF record.
..signature:class VcfRecord
..description:
We store most information as strings and without structure since the VCF format's definition is quite loose.
We plan to provide classes for structured access to these strings later.
..remarks:
Although all positions in the VCF text format are 1-based, they are stored 0-based in the @Class.VcfRecord@.
Positions in string members are stored verbatim in the @Class.VcfRecord@'s members, e.g. 1-based.
..remarks:
Invalid qualities are stored as a float $NaN$ (not a number).
To test a float quality $q$ for being $NaN$, test for $q != q$.
Only $NaN$ has the property that $NaN != NaN$.
..include:seqan/vcf_io.h

.Memfunc.VcfRecord#VcfRecord
..class:Class.VcfRecord
..signature:VcfRecord::VcfRecord()
..summary:Default constructor.

.Memfunc.VcfRecord#MISSING_QUAL
..class:Class.VcfRecord
..signature:float VcfRecord::MISSING_QUAL()
..summary:Return IEEE $NaN$ float value.

.Memvar.VcfRecord#INVALID_REFID
..class:Class.VcfRecord
..summary:Static member as marker for invalid reference ($__int32$)

.Memvar.VcfRecord#INVALID_POS
..class:Class.VcfRecord
..summary:Static member as marker for invalid position ($__int32$)

.Memvar.VcfRecord#rID
..class:Class.VcfRecord
..summary:Static member as marker for invalid reference ($__int32$)

.Memvar.VcfRecord#beginPos
..class:Class.VcfRecord
..summary:Position of the VCF record ($__int32$).

.Memvar.VcfRecord#id
..class:Class.VcfRecord
..summary:Textual identifier of the variant (@Shortcut.CharString@).

.Memvar.VcfRecord#ref
..class:Class.VcfRecord
..summary:Bases in the reference (@Shortcut.CharString@).

.Memvar.VcfRecord#alt
..class:Class.VcfRecord
..summary:Alternative bases in the variants, comma-separated if multiple (@Shortcut.CharString@).

.Memvar.VcfRecord#qual
..class:Class.VcfRecord
..summary:Quality, $NaN$ if invalid ($float$).

.Memvar.VcfRecord#filter
..class:Class.VcfRecord
..summary:Value of FILTER field, empty if "." in VCF file (@Shortcut.CharString@).

.Memvar.VcfRecord#info
..class:Class.VcfRecord
..summary:Value of the INFO field, empty if "." in VCF file (@Shortcut.CharString@).

.Memvar.VcfRecord#format
..class:Class.VcfRecord
..summary:Value of the VCF FORMAT field, empty if "." in VCF file (@Shortcut.CharString@).

.Memvar.VcfRecord#genotypeInfos
..class:Class.VcfRecord
..summary:Genotype information, as in VCF file (@Class.StringSet@<@Shortcut.CharString@>).
*/

class VcfRecord
{
public:
    // Constant for invalid reference id.
    static const __int32 INVALID_REFID = -1;
    // Constant for invalid position.
    static const __int32 INVALID_POS = -1;

    // Numeric id of the reference sequence.
    __int32 rID;
    // Position on the reference.
    __int32 beginPos;
    // Textual identifier of the variant.
    CharString id;
    // Bases in the reference.
    CharString ref;
    // Bases in the alternatives, comma-separated.
    CharString alt;
    // Quality
    float qual;
    // Value of FILTER field.
    CharString filter;
    // Value of INFO field.
    CharString info;
    // Value of FORMAT field.
    CharString format;
    // The genotype infos.
    StringSet<CharString> genotypeInfos;

    // Default constructor.
    VcfRecord() : rID(INVALID_REFID), beginPos(INVALID_POS), qual(MISSING_QUAL())
    {}

    // Actually this is IEEE NaN.
    static float MISSING_QUAL()
    {
        union {
            __uint32 u;
            float f;
        } tmp;
        tmp.u = 0x7F800001;
        return tmp.f;
    }
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
.Function.VcfRecord#clear
..cat:VCF I/O
..class:Class.VcfRecord
..summary:Clear a @Class.VcfRecord@.
..signature:void clear(record)
..param.record:The @Class.VcfRecord@ to clear.
...type:Class.VcfRecord
..include:seqan/vcf_io.h
*/

inline void clear(VcfRecord & record)
{
    clear(record.id);
    clear(record.ref);
    clear(record.alt);
    clear(record.filter);
    clear(record.info);
    clear(record.format);
    clear(record.genotypeInfos);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_RECORD_H_
