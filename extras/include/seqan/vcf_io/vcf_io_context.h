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

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_IO_CONTEXT_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_IO_CONTEXT_H_

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
.Class.VcfIOContext
..cat:VCF I/O
..signature:class VcfIOContext
..summary:The I/O context to use for VCF I/O.
..description:
VcfIOContext objects store the names of (and provide a cache for) reference and sample names.
@Class.StringSet@ of @Shortcut.CharString@ are used for the name stores.
..include:seqan/vcf_io.h

.Memfunc.VcfIOContext#VcfIOContext
..class:Class.VcfIOContext
..signature:VcfIOContext()
..signature:VcfIOContext(sequenceNames, sampleNames)
..summary:Constructor.
..remarks:Default constructor or construction with references to sequence and sample names.

.Memvar.VcfIOContext#sequenceNames
..class:Class.VcfIOContext
..summary:Names of the reference sequences.

.Memvar.VcfIOContext#sequenceNamesCache
..class:Class.VcfIOContext
..summary:Name store cache for of the reference names.

.Memvar.VcfIOContext#sampleNames
..class:Class.VcfIOContext
..summary:Names of the samples.

.Memvar.VcfIOContext#sampleNamesCache
..class:Class.VcfIOContext
..summary:Name store cache for the sample names.
*/

class VcfIOContext
{
public:
    typedef StringSet<CharString> TNameStore_;
    typedef NameStoreCache<TNameStore_> TNameStoreCache_;

    // Pointer to the sequence names store.
    TNameStore_ * sequenceNames;
    // Cache for the sequence name lookup.
    NameStoreCache<TNameStore_> sequenceNamesCache;
    // Pointer to the sample name store.
    TNameStore_ * sampleNames;
    // Cache for the sample name lookup.
    NameStoreCache<TNameStore_> sampleNamesCache;

    // Default constructor.
    VcfIOContext() :
            sequenceNames(), sequenceNamesCache(*sequenceNames),
            sampleNames(), sampleNamesCache(*sampleNames)
    {}

    // Construct directly with references to stores.
    VcfIOContext(TNameStore_ & sequenceNames, TNameStore_ & sampleNames) :
            sequenceNames(&sequenceNames),
            sequenceNamesCache(sequenceNames),
            sampleNames(&sampleNames),
            sampleNamesCache(sequenceNames)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_IO_VCF_IO_CONTEXT_H_
