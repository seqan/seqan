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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Smart file for reading/writing files in Vcf format.
// ==========================================================================
// TODO(weese:) add Bcf I/O and integrate it

#ifndef SEQAN_VCF_IO_VCF_FILE_H_
#define SEQAN_VCF_IO_VCF_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

typedef SmartFile<Vcf, Input>   VcfFileIn;
typedef SmartFile<Vcf, Output>  VcfFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct SmartFileContext<SmartFile<Vcf, TDirection, TSpec>, TStorageSpec>
{
    typedef StringSet<CharString>                                   TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    typedef VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<SmartFile<Vcf, TDirection, TSpec> >
{
    typedef Vcf Type;
};

// ----------------------------------------------------------------------------
// Function readRecord(); VcfRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(VcfHeader & record, SmartFile<Vcf, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

template <typename TSpec>
inline void
readRecord(VcfRecord & record, SmartFile<Vcf, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); VcfRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(SmartFile<Vcf, Output, TSpec> & file, VcfHeader & record)
{
    write(file.iter, record, context(file), file.format);
}

template <typename TSpec>
inline void
writeRecord(SmartFile<Vcf, Output, TSpec> & file, VcfRecord & record)
{
    write(file.iter, record, context(file), file.format);
}

}  // namespace seqan

#endif // SEQAN_VCF_IO_VCF_FILE_H_
