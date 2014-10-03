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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Smart file for reading/writing files in ROI format.
// ==========================================================================

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_ROI_IO_ROI_STREAM_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_ROI_IO_ROI_STREAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Typedefs
// ============================================================================

typedef SmartFile<Roi, Input>   RoiFileIn;
typedef SmartFile<Roi, Output>  RoiFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct SmartFileContext<SmartFile<Roi, TDirection, TSpec>, TStorageSpec>
{
    typedef RoiIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<SmartFile<Roi, TDirection, TSpec> >
{
    typedef Roi Type;
};

// ----------------------------------------------------------------------------
// Function readRecord(); RoiRecord
// ----------------------------------------------------------------------------

// convient RoiFile variant
template <typename TSpec>
inline void
readRecord(RoiRecord & record, SmartFile<Roi, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(); RoiHeader
// ----------------------------------------------------------------------------

// convient RoiFile variant
template <typename TSpec>
inline void
readRecord(RoiHeader & record, SmartFile<Roi, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); RoiRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(SmartFile<Roi, Output, TSpec> & file, RoiRecord const & record)
{
    writeRecord(file.iter, record, file.format);
}

// ----------------------------------------------------------------------------
// Function write(); RoiRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(SmartFile<Roi, Output, TSpec> & file, RoiHeader const & record)
{
    writeRecord(file.iter, record, file.format);
}

}  // namespace seqan

#endif // SEQAN_EXTRAS_INCLUDE_SEQAN_ROI_IO_ROI_STREAM_H_

