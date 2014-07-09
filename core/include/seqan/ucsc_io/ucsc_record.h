
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
// Author: Jochen Singer<jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_UCSC_RECORD_H_
#define CORE_INCLUDE_SEQAN_UCSC_RECORD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

class UcscRecord;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class UcscRecord
{
public:

    CharString      transName;
    CharString      contigName;
    __int64         cdsBegin;
    __int64         cdsEnd;
    String<__int64> exonBegin;
    String<__int64> exonEnd;
    CharString      proteinName;

    __uint64        annotationBeginPos;
    __uint64        annotationEndPos;

    enum {KNOWN_GENE, KNOWN_ISOFORMS} format;
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

/*!
 * @fn BamAlignmentRecord#clear
 * @brief Clear BamAlignmentRecord.
 *
 * @signature void clear(record);
 *
 * @param record The BamAlignmentRecord to clear.
 *
 * Clears all strings and resets it to default initialization state.
 */

///.Function.clear.param.object.type:Class.BamAlignmentRecord
///.Function.clear.class:Class.BamAlignmentRecord

inline void clear(UcscRecord & record)
{
    clear(record.transName);
    clear(record.contigName);
    record.cdsBegin = -1;
    record.cdsEnd = -1;
    clear(record.exonBegin);
    clear(record.exonEnd);
    clear(record.proteinName);

    record.annotationBeginPos = -1;
    record.annotationEndPos = -1;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
