// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================

// TODO(holtgrew): Parse more than just the key/value pair.

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_HEADER_RECORD_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_HEADER_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class rnaHeaderRecord
// ----------------------------------------------------------------------------

/*!
 * @class rnaHeaderRecord
 * @headerfile <seqan/rna_io.h>
 * @brief Store key/value pair for rna header records.
 *
 * @signature class rnaHeaderRecord;
 *
 * @var CharString rnaHeaderRecord::key;
 * @brief Key of the header record.
 *
 * @var CharString rnaHeaderRecord::value;
 * @brief Value of the header record.
 */

/*!
 * @fn rnaHeaderRecord::rnaHeaderRecord
 * @brief Constructor
 *
 * @signature rnaHeaderRecord::rnaHeaderRecord();
 * @signature rnaHeaderRecord::rnaHeaderRecord(key, value);
 *
 * @param[in] key   Key of the header record, @link CharString @endlink.
 * @param[in] value Key of the header record, @link CharString @endlink.
 */

/*!
 * @fn rnaHeaderRecord#clear
 *
 * @brief Clear a rnaHeaderRecord.
 * @signature void clear(record);
 *
 * @param[in,out] record The rnaHeaderRecord to clear.
 */

class RnaHeaderRecord
{
public:
    // Record's key.
    CharString key;
    // Record's value.
    CharString value;

    // Default constructor.
    RnaHeaderRecord()
    {}

    // Construct directly with key/value.
    RnaHeaderRecord(CharString const & key, CharString const & value) :
            key(key), value(value)
    {}
};

// ============================================================================
// Functions
// ============================================================================

inline void clear(RnaHeaderRecord & record)
{
    clear(record.key);
    clear(record.value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_HEADER_RECORD_H_
