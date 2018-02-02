// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Authors: Joerg Winkler <j.winkler@fu-berlin.de>
//          Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_HEADER_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_HEADER_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class rnaHeader
// ----------------------------------------------------------------------------

/*!
 * @class RnaHeader
 * @headerfile <seqan/rna_io.h>
 * @brief A container for labels of several RNA structures.
 *
 * @signature class RnaHeader;
 *
 * The container stores all kinds of data that can be obtained by reading RNA structure file headers.
 */
class RnaHeader
{
public:
    /*!
     * @var CharString RnaHeader::description
     * @brief Free text for file description.
     */
    CharString description;

    /*!
     * @var StringSet<CharString> RnaHeader::seqLabels
     * @brief List of sequence names.
     */
    StringSet<CharString> seqLabels;

    /*!
     * @var StringSet<CharString> RnaHeader::fixLabels
     * @brief List of fixed structure computation methods.
     */
    StringSet<CharString> fixLabels;

    /*!
     * @var StringSet<CharString> RnaHeader::bppLabels
     * @brief List of base pair probability matrix computation methods.
     */
    StringSet<CharString> bppLabels;

    /*!
     * @var StringSet<CharString> RnaHeader::typeLabels
     * @brief List of types of biological validated data.
     */
    StringSet<CharString> typeLabels;

    /*!
     * @fn RnaHeader::RnaHeader
     * @brief The constructor.
     * @signature RnaHeader::RnaHeader()
     */
    RnaHeader() : description("") {}
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

inline void clear(RnaHeader & header)
{
    clear(header.description);
    clear(header.seqLabels);
    clear(header.fixLabels);
    clear(header.bppLabels);
    clear(header.typeLabels);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_HEADER_H_
