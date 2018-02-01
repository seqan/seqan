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
// Authors: Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// Class RnaIOContext, accessor functions.
// ==========================================================================

#ifndef INCLUDE_SEQAN_RNA_IO_RNA_IO_CONTEXT_H_
#define INCLUDE_SEQAN_RNA_IO_RNA_IO_CONTEXT_H_

#include <seqan/store.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class RnaIOContext
 * @headerfile <seqan/rna_io.h>
 * @brief File context for sharing information between @link RnaRecord @endlink and @link RnaHeader @endlink.
 *
 * @signature class RnaIOCOntext;
 *
 * As the labels and identifiers defined in the header are needed for parsing the record data,
 * this container stores all the shared information.
 */
class RnaIOContext
{
public:
    /*!
     * @var StringSet<CharString> RnaIOContext::seqLabels
     * @brief Descriptions for the sequences (value of S.. field).
     */
    StringSet<CharString> seqLabels;

    /*!
     * @var StringSet<CharString> RnaIOContext::fixLabels
     * @brief Descriptions for the fixed structure graphs (value of F.. field).
     */
    StringSet<CharString> fixLabels;

    /*!
     * @var StringSet<CharString> RnaIOContext::bppLabels
     * @brief Descriptions for the base pair probability structure graphs (value of M.. field).
     */
    StringSet<CharString> bppLabels;

    /*!
     * @var StringSet<CharString> RnaIOContext::seqIdent
     * @brief Identifiers for the sequences (starting with S).
     */
    StringSet<CharString> seqIdent;

    /*!
     * @var StringSet<CharString> RnaIOContext::fixIdent
     * @brief Identifiers for the fixed structure graphs (starting with F).
     */
    StringSet<CharString> fixIdent;

    /*!
     * @var StringSet<CharString> RnaIOContext::bppIdent
     * @brief Identifiers for the base pair probability structure graphs (starting with M).
     */
    StringSet<CharString> bppIdent;

    /*!
     * @var StringSet<CharString> RnaIOContext::typIdent
     * @brief Identifiers for the types of biological data.
     */
    StringSet<CharString> typIdent;
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

inline void clear(RnaIOContext & context)
{
    clear(context.seqLabels);
    clear(context.fixLabels);
    clear(context.bppLabels);
    clear(context.seqIdent);
    clear(context.fixIdent);
    clear(context.bppIdent);
    clear(context.typIdent);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_RNA_IO_RNA_IO_CONTEXT_H_
