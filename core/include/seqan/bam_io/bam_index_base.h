// ==========================================================================
//                              bam_index_base.h
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

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BASE_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class BamIndex
 * @headerfile <seqan/bam_io.h>
 *
 * @brief Access to BAM indices.
 *
 * @signature template <typename TSpec>
 *            class BamIndex;
 *
 * @section Remarks
 *
 * This is an abstract class, don't use it itself but its specializations.
 */

/**
.Class.BamIndex
..cat:BAM I/O
..summary:Access to BAM Indices.
..signature:BamIndex<TSpec>
..param.TSpec:Tag to specialize index.
..remarks:This is an abstract class, don't use it itself but its specializations.
..include:seqan/bam_io.h
*/

template <typename TSpec>
class BamIndex;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function BamIndex#jumpToRegion
// ----------------------------------------------------------------------------

/*!
 * @fn BamIndex#jumpToRegion
 * @brief Seek in BAM BGZF stream using an index.
 *
 * You provide a region <tt>[pos, posEnd)</tt> on the reference <tt>refID</tt> that you want to jump to and the function
 * jumps to the first alignment in this region, if any.
 *
 * @signature bool jumpToRegion(stream, hasAlignments, bamIOContext, refID, pos, posEnd, index);
 *
 * @param[in,out] stream        The @link BgzfStream @endlink to jump with.
 * @param[out]    hasAlignments A <tt>bool</tt> that is set true if the region <tt>[pos, posEnd)</tt> has any
 *                              alignments.
 * @param[in,out] bamIOContext  The @link BamIOContext @endlink to use for jumping.
 * @param[in]     refID         The reference id to jump to (<tt>__int32</tt>).
 * @param[in]     pos           The begin of the region to jump to.
 * @param[in]     posEnd        The end of the region to jump to.
 * @param[in]     index         The @link BamIndex @endlink to use for the jumping.
 *
 * @return bool true if seeking was successful, false if not.
 *
 * @section Remarks
 *
 * This function fails if <tt>refID</tt>/<tt>pos</tt> are invalid.
 */

// ----------------------------------------------------------------------------
// Function jumpToOrphans
// ----------------------------------------------------------------------------

/*!
 * @fn BamIndex#jumpToOrphans
 * @brief Seek to orphans block in BAM BGZF stream using an index.
 *
 * @signature bool jumpToOrphans(stream, hasAlignments, bamIOContext, index);
 *
 * @param[in,out] stream         The @link BgzfStream @endlink object to jump with.
 * @param[out]    hasAlignments  A <tt>bool</tt> that is set to true if there are any orphans.
 * @param[in,out] bamIOContext   The @link BamIOContext @endlink to use for the state.
 * @param[in]     index          The index to use for jumping.
 */

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BASE_H_
