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

#ifndef EXTRAS_INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNMENT_OPTIONS_H_
#define EXTRAS_INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNMENT_OPTIONS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConsensusAlignmentOptions
// ----------------------------------------------------------------------------

/*!
 * @class ConsensusAlignmentOptions
 * @headerfile <seqan/consensus.h>
 * @brief Configuration for @link consensusAlignment2 @endlink
 *
 * @signature struct ConsensusAlignmentOptions;
 *
 * @var unsigned ConsensusAlignmentOptions::INVALID;
 * @brief Static member variable with a marker for invalid contig ids.
 *
 * @var unsigned ConsensusAlignmentOptions::contigID;
 * @brief The id of the contig to compute the consensus for, defaults to <tt>INVALID</tt>.
 *
 * Set to <tt>INVALID</tt> for all.
 *
 * @var bool ConsensusAlignmentOptions::useContigID;
 * @brief Whether or not to use the value of <tt>contigID</tt>.
 *
 * When this variable is set to <tt>false</tt> then the value of <tt>contigID</tt> is ignored and treated
 * as if it was <tt>INVALID</tt>.
 *
 * @var bool ConsensusAlignmentOptions::usePositions;
 * @brief Whether or not to use positions of the in the @link FragmentStore::alignedReadStore @endlink.
 *
 * When set to <tt>false</tt>, then the @link FragmentStore::alignedReadStore @endlink will be cleared and filled
 * with new entries, one for each read.
 *
 * @var unsigned ConsensusAlignmentOptions::posDelta;
 * @brief Positions are considered with an environment of <tt>posDelta</tt>.
 *
 * When positions are not used then this value is ignored.
 *
 * See @link consensusAlignment2 @endlink for more details.
 *
 * @var bool ConsensusAlignmentOptions::runRealignment;
 * @brief Perform a realignment using standard parameters, depending on the sequence length.
 *
 * Defaults to <tt>true</tt>.
 */

struct ConsensusAlignmentOptions
{
    static const unsigned INVALID = (unsigned)-1;

    ConsensusAlignmentOptions() :
            contigID(INVALID), useContigID(true), usePositions(true), posDelta(30),
            runRealignment(true), verbosity(0), overlapMaxErrorRate(5), overlapMinLength(20),
            overlapMinCount(3), kMerSize(20), kMerMaxOcc(200)
    {}

    unsigned contigID;
    bool useContigID;
    bool usePositions;
    unsigned posDelta;  // TODO(holtgrew): Rename to overlapWindowSize
    bool runRealignment;

    int verbosity;

    int overlapMaxErrorRate;
    int overlapMinLength;
    int overlapMinCount;
    int kMerSize;
    int kMerMaxOcc;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNMENT_OPTIONS_H_
