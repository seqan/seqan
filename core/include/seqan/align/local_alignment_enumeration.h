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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class LocalAlignmentEnumerator
// ----------------------------------------------------------------------------

template <typename TScore, typename TSpec>
class LocalAlignmentEnumerator;

/**
.Class.LocalAlignmentEnumerator
..cat:Alignments
..summary:Enumerate local alignments using the Waterman-Eggert algorithm.
..signature:LocalAlignmentEnumerator<TScore, TSpec>
..param.TScore:The @Class.Score@ type to use.
...type:Class.Score
..param.TSpec:Specialization tag.
..example.text:See the specializations for usage examples.
..cite:Waterman MS, Eggert M: A new algorithm for best subsequence alignments with application to tRNA-rRNA comparisons. J Mol Biol 1987, 197(4):723-728.
..include:seqan/align.h

.Spec.Unbanded LocalAlignmentEnumerator
..cat:Alignments
..general:Class.LocalAlignmentEnumerator
..summary:Unbanded enumeration of local alignments using the Waterman-Eggert algorithm.
..signature:LocalAlignmentEnumerator<TScore, Unbanded>
..example.text:Enumerate all alignments into an @Class.Align@ object.
..example.code:
SimpleScore scoringScheme(2, -1, -1, -2);
LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, 5);

Dna5String seqH = "CGAGAGAGACCGAGA";
Dna5String seqV = "TTCTGAGATCCGTTTTT";

Align<Dna5String> align;
resize(rows(align), 2);
assignSource(row(align), 0, seqH);
assignSource(row(align), 1, seqV);

int i = 0;
while (nextLocalAlignment(align, enumerator))
{
    std::cout << i << "-th alignment:\n";
    std::cout << align << "\n\n";
    std::cout << "score == " << getScore(enumerator) << "\n";
}
..include:seqan/align.h

.Memfunc.Unbanded LocalAlignmentEnumerator#LocalAlignmentEnumerator
..class:Spec.Unbanded LocalAlignmentEnumerator
..summary:Constructor
..signature:LocalAlignmentEnumerator(score, [cutoff])
..param.score:The scoring scheme to use for the alignments.
...type:Class.Score
..param.cutoff:Alignments with scores < $cutoff$ will be discarded.
...default:0
...type:nolink:$int$

.Spec.Banded LocalAlignmentEnumerator
..cat:Alignments
..general:Class.LocalAlignmentEnumerator
..signature:LocalAlignmentEnumerator<TScore, Banded>
..summary:Banded enumeration of local alignments using the Waterman-Eggert algorithm.
..example.text:Enumerate all alignments in the band between $-3$ and $0$ into an @Class.Align@ object.
..example.code:
SimpleScore scoringScheme(2, -1, -1, -2);
LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, 5, -3, 0);

Dna5String seqH = "CGAGAGAGACCGAGA";
Dna5String seqV = "TTCTGAGATCCGTTTTT";

Align<Dna5String> align;
resize(rows(align), 2);
assignSource(row(align), 0, seqH);
assignSource(row(align), 1, seqV);

int i = 0;
while (nextLocalAlignment(align, enumerator))
{
    std::cout << i << "-th alignment:\n";
    std::cout << align << "\n\n";
    std::cout << "score == " << getScore(enumerator) << "\n";
}
..include:seqan/align.h

.Memfunc.Banded LocalAlignmentEnumerator#LocalAlignmentEnumerator
..class:Spec.Banded LocalAlignmentEnumerator
..summary:Constructor
..signature:LocalAlignmentEnumerator(score, upperDiag, lowerDiag, [cutoff])
..param.score:The scoring scheme to use for the alignments.
...type:Class.Score
..param.upperDiag:Upper diagonal of the band.
...type:nolink:$int$
..param.lowerDiag:Lower diagonal of the band.
...type:nolink:$int$
..param.cutoff:Alignments with scores < $cutoff$ will be discarded.
...type:nolink:$int$
...default:0
*/

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getScore()
// ----------------------------------------------------------------------------

/**
.Function.LocalAlignmentEnumerator#getScore
..cat:Alignments
..summary:Compute next suboptimal local alignment.
..signature:getScore(enumerator)
..param.enumerator:The local alignment enumerator to use.
...type:Class.LocalAlignmentEnumerator
..returns:
The score of the previously computed alignment.
(Type: @Metafunction.Value@ of $enumerator$'s class.)
..include:seqan/align.h
 */

// ----------------------------------------------------------------------------
// Function nextLocalAlignment()
// ----------------------------------------------------------------------------

/**
.Function.nextLocalAlignment
..cat:Alignments
..summary:Compute next suboptimal local alignment.
..signature:nextLocalAlignment(align,        enumerator)
..signature:nextLocalAlignment(gapsH, gapsV, enumerator)
..param.align:The @Class.Align@ object to use for the alignment representation.
...type:Class.Align
..param.gapsH:The @Class.Gaps@ object to use for the horizontal sequence in the alignment matrix.
...type:Class.Gaps
..param.gapsV:The @Class.Gaps@ object to use for the vertical sequence in the alignment matrix.
...type:Class.Gaps
..param.enumerator:The @Class.LocalAlignmentEnumerator@ object to use.
...type:Class.LocalAlignmentEnumerator
..returns:$true$ if another suboptimal alignment above the given threshold was found, $false$ otherwise.
...type:nolink:$bool$
..include:seqan/align.h
*/

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_
