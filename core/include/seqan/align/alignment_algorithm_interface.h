// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Tags and interface documentation for the alignment algorithms.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_ALGORITHM_INTERFACE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_ALGORITHM_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function globalAlignment()
// ----------------------------------------------------------------------------

/**
.Function.globalAlignment
..summary:Computes the best global pairwise alignment.
..cat:Alignments
..signature:globalAlignment(align,          scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignment(gapsH, gapsV,   scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignment(frags, strings, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignment(alignmentGraph, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..param.align:
An @Class.Align@ object that stores the alignment.
The number of rows must be 2 and the sequences must have already been set.
$row(align, 0)$ is the horizontal one in the alignment matrix alignment, $row(align, 1)$ is the vertical one.
...type:Class.Align
..param.gapsH:Horizontal gapped sequenc in alignment matrix.
...type:Class.Gaps
..param.gapsV:Vertical gapped sequenc in alignment matrix.
...type:Class.Gaps
..param.frags:String of @Class.Fragment@ objects.  The sequence with id $0$ is the horizontal one, the sequence with id $1$ is the vertical one.
..param.alignmentGraph:
@Spec.Alignment Graph@ object to store the alignment in.
...type:Spec.Alignment Graph
...remarks:The underlying @Class.StringSet@ must be an @Spec.Owner|Owner StringSet@.
..param.strings:A @Class.StringSet@ containing two sequences.
...type:Class.StringSet
..param.scoringScheme:
The scoring scheme to use for the alignment.
Note that the user is responsible for ensuring that the scoring scheme is compatible with $algorithmTag$.
...type:Class.Score
..param.alignConfig:The @Class.AlignConfig@ to use for the alignment.
...type:Class.AlignConfig
..param.lowerDiag:Optional lower diagonal.
...type:nolink:$int$
..param.upperDiag:Optional upper diagonal.
...type:nolink:$int$
..param.algorithmTag:The Tag for picking the alignment algorithm.
...type:Tag.Pairwise Global Alignment Algorithms.tag.Gotoh
...type:Tag.Pairwise Global Alignment Algorithms.tag.NeedlemanWunsch
...type:Tag.Pairwise Global Alignment Algorithms.tag.Hirschberg
...type:Tag.Pairwise Global Alignment Algorithms.tag.MyersHirschberg
..returns:An integer with the alignment score, as given by the @Metafunction.Value@ metafunction of the @Class.Score@ type.
..remarks:
There exist multiple overloads for this function with four configuration dimensions.
..remarks:
First, you can select whether begin and end gaps are free in either sequence using $alignConfig$.
..remarks:
Second, you can select the type of the target storing the alignment.
This can be either an @Class.Align@ object, two @Class.Gaps@ objects, a @Spec.Alignment Graph@, or a string of @Class.Fragment@ objects.
@Class.Align@ objects provide an interface to tabular alignments with the restrction of all rows having the same type.
Using two @Class.Gaps@ objects has the advantage that you an align sequences with different types, for example @Shortcut.DnaString@ and @Shortcut.Dna5String@.
@Spec.Alignment Graph|Alignment Graphs@ provide a graph-based representation of segment-based collinear alignments.
Using @Class.Fragment@ strings is useful for collecting many pairwise alignments, for example in the construction of @Spec.Alignment Graph|Alignment Graphs@ for multiple-sequence alignments (MSA).
..remarks:
Third, you can optionally give a band for the alignment using $lowerDiag$ and $upperDiag$.
The center diagonal has index $0$, the $i$th diagonal below has index $-i$, the $i$th above has index $i$.
..remarks:
Fourth, you can select the algorithm to use with $algorithmTag$.
This can be one of @Tag.Pairwise Global Alignment Algorithms.value.NeedlemanWunsch@ and @Tag.Pairwise Global Alignment Algorithms.value.Gotoh@.
The Needleman-Wunsch algorithm supports scoring schemes with linear gap costs only while Gotoh's algorithm also allows affine gap costs.
..remarks:
The available alignment algorithms all have some restrictions.
Gotoh's algorithm can handle arbitrary substitution and affine gap scores.
Needleman-Wunsch is limited to linear gap scores.
The implementation of Hirschberg's algorithm is further limited that it does not support $alignConfig$ objects or banding.
The implementation of the Myers-Hirschberg algorithm further limits this to only support edit distance (as scores, matches are scored with 0, mismatches are scored with -1).
..remarks:
The examples below show some common use cases.
..example.text:Global alignment of two sequences using an @Class.Align@ object and the Needleman-Wunsch algorithm.
..example.code:
Dna5String seqH = "CGATT";
Dna5String seqV = "CGAAATT";

Align<Dna5String> align;
resize(rows(align), 2);
assignSource(row(align, 0), seqH);
assignSource(row(align, 0), seqV);
Score<int, Simple> scoringScheme(2, -1, -2);
AlignConfig<> alignConfig;

int result = globalAlignment(align, scoringScheme, alignConfig,
                             NeedlemanWunsch());
..example.text:Global banded alignment of two sequences using two @Class.Gaps@ objects and the Gotoh algorithm.
..example.code:
Dna5String seqH = "CGATT";
Gaps<Dna5String, ArrayGaps> gapsH(seqH);
DnaString seqV = "CGAAATT";
Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);

Score<int, Simple> scoringScheme(5, -3, -1, -5);
AlignConfig<> alignConfig;

int result = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, -2, 2);
..see:Function.localAlignment
..see:Function.globalAlignmentScore
..include:seqan/align.h
..wiki:Tutorial/Alignments
..cite:Needleman SB, Wunsch CD: A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol 1970, 48(3): 443-53.
..cite:Gotoh O: An improved algorithm for matching biological sequences. J Mol Biol 1982, 162(3):705-8
.
*/

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()
// ----------------------------------------------------------------------------

/**
.Function.globalAlignmentScore
..summary:Computes the best global pairwise alignment score.
..cat:Alignments
..signature:globalAlignmentScore(seqH, seqV, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignmentScore(strings,    scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignmentScore(seqH, seqV, {MyersBitVector | MyersHirschberg})
..signature:globalAlignmentScore(strings,    {MyersBitVector | MyersHirschberg})
..param.seqH:Horizontal gapped sequence in alignment matrix.
...type:Class.String
..param.seqV:Vertical gapped sequence in alignment matrix.
...type:Class.String
..param.strings:A @Class.StringSet@ containing two sequences.
...type:Class.StringSet
..param.scoringScheme:
The scoring scheme to use for the alignment.
Note that the user is responsible for ensuring that the scoring scheme is compatible with $algorithmTag$.
...type:Class.Score
..param.alignConfig:The @Class.AlignConfig@ to use for the alignment.
...type:Class.AlignConfig
..param.lowerDiag:Optional lower diagonal.
...type:nolink:$int$
..param.upperDiag:Optional upper diagonal.
...type:nolink:$int$
..param.algorithmTag:The Tag for picking the alignment algorithm.
...type:Tag.Pairwise Global Alignment Algorithms.tag.Gotoh
...type:Tag.Pairwise Global Alignment Algorithms.tag.NeedlemanWunsch
...type:Tag.Pairwise Global Alignment Algorithms.tag.Hirschberg
...type:Tag.Pairwise Global Alignment Algorithms.tag.MyersHirschberg
...type:Tag.Pairwise Global Alignment Algorithms.tag.MyersBitVector
..returns:An integer with the alignment score, as given by the @Metafunction.Value@ metafunction of the @Class.Score@ type.
..remarks:
This function does not perform the (linear time) traceback step after the (mostly quadratic time) dynamic programming step.
Note that Myers' bit-vector algorithm does not compute an alignment (only in the Myers-Hirschberg variant) but scores can be computed using $globalAlignmentScore$.
..remarks:
The same limitations to algorithms as in @Function.globalAlignment@ apply.
Furthermore, the $MyersBitVector$ and $MyersHirschberg$ variants can only be used without any other parameter.
..see:Function.globalAlignment
*/

// ----------------------------------------------------------------------------
// Function localAlignment()
// ----------------------------------------------------------------------------

/**
.Function.localAlignment
..summary:Computes the best pairwise local alignment using the Smith-Waterman algorithm.
..cat:Alignments
..signature:localAlignment(align,          scoringScheme, [lowerDiag, upperDiag])
..signature:localAlignment(gapsH, gapsV,   scoringScheme, [lowerDiag, upperDiag])
..signature:localAlignment(fragmentString, scoringScheme, [lowerDiag, upperDiag])
..param.align:
An @Class.Align@ object that stores the alignment.
The number of rows must be 2 and the sequences must have already been set.
$align[0]$ is the horizontal one in the alignment matrix alignment, $align[1]$ is the vertical one.
...type:Class.Align
..param.gapsH:Horizontal gapped sequenc in alignment matrix.
...type:Class.Gaps
..param.gapsV:Vertical gapped sequenc in alignment matrix.
...type:Class.Gaps
..param.fragmentString
String of @Class.Fragment@ objects.
The sequence with id $0$ is the horizontal one, the sequence with id $1$ is the vertical one.
..param.scoringScheme:
The scoring scheme to use for the alignment.
Note that the user is responsible for ensuring that the scoring scheme is compatible with $algorithmTag$.
...type:Class.Score
..param.lowerDiag:Optional lower diagonal.
...type:nolink:$int$
..param.upperDiag:Optional upper diagonal.
...type:nolink:$int$
..returns:An integer with the alignment score, as given by the @Metafunction.Value@ metafunction of the @Class.Score@ type.
..remarks:The Waterman-Eggert algorithm (local alignment with declumping) is available through the @Class.LocalAlignmentEnumerator@ class.
..remarks:
When using @Class.Gaps@ and @Class.Align@ objects, only parts (i.e. one infix) of each sequence will be aligned.
This will be presented to the user by setting the clipping begin and end position of the gaps (the rows in the case of @Class.Align@ objects).
When using @Class.Fragment@ strings, these parts of the sequences will not appear in any fragment.
..remarks:
There exist multiple overloads for this function with two configuration dimensions.
..remarks:
First, you can select the type of the target storing the alignment.
This can be either an @Class.Align@ object, two @Class.Gaps@ objects, or a string of @Class.Fragment@ objects.
@Class.Align@ objects provide an interface to tabular alignments with the restrction of all rows having the same type.
Using two @Class.Gaps@ objects has the advantage that you an align sequences with different types, for example @Shortcut.DnaString@ and @Shortcut.Dna5String@.
Using @Class.Fragment@ strings is useful for collecting many pairwise alignments, for example in the construction of @Spec.Alignment Graph|Alignment Graphs@ for multiple-sequence alignments (MSA).
..remarks:
Second, you can optionally give a band for the alignment using $lowerDiag$ and $upperDiag$.
The center diagonal has index $0$, the $i$th diagonal below has index $-i$, the $i$th above has index $i$.
..remarks:
The examples below show some common use cases.
..example.text:Local alignment of two sequences using an @Class.Align@ object.
..example.code:
Dna5String seqH = "CGATT";
Dna5String seqV = "CGAAATT";

Align<Dna5String> align;
resize(rows(align), 2);
assignSource(row(align, 0), seqH);
assignSource(row(align, 0), seqV);
Score<int, Simple> scoringScheme(2, -1, -2);

int result = localAlignment(align, scoringScheme);
..example.text:Local banded alignment of two sequences using two @Class.Gaps@ objects.
..example.code:
Dna5String seqH = "CGATT";
Gaps<Dna5String, ArrayGaps> gapsH(seqH);
DnaString seqV = "CGAAATT";
Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);

Score<int, Simple> scoringScheme(5, -3, -1, -5);

int result = localAlignment(gapsH, gapsV, scoringScheme, -2, 2);
..see:Function.globalAlignment
..see:Class.LocalAlignmentEnumerator
..include:seqan/align.h
..wiki:Tutorial/Alignments
..cite:Smith TF, Waterman, MS: Identification of Common Molecular Subsequences. J Mol Biol 1981, 147(1):195-7.
.
*/

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_ALGORITHM_INTERFACE_H_
