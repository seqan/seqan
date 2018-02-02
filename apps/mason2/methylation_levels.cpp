// ==========================================================================
//                         Mason - A Read Simulator
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "methylation_levels.h"

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_fixVariationLevels()
// ----------------------------------------------------------------------------

void fixVariationLevels(MethylationLevels & levels,
                        TRng & rng,
                        seqan::Dna5String const & contig,
                        seqan::String<std::pair<int, bool> > const & varPoints,
                        MethylationLevelSimulatorOptions const & options)
{
    MethylationLevelSimulator methSim(rng, options);
    seqan::Shape<seqan::Dna5> shape2, shape3;
    resize(shape2, 2);
    resize(shape3, 3);

    for (unsigned i = 0; i < length(varPoints); ++i)
    {
        int pos = varPoints[i].first;
        if (varPoints[i].second)  // is SNP
        {
            if (pos > 2)
            {
                levels.forward[pos - 2] = levels.reverse[pos - 2] = '!';
                methSim.handleTwoMer(levels, pos - 2, hash(shape2, iter(contig, pos - 2, seqan::Standard())));
                methSim.handleThreeMer(levels, pos - 2, hash(shape3, iter(contig, pos - 2, seqan::Standard())));
            }
            if (pos > 1)
            {
                levels.forward[pos - 1] = levels.reverse[pos - 1] = '!';
                methSim.handleTwoMer(levels, pos - 1, hash(shape2, iter(contig, pos - 1, seqan::Standard())));
            }
            levels.forward[pos] = levels.reverse[pos] = '!';
            if (pos + 1 < (int)length(contig))
            {
                levels.forward[pos + 1] = levels.reverse[pos + 1] = '!';
                methSim.handleTwoMer(levels, pos, hash(shape2, iter(contig, pos, seqan::Standard())));
            }
            if (pos + 2 < (int)length(contig))
            {
                levels.forward[pos + 2] = levels.reverse[pos + 2] = '!';
                methSim.handleTwoMer(levels, pos + 1, hash(shape2, iter(contig, pos + 1, seqan::Standard())));
                methSim.handleThreeMer(levels, pos, hash(shape3, iter(contig, pos, seqan::Standard())));
            }
        }
        else  // is no SNP but breakpoint
        {
            // TODO(holtgrew): Double-check for correctness, might recompute too much around breakpoints.
            if (pos > 2)
            {
                levels.forward[pos - 2] = levels.reverse[pos - 2] = '!';
                methSim.handleTwoMer(levels, pos - 2, hash(shape2, iter(contig, pos - 2, seqan::Standard())));
                methSim.handleThreeMer(levels, pos - 2, hash(shape3, iter(contig, pos - 2, seqan::Standard())));
            }
            if (pos > 1)
            {
                levels.forward[pos - 1] = levels.reverse[pos - 1] = '!';
                methSim.handleTwoMer(levels, pos - 1, hash(shape2, iter(contig, pos - 1, seqan::Standard())));
            }
            levels.forward[pos] = levels.reverse[pos] = '!';
            if (pos + 1 < (int)length(contig))
            {
                methSim.handleTwoMer(levels, pos, hash(shape2, iter(contig, pos, seqan::Standard())));
                levels.forward[pos + 1] = levels.reverse[pos + 1] = '!';
            }
            if (pos + 2 < (int)length(contig))
            {
                methSim.handleTwoMer(levels, pos + 1, hash(shape2, iter(contig, pos + 1, seqan::Standard())));
                methSim.handleThreeMer(levels, pos, hash(shape3, iter(contig, pos, seqan::Standard())));
            }
        }
    }
}
