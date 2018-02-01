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

#ifndef APPS_MASON2_METHYLATION_LEVELS_H_
#define APPS_MASON2_METHYLATION_LEVELS_H_

#include <seqan/index.h>  // for Shape<>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>  // for the journal
#include <seqan/random.h> // for beta distribution

#include "mason_types.h"
#include "mason_options.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class MethylationLevels
// --------------------------------------------------------------------------

// Stores methylation levels separately for forward and reverse strand.

struct MethylationLevels
{
    // Forward and reverse levels, encoded as round(level / 0.125) + 33.
    seqan::CharString forward, reverse;

    void resize(unsigned len)
    {
        seqan::resize(forward, len, '!');
        seqan::resize(reverse, len, '!');
    }

    void clear()
    {
        seqan::clear(forward);
        seqan::clear(reverse);
    }

    // Translate character in forward/reverse to level (0..80).
    inline char charToLevel(char c) const
    {
        if (c < '>')  // '>' cannot be used as value
            return c - 33;
        else
            return c - 34;
    }

    // Translate level (0..80) to character in forward/reverse.
    inline char levelToChar(char c) const
    {
        if (c + '!' < '>')
            return c + 33;
        else
            return c + 34;
    }

    // Returns methylation level for forward strand at position i.
    inline float levelF(unsigned i) const
    {
        return (charToLevel(forward[i]) * 0.0125);
    }

    // Sets methylation level for forward strand at position i.
    inline void setLevelF(unsigned i, float level)
    {
        SEQAN_ASSERT_GEQ(level, 0.0);
        SEQAN_ASSERT_LEQ(level, 1.0);
        // std::cerr << "forward[i] = " << levelToChar(round(level / 0.0125)) << " ~ " << (level / 0.0125) << " ~ " << level << "\n";
        char c = levelToChar(static_cast<char>(round(level / 0.0125)));
        SEQAN_ASSERT_NEQ((int)c, (int)'>');
        forward[i] = std::max(forward[i], c);
    }

    // Returns methylation level for reverse strand at position i.
    inline float levelR(unsigned i) const
    {
        return (charToLevel(reverse[i]) * 0.0125);
    }

    // Sets methylation level for reverse strand at position i if level is larger than current.
    inline void setLevelR(unsigned i, float level)
    {
        SEQAN_ASSERT_GEQ(level, 0.0);
        SEQAN_ASSERT_LEQ(level, 1.0);
        // std::cerr << "reverse[" << i << "] = " << levelToChar(round(level / 0.0125)) << " ~ " << (level / 0.0125) << " ~ " << level << "\n";
        char c = levelToChar(static_cast<char>(round(level / 0.0125)));
        SEQAN_ASSERT_NEQ((int)c, (int)'>');
        reverse[i] = std::max(reverse[i], c);
    }
};

inline
void swap(MethylationLevels & lhs, MethylationLevels & rhs)
{
    swap(lhs.forward, rhs.forward);
    swap(lhs.reverse, rhs.reverse);
}

// --------------------------------------------------------------------------
// Class MethylationLevelSimulator
// --------------------------------------------------------------------------

// Simulate methylation levels for a Dna sequence/contig on forward and reverse strand.

class MethylationLevelSimulator
{
public:
    // Options for the mu/sigma values.
    MethylationLevelSimulatorOptions const & options;

    // Random number generator to use.
    TRng & rng;

    // Beta probability density functions for level generation.
    seqan::BetaDistribution<double> pdfCG, pdfCHG, pdfCHH; // TODO::(smehringer) change to beta distribution!

    MethylationLevelSimulator(TRng & rng, MethylationLevelSimulatorOptions const & options) :
            options(options), rng(rng),
            pdfCG(seqan::cvtBetaDistParam(options.methMuCG,options.methSigmaCG)),
            pdfCHG(seqan::cvtBetaDistParam(options.methMuCHG, options.methSigmaCHG)),
            pdfCHH(seqan::cvtBetaDistParam(options.methMuCHH, options.methSigmaCHH))
    {}

    // Simulate methylation levels for the sequence in contig.  The results are stored in levels.
    void run(MethylationLevels & levels, seqan::Dna5String const & contig)
    {
        levels.resize(length(contig));

        typedef seqan::Iterator<seqan::Dna5String const>::Type TContigIter;
        TContigIter it = begin(contig, seqan::Standard());
        TContigIter itEnd = end(contig, seqan::Standard()) - 3;

        // We will go over the contig with hashes to search for patterns efficiently.
        seqan::Shape<seqan::Dna5> shape2, shape3;
        if (levels.forward[0] != '!') SEQAN_ASSERT_EQ_MSG(contig[0], 'C', "pos = %d", 0);
        if (levels.reverse[0] != '!') SEQAN_ASSERT_EQ_MSG(contig[0], 'G', "pos = %d", 0);
        if (length(contig) >= 2u)
        {
            resize(shape2, 2);
            hash(shape2, it);
            handleTwoMer(levels, 0, value(shape2));
            if (levels.forward[1] != '!') SEQAN_ASSERT_EQ_MSG(contig[1], 'C', "pos = %d", 1);
            if (levels.reverse[1] != '!') SEQAN_ASSERT_EQ_MSG(contig[1], 'G', "pos = %d", 1);
        }
        if (length(contig) >= 3u)
        {
            resize(shape3, 3);
            hash(shape3, it);
            handleThreeMer(levels, 0, value(shape3));
            if (levels.forward[2] != '!') SEQAN_ASSERT_EQ_MSG(contig[2], 'C', "pos = %d", 2);
            if (levels.reverse[2] != '!') SEQAN_ASSERT_EQ_MSG(contig[2], 'G', "pos = %d", 2);
        }
        ++it;
        unsigned pos = 1;
        for (; (pos + 3 < length(contig)) && (it != itEnd); ++it, ++pos)
        {
            hashNext(shape2, it);
            hashNext(shape3, it);
            handleTwoMer(levels, pos, value(shape2));
            handleThreeMer(levels, pos, value(shape3));
            if (levels.forward[pos] != '!') SEQAN_ASSERT_EQ_MSG(contig[pos], 'C', "pos = %u", pos);
            if (levels.reverse[pos] != '!') SEQAN_ASSERT_EQ_MSG(contig[pos], 'G', "pos = %u", pos);
        }
        if (pos + 1 < length(contig))
        {
            hashNext(shape2, it++);
            handleTwoMer(levels, pos++, value(shape2));
        }
    }

    // Handle 3mer, forward case.
    void handleThreeMer(MethylationLevels & levels, unsigned pos, unsigned hashValue)
    {
        // seqan::Dna5String dbg, dbg2;
        // unhash(dbg, hashValue, 3);
        // dbg2 = dbg;
        // reverse(dbg2);
        switch (hashValue)
        {
            case 27:    // CAG
            case 42:    // CTG
                // std::cerr << "CHG fw    \t" << dbg << "\t" << dbg2 << "\t" << pos << "\n";
                levels.setLevelF(pos, pdfCHG(rng));
                levels.setLevelR(pos + 2, pdfCHG(rng));
                break;
            case 32:  // CCG
                levels.setLevelF(pos, pdfCHG(rng));
                levels.setLevelR(pos + 2, pdfCG(rng)); 
                break;
            case 25:    // CAA
            case 26:    // CAC
            case 28:    // CAT
            case 30:    // CCA
            case 31:    // CCC
            case 33:    // CCT
            case 40:    // CTA
            case 41:    // CTC
            case 43:    // CTT
                // std::cerr << "CHH       \t" << dbg << "\t" << dbg2 << "\t" << pos << "\n";
                levels.setLevelF(pos, pdfCHH(rng));
                break;
            case 2:     // AAG
            case 12:    // AGG
            case 17:    // ATG
            case 52:    // GAG
            case 62:    // GGG
            case 67:    // GTG
            case 77:    // TAG
            case 87:    // TGG
            case 92:    // TTG
                levels.setLevelR(pos + 2, pdfCHH(rng));
                break;

            default:
                // nop
                break;
        }
    }

    // Handle 2mer.
    void handleTwoMer(MethylationLevels & levels, unsigned pos, unsigned hashValue)
    {
        if (hashValue == 7)  // CpG forward (symmetric, also reverse)
        {
            // seqan::Dna5String dbg;
            // unhash(dbg, hashValue, 2);
            // std::cerr << "CpG     \t" << dbg << "\n";
            levels.setLevelF(pos, pdfCG(rng));
            levels.setLevelR(pos + 1, pdfCG(rng));
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_fixVariationLevels()
// ----------------------------------------------------------------------------

// Fix variation levels on contig given the points (second == true -> SNP, second == false -> breakpoint).

void fixVariationLevels(MethylationLevels & levels,
                        TRng & rng,
                        seqan::Dna5String const & contig,
                        seqan::String<std::pair<int, bool> > const & varPoints,
                        MethylationLevelSimulatorOptions const & options);

#endif  // #ifndef APPS_MASON2_METHYLATION_LEVELS_H_
