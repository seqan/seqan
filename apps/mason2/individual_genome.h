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
// Representation of an individual genome with variants.
//
// The genome is split into its contigs (aka chromosomes).  Each chromosome
// can then have a number of haplotypes.
// ==========================================================================

// TODO(holtgrew): Extend this with C+G biases in sequencing.
// TODO(holtgrew): Extend this with bisulphite models.

#ifndef APPS_MASON2_INDIVIDUAL_GENOME_H_
#define APPS_MASON2_INDIVIDUAL_GENOME_H_

#include <memory>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class SmallVariation
// ----------------------------------------------------------------------------

// Represent a small variation (SNP/insert/deletion) in a genome.

struct SmallVariation
{
    // The type of the variation represented by the record.
    enum VariationType
    {
        NONE,
        SNP,
        INSERTION,
        DELETION
    };

    // The idx of the contig the variation is on.  -1 if not set.
    int rId;
    // The position of the variation in the sequence.  -1 if not set.
    int pos;
    // The type of the variant.
    VariationType vType;
    // The length of the variant.  For SNPs, this is 1, otherwise the length of the indel.
    int len;
    // The sequence that is to be inserted here, if any.
    seqan::Char5String seq;

    SmallVariation() : rId(-1), pos(-1), vType(NONE), len(0)
    {}
};

// ----------------------------------------------------------------------------
// Class GenomeVariantOptions
// ----------------------------------------------------------------------------

// Options for creating variants.

struct GenomeVariantOptions
{
    // The number of haplotypes to simulate.
    unsigned numHaplotypes;
    // The per-base probability for having a SNP at a given base.
    double snpRate;
    // The per-base probability of creating an insert at a given base.
    double indelRate;
    // The smallest indel size to generate.
    unsigned indelRangeMax;
    // The larges indel size to generate.
    unsigned indelRangeMin;
    // Whether or not to insert Ns into variants.
    bool noN;

    GenomeVariantOptions() :
            numHaplotypes(1), snpRate(0.001), indelRate(0.001), indelRangeMin(1), indelRangeMax(6),
            noN(false)
    {}
};

// ----------------------------------------------------------------------------
// Class GenomeVariantManager
// ----------------------------------------------------------------------------

// Manages the genome variants.
//
// This class handles the broking of contigs and variants thereof.

class GenomeVariantManager
{
public:
    // Type for the variations for one haplotype of one chromosome.
    typedef seqan::String<SmallVariation> THaplotypeVariations;
    // Type for storing the haplotypes for one contig.
    typedef seqan::String<THaplotypeVariations> TContigHaploVariations;
    // Type for storing the variation information for the whole genome.
    typedef seqan::String<TContigHaploVariations> TGenomeHaploVariations;

    // The variations for the whole genome.
    TGenomeHaploVariations _genomeHaploVariations;
    // The number of contigs.
    unsigned _numContigs;

    GenomeVariantManager(unsigned numContigs) : _numContigs(numContigs)
    {}

    // Generate the haplotype variations for one contig with the given length.
    //
    // The variations will be stored internally.
    void generateContigHaploVariations(unsigned contigNo, unsigned len, unsigned numHaplotypes,
                                       GenomeVariantOptions const & options);

    // Get the haplotype variations for one contig.
    void getContigHaploVariations(TContigHaploVariations & result, unsigned contigNo);

    // Get the number of contigs.
    unsigned numContigs() const;

    // Get the number of haplotypes for the given contig.
    unsigned numHaplotypes(unsigned contigNo) const;

    // Apply the variations for the given contig and haplotype to the given Journaled String.
    void applyVariations(seqan::String<seqan::Dna5, seqan::Journaled<> > & haplotype,
                         unsinged contigNo, unsigned haplotypeNo) const;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef APPS_MASON2_INDIVIDUAL_GENOME_H_
