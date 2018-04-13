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

#ifndef APPS_MASON2_VCF_MATERIALIZATION_H_
#define APPS_MASON2_VCF_MATERIALIZATION_H_

#include <stdexcept>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>

#include "genomic_variants.h"
#include "methylation_levels.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfMaterializer
// ----------------------------------------------------------------------------

// Allows the contig- and haplotype-wise construction of haplotypes stored in VCF files.

class VcfMaterializer
{
public:
    // The random number generator to use for methylation simulation, if any.
    TRng & rng;
    // Options for the methylation simulation.
    MethylationLevelSimulatorOptions const * methOptions;

    // The PositionMap is built for each contig to map between large variants, small variants, and original coordinate
    // system.
    PositionMap posMap;

    // ------------------------------------------------------------------------
    // Paths
    // ------------------------------------------------------------------------

    // Path to reference file.
    seqan::CharString fastaFileName;
    // Path to VCF file.
    seqan::CharString vcfFileName;
    // Path to methylation FASTA file.
    seqan::CharString methFastaFileName;

    // ------------------------------------------------------------------------
    // State for position in reference
    // ------------------------------------------------------------------------

    // The index of the current contig.
    int currRID;
    // The index of the current haplotype of the contig.
    int nextHaplotype;
    // The number of haplotypes (set after call to init()).
    int numHaplotypes;
    // The variants for the current contig.
    Variants contigVariants;
    // Current contig reference sequence.
    seqan::Dna5String contigSeq;
    // Current methylation levels.
    MethylationLevels currentLevels;

    // ------------------------------------------------------------------------
    // File Input
    // ------------------------------------------------------------------------

    // The FAI Index to load the reference sequence from.
    seqan::FaiIndex faiIndex;
    // The FAI Index to load the methylation sequences from.
    seqan::FaiIndex methFaiIndex;
    // The VCF stream to load from and the VCF heade.r
    seqan::VcfFileIn vcfFileIn;
    seqan::VcfHeader vcfHeader;
    // The current VCF record.  rID == INVALID_REFID if invalid, used for termination.
    seqan::VcfRecord vcfRecord;

    VcfMaterializer(TRng & rng) : rng(rng), currRID(-1), nextHaplotype(0), numHaplotypes(0)
    {}

    // If you give methFastaFileName, then you also have to set methOptions.
    //
    // The methylation simulation assumes that there is an methylation options object.
    VcfMaterializer(TRng & rng,
                    char const * fastaFileName,
                    char const * vcfFileName,
                    char const * methFastaFileName = "",
                    MethylationLevelSimulatorOptions const * methOptions = 0) :
            rng(rng), methOptions(methOptions), fastaFileName(fastaFileName), vcfFileName(vcfFileName),
            methFastaFileName(methFastaFileName), currRID(-1), nextHaplotype(0), numHaplotypes(0)
    {}

    // Call to open all files.
    //
    // Throws: MasonIOException
    void init();

    // Materialize next contig.
    //
    // Write sequence to seq, reference id to rID, haplotype to haplotype.  Returns true if the materialization could be
    // done and false if there are no more contigs to materialize.
    //
    // Call init() before calling materializeNext().
    //
    // Throws: MasonIOException
    bool materializeNext(seqan::Dna5String & seq,
                         std::vector<SmallVarInfo> & varInfos,
                         std::vector<std::pair<int, int> > & breakpoints,
                         int & rID, int & haplotype);   

    // Similar to the one above but loads methylation levels into levels.  Can only work if methFastFileName is not
    // empty.
    bool materializeNext(seqan::Dna5String & seq,
                         MethylationLevels & levels,
                         std::vector<SmallVarInfo> & varInfos,
                         std::vector<std::pair<int, int> > & breakpoints,
                         int & rID, int & haplotype);

private:

    bool _materializeNext(seqan::Dna5String & seq, MethylationLevels * levels,
                          std::vector<SmallVarInfo> & varInfos,
                          std::vector<std::pair<int, int> > & breakpoints,
                          int & rID, int & haplotype);

    // Load variants of next contig into variants.
    int _loadVariantsForContig(Variants & variants, int rID);

    // Append VCF record to variants.
    void _appendToVariants(Variants & variants, seqan::VcfRecord const & vcfRecord);

    // Append chunk of 6 BND records to variants.
    void _appendToVariantsBnd(Variants & variants, std::vector<seqan::VcfRecord> const & vcfRecords);

    // Load the levels for the contig with the given rID.
    void _loadLevels(int rID);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef APPS_MASON2_VCF_MATERIALIZATION_H_
