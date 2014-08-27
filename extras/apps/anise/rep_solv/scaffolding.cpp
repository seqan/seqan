// ==========================================================================
//                                   ANISE
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

#include "scaffolding.h"

#include <seqan/basic.h>

#include "asm/overlapper.h"

#include "rep_solv/contig_graph.h"
#include "rep_solv/options.h"

#include "scaffolder/gpm.h"
#include "scaffolder/gpm_options.h"
#include "scaffolder/mate_link.h"
#include "scaffolder/overlap_resolution.h"
#include "scaffolder/scaffolding_result.h"

namespace rep_solv {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ScaffoldBuilder
// ----------------------------------------------------------------------------

class ScaffoldBuilder
{
public:
    ScaffoldBuilder(seqan::StringSet<seqan::Dna5String> & out,
                    ContigGraph const & cg,
                    seqan::StringSet<seqan::Dna5String> const & contigs,
                    Options const & options) :
            out(out), cg(cg), contigs(contigs), options(options), ovlOptions(buildOvlOptions()),
            overlapper(ovlOptions), overlapExists(mkOverlapExists()), computeOverlap(mkComputeOverlap())
    {}

    void run();

private:

    // Create overlapper options from options.
    assembler::OverlapperOptions buildOvlOptions() const;
    // Compute overlap, result.errors == Overlap::INVALID incase of non-existent overlap.
    assembler::Overlap doComputeOverlap(unsigned seq0, unsigned seq1, unsigned ovlLen, unsigned bandwidth) const;
    // Create function for overlap existence checking.
    scaffolder::TOverlapExistsFunc mkOverlapExists() const;
    // Create function for overlap computation.
    scaffolder::TComputeOverlapFunc mkComputeOverlap() const;

    // Input / Output

    seqan::StringSet<seqan::Dna5String> & out;
    ContigGraph const & cg;
    seqan::StringSet<seqan::Dna5String> const & contigs;
    Options const & options;

    // State / Members

    // Overlapper and overlapper options thereof.
    assembler::OverlapperOptions ovlOptions;
    assembler::Overlapper overlapper;

    // Contig verification functions for scaffolder.
    scaffolder::TOverlapExistsFunc overlapExists;
    scaffolder::TComputeOverlapFunc computeOverlap;
};

assembler::Overlap ScaffoldBuilder::doComputeOverlap(
        unsigned seq0, unsigned seq1, unsigned ovlLen, unsigned bandwidth) const
{
    // Subsequently ignored overlap and fragment string.
    assembler::Overlap ovl;
    assembler::TFragments frags;
    // Compute diagonal from bandwidth.
    int diag = (int)length(contigs[seq0]) - ovlLen;
    // Compute overlap.
    overlapper.computeOverlap(ovl, frags, contigs[seq0], contigs[seq1], seq0, seq1, diag, bandwidth);
    return ovl;
}

assembler::OverlapperOptions ScaffoldBuilder::buildOvlOptions() const
{
    assembler::OverlapperOptions result;

    result.overlapErrorRate = options.maxContigOverlapErrorRate;
    result.overlapMinLength = options.contigOverlapMinLen;
    result.logging = (options.verbosity >= 3);

    return result;
}

scaffolder::TOverlapExistsFunc ScaffoldBuilder::mkOverlapExists() const
{
    return [&](unsigned seq0, unsigned seq1, unsigned ovlLen, unsigned bandwidth) {
        return (doComputeOverlap(seq0, seq1, ovlLen, bandwidth).errors !=
                assembler::Overlap::INVALID);
    };
}

// Create function for overlap computation.
scaffolder::TComputeOverlapFunc ScaffoldBuilder::mkComputeOverlap() const
{
    return [&](unsigned seq0, unsigned seq1, unsigned ovlLen, unsigned bandwidth) {
        return doComputeOverlap(seq0, seq1, ovlLen, bandwidth);
    };
}

void ScaffoldBuilder::run()
{
    // Links have already been transitively reduced and bundled in the caller (including filtering against conflicts and
    // so on).  We just need to call the greedy path merging here and then apply the result and write out to out.
    std::vector<scaffolder::ContigEdgeLabel> contigInfos;  // collect contig infos (lengths)
    for (auto const & node : cg.node)
        contigInfos.push_back(scaffolder::ContigEdgeLabel(cg.contig[node].length));
    std::vector<scaffolder::MateLink> mateLinks;  // collect mate links
    for (lemon::SmartGraph::EdgeIt edge(cg.graph); edge != lemon::INVALID; ++edge)
    {
        auto const & label = cg.link[edge];
        if (label.leftID == EdgeLabel::INVALID || label.rightID == EdgeLabel::INVALID ||
            label.leftID == EdgeLabel::SOURCE  || label.rightID == EdgeLabel::SOURCE ||
            label.leftID == EdgeLabel::TARGET  || label.rightID == EdgeLabel::TARGET)
            continue;  // skip s-/t-links.
        auto dist = label.uniqueLinks ? label.uniqueLinkDistance : label.overlapDistance;
        mateLinks.push_back(scaffolder::MateLink(
                label.leftID, label.rightID,
                scaffolder::MateEdgeLabel(dist.mean, dist.sd, label.duplicateLinks, label.uniqueLinks)));
        if (options.verbosity >= 3)
            std::cerr << "MATE LINK FOR SCAFFOLDING\t" << mateLinks.back() << "\n";
    }

    // Actually call greedy path merging.
    scaffolder::PathMergingOptions pmOptions;
    pmOptions.verbosity = options.verbosity;
    scaffolder::ScaffoldingResult scaffoldingResult;
    greedyPathMerging(scaffoldingResult, mateLinks, contigInfos, overlapExists, computeOverlap, pmOptions);

    // Compute resulting scaffolds.
    resolveOverlaps(out, scaffoldingResult, contigs);
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function buildScaffold()
// ----------------------------------------------------------------------------

void buildScaffold(seqan::StringSet<seqan::Dna5String> & out,
                   ContigGraph const & cg,
                   seqan::StringSet<seqan::Dna5String> const & contigs,
                   Options const & options)
{
    ScaffoldBuilder helper(out, cg, contigs, options);
    helper.run();
}

}  // namespace rep_solv
