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

#include "genomic_variants.h"

std::ostream & operator<<(std::ostream & out, SnpRecord const & record)
{
    out << "SnpRecord(" << record.haplotype << ", " << record.rId << ", " << record.pos
        << ", " << record.to << ")";
    return out;
}

std::ostream & operator<<(std::ostream & out, SmallIndelRecord const & record)
{
    out << "SnpRecord(" << record.haplotype << ", " << record.rId << ", " << record.pos
        << ", " << record.size << ", " << record.seq << ")";
    return out;
}

int StructuralVariantRecord::endPosition() const
{
    if (pos == -1)
        return std::numeric_limits<int>::max();

    switch (kind)
    {
        case INDEL:
            if (size > 0)
                return pos;
            else
                return pos + size;
        case INVERSION:
            return pos + size;
        case TRANSLOCATION:
            return targetPos;
        case DUPLICATION:
            return targetPos;
        default:
            return -1;
    }
}

bool StructuralVariantRecord::nearBreakend(int query) const
{
    if (pos == -1)
        return false;  // invalid/sentinel has no breakends

    switch (kind)
    {
        case INDEL:
            if (size > 0)
                return (query == pos || query == pos + 1);
            else
                return (query == pos || query == pos + 1 ||
                        query == pos - size || query == pos - size + 1);
        case INVERSION:
            return (query == pos || query == pos + 1 ||
                    query == pos + size || query == pos + size + 1);
        case TRANSLOCATION:
            return (query == pos || query == pos + 1 ||
                    query == targetPos || query == targetPos + 1);
        case DUPLICATION:
            return (query == pos || query == pos + 1 ||
                    query == targetPos || query == targetPos + 1);
        default:
            return false;
    }
}

std::ostream & operator<<(std::ostream & out, StructuralVariantRecord const & record)
{
    char const * kind;
    switch (record.kind)
    {
        case StructuralVariantRecord::TRANSLOCATION:
            kind = "TRANSLOCATION";
            break;
        case StructuralVariantRecord::INDEL:
            kind = "INDEL";
            break;
        case StructuralVariantRecord::INVERSION:
            kind = "INVERSION";
            break;
        case StructuralVariantRecord::DUPLICATION:
            kind = "DUPLICATION";
            break;
        default:
            kind = "INVALID";
            break;
    }
    out << "StructuralVariantRecord(kind=" << kind << ", haplotype=" << record.haplotype
        << ", rId=" << record.rId << ", pos=" << record.pos << ", size=" << record.size
        << ", targetRId=" << record.targetRId << ", targetPos=" << record.targetPos
        << ", seq=\"" << record.seq << "\")";
    return out;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_runImpl()
// ----------------------------------------------------------------------------

int VariantMaterializer::_runImpl(
        seqan::Dna5String * resultSeq,
        PositionMap * posMap,
        MethylationLevels * resultLvls,
        std::vector<SmallVarInfo> & varInfos,
        std::vector<std::pair<int, int> > & breakpoints,
        seqan::Dna5String const * ref,
        MethylationLevels const * refLvls,
        int haplotypeId)
{
    breakpoints.clear();
    clear(*resultSeq);
    if (resultLvls)
        resultLvls->clear();

    // Apply small variants.  We get a sequence with the small variants and a journal of the difference to contig.
    seqan::Dna5String seqSmallVariants;
    TJournalEntries journal;
    MethylationLevels levelsSmallVariants;  // only used if revLevels != 0
    MethylationLevels * smallLvlsPtr = refLvls ? &levelsSmallVariants : 0;
    std::vector<SmallVarInfo> smallVarInfos;
    if (_materializeSmallVariants(seqSmallVariants, journal, smallLvlsPtr, smallVarInfos, *ref, *variants,
                                  refLvls, haplotypeId) != 0)
        return 1;

    // Build position map for large variant -> small variant and small variant <-> reference position mapping.
    posMap->reinit(journal);  // build mapping from small variant to reference positions

    // Apply structural variants and build the interval tree of posMap
    if (_materializeLargeVariants(*resultSeq, resultLvls, varInfos, breakpoints, *posMap, journal, seqSmallVariants,
                                  smallVarInfos, *variants, smallLvlsPtr, haplotypeId) != 0)
        return 1;

    // Sort resulting variant infos.
    std::sort(varInfos.begin(), varInfos.end());

    // std::sort(breakpoints.begin(), breakpoints.end());
    // std::cout << "SV Breakpoints\n";
    // for (unsigned i = 0; i < breakpoints.size(); ++i)
    //     std::cout << "  " << breakpoints[i] << "\n";

    // Copy out SV breakpoints.
    posMap->svBreakpoints.insert(breakpoints.begin(), breakpoints.end());

    return 0;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_materializeSmallVariants()
// ----------------------------------------------------------------------------

int VariantMaterializer::_materializeSmallVariants(
        seqan::Dna5String & seq,
        TJournalEntries & journal,
        MethylationLevels * levelsSmallVariants,
        std::vector<SmallVarInfo> & smallVarInfos,
        seqan::Dna5String const & contig,
        Variants const & variants,
        MethylationLevels const * levels,
        int hId)
{
    if (methSimOptions)
    {
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levelsSmallVariants != 0));
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levels != 0));
    }

    // Clear journal and output methylation levels.
    reinit(journal, length(contig));
    if (levelsSmallVariants)
        levelsSmallVariants->clear();
    // Store variation points with a flag whether it is a SNP (true) or a breakpoint (false).
    seqan::String<std::pair<int, bool> > varPoints;

    // Fors this, we have to iterate in parallel over SNP and small indel records.
    //
    // Current index in snp/small indel array.
    unsigned snpsIdx = 0;
    unsigned smallIndelIdx = 0;
    // Current SNP record, default to sentinel.
    SnpRecord snpRecord;
    snpRecord.rId = std::numeric_limits<int>::max();
    if (snpsIdx < length(variants.snps))
        snpRecord = variants.snps[snpsIdx++];
    // Current small indel record, default to sentinel.
    SmallIndelRecord smallIndelRecord;
    smallIndelRecord.rId = std::numeric_limits<int>::max();
    if (smallIndelIdx < length(variants.smallIndels))
        smallIndelRecord = variants.smallIndels[smallIndelIdx++];
    // Track last position from contig appended to seq so far.
    int lastPos = 0;
    if (verbosity >= 3)
        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

    // TODO(holtgrew): Extract contig building into their own functions.
    if (verbosity >= 2)
        std::cerr << "building output\n";
    while (snpRecord.rId != std::numeric_limits<int>::max() || smallIndelRecord.rId != std::numeric_limits<int>::max())
    {
        // TODO(holtgrew): Extract SNP and small indel handling into their own functions.
        if (snpRecord.getPos() < smallIndelRecord.getPos())  // process SNP records
        {
            if (snpRecord.haplotype == hId)  // Ignore all but the current contig.
            {
                if (verbosity >= 3)
                    std::cerr << "append(seq, infix(contig, " << lastPos << ", " << snpRecord.pos << ") " << __LINE__ << "\n";
                // Append interim sequence and methylation levels->
                append(seq, infix(contig, lastPos, snpRecord.pos));
                if (methSimOptions && methSimOptions->simulateMethylationLevels)
                {
                    append(levelsSmallVariants->forward, infix(levels->forward, lastPos, snpRecord.pos + 1));
                    append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, snpRecord.pos + 1));
                    appendValue(varPoints, std::make_pair((int)length(seq), true));      // variation points before/after SNP
                    appendValue(varPoints, std::make_pair((int)length(seq) + 1, true));
                }

                SEQAN_ASSERT_GEQ(snpRecord.pos, lastPos);
                if (verbosity >= 3)
                    std::cerr << "appendValue(seq, " << snpRecord.to << "')\n";
                appendValue(seq, snpRecord.to);
                lastPos = snpRecord.pos + 1;
                if (verbosity >= 3)
                    std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

                // Register SNP as small variant info.
                smallVarInfos.push_back(SmallVarInfo(SmallVarInfo::SNP, length(seq) - 1, 1));
            }

            if (snpsIdx >= length(variants.snps))
                snpRecord.rId = std::numeric_limits<int>::max();
            else
                snpRecord = variants.snps[snpsIdx++];
        }
        else
        {
            if (smallIndelRecord.haplotype == hId)  // Ignore all but the current contig.
            {
                if (smallIndelRecord.size > 0)
                {
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") "
                                  << __LINE__ << "\n";

                    // Simulate methylation levels for insertion.
                    MethylationLevels lvls;
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        MethylationLevelSimulator methSim(*rng, *methSimOptions);
                        methSim.run(lvls, smallIndelRecord.seq);
                    }

                    // Append interim sequence and methylation levels->
                    append(seq, infix(contig, lastPos, smallIndelRecord.pos));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, smallIndelRecord.pos));
                        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, smallIndelRecord.pos));
                        appendValue(varPoints, std::make_pair((int)length(seq), false));  // variation point before insertion
                    }

                    SEQAN_ASSERT_GEQ(smallIndelRecord.pos, lastPos);
                    if (verbosity >= 3)
                        std::cerr << "append(seq, \"" << smallIndelRecord.seq << "\") " << __LINE__ << "\n";
                    // Register insertion as small variant info.
                    for (unsigned i = 0; i < length(smallIndelRecord.seq); ++i)
                        smallVarInfos.push_back(SmallVarInfo(SmallVarInfo::INS, length(seq) + i, 1));
                    // Append novel sequence and methylation levels->
                    append(seq, smallIndelRecord.seq);
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        append(levelsSmallVariants->forward, lvls.forward);
                        append(levelsSmallVariants->reverse, lvls.reverse);
                        appendValue(varPoints, std::make_pair((int)length(seq), false));  // variation point after insertion
                    }
                    lastPos = smallIndelRecord.pos;
                    recordInsertion(journal, hostToVirtualPosition(journal, smallIndelRecord.pos),
                                    0, smallIndelRecord.size);
                    if (verbosity >= 3)
                        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";
                }
                else  // deletion
                {
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << smallIndelRecord.pos << ") " << __LINE__ << "\n";
                    // Append interim sequence and methylation levels->
                    append(seq, infix(contig, lastPos, smallIndelRecord.pos));  // interim chars
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair((int)length(seq), false));  // variation point at deletion
                        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, smallIndelRecord.pos));
                        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, smallIndelRecord.pos));
                    }

                    lastPos = smallIndelRecord.pos - smallIndelRecord.size;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));
                    recordErase(journal,
                                hostToVirtualPosition(journal, smallIndelRecord.pos),
                                hostToVirtualPosition(journal, smallIndelRecord.pos - smallIndelRecord.size));
                    if (verbosity >= 3)
                        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

                    // Register deletion as small variant info.
                    smallVarInfos.push_back(SmallVarInfo(SmallVarInfo::DEL, length(seq), -smallIndelRecord.size));
                }
            }

            if (smallIndelIdx >= length(variants.smallIndels))
                smallIndelRecord.rId = std::numeric_limits<int>::max();
            else
                smallIndelRecord = variants.smallIndels[smallIndelIdx++];
        }
    }
    // Insert remaining characters.
    if (verbosity >= 3)
        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ")\n";
    append(seq, infix(contig, lastPos, length(contig)));

    if (methSimOptions && methSimOptions->simulateMethylationLevels)
    {
        append(levelsSmallVariants->forward, infix(levels->forward, lastPos, length(contig)));
        append(levelsSmallVariants->reverse, infix(levels->reverse, lastPos, length(contig)));

        SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants->forward));
        SEQAN_ASSERT_EQ(length(seq), length(levelsSmallVariants->reverse));

        fixVariationLevels(*levelsSmallVariants, *rng, seq, varPoints, *methSimOptions);
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function VariantMaterializer::_materializeLargeVariants()
// ----------------------------------------------------------------------------

int VariantMaterializer::_materializeLargeVariants(
        seqan::Dna5String & seq,
        MethylationLevels * levelsLargeVariants,
        std::vector<SmallVarInfo> & varInfos,
        std::vector<std::pair<int, int> > & breakpoints,
        PositionMap & positionMap,
        TJournalEntries const & journal,
        seqan::Dna5String const & contig,
        std::vector<SmallVarInfo> const & smallVarInfos,
        Variants const & variants,
        MethylationLevels const * levels,
        int hId)
{
    if (methSimOptions)
    {
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levelsLargeVariants != 0));
        SEQAN_ASSERT_EQ(methSimOptions->simulateMethylationLevels, (levels != 0));
    }

    // We will record all intervals for the positionMap.svIntervalTree in this String.
    seqan::String<GenomicInterval> intervals;

    // Clear output methylation levels->
    if (levelsLargeVariants)
        levelsLargeVariants->clear();
    // Store variation points.  We reuse the fixVariationLevels() function from small indel/snp simulation and thus
    // have to store a bool that is always set to false.
    seqan::String<std::pair<int, bool> > varPoints;

    // Track last position from contig appended to seq so far.
    int lastPos = 0;
    if (verbosity >= 3)
        std::cerr << __LINE__ << "\tlastPos == " << lastPos << "\n";

    // Pointer to the current small variant to write out translated to varInfo.
    std::vector<SmallVarInfo>::const_iterator itSmallVar = smallVarInfos.begin();

    // Number of bytes written out so far/current position in variant.
    unsigned currentPos = 0;

    for (unsigned i = 0; i < length(variants.svRecords); ++i)
    {
        if (variants.svRecords[i].haplotype != hId)  // Ignore all but the current contig.
            continue;
        // We obtain a copy of the current SV record since we translate its positions below.
        StructuralVariantRecord svRecord = variants.svRecords[i];

        // Translate positions and lengths of SV record.
        if (verbosity >= 2)
            std::cerr << "  Translating SvRecord\n  " << svRecord << '\n';
        svRecord.pos = hostToVirtualPosition(journal, svRecord.pos);
        SEQAN_ASSERT_LT(svRecord.pos, (int)length(contig));
        // We do not need to adjust the sizes for insertions.
        if (svRecord.kind != StructuralVariantRecord::INDEL || svRecord.size < 0)
            svRecord.size = hostToVirtualPosition(journal, svRecord.pos + svRecord.size) -
                    hostToVirtualPosition(journal, svRecord.pos);
        if (svRecord.targetPos != -1)
            svRecord.targetPos = hostToVirtualPosition(journal, svRecord.targetPos);
        if (verbosity >= 2)
            std::cerr << "  => " << svRecord << '\n';

        // Copy out small variant infos for interim chars.
        for (; itSmallVar != smallVarInfos.end() && itSmallVar->pos < svRecord.pos; ++itSmallVar)
        {
            int offset = (int)currentPos - lastPos;
            varInfos.push_back(*itSmallVar);
            varInfos.back().pos += offset;
        }

        // Copy from contig to seq with SVs.
        if (verbosity >= 3)
            std::cerr << "lastPos == " << lastPos << "\n";
        append(seq, infix(contig, lastPos, svRecord.pos));  // interim chars
        if (methSimOptions && methSimOptions->simulateMethylationLevels)
        {
            append(levelsLargeVariants->forward, infix(levels->forward, lastPos, svRecord.pos));
            append(levelsLargeVariants->reverse, infix(levels->reverse, lastPos, svRecord.pos));
            appendValue(varPoints, std::make_pair((int)length(seq), false));
        }
        if (currentPos != length(seq))
            appendValue(intervals, GenomicInterval(currentPos, length(seq), lastPos, svRecord.pos,
                                                   '+', GenomicInterval::NORMAL));
        currentPos = length(seq);
        if (verbosity >= 3)
            std::cerr << "append(seq, infix(contig, " << lastPos << ", " << svRecord.pos << ") " << __LINE__ << " (interim)\n";
        switch (svRecord.kind)
        {
            case StructuralVariantRecord::INDEL:
                {
                    if (svRecord.size > 0)  // insertion
                    {
                        SEQAN_ASSERT_EQ((int)length(svRecord.seq), svRecord.size);

                        // Simulate methylation levels for insertion.
                        MethylationLevels lvls;
                        if (methSimOptions && methSimOptions->simulateMethylationLevels)
                        {
                            MethylationLevelSimulator methSim(*rng, *methSimOptions);
                            methSim.run(lvls, svRecord.seq);
                        }

                        // Append novel sequence and methylation levels.
                        append(seq, svRecord.seq);
                        if (methSimOptions && methSimOptions->simulateMethylationLevels)
                        {
                            append(levelsLargeVariants->forward, lvls.forward);
                            append(levelsLargeVariants->reverse, lvls.reverse);
                            appendValue(varPoints, std::make_pair((int)length(seq), false));  // variation point after insertion
                        }
                        if (currentPos != length(seq))
                            appendValue(intervals, GenomicInterval(currentPos, length(seq), -1, -1,
                                                                   '+', GenomicInterval::INSERTED));
                        if (verbosity >= 3)
                            std::cerr << "append(seq, svRecord.seq (length == " << length(svRecord.seq) << ") " << __LINE__ << " (insertion)\n";
                        lastPos = svRecord.pos;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        // Copy out breakpoints.
                        breakpoints.push_back(std::make_pair(currentPos, variants.posToIdx(Variants::SV, i)));
                        breakpoints.push_back(std::make_pair((int)length(seq), variants.posToIdx(Variants::SV, i)));

                        currentPos = length(seq);
                    }
                    else  // deletion
                    {
                        lastPos = svRecord.pos - svRecord.size;
                        SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                        // Copy out breakpoint.
                        breakpoints.push_back(std::make_pair(currentPos, variants.posToIdx(Variants::SV, i)));
                    }
                }
                break;
            case StructuralVariantRecord::INVERSION:
                {
                    unsigned oldLen = length(seq);
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair((int)length(seq), false));  // variation point at deletion
                        append(levelsLargeVariants->forward, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                        reverse(infix(levelsLargeVariants->forward, oldLen, length(levelsLargeVariants->forward)));
                        append(levelsLargeVariants->reverse, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        reverse(infix(levelsLargeVariants->reverse, oldLen, length(levelsLargeVariants->reverse)));
                    }
                    if (currentPos != length(seq))
                        appendValue(intervals, GenomicInterval(currentPos, length(seq), svRecord.pos, svRecord.pos + svRecord.size,
                                                               '-', GenomicInterval::INVERTED));

                    // Copy out small variant infos for inversion.
                    for (; itSmallVar != smallVarInfos.end() && itSmallVar->pos < svRecord.pos + svRecord.size; ++itSmallVar)
                    {
                        varInfos.push_back(*itSmallVar);
                        varInfos.back().pos = currentPos + svRecord.size - (varInfos.back().pos - lastPos);
                    }

                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (inversion)\n";
                    reverseComplement(infix(seq, oldLen, length(seq)));
                    lastPos = svRecord.pos + svRecord.size;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out breakpoints.
                    breakpoints.push_back(std::make_pair(currentPos, variants.posToIdx(Variants::SV, i)));
                    breakpoints.push_back(std::make_pair((int)length(seq), variants.posToIdx(Variants::SV, i)));

                    currentPos = length(seq);
                }
                break;
            case StructuralVariantRecord::TRANSLOCATION:
                {
                    SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                    append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair((int)length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                    }
                    if (currentPos != length(seq))
                        appendValue(intervals, GenomicInterval(currentPos, length(seq), svRecord.pos + svRecord.size, svRecord.targetPos,
                                                               '+', GenomicInterval::NORMAL));
                    unsigned tmpCurrentPos = length(seq);
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair((int)length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    if (tmpCurrentPos != length(seq))
                        appendValue(intervals, GenomicInterval(tmpCurrentPos, length(seq), svRecord.pos, svRecord.pos + svRecord.size,
                                                               '+', GenomicInterval::NORMAL));
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << " (translocation)\n"
                                  << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                    lastPos = svRecord.targetPos;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Copy out small variant infos for translocation, shift left to right and righ to left but keep
                    // center intact.
                    for (; itSmallVar != smallVarInfos.end() && itSmallVar->pos < svRecord.pos; ++itSmallVar)
                    {
                        int offset = (int)currentPos - lastPos;
                        varInfos.push_back(*itSmallVar);
                        varInfos.back().pos += offset;

                        int bpLeft = svRecord.pos + svRecord.size;
                        int bpRight = svRecord.targetPos;
                        if (itSmallVar->pos < bpLeft)
                            varInfos.back().pos -= (svRecord.targetPos - svRecord.pos);
                        else if (itSmallVar->pos >= bpRight)
                            varInfos.back().pos += (svRecord.targetPos - svRecord.pos);
                    }

                    // Copy out breakpoints.
                    breakpoints.push_back(std::make_pair(currentPos, variants.posToIdx(Variants::SV, i)));
                    breakpoints.push_back(std::make_pair(currentPos + svRecord.targetPos - svRecord.pos - svRecord.size, variants.posToIdx(Variants::SV, i)));
                    breakpoints.push_back(std::make_pair((int)length(seq), variants.posToIdx(Variants::SV, i)));

                    currentPos = length(seq);
                }
                break;
            case StructuralVariantRecord::DUPLICATION:
                {
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    SEQAN_ASSERT_GEQ(svRecord.targetPos, svRecord.pos + svRecord.size);
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)  // first copy
                    {
                        appendValue(varPoints, std::make_pair((int)length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    if (currentPos != length(seq))
                        appendValue(intervals, GenomicInterval(currentPos, length(seq), svRecord.pos, svRecord.pos + svRecord.size,
                                                               '+', GenomicInterval::DUPLICATED));
                    unsigned tmpCurrentPos = length(seq);
                    append(seq, infix(contig, svRecord.pos + svRecord.size, svRecord.targetPos));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)
                    {
                        appendValue(varPoints, std::make_pair((int)length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos + svRecord.size, svRecord.targetPos));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos + svRecord.size, svRecord.targetPos));
                    }
                    if (tmpCurrentPos != length(seq))
                        appendValue(intervals, GenomicInterval(tmpCurrentPos, length(seq), svRecord.pos + svRecord.size, svRecord.targetPos,
                                                               '+', GenomicInterval::NORMAL));
                    tmpCurrentPos = length(seq);
                    append(seq, infix(contig, svRecord.pos, svRecord.pos + svRecord.size));
                    if (methSimOptions && methSimOptions->simulateMethylationLevels)  // second copy
                    {
                        appendValue(varPoints, std::make_pair((int)length(seq), false));
                        append(levelsLargeVariants->forward, infix(levels->forward, svRecord.pos, svRecord.pos + svRecord.size));
                        append(levelsLargeVariants->reverse, infix(levels->reverse, svRecord.pos, svRecord.pos + svRecord.size));
                    }
                    if (tmpCurrentPos != length(seq))
                        appendValue(intervals, GenomicInterval(tmpCurrentPos, length(seq), svRecord.pos, svRecord.pos + svRecord.size,
                                                               '+', GenomicInterval::NORMAL));
                    if (verbosity >= 3)
                        std::cerr << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << " (duplication)\n"
                                  << "append(seq, infix(contig, " << svRecord.pos + svRecord.size << ", " << svRecord.targetPos << ") " << __LINE__ << "\n"
                                  << "append(seq, infix(contig, " << svRecord.pos << ", " << svRecord.pos + svRecord.size << ") " << __LINE__ << "\n";
                    lastPos = svRecord.targetPos;
                    SEQAN_ASSERT_LT(lastPos, (int)length(contig));

                    // Write out small variant infos for duplication.
                    for (; itSmallVar != smallVarInfos.end() && itSmallVar->pos < svRecord.pos + svRecord.size; ++itSmallVar)
                    {
                        int offset = (int)currentPos - lastPos;
                        varInfos.push_back(*itSmallVar);
                        varInfos.back().pos += offset;

                        if (itSmallVar->pos < svRecord.pos + svRecord.size)
                        {
                            varInfos.push_back(*itSmallVar);
                            varInfos.back().pos += (svRecord.targetPos - svRecord.pos);
                        }
                    }

                    // Copy out breakpoints.
                    breakpoints.push_back(std::make_pair(currentPos, variants.posToIdx(Variants::SV, i)));
                    breakpoints.push_back(std::make_pair(currentPos + svRecord.pos + svRecord.size - svRecord.pos, variants.posToIdx(Variants::SV, i)));
                    breakpoints.push_back(std::make_pair(currentPos + svRecord.pos + svRecord.size - svRecord.pos + svRecord.targetPos - (svRecord.pos + svRecord.size), variants.posToIdx(Variants::SV, i)));
                    breakpoints.push_back(std::make_pair((int)length(seq), variants.posToIdx(Variants::SV, i)));

                    currentPos = length(seq);
                }
                break;
            default:
                return 1;
        }
    }
    if (verbosity >= 3)
        std::cerr << "append(seq, infix(contig, " << lastPos << ", " << length(contig) << ") "
                  << __LINE__ << " (last interim)\n";
    append(seq, infix(contig, lastPos, length(contig)));
    if (methSimOptions && methSimOptions->simulateMethylationLevels)
    {
        append(levelsLargeVariants->forward, infix(levels->forward, lastPos, length(contig)));
        append(levelsLargeVariants->reverse, infix(levels->reverse, lastPos, length(contig)));

        SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants->forward));
        SEQAN_ASSERT_EQ(length(seq), length(levelsLargeVariants->reverse));

        fixVariationLevels(*levelsLargeVariants, *rng, seq, varPoints, *methSimOptions);
    }
    if (currentPos != length(seq))
        appendValue(intervals, GenomicInterval(currentPos, length(seq), lastPos, length(contig),
                                               '+', GenomicInterval::NORMAL));

    // Copy out small variant infos for trailing characters.
    for (; itSmallVar != smallVarInfos.end(); ++itSmallVar)
    {
        int offset = (int)currentPos - lastPos;
        varInfos.push_back(*itSmallVar);
        varInfos.back().pos += offset;
    }

    // Build the interval trees of the positionMap.
    seqan::String<PositionMap::TInterval> svIntervals, svIntervalsSTL;
    for (unsigned i = 0; i < length(intervals); ++i)
        appendValue(svIntervals, PositionMap::TInterval(
                intervals[i].svBeginPos, intervals[i].svEndPos, intervals[i]));
    for (unsigned i = 0; i < length(intervals); ++i)
        if (intervals[i].smallVarBeginPos != -1)  // ignore insertions
            appendValue(svIntervalsSTL, PositionMap::TInterval(
                    intervals[i].smallVarBeginPos, intervals[i].smallVarEndPos, intervals[i]));
    createIntervalTree(positionMap.svIntervalTree, svIntervals);
    createIntervalTree(positionMap.svIntervalTreeSTL, svIntervalsSTL);

    return 0;
}

// --------------------------------------------------------------------------
// Function PositionMap::overlapsWithBreakpoint()
// --------------------------------------------------------------------------

bool PositionMap::overlapsWithBreakpoint(int svBeginPos, int svEndPos) const
{
    std::set<std::pair<int, int> >::const_iterator it = svBreakpoints.upper_bound(std::make_pair(svBeginPos, std::numeric_limits<int>::max()));
    return (it != svBreakpoints.end() && it->first < svEndPos);
}

// --------------------------------------------------------------------------
// Function PositionMap::getGenomicInterval()
// --------------------------------------------------------------------------

GenomicInterval PositionMap::getGenomicInterval(int svPos) const
{
    seqan::String<GenomicInterval> intervals;
    findIntervals(intervals, svIntervalTree, svPos);
    SEQAN_ASSERT_EQ(length(intervals), 1u);
    return intervals[0];
}

// --------------------------------------------------------------------------
// PositionMap::getGenomicIntervalSmallVarPos()
// --------------------------------------------------------------------------

// Returns the GenomicInterval on the sequence using a position on the small var reference.
GenomicInterval PositionMap::getGenomicIntervalSmallVarPos(int smallVarPos) const
{
    seqan::String<GenomicInterval> intervals;
    findIntervals(intervals, svIntervalTreeSTL, smallVarPos);
    return intervals[0];
}

// --------------------------------------------------------------------------
// Function PositionMap::toSmallVarInterval()
// --------------------------------------------------------------------------

std::pair<int, int> PositionMap::toSmallVarInterval(int svBeginPos, int svEndPos) const
{
    SEQAN_ASSERT(!overlapsWithBreakpoint(svBeginPos, svEndPos));
    GenomicInterval gi = getGenomicInterval(svBeginPos);
    if (gi.kind == GenomicInterval::INSERTED)
    {
        // novel sequence, cannot be projected
        return std::make_pair(-1, -1);
    }
    if (gi.kind != GenomicInterval::INVERTED)
    {
        // forward
        return std::make_pair(gi.smallVarBeginPos + (svBeginPos - gi.svBeginPos),
                              gi.smallVarBeginPos + (svEndPos - gi.svBeginPos));
    }
    else
    {
        // reverse
        return std::make_pair(gi.smallVarBeginPos + (gi.svEndPos - svBeginPos),
                              gi.smallVarBeginPos + (gi.svEndPos - svEndPos));
    }

    return std::make_pair(-1, -1);  // cannot reach here
}

// --------------------------------------------------------------------------
// Function PositionMap::smallVarToLargeVarInterval()
// --------------------------------------------------------------------------

std::pair<int, int> PositionMap::smallVarToLargeVarInterval(int beginPos, int endPos) const
{
    // SEQAN_ASSERT(!overlapsWithBreakpoint(svBeginPos, svEndPos));
    GenomicInterval gi = getGenomicIntervalSmallVarPos(beginPos);
    SEQAN_ASSERT_NEQ(gi.kind, GenomicInterval::INSERTED);

    if (gi.kind != GenomicInterval::INVERTED)
    {
        // forward
        return std::make_pair(gi.svBeginPos + (beginPos - gi.smallVarBeginPos),
                              gi.svBeginPos + (endPos - gi.smallVarBeginPos));
    }
    else
    {
        // reverse
        return std::make_pair(gi.svBeginPos + (gi.svEndPos - beginPos),
                              gi.svBeginPos + (gi.svEndPos - endPos));
    }

    return std::make_pair(-1, -1);  // cannot reach here
}

// --------------------------------------------------------------------------
// Function PositionMap::orignalToSmallVarInterval()
// --------------------------------------------------------------------------

std::pair<int, int> PositionMap::originalToSmallVarInterval(int beginPos, int endPos) const
{
    // TODO(holtgrew): Project to the left at the end.

    // Get anchor gaps objects from anchors.
    TGaps refGaps(seqan::Nothing(), refGapAnchors);
    TGaps smallVarGaps(seqan::Nothing(), smallVarGapAnchors);

    // Translate begin and end position.
    int beginPos2 = toViewPosition(refGaps, beginPos);
    int beginPos3 = toSourcePosition(smallVarGaps, beginPos2);
    int endPos2 = toViewPosition(refGaps, endPos);
    int endPos3 = toSourcePosition(smallVarGaps, endPos2);
    return std::make_pair(beginPos3, endPos3);
}

// --------------------------------------------------------------------------
// Function PositionMap::toOriginalInterval()
// --------------------------------------------------------------------------

std::pair<int, int> PositionMap::toOriginalInterval(int smallVarBeginPos, int smallVarEndPos) const
{
    // TODO(holtgrew): Project to the left at the end.

    // Get anchor gaps objects from anchors.
    TGaps refGaps(seqan::Nothing(), refGapAnchors);
    TGaps smallVarGaps(seqan::Nothing(), smallVarGapAnchors);

    // Translate begin and end position.
    int smallVarBeginPos2 = toViewPosition(smallVarGaps, smallVarBeginPos);
    int smallVarBeginPos3 = toSourcePosition(refGaps, smallVarBeginPos2);
    int smallVarEndPos2 = toViewPosition(smallVarGaps, smallVarEndPos);
    int smallVarEndPos3 = toSourcePosition(refGaps, smallVarEndPos2);
    return std::make_pair(smallVarBeginPos3, smallVarEndPos3);
}

// --------------------------------------------------------------------------
// Function PositionMap::reinit()
// --------------------------------------------------------------------------

void PositionMap::reinit(TJournalEntries const & journal)
{
    // Reset the interval tree and breakpoints.
    // TODO(holtgrew): Better API support for IntervalTree?
    svIntervalTree = TIntervalTree();
    svIntervalTreeSTL = TIntervalTree();
    svBreakpoints.clear();
    clear(refGapAnchors);
    clear(smallVarGapAnchors);

    // Convert the journal to two gaps.
    //
    // Get anchor gaps objects from anchors.
    typedef seqan::Iterator<TGaps, seqan::Standard>::Type TGapsIter;
    TGaps refGaps(seqan::Nothing(), refGapAnchors);
    TGapsIter itRef = begin(refGaps, seqan::Standard());
    TGaps smallVarGaps(seqan::Nothing(), smallVarGapAnchors);
    TGapsIter itVar = begin(smallVarGaps, seqan::Standard());

    // Iterate over the journal.
    typedef seqan::Iterator<TJournalEntries const, seqan::Standard>::Type TJournalEntriesIt;
    TJournalEntriesIt it = begin(journal, seqan::Standard());
    SEQAN_ASSERT_NEQ(it->segmentSource, seqan::SOURCE_NULL);
    SEQAN_ASSERT_EQ(it->virtualPosition, 0u);

    unsigned lastRefPos = std::numeric_limits<unsigned>::max();  // Previous position from reference.
    for (; it != end(journal, seqan::Standard()); ++it)
    {
        // std::cerr << *it << "\n";
        SEQAN_ASSERT_NEQ(it->segmentSource, seqan::SOURCE_NULL);
        if (it->segmentSource == seqan::SOURCE_ORIGINAL)
        {
            if (lastRefPos == std::numeric_limits<unsigned>::max())
            {
                if (it->physicalPosition != 0)
                {
                    insertGaps(itRef, it->physicalPosition);
                    itRef += it->physicalPosition;
                    itVar += it->physicalPosition;
                    lastRefPos = it->physicalPosition + it->length;
                    // std::cerr << "INSERT REF GAPS\t" << it->physicalPosition << "\n";
                }
                itRef += it->length;
                itVar += it->length;
                // std::cerr << "FORWARD\t" << it->length << "\n";
            }
            else
            {
                if (it->physicalPosition != lastRefPos)
                {
                    int len = it->physicalPosition - lastRefPos;
                    insertGaps(itVar, len);
                    // std::cerr << "INSERT VAR GAPS\t" << len << "\n";
                    itRef += len;
                    itVar += len;
                    // std::cerr << "FORWARD\t" << len << "\n";
                }
                itRef += it->length;
                itVar += it->length;
                // std::cerr << "2 FORWARD\t" << it->length << "\n";
            }
            lastRefPos = it->physicalPosition + it->length;
        }
        else
        {
            insertGaps(itRef, it->length);
            // std::cerr << "INSERT REF GAPS\t" << it->length << "\n";
            itRef += it->length;
            itVar += it->length;
            // std::cerr << "FORWARD\t" << it->length << "\n";
        }
    }

    // std::cerr << "--> done\n";

    // typedef seqan::Gaps<seqan::CharString, seqan::AnchorGaps<TGapAnchors> > TGaps2;
    // seqan::CharString seqH;
    // seqan::CharString seqV;
    // for (unsigned i = 0; i < 1000; ++i)
    // {
    //     appendValue(seqH, 'X');
    //     appendValue(seqV, 'X');
    // }
    // TGaps2 gapsH(seqH, refGapAnchors);
    // TGaps2 gapsV(seqV, smallVarGapAnchors);

    // std::cerr << "REF\t" << gapsH << "\n"
    //           << "VAR\t" << gapsV << "\n";
}
