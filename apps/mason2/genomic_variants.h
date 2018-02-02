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
// Data structures for storing genomic variants and routines to work on them.
//
// Variants are split into SNPs, small indels, and large structural variants.
// There is a data structure for each class of variant.  Doing such a
// separation allows for almost simple translation between the reference
// coordinate system and the one for the genome including all variants.
//
// There are routines provided for
// ==========================================================================

// TODO(holtgrew): Could put the coordinate system into types. This would add static checking on coordinate system.

#ifndef APPS_MASON2_GENOMIC_VARIANTS_H_
#define APPS_MASON2_GENOMIC_VARIANTS_H_

#include <seqan/align.h>
#include <seqan/misc/interval_tree.h>

#include "methylation_levels.h"

// ============================================================================
// Forwards
// ============================================================================

typedef seqan::JournalEntries<seqan::JournalEntry<unsigned, int>, seqan::SortedArray> TJournalEntries;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class SnpRecord
// --------------------------------------------------------------------------

// Represents a SNP in one haplotype.
//
// All coordinates are given as source coordinates.

struct SnpRecord
{
    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // The target nucleotide.
    seqan::Dna5 to;

    SnpRecord() : rId(-1), pos(-1), haplotype(-1), to('\0')
    {}

    SnpRecord(int haplotype, int rId, int pos, char to) :
            rId(rId), pos(pos), haplotype(haplotype), to(to)
    {}

    bool operator<(SnpRecord const & other) const
    {
        if (rId < other.rId || (rId == other.rId && pos < other.pos) ||
            (rId == other.rId && pos == other.pos && haplotype < other.haplotype))
            return true;
        return false;
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }
};

std::ostream & operator<<(std::ostream & out, SnpRecord const & record);

// --------------------------------------------------------------------------
// Class SmallIndelRecord
// --------------------------------------------------------------------------

// Represents a small indel.
//
// All coordinates are given as source coordinates.

struct SmallIndelRecord
{
    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // The size of the indel, negative numbers for deletions, positive numbers for insertions.
    int size;

    // The inserted sequence if any.
    seqan::CharString seq;

    SmallIndelRecord() : rId(-1), pos(-1), haplotype(-1), size(0)
    {}

    SmallIndelRecord(int haplotype, int rId, int pos, int size, seqan::CharString const & seq) :
            rId(rId), pos(pos), haplotype(haplotype), size(size), seq(seq)
    {}

    bool operator<(SmallIndelRecord const & other) const
    {
        return getPos() < other.getPos();
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }
};

std::ostream & operator<<(std::ostream & out, SmallIndelRecord const & record);

// --------------------------------------------------------------------------
// Class SmallVarInfo
// --------------------------------------------------------------------------

// Information about small variants.

struct SmallVarInfo
{
    enum Kind { SNP, INS, DEL };

    Kind kind;
    int pos;
    int count;

    SmallVarInfo() : kind(SNP), pos(-1), count(0) {}
    SmallVarInfo(Kind kind, int pos, int count) : kind(kind), pos(pos), count(count) {}

    bool operator<(SmallVarInfo const & other) const
    { return (pos < other.pos || (pos == other.pos && kind < other.kind)); }
};

// --------------------------------------------------------------------------
// Class StructuralVariantRecord
// --------------------------------------------------------------------------

// Store structural variant information.
//
// All coordinates are given as source coordinates.

struct StructuralVariantRecord
{
    enum Kind
    {
        INVALID,
        INDEL,
        INVERSION,
        TRANSLOCATION,
        DUPLICATION
    };

    // The kind of the SV.
    Kind kind;

    // Reference id and position on the reference.
    int rId;
    int pos;

    // The haplotype that this variation belongs to.
    int haplotype;

    // In case of indel, negative numbers give deletions, positionve numbers give insertions.  In case of inversions,
    // the length of the inverted sequence.  In case of translocation the length of the translocated sequence, and in
    // the case the length of the duplicated sequence.
    int size;

    // Target reference id and position on this reference.  Currently, we only have SVs on the same chromosome.  Used in
    // case of translocation and duplication.
    int targetRId;
    int targetPos;

    // The inserted sequence if any.
    seqan::CharString seq;

    StructuralVariantRecord() :
            kind(INVALID), rId(-1), pos(-1), haplotype(-1), size(0), targetRId(-1), targetPos(-1)
    {}

    StructuralVariantRecord(Kind kind, int haplotype, int rId, int pos, int size,
                            int targetRId = -1, int targetPos = -1) :
            kind(kind), rId(rId), pos(pos), haplotype(haplotype), size(size), targetRId(targetRId),
            targetPos(targetPos)
    {}

    bool operator<(StructuralVariantRecord const & other) const
    {
        return getPos() < other.getPos();
    }

    std::pair<int, int> getPos() const
    {
        return std::make_pair(rId, pos);
    }

    // Returns true if query is within one base of the breakend.
    bool nearBreakend(int query) const;

    // Returns end position of SV.
    //
    // Return max value of int in case it is marked as invalid.
    int endPosition() const;

    // Return true if this SV overlaps with other within one bp.
    bool overlapsWith(StructuralVariantRecord const & other) const
    {
        return (other.pos <= endPosition()) && (pos <= other.endPosition());
    }
};

std::ostream & operator<<(std::ostream & out, StructuralVariantRecord const & record);

// --------------------------------------------------------------------------
// Class Variants
// --------------------------------------------------------------------------

// Contains the simulated variants.

struct Variants
{
    // SNP records.
    seqan::String<SnpRecord> snps;

    // Small indel records.
    seqan::String<SmallIndelRecord> smallIndels;

    // Structural variation record.
    seqan::String<StructuralVariantRecord> svRecords;

    // Names of the snps, smallIndels, SVs, as in the VCF file.
    seqan::StringSet<seqan::CharString> snpIDs;
    seqan::StringSet<seqan::CharString> smallIndelIDs;
    seqan::StringSet<seqan::CharString> svIDs;

    void clear()
    {
        seqan::clear(snps);
        seqan::clear(smallIndels);
        seqan::clear(svRecords);

        seqan::clear(snpIDs);
        seqan::clear(smallIndelIDs);
        seqan::clear(svIDs);
    }

    enum IndexKind
    {
        SNP,
        SMALL_INDEL,
        SV
    };

    // Return the name of a variant from its index.
    seqan::CharString getVariantName(int idx) const
    {
        std::pair<IndexKind, int> res = resolveIdx(idx);
        if (res.first == SNP)
        {
            if (res.second >= (int)length(snpIDs))
                return ".";
            return snpIDs[res.second];
        }
        else if (res.first == SMALL_INDEL)
        {
            if (res.second >= (int)length(smallIndelIDs))
                return ".";
            return smallIndelIDs[res.second];
        }
        else
        {
            if (res.second >= (int)length(svIDs))
                return ".";
            return svIDs[res.second];
        }
    }

    // Translate number of variant with type to ids.
    int posToIdx(IndexKind kind, int pos) const
    {
        switch (kind)
        {
            case SNP:
                return pos;
            case SMALL_INDEL:
                return pos + length(snps);
            case SV:
                return pos + length(snps) + length(smallIndels);
            default:
                SEQAN_FAIL("Cannot reach here.");
        }

        SEQAN_FAIL("Invalid kind!");
        return -1;
    }

    // Each variant has a numeric index.  Indices go up to (length(snps) + length(smallIndels) + length(svRecords)) - 1
    // and are "stacked".
    std::pair<IndexKind, int> resolveIdx(int idx) const
    {
        if (idx < (int)length(snps))
            return std::make_pair(SNP, idx);
        else if (idx < (int)(length(snps) + (int)length(smallIndels)))
            return std::make_pair(SMALL_INDEL, (int)(idx - length(snps)));
        else if (idx < (int)(length(snps) + (int)(length(smallIndels) + length(svRecords))))
            return std::make_pair(SV, (int)(idx - length(snps) - length(smallIndels)));
        else
            SEQAN_FAIL("Invalid idx!");
        return std::make_pair(SNP, -1);
    }
};

// --------------------------------------------------------------------------
// Class GenomicInterval
// --------------------------------------------------------------------------

// We annotate intervals (in an interval tree) from the genome with structural variants with intervals in the sequence
// with small variants.  These intervals are stored as GenomicInterval objects.

struct GenomicInterval
{
    // The kind of the interval.  Can be inserted, inverted, duplicated, or anything else (normal = translocation,
    // untouched.
    enum Kind
    {
        NORMAL,
        INSERTED,
        INVERTED,
        DUPLICATED
    };

    // Begin and end position on sequence with SVs.
    int svBeginPos;
    int svEndPos;

    // Begin and end position on sequence without SVs.  Set to [-1, -1) in case of insertion.
    int smallVarBeginPos;
    int smallVarEndPos;

    // The strand of the interval '+'/'-'.
    char strand;

    // The kind of the interval.
    Kind kind;

    explicit
    GenomicInterval(int svBeginPos = -1, int svEndPos = -1,
                    int smallVarBeginPos = -1, int smallVarEndPos = -1,
                    char strand = '.', Kind kind = NORMAL) :
            svBeginPos(svBeginPos), svEndPos(svEndPos), smallVarBeginPos(smallVarBeginPos),
            smallVarEndPos(smallVarEndPos), strand(strand), kind(kind)
    {}

    bool operator==(GenomicInterval const & other) const
    {
        return (svBeginPos == other.svBeginPos) &&
                (svEndPos == other.svEndPos) &&
                (smallVarBeginPos == other.smallVarBeginPos) &&
                (smallVarEndPos == other.smallVarEndPos) &&
                (strand == other.strand);
    }

    bool operator!=(GenomicInterval const & other) const
    {
        return !(*this == other);
    }
};

inline bool cmpIntervalSmallVarPos(GenomicInterval const & lhs,
                                   GenomicInterval const & rhs)
{
    return std::make_pair(lhs.smallVarBeginPos, lhs.smallVarEndPos) < std::make_pair(rhs.smallVarBeginPos, rhs.smallVarEndPos);
}

// --------------------------------------------------------------------------
// Class PositionMap
// --------------------------------------------------------------------------

// Used for translating between the sequence with SVs, small variants, and original sequence.

class PositionMap
{
public:
    typedef int TValue;
    typedef GenomicInterval TCargo;
    typedef seqan::IntervalAndCargo<TValue, TCargo> TInterval;
    typedef seqan::IntervalTree<TValue, TCargo> TIntervalTree;

    typedef seqan::String<seqan::GapAnchor<int> > TGapAnchors;
    typedef seqan::Gaps<seqan::Nothing, seqan::AnchorGaps<TGapAnchors> > TGaps;

    // TODO(holtgrew): We need a function *in this class* that builds the large variants data strutures for better encapsulation!

    // Gap anchors for gaps for translating between original and small variant coordinate system.
    TGapAnchors refGapAnchors, smallVarGapAnchors;
    // The mapping from the genome with large variants to the one with small variants.
    TIntervalTree svIntervalTree;
    // The mapping from the genome with small variants to the one with large variants.
    TIntervalTree svIntervalTreeSTL;  // small-to-large
    // The breakpoints (as [(point, idx)] on the sequence with variants.
    std::set<std::pair<int, int> > svBreakpoints;

    // Returns true if the interval on the sequence with structural variants overlaps with a breakpoint.
    bool overlapsWithBreakpoint(int svBeginPos, int svEndPos) const;

    // Returns the GenomicInterval on the sequence with small variants for the given position on the sequence with SVs.
    GenomicInterval getGenomicInterval(int svPos) const;

    // Returns the GenomicInterval on the sequence using a position on the small var reference.
    GenomicInterval getGenomicIntervalSmallVarPos(int smallVarPos) const;

    // Translates an interval on the sequence with small variants to an interval on the original sequence.  The
    // translation is done in such a way that when the begin position is in an insertion, it is projected to the right
    // of the gap and if the end position is in an insertion, it is projected to the left of the gap.
    std::pair<int, int> toOriginalInterval(int smallVarBeginPos, int smallVarEndPos) const;

    // Translates an interval on the sequence with SVs to an interval on the sequence with small variants.
    //
    // Returns (-1, -1) if the interval lies in an insertion.
    //
    // Returns (a, b), a > b if on the reverse strand
    //
    // The interval must not overlap with a breakpoint.
    std::pair<int, int> toSmallVarInterval(int svBeginPos, int svEndPos) const;

    // Translate the interval on the original sequence into coordinates with small variants.
    std::pair<int, int> originalToSmallVarInterval(int beginPos, int endPos) const;

    // Translate the interval on the reference with small variants to the one with large variants.
    std::pair<int, int> smallVarToLargeVarInterval(int beginPos, int endPos) const;

    // Reset the PositionMap with the length of the original sequence.
    void reinit(TJournalEntries const & journal);
};

// --------------------------------------------------------------------------
// Class VariantMaterializer
// --------------------------------------------------------------------------

// Materialize variants stored in a Variants object.
//
// Note that the class assumes that all variants come from the same contig and haplotype.

// TODO(holtgrew): Rename to ContigMaterializer.

class VariantMaterializer
{
public:
    // The random number generator to use for methylation levels.
    TRng * rng;
    // The Variants to materialize for.
    Variants const * variants;
    // Options for the methylation level simulator.  Methylation simulation is required for fixing methylation levels.
    MethylationLevelSimulatorOptions const * methSimOptions;

    // Verbosity.
    int verbosity;

    VariantMaterializer() : rng(), variants(), methSimOptions(), verbosity(1)
    {}

    VariantMaterializer(TRng & rng, Variants const & variants) :
            rng(&rng), variants(&variants), methSimOptions(), verbosity(1)
    {}

    VariantMaterializer(TRng & rng, Variants const & variants, MethylationLevelSimulatorOptions const & methSimOptions) :
            rng(&rng), variants(&variants), methSimOptions(&methSimOptions), verbosity(1)
    {}

    // Materialize the variants from the haplotype with the given id in *variants to result given the reference sequence refSeq.
    //
    // Breakpoints is a vector of (point, id) where point is a point on the materialized contig with variants and id is
    // an integer index into variants.  See Variants::resolveIdx() for more information.
    int run(seqan::Dna5String & resultSeq,
            PositionMap & posMap,
            std::vector<SmallVarInfo> & smallVars,
            std::vector<std::pair<int, int> > & breakpoints,
            seqan::Dna5String const & refSeq,
            int haplotypeId)
    {
        return _runImpl(&resultSeq, &posMap, 0, smallVars, breakpoints, &refSeq, 0, haplotypeId);
    }

    // Same as the run() above, but including reference levels.
    int run(seqan::Dna5String & resultSeq,
            PositionMap & posMap,
            MethylationLevels & resultLvls,
            std::vector<SmallVarInfo> & smallVars,
            std::vector<std::pair<int, int> > & breakpoints,
            seqan::Dna5String const & refSeq,
            MethylationLevels const & refLvls,
            int haplotypeId)
    {
        return _runImpl(&resultSeq, &posMap, &resultLvls, smallVars, breakpoints, &refSeq, &refLvls, haplotypeId);
    }

    // Implementation of the materialization, uses pointers instead of references for deciding whether materializing
    // levels or not.
    int _runImpl(seqan::Dna5String * resultSeq,
                 PositionMap * posMap,
                 MethylationLevels * resultLvls,
                 std::vector<SmallVarInfo> & smallVars,
                 std::vector<std::pair<int, int> > & breakpoints,
                 seqan::Dna5String const * ref,
                 MethylationLevels const * refLvls,
                 int haplotypeId);

    // Materialization of the small variants.
    //
    // Levels passed as NULL if not given.
    int _materializeSmallVariants(seqan::Dna5String & seq,
                                  TJournalEntries & journal,
                                  MethylationLevels * levelsSmallVariants,
                                  std::vector<SmallVarInfo> & smallVarInfos,
                                  seqan::Dna5String const & contig,
                                  Variants const & variants,
                                  MethylationLevels const * levels,
                                  int hId);

    // Materialization of the large variants.
    //
    // Levels passed as NULL if not given.
    int _materializeLargeVariants(seqan::Dna5String & seq,
                                  MethylationLevels * levelsLargeVariants,
                                  std::vector<SmallVarInfo> & varInfos,
                                  std::vector<std::pair<int, int> > & breakpoints,
                                  PositionMap & positionMap,
                                  TJournalEntries const & journal,
                                  seqan::Dna5String const & contig,
                                  std::vector<SmallVarInfo> const & smallVarInfos,
                                  Variants const & variants,
                                  MethylationLevels const * levels,
                                  int hId);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef APPS_MASON2_GENOMIC_VARIANTS_H_
