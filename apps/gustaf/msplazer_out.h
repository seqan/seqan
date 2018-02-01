// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011-2018, Kathrin Trappe, FU Berlin
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
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_APPS_GUSTAF_MSPLAZER_OUT_H_
#define SEQAN_APPS_GUSTAF_MSPLAZER_OUT_H_

#include <iostream>
#include <fstream>
#include <seqan/gff_io.h>
#include <seqan/vcf_io.h>
#include <seqan/file.h>

#include "../stellar/stellar.h"
#include "msplazer.h"
#include "gustaf_matepairs.h"

using namespace seqan;

struct DotDrawingMSplazer_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazer;

struct DotDrawingMSplazerBestChain_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazerBestchain;


template <typename TGraph, typename TVertexDescriptor, typename TScoreAlloc, typename TMatch, // typename TFile, 
          typename TBreakpoint, typename TPos, typename TMatchAlloc, typename TID>
// typename TBreakpointAlloc, typename TMatchAlloc> // Requires Value<SparsePropertyMap> specialisation in msplazer.h
void
write(std::ostream & out,
      // MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, TBreakpointAlloc, // Requires Value<SparsePropertyMap> specialisation in msplazer.h
      MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, SparsePropertyMap<TBreakpoint, TPos>,
                    TMatchAlloc> const & msplazerchain,
      TMatch const & queryMatches,
      unsigned const & queryLength,
      TID const & queryID,
      DotDrawingMSplazer const &)
{
    // IOREV _doc_ _batchreading_
    // typedef typename Value<TBreakpointAlloc>::Type TBreakpoint; // Requires Value<SparsePropertyMap> specialisation in msplazer.h
    typedef typename TBreakpoint::TId TId;

    // _writeGraphType(file,g,DotDrawing());
    out << "digraph G {\n";
    out << '\n';
    out << "/* Graph Attributes */\n";
    out << "graph [rankdir = LR, clusterrank = local];\n";
    out << '\n';
    out << "/* Node Attributes */\n";
    out << "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n";
    out << '\n';
    out << "/* Edge Attributes */\n";
    out << "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n";
    out << '\n';

    out << "/* Nodes */\n";
    typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
    TConstIter it(msplazerchain.graph);
    unsigned i = 0;
    bool atEndV = false;
    for (; !atEnd(it); ++it)
    {
        TId sId;
        if (i < length(queryMatches))
        {
            out << *it;
            out << " [label = \"";
            out << "chr: ";
            _getShortId(sId, queryMatches[i].id);
            out << sId;
            out << "\\n db: ";
            out << queryMatches[i].begin1 + 1;
            out << "...";
            out << queryMatches[i].end1;
            out << "  ";
            out << (queryMatches[i].orientation ? '+' : '-');
            out << "\\n read: ";
            out << queryMatches[i].begin2 + 1;
            out << "...";
            out << queryMatches[i].end2;
            out << "\\n";
            out << i;
            if (_isLeftMate(queryMatches[i], msplazerchain.mateJoinPosition))
                out << "\", fontcolor= green];\n";
            else
                out << "\", fontcolor= blue];\n";
        }
        else if (!atEndV)
        {
            out << *it;
            out << " [label = \"start";
            out << "\"];\n";
            atEndV = true;
        }
        else if (atEndV)
        {
            out << *it;
            out << " [label = \"end";
            out << "\\n";
            out << queryLength;
            out << "\"];\n";
        }
        else
            std::cerr << "in writing dot: default vertex??" << std::endl;
        ++i;
    }

    out << '\n';
    out << "/* Edges */\n";
    typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
    TConstEdIter itEd(msplazerchain.graph);
    for (; !atEnd(itEd); ++itEd)
    {
        TVertexDescriptor sc = sourceVertex(itEd);
        TVertexDescriptor tr = targetVertex(itEd);
        out << sc;
        _writeEdgeType(out, msplazerchain.graph, DotDrawing());
        out << tr;
        out << " [label = \"";
        out << getCargo(*itEd);
        TBreakpoint bp;
        bool foundBP = getProperty(msplazerchain.breakpoints, value(itEd), bp);
        if (foundBP)
            out << "*";
        out << "\"];\n";
    }
        out << '\n'
        << "overlap=false\n label = \"readID: " << queryID << "\"fontsize=10;\n"
        << "}\n";
}

// DotWriting call for read graphs
template <typename TSequence, typename TId, typename TMSplazerChain>
void _writeDotfiles(StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & stellarMatches,
                    StringSet<TSequence> const & queries,
                    StringSet<TId> const & queryIDs,
                    String<TMSplazerChain> & queryChains,
                    MSplazerOptions const & msplazerOptions)
{
    for (unsigned i = 0; i < length(queryChains); ++i)
    {
        if (!queryChains[i].isEmpty)
        {
            std::string fn = "read" + toString(i + 1) + '_' + toString(msplazerOptions.jobName) + ".dot";
            std::ofstream f( toCString("read" + toString(i + 1) + '_' + toString(msplazerOptions.jobName) + ".dot"));
            if (!f.good())
                std::cerr << "Error while writing dot files!" << std::endl;
            write(f, queryChains[i], stellarMatches[i].matches, length(queries[i]), queryIDs[i], DotDrawingMSplazer());
            f.close();
        }
    }
}

template <typename TBreakpoint>
void _setGffRecordType(GffRecord & record, TBreakpoint & bp)
{
    switch(bp.svtype)
    {
        case TBreakpoint::INSERTION:
        record.type = "insertion";
        break;
        case TBreakpoint::DELETION:
        record.type = "deletion";
        break;
        case TBreakpoint::INVERSION:
        record.type = "inversion";
        break;
        case TBreakpoint::SEQAN_TANDEM:
        record.type = "tandem";
        break;
        case TBreakpoint::DISPDUPLICATION:
        record.type = "duplication";
        break;
        case TBreakpoint::INTERTRANSLOCATION:
        record.type = "inter-chr-translocation";
        break;
        case TBreakpoint::TRANSLOCATION:
        record.type = "translocation";
        break;
        case TBreakpoint::BREAKEND:
        record.type = "breakend";
        break;
        default:
        record.type = "invalid";
    }

}

template <typename TBreakpoint>
inline void _fillGffRecordDuplication(GffRecord & record, TBreakpoint & bp, unsigned id)
{
    typedef typename TBreakpoint::TId TId;
    TId sId;
    typedef typename TBreakpoint::TPos TPos;
    _getShortId(sId, bp.startSeqId);
    record.ref = sId;
    record.source = "GUSTAF";
    record.type = "duplication";
    TPos begin, end, target = std::numeric_limits<unsigned>::max();
    // Using set function for VCF duplication to set positions
    _setVcfRecordDuplicationPos(bp, begin, end, target);
    if (begin > end)
        std::swap(begin, end);
    record.beginPos = begin;
    record.endPos = end;
    /*
    std::cerr << "#####################################################" << std::endl;
    std::cerr << bp << std::endl;
    std::cerr << "begin " << begin << " end " << end << " target " << std::cerr;
    std::cerr << "record end " << record.endPos << std::endl;
    std::cerr << "#####################################################" << std::endl;
    */
    record.strand = '+';
    appendValue(record.tagNames, "ID");
    appendValue(record.tagValues, toString(id));
    if (target != std::numeric_limits<unsigned>::max())
    {
        appendValue(record.tagNames, "size");
        appendValue(record.tagValues, toString((end - begin)));
        appendValue(record.tagNames, "targetPos");
        appendValue(record.tagValues, toString(target));
    }
    else
    {
        std::stringstream dpos;
        dpos << (begin + 1) << "|" << end;
        appendValue(record.tagNames, "size");
        appendValue(record.tagValues, "imprecise");
        appendValue(record.tagNames, "targetPos");
        appendValue(record.tagValues, dpos.str());
    }
    appendValue(record.tagNames, "support");
    appendValue(record.tagValues, toString(bp.support));
    appendValue(record.tagNames, "supportIds");

    std::stringstream s;
    for (unsigned i = 0; i < length(bp.supportIds); ++i)
    {
        s << bp.supportIds[i];
        s << ',';
    }

    appendValue(record.tagValues, s.str());
}

template <typename TBreakpoint>
inline void _fillGffRecord(GffRecord & record, TBreakpoint & bp, unsigned id)
{
    typedef typename TBreakpoint::TId TId;
    TId sId;
    typedef typename TBreakpoint::TPos TPos;
    _getShortId(sId, bp.startSeqId);
    record.ref = sId;
    record.source = "GUSTAF";
    _setGffRecordType(record, bp);
//    record.type = bp.svtype;
    record.beginPos = bp.startSeqPos;
    if (bp.svtype == 2 || bp.svtype == 3) // 2=deletion;3=inversion;
        record.endPos = bp.endSeqPos;
    else
        record.endPos = bp.startSeqPos + 1;
    if (bp.startSeqStrand)
        record.strand = '+';
    else
        record.strand = '-';

    appendValue(record.tagNames, "ID");
    appendValue(record.tagValues, toString(id));
    if (bp.svtype == 1) // 1=insertion
    {
        appendValue(record.tagNames, "size");
        appendValue(record.tagValues, "-" + toString(length(bp.insertionSeq)));
        appendValue(record.tagNames, "seq");
        appendValue(record.tagValues, bp.insertionSeq);
    }
    else if (bp.svtype == 2 || bp.svtype == 8) // 2=deletion;8=breakend
    {
        appendValue(record.tagNames, "size");
        appendValue(record.tagValues, toString(static_cast<TPos>(bp.endSeqPos - bp.startSeqPos)));
    }
    else if (bp.svtype == 7) // 7=translocation
    {
        appendValue(record.tagNames, "middlePos");
        appendValue(record.tagValues, toString(bp.dupMiddlePos));
        appendValue(record.tagNames, "endPos");
        appendValue(record.tagValues, toString(bp.endSeqPos));
    }
    else
    {
        appendValue(record.tagNames, "endChr");
        _getShortId(sId, bp.endSeqId);
        appendValue(record.tagValues, sId);
        appendValue(record.tagNames, "endPos");
        appendValue(record.tagValues, toString(bp.endSeqPos));
        appendValue(record.tagNames, "endStrand");
        if (bp.endSeqStrand)
            appendValue(record.tagValues, "+");
        else
            appendValue(record.tagValues, "-");
    }
    appendValue(record.tagNames, "support");
    appendValue(record.tagValues, toString(bp.support));
    appendValue(record.tagNames, "supportIds");

    std::stringstream s;
    for (unsigned i = 0; i < length(bp.supportIds); ++i)
    {
        s << bp.supportIds[i];
        s << ',';
    }

    appendValue(record.tagValues, s.str());
    appendValue(record.tagNames, "breakpoint");
    appendValue(record.tagValues, toString(bp.readStartPos + 1));
}

// Breakpoint writing call
template <typename TBreakpoint>
bool _writeGlobalBreakpoints(String<TBreakpoint> & globalBreakpoints,
                             MSplazerOptions const & msplazerOptions, Gff /*tag*/)
{
    // Creating GFF record and writing to gff file
    std::string fn_gff = toString(msplazerOptions.gffOutFile);
    seqan::GffFileOut gffOut;
    if (!open(gffOut, fn_gff.c_str()))
    {
        std::cerr << "Error: Could not open file\n";
        return 1;
    }
    if (length(globalBreakpoints) == 0)
    {
        std::cerr << "Empty output list, skip writing gff" << std::endl;
        return 1;
    }

    GffRecord gff_record;
    for (unsigned i = 0; i < length(globalBreakpoints); ++i)
    {
        TBreakpoint & tempBP = globalBreakpoints[i];
        // Added support check in 2phase breakpoint extraction adaption (but only for bp, not be...)
        //if (tempBP.svtype != 0) // 0=invalid
        if (tempBP.svtype != 0 && tempBP.support >= msplazerOptions.support) // 0=invalid
        {
            if (tempBP.svtype == TBreakpoint::DISPDUPLICATION && tempBP.translSuppStartPos && tempBP.translSuppEndPos)
                tempBP.svtype = TBreakpoint::TRANSLOCATION;
            // Fill record
            if (tempBP.svtype == TBreakpoint::DISPDUPLICATION)
                _fillGffRecordDuplication(gff_record, tempBP, i);
            else
                _fillGffRecord(gff_record, tempBP, i);
            // Write record
            try
            {
                writeRecord(gffOut, gff_record);
            }
            catch (seqan::IOError const & ioErr)
            {
                std::cerr << "Error while writing breakpoint gff record (" << ioErr.what() << ")!\n";
            }
            clear(gff_record);
        }
    }

    std::cout << " completed gff writing " << msplazerOptions.gffOutFile << std::endl;
    return 0;
}

// Functions to fill the VCF entries for the specific SV type
// Begin position is always position before the variant, end position behind the variant.
// Target position is position before the event.
template <typename TBreakpoint, typename TSequence>
inline void _fillVcfRecordInsertion(VcfRecord & record, TBreakpoint & bp, TSequence & ref, unsigned id)
{
    record.rID = id;
    record.beginPos = bp.startSeqPos - 1; // In VCF position before event
    record.filter = "PASS";
    std::stringstream ss;
    ss << "SVTYPE=INS";
    SEQAN_ASSERT_GEQ_MSG(bp.endSeqPos, bp.startSeqPos, "Insertion end position smaller than begin position!");
    ss << ";SVLEN=" << length(bp.insertionSeq);
    //if (bp.similar != std::numeric_limits<unsigned>::max())
        ss << ";BM=" << bp.similar;
    ss << ";DP=" << bp.support;
    record.info = ss.str();

    // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
    // deletion of length k.
    if (bp.startSeqPos != 0)
        appendValue(record.ref, ref[bp.startSeqPos - 1]); // position/base before SV except for position 1, then it is the base AFTER the event!
    else
        appendValue(record.ref, ref[length(bp.insertionSeq)]); // position/base before SV except for position 1, then it is the base AFTER the event!

    // Compute ALT columns and a map to the ALT.
    if (length(bp.insertionSeq) < 20)
    {
        appendValue(record.alt, record.ref[0]);
        append(record.alt, bp.insertionSeq);
    } else
        append(record.alt, "<INS>");

    // Create genotype infos.
    appendValue(record.genotypeInfos, "1");
}

template <typename TBreakpoint, typename TSequence>
inline void _fillVcfRecordDeletion(VcfRecord & record, TBreakpoint & bp, TSequence & ref, unsigned id)
{
    record.rID = id;
    record.beginPos = bp.startSeqPos - 1;
    record.filter = "PASS";
    std::stringstream ss;
    ss << "SVTYPE=DEL";
    ss << ";END=" << bp.endSeqPos + 1 - 1; // 1-base adjustment, -1 bc endPos is behind last variant position
    SEQAN_ASSERT_GEQ_MSG(bp.endSeqPos, bp.startSeqPos, "Deletion end position smaller than begin position!");
    ss << ";SVLEN=-" << bp.endSeqPos-bp.startSeqPos;
    //if (bp.similar != std::numeric_limits<unsigned>::max())
        ss << ";BM=" << bp.similar;
    ss << ";DP=" << bp.support;
    record.info = ss.str();

    // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
    // deletion of length k.
    if (bp.startSeqPos != 0)
        appendValue(record.ref, ref[bp.startSeqPos - 1]);// position/base before SV except for position 1, then it is the base AFTER the event!
    else
        appendValue(record.ref, ref[bp.endSeqPos]); // position/base before SV except for position 1, then it is the base AFTER the event!


    // Compute ALT columns and a map to the ALT.
    if ((bp.endSeqPos - bp.startSeqPos) < 21)
    {
        appendValue(record.alt, record.ref[0]); // std::max(0,
        append(record.ref, infix(ref, bp.startSeqPos, bp.endSeqPos)); // Deletions on reverse strand??? correct positions??
    } else
        append(record.alt, "<DEL>");


    // Create genotype infos.
    appendValue(record.genotypeInfos, "1");
}

template <typename TBreakpoint, typename TSequence>
inline void _fillVcfRecordInversion(VcfRecord & record, TBreakpoint & bp, TSequence & ref, unsigned id)
{
    record.rID = id;
    record.beginPos = bp.startSeqPos - 1;
    record.filter = "PASS";
    std::stringstream ss;
    ss << "SVTYPE=INV";
    ss << ";END=" << bp.endSeqPos + 1 - 1; // 1-base adjustment, -1 bc endPos is behind last variant position
    SEQAN_ASSERT_GEQ_MSG(bp.endSeqPos, bp.startSeqPos, "Inversion end position smaller than begin position!");
    ss << ";SVLEN=" << bp.endSeqPos-bp.startSeqPos;
    //if (bp.similar != std::numeric_limits<unsigned>::max())
        ss << ";BM=" << bp.similar;
    ss << ";DP=" << bp.support;
    record.info = ss.str();

    // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
    // deletion of length k.
    if (bp.startSeqPos != 0)
        appendValue(record.ref, ref[bp.startSeqPos - 1]);
    else
        appendValue(record.ref, ref[bp.endSeqPos]);

    // Compute ALT columns and a map to the ALT.
    append(record.alt, "<INV>");

    // Create genotype infos.
    appendValue(record.genotypeInfos, "1");
}

template <typename TBreakpoint, typename TSequence>
inline void _fillVcfRecordTandem(VcfRecord & record, TBreakpoint & bp, TSequence & ref, unsigned id)
{
    record.rID = id;
    record.beginPos = bp.startSeqPos - 1;
    record.filter = "PASS";
    std::stringstream ss;
    ss << "SVTYPE=DUP";
    ss << ";END=" << bp.endSeqPos + 1 - 1; // 1-base adjustment, -1 bc endPos is behind last variant position
    SEQAN_ASSERT_GEQ_MSG(bp.endSeqPos, bp.startSeqPos, "Tandem duplication end position smaller than begin position!");
    ss << ";SVLEN=" << bp.endSeqPos-bp.startSeqPos - 1; // -1 bc positions are flanking the variant region
    //if (bp.similar != std::numeric_limits<unsigned>::max())
        ss << ";BM=" << bp.similar;
    ss << ";DP=" << bp.support;
    record.info = ss.str();

    // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
    // deletion of length k.
    if (bp.startSeqPos != 0)
        appendValue(record.ref, ref[bp.startSeqPos - 1]);
    else
        appendValue(record.ref, ref[bp.endSeqPos]);

    // Compute ALT columns and a map to the ALT.
    append(record.alt, "<DUP:TANDEM>");

    // Create genotype infos.
    appendValue(record.genotypeInfos, "1");
}

// This function determines begin, end and target position for duplication output depending on whether the target
// position is known and if it is upstream or downstream of the duplicated region. Since the positions are directly
// passed to the output stream and not via the vcf record, we add +1 here for the last position.
template <typename TBreakpoint, typename TPos>
inline bool _setVcfRecordDuplicationPos(TBreakpoint & bp, TPos & begin, TPos & end, TPos & target)
{
    if (bp.dupMiddlePos != std::numeric_limits<unsigned>::max())
    {
	// Downstream duplication dup(middlePos, endPos, startPos)
        if (bp.dupTargetPos == bp.startSeqPos)
        {
            begin = bp.dupMiddlePos;
            end = bp.endSeqPos;
            target = bp.startSeqPos;
            return true;
        }
	// Upstream duplication dup(startPos, middlePos, endPos)
        begin = bp.startSeqPos;
        end = bp.dupMiddlePos;
        target = bp.endSeqPos;
        return true;
    }
    begin = bp.startSeqPos;
    end = bp.endSeqPos;
    return false;
}
template <typename TBreakpoint, typename TSequence>
inline void _fillVcfRecordDuplication(VcfRecord & record, TBreakpoint & bp, TSequence & ref, unsigned id)
{
    typedef typename TBreakpoint::TPos TPos;
    TPos begin, end, target = std::numeric_limits<unsigned>::max();
    std::stringstream ss;
    if (!_setVcfRecordDuplicationPos(bp, begin, end, target))
        ss << "IMPRECISE;";

    if (begin > end)
        std::swap(begin, end);
    record.rID = id;
    record.beginPos = begin - 1; // Position before event
    record.filter = "PASS";
    ss << "SVTYPE=DUP";
    ss << ";END=" << end - 1 + 1; // -1 to set end as last position of variant, +1 for 1-base adjustment
    SEQAN_ASSERT_GEQ_MSG(bp.endSeqPos, bp.startSeqPos, "Duplication end position smaller than begin position!");
    if (target != std::numeric_limits<unsigned>::max())
    {
        ss << ";SVLEN=" << end-begin;
        ss << ";TARGETPOS=" << target - 1 + 1; // -1 to set target as position before insertion, +1 for 1-base adjustment
    }
    //if (bp.similar != std::numeric_limits<unsigned>::max())
        ss << ";BM=" << bp.similar;
    ss << ";DP=" << bp.support;
    record.info = ss.str();

    // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
    // deletion of length k.
    if (bp.startSeqPos != 0)
        appendValue(record.ref, ref[begin-1]);
    else
        appendValue(record.ref, ref[end]);

    // Compute ALT columns and a map to the ALT.
    append(record.alt, "<DUP>");

    // Create genotype infos.
    appendValue(record.genotypeInfos, "1");
}

template <typename TBreakpoint, typename TSequence>
inline void _fillVcfRecordBreakend(VcfRecord & record, TBreakpoint & bp, TSequence & ref, unsigned id)
{
    record.rID = id;
    record.beginPos = bp.startSeqPos - 1;
    record.filter = "PASS";
    std::stringstream ss;
    // ss << "IMPRECISE;";
    ss << "SVTYPE=BND";
    // ss << ";CIPOS=-5,5";
    //if (bp.similar != std::numeric_limits<unsigned>::max())
        //ss << ";BM=" << bp.similar;
    ss << ";DP=" << bp.support;
    record.info = ss.str();

    // Compute the number of bases in the REF column (1 in case of insertion and (k + 1) in the case of a
    // deletion of length k.
    // Breakend unlikely to be at position 0, and start equals end anyway...
    appendValue(record.ref, ref[bp.startSeqPos - 1]);

    // Compute ALT columns and a map to the ALT.
    // Storing if breakend is left of (0) or right of (1) breakpoint position
    if (bp.breakend) // breakend=1; right end breakend
    {
        appendValue(record.alt, record.ref[0]);
        appendValue(record.alt, '.');
    } else
    {
        appendValue(record.alt, '.');
        appendValue(record.alt, record.ref[0]);
    }

    // Create genotype infos.
    appendValue(record.genotypeInfos, "1");
}

template <typename TBreakpoint, typename TSequence>
inline bool _fillVcfRecord(VcfRecord & record, TBreakpoint & bp, TSequence & ref, unsigned id)
{
    switch(bp.svtype)
    {
        case TBreakpoint::INSERTION:
        _fillVcfRecordInsertion(record, bp, ref, id);
        return 1;
        case TBreakpoint::DELETION:
        _fillVcfRecordDeletion(record, bp, ref, id);
        return 1;
        case TBreakpoint::INVERSION:
        _fillVcfRecordInversion(record, bp, ref, id);
        return 1;
        case TBreakpoint::SEQAN_TANDEM:
        _fillVcfRecordTandem(record, bp, ref, id);
        return 1;
        case TBreakpoint::DISPDUPLICATION:
        _fillVcfRecordDuplication(record, bp, ref, id);
        return 1;
        // Translocation are handled seperately bc there is more than 1 record
        case TBreakpoint::INTERTRANSLOCATION:
        return 1;
        case TBreakpoint::TRANSLOCATION:
        return 1;
        case TBreakpoint::BREAKEND:
        _fillVcfRecordBreakend(record, bp, ref, id);
        return 1;

        default:
        return 0;
    }
    return 0;
}

template <typename TBreakpoint, typename TSequence>
inline bool _writeVcfTranslocation(VcfFileOut & vcfOut, TBreakpoint & bp, TSequence & ref, TSequence & ref2, unsigned id, unsigned id2, unsigned bp_id)
{
    // Translocation is encoded with 6 breakend entries (or 4 in the case where the "middle" position
    // of the translocation is not known)
    // The six entries mark new adjacencies induced by the translocation event:
    // 1st entry: position before first split (equals startSeqPosition - 1 of breakpoint) ALT: dupMiddlePos
    // 2nd: pos. after first split (startSeqPos) ALT: endSeqPos - 1
    // 3rd: pos. before second split (dupMiddlePos - 1) ALT: endSeqPos
    // 4th: pos. after second split (dupMiddlePos) ALT: startSeqPos - 1
    // 5th: pos. before third split (endSeqPos - 1) ALT: startSeqPos
    // 6th: pos. after third split (endSeqPos) ALT: dupMiddlePos - 1
    // Note that the positions have to be adjusted for 0-1 base change, only the record.beginPos is beeing adjusted
    // automatically
    typedef typename TBreakpoint::TId TId;
    TId sId;
    VcfRecord record;
    // Record values shared by (almost) all entries
    record.rID = id;
    record.filter = "PASS";
    std::stringstream ss;
    ss << "SVTYPE=BND";
    ss << ";EVENT=Trans" << bp_id;
    //if (bp.similar != std::numeric_limits<unsigned>::max())
        ss << ";BM=" << bp.similar;
    ss << ";DP=" << bp.support;
    record.info = ss.str();
    // Create genotype infos.
    appendValue(record.genotypeInfos, "1");

    // 1st entry
    record.id = "BND_" + toString(bp_id) + "_1";
    // No +1 adjustment
    // record.beginPos = (bp.startSeqPos > 0) ? bp.startSeqPos : 0;
    record.beginPos = bp.startSeqPos - 1;

    // Positions encoded in REF and ALT coloumn have to be already 1-based
    // whereas record.beginPos remains 0-based and is adjusted within writeVcfRecord() function
    // Compute base for REF and ALT columns
    if (bp.startSeqPos != 0)
    {
        appendValue(record.ref, ref[bp.startSeqPos - 1]);
        appendValue(record.alt, ref[bp.startSeqPos - 1]);
    } else
    {
        appendValue(record.ref, 'N');
        appendValue(record.alt, '.');
    }

    // Compute ALT columns
    if (bp.dupMiddlePos != std::numeric_limits<unsigned>::max())
    {
        _getShortId(sId, bp.midPosId);
        std::stringstream alt1;
        alt1 << "[";
        alt1 << sId; // alt1 << bp.midPosId;
        alt1 << ':';
        alt1 << bp.dupMiddlePos + 1; // 1-based adjustment
        alt1 << "[";
        append(record.alt, alt1.str());
    } else
        appendValue(record.alt, '.');

    // Write record and clear REF and ALT values
    try
    {
        writeRecord(vcfOut, record);
    }
    catch (seqan::IOError const & ioErr)
    {
        std::cerr << "Error while writing breakpoint translocation entry 1 vcf record!" << std::endl;
        return 1;
    }
    clear(record.ref);
    clear(record.alt);

    // 2nd entry
    // 2nd: pos. after first split (startSeqPos) ALT: endSeqPos - 1
    record.id = "BND_" + toString(bp_id) + "_2";
    // No 1-base adjustment
    record.beginPos = bp.startSeqPos;

    // Compute base for REF and ALT columns
    appendValue(record.ref, ref[bp.startSeqPos]);

    // Compute ALT columns
    if (bp.startSeqStrand != bp.endSeqStrand)
    {
        std::stringstream alt1;
        alt1 << "[";
        alt1 << bp.endSeqId;
        alt1 << ':';
        // -1 for position before split, -1 for end position adjustment, +1 for 1-base adjustment
        alt1 << bp.endSeqPos - 1 + 1; // 1-based adjustment, writing explicitely to keep overview
        alt1 << "[";
        append(record.alt, alt1.str());
    } else
    {
        std::stringstream alt1;
        alt1 << "]";
        alt1 << bp.endSeqId;
        alt1 << ':';
        alt1 << bp.endSeqPos - 1 + 1; // 1-based adjustment
        alt1 << "]";
        append(record.alt, alt1.str());
    }
    appendValue(record.alt, ref[bp.startSeqPos]);

    // Write record and clear REF and ALT values
    try
    {
        writeRecord(vcfOut, record);
    }
    catch (seqan::IOError const & ioErr)
    {
        std::cerr << "Error while writing breakpoint translocation entry 2 vcf record!" << std::endl;
        return 1;
    }
    clear(record.ref);
    clear(record.alt);

    // 3rd and 4th entry only exist if the middle (bp.dupMiddlePos) of translocation is known
    if (bp.dupMiddlePos != std::numeric_limits<unsigned>::max())
    {
        // 3rd entry
        // 3rd: pos. before second split (dupMiddlePos - 1) ALT: endSeqPos
        record.id = "BND_" + toString(bp_id) + "_3";
        // No 1-base adjustment
        record.beginPos = bp.dupMiddlePos - 1;

        // Compute base for REF and ALT columns
        appendValue(record.ref, ref[bp.dupMiddlePos - 1]);
        appendValue(record.alt, ref[bp.dupMiddlePos - 1]);

        // Compute ALT columns
        if (bp.startSeqStrand != bp.endSeqStrand)
        {
            std::stringstream alt1;
            alt1 << "]";
            alt1 << bp.endSeqId;
            alt1 << ':';
            alt1 << bp.endSeqPos + 1; // 1-based adjustment
            alt1 << "]";
            append(record.alt, alt1.str());
        } else
        {
            std::stringstream alt1;
            alt1 << "[";
            alt1 << bp.endSeqId;
            alt1 << ':';
            alt1 << bp.endSeqPos + 1; // 1-based adjustment
            alt1 << "[";
            append(record.alt, alt1.str());
        }

        // Write record and clear REF and ALT values
        try
        {
            writeRecord(vcfOut, record);
        }
        catch (seqan::IOError const & ioErr)
        {
            std::cerr << "Error while writing breakpoint translocation entry 3 vcf record!" << std::endl;
            return 1;
        }
        clear(record.ref);
        clear(record.alt);

        // 4th entry
        // 4th: pos. after second split (dupMiddlePos) ALT: startSeqPos - 1
        record.id = "BND_" + toString(bp_id) + "_4";
        record.beginPos = bp.dupMiddlePos;

        // Compute ALT columns
        std::stringstream alt4;
        alt4 << "]";
        alt4 << sId; // alt1 << bp.midPosId;
        alt4 << ':';
        if (bp.startSeqPos > 0)
            alt4 << bp.startSeqPos - 1 + 1; // 1-based adjustment
        else 
            alt4 << 0; // NO(!) 1-based adjustment
        alt4 << "]";
        append(record.alt, alt4.str());

        // Compute base for REF and ALT columns
        appendValue(record.ref, ref[bp.dupMiddlePos]);
        appendValue(record.alt, ref[bp.dupMiddlePos]);

        // Write record and clear REF and ALT values
        try
        {
            writeRecord(vcfOut, record);
        }
        catch (seqan::IOError const & ioErr)
        {
            std::cerr << "Error while writing breakpoint translocation entry 4 vcf record!" << std::endl;
            return 1;
        }
        clear(record.ref);
        clear(record.alt);

    }

    record.rID = id2;
    // 5th entry
    // 5th: pos. before third split (endSeqPos - 1) ALT: startSeqPos
    record.id = "BND_" + toString(bp_id) + "_5";
    record.beginPos = bp.endSeqPos - 1;

    // Compute base for REF and ALT columns
    appendValue(record.ref, ref2[bp.endSeqPos - 1]);
    appendValue(record.alt, ref2[bp.endSeqPos - 1]);

    // Compute ALT columns
    if (bp.startSeqStrand != bp.endSeqStrand)
    {
        std::stringstream alt1;
        alt1 << "]";
        alt1 << bp.startSeqId;
        alt1 << ':';
        alt1 << bp.startSeqPos + 1; // 1-based adjustment
        alt1 << "]";
        append(record.alt, alt1.str());
    } else
    {
        std::stringstream alt1;
        alt1 << "[";
        alt1 << bp.startSeqId;
        alt1 << ':';
        alt1 << bp.startSeqPos + 1; // 1-based adjustment
        alt1 << "[";
        append(record.alt, alt1.str());
    }

    // Write record and clear REF and ALT values
    try
    {
        writeRecord(vcfOut, record);
    }
    catch (seqan::IOError const & ioErr)
    {
        std::cerr << "Error while writing breakpoint translocation entry 5 vcf record!" << std::endl;
        return 1;
    }
    clear(record.ref);
    clear(record.alt);

    // 6th entry
    // 6th: pos. after third split (endSeqPos) ALT: dupMiddlePos - 1
    record.id = "BND_" + toString(bp_id) + "_6";
    record.beginPos = bp.endSeqPos;

    // Compute ALT columns
    if (bp.dupMiddlePos != std::numeric_limits<unsigned>::max())
    {
        std::stringstream alt1;
        if (bp.endSeqStrand != bp.midPosStrand)
            alt1 << "[";
        else
            alt1 << "]";
        alt1 << sId; // alt1 << bp.midPosId;
        alt1 << ':';
        alt1 << bp.dupMiddlePos - 1 + 1; // 1-based adjustment
        if (bp.endSeqStrand != bp.midPosStrand)
            alt1 << "[";
        else
            alt1 << "]";
        append(record.alt, alt1.str());
    } else
        appendValue(record.alt, '.');

    // Compute base for REF and ALT columns
    if (bp.endSeqPos < length(ref2) - 1)
    {
        appendValue(record.ref, ref2[bp.endSeqPos]);
        appendValue(record.alt, ref2[bp.endSeqPos]);
    } else
    {
        appendValue(record.ref, 'N');
        appendValue(record.alt, '.');
    }

    // Write record and clear REF and ALT values
    try
    {
        writeRecord(vcfOut, record);
    }
    catch (seqan::IOError const & ioErr)
    {
        std::cerr << "Error while writing breakpoint translocation entry 6 vcf record!" << std::endl;
        return 1;
    }
    clear(record.ref);
    clear(record.alt);

    return 0;
}

template <typename TSequence, typename TId>
void _fillVcfHeader(seqan::VcfHeader & vcfHeader,
                    seqan::VcfFileOut & vcfFileOut,
                    StringSet<TSequence> & databases,
                    StringSet<TId> & databaseIDs,
                    MSplazerOptions const & msplOpt)
{
    // Header record entries
    appendValue(vcfHeader, VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(vcfHeader, VcfHeaderRecord("source", "GUSTAF"));
    appendValue(vcfHeader, VcfHeaderRecord("reference", msplOpt.databaseFile));
    for (unsigned i = 0; i < length(msplOpt.queryFile); ++i)
        appendValue(vcfHeader, VcfHeaderRecord("reads", msplOpt.queryFile[i]));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
            "INFO", "<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
            "INFO", "<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record (for DUP last position of duplicated sequence)\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
            "INFO", "<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
            "INFO", "<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"));
    // appendValue(vcfHeader, seqan::VcfHeaderRecord(
    //            "INFO", "<ID=BND,Description=\"Breakend\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
            "INFO", "<ID=EVENT,Number=1,Type=String,Description=\"Event identifier for breakends.\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
            "INFO", "<ID=TARGETPOS,Number=1,Type=String,Description=\"Target position for duplications (position before inserted sequence).\">"));
    appendValue(vcfHeader, VcfHeaderRecord(
            "INFO", "<ID=BM,Number=1,Type=Integer,Description=\"Breakpoint Mate: similar breakpoints have same ID\">"));
    appendValue(vcfHeader, VcfHeaderRecord(
            "INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Number of Supporting Reads/Read Depth for Variant\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
            "ALT", "<ID=INV,Description=\"Inversion\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
                "ALT", "<ID=DUP,Description=\"Duplication\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord(
                "ALT", "<ID=DUP:TANDEM,Description=\"Tandem Duplication\">"));

    // Fill sequence names
    for (unsigned i = 0; i < length(databaseIDs); ++i)
    {
        seqan::CharString contigStr = "<ID=";
        append(contigStr, databaseIDs[i]);
        append(contigStr, ",length=");
        std::stringstream ss;
        ss << length(databases[i]);
        append(contigStr, ss.str());
        append(contigStr, ">");
        appendValue(vcfHeader, seqan::VcfHeaderRecord("contig", contigStr));

        appendName(contigNamesCache(context(vcfFileOut)), databaseIDs[i]);
    }
}

template <typename TId>
int32_t _getrID(StringSet<TId> & databaseIDs, TId dbID)
{
    for (unsigned i = 0; i < length(databaseIDs); ++i)
    {
        TId sID;
	_getShortId(sID, databaseIDs[i]);
        if (sID == dbID)
            return static_cast<int32_t>(i);
    }
    return std::numeric_limits<int>::max();
}

// Breakpoint writing call
template <typename TBreakpoint, typename TSequence, typename TId>
bool _writeGlobalBreakpoints(String<TBreakpoint> & globalBreakpoints,
                             StringSet<TSequence> & databases,
                             StringSet<TId> & databaseIDs,
                             MSplazerOptions const & msplazerOptions,
                             Vcf /*tag*/)
{
    // Creating GFF record and writing to gff file
    std::string fn_vcf = toString(msplazerOptions.vcfOutFile);
    seqan::VcfFileOut vcfOut;
    if (length(globalBreakpoints) == 0)
    {
        std::cerr << "Empty output list, skip writing vcf" << std::endl;
        return 1;
    }

    if (!open(vcfOut, fn_vcf.c_str()))
    {
        std::cerr << "Error while opening vcf breakpoint file!" << std::endl;
        return false;
    }
    seqan::VcfHeader vcfHeader;
    _fillVcfHeader(vcfHeader, vcfOut, databases, databaseIDs, msplazerOptions);
    try
    {
        writeHeader(vcfOut, vcfHeader);
    }
    catch (seqan::IOError const & ioErr)
    {
        std::cerr << "Error writing VCF record.\n";
        return false;
    }

    VcfRecord vcf_record;
    int32_t id = std::numeric_limits<int>::max();
    for (unsigned i = 0; i < length(globalBreakpoints); ++i)
    {
        TBreakpoint & bp = globalBreakpoints[i];
        if (bp.svtype != 0 && bp.support >= msplazerOptions.support) // 0=invalid
        {

            if (bp.svtype == TBreakpoint::DISPDUPLICATION && bp.translSuppStartPos && bp.translSuppEndPos)
                bp.svtype = TBreakpoint::TRANSLOCATION;
            id = _getrID(databaseIDs, bp.startSeqId);
            if (bp.svtype != 6 && bp.svtype != 7) // 6=inter-chr-translocation; 7=translocation
            {
                // Fill and write record
                if (_fillVcfRecord(vcf_record, bp, databases[id], id))
                    try
                    {
                        writeRecord(vcfOut, vcf_record);
                    }
                    catch (seqan::IOError const & ioErr)
                    {
                        std::cerr << "Error while writing breakpoint vcf record!" << std::endl;
                        return false;
                    }
                clear(vcf_record);
            } else
            {
                // extra write function because we have to write 6 records here instead of 1
                int32_t id2 = std::numeric_limits<int>::max();
                id2 = _getrID(databaseIDs, bp.endSeqId);
                if (_writeVcfTranslocation(vcfOut, bp, databases[id], databases[id2], id, id2, i))
                        std::cerr << "Error while writing breakpoint translocation vcf record!" << std::endl;

            }
        }
    }
    std::cout << " completed vcf writing " << msplazerOptions.vcfOutFile << std::endl;
    return 0;
}


// /////////////////////////////////////////////////////////////////////////////
// Writes parameters from options object to std::cout
template <typename TOptions>
void
_writeParams(TOptions & options)
{
// IOREV _notio_
    // Output user specified parameters
    // if (options.outDir != "")
    //  std::cout << "Output directory        : " << options.outDir << std::endl;
    if (options.jobName != "")
        std::cout << "Job name        : " << options.jobName << std::endl;
    std::cout << "Thresholds:" << std::endl;
    std::cout << "  overlap threshold (oth)          : " << options.simThresh << std::endl;
    std::cout << "  gap threshold (gth)              : " << options.gapThresh << std::endl;
    std::cout << "  initial gap threshold (ith)      : " << options.initGapThresh << std::endl;
    std::cout << "Penalties:" << std::endl;
    std::cout << "  translocation penalty (tp)       : " << options.diffDBPen << std::endl;
    std::cout << "  inversion penalty (ip)           : " << options.diffStrandPen << std::endl;
    std::cout << "  order penalty (op)               : " << options.diffOrderPen << std::endl;

    std::cout << "  required read support (st)       : " << options.support << std::endl;
}

// ////////////////////
// SAM Output
// TODO(ktrappe): Use SamBamIO, see tutorial
template <typename TSequence, typename TId>
void _writeSamFile(CharString outfile,
                   StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & stellarMatches,
                   StringSet<TSequence> & databases,
                   StringSet<CharString> databaseIDs,
                   StringSet<CharString> queryIDs)
{

    // CharString samOutFile2 = "msplazerStellarOut.sam";
    std::ofstream file2;
    file2.open(toCString(outfile), ::std::ios_base::out | ::std::ios_base::app);
    if (!file2.is_open())
    {
        std::cerr << "Could not open output file " << outfile << std::endl;
    }
    else
    {

        typedef StellarMatch<TSequence, TId> TMatch;
        typedef typename Iterator<String<TMatch> >::Type TIterator;
        typedef typename Row<typename TMatch::TAlign>::Type TMatchRow;
        typedef typename Size<TMatchRow>::Type TRowSize;
        TRowSize matchNumber, alignLen, alignDistance;

        //  std::cerr << "length: " << length(stellarMatches[0].matches);


        file2 << "@HD\tVN:1.0\tSO:coordinate\n";
        for (unsigned i = 0; i < length(databaseIDs); ++i)
            file2 << "@SQ\tSN:" << databaseIDs[i] << "\tLN:" << length(databases[i]) << "\n";
        // TODO(ktrappe): edit and flag CIGAR with clipped pos, mind strang direction
        for (unsigned i = 0; i < length(stellarMatches); ++i)
        {

            TIterator itStellarMatches = begin(stellarMatches[i].matches);
            TIterator itEndStellarMatches = end(stellarMatches[i].matches);

            if (itStellarMatches != itEndStellarMatches)
            {
                stellarMatches[i].lengthAdjustment = _computeLengthAdjustment(length(source((*itStellarMatches).row1)),
                                                                              length(source((*itStellarMatches).row2)));
            }

            //  std::cerr << "test itStellarMatches: " << *itStellarMatches
            // << " itEndStellarMatches: " << *itEndStellarMatches << std::endl;
            while (itStellarMatches < itEndStellarMatches)
            {

                // for(unsigned stellarmatches = 0; stellarmatches < length(matches[i].matches); ++stellarmatches){
                if ((*itStellarMatches).orientation)
                {
                    std::stringstream cigar, mutations;
                    _getCigarLine((*itStellarMatches).row1, (*itStellarMatches).row2, cigar, mutations);
                    //                                std::cout << "cigar string: " << cigar.str() << std::endl;
//                                    std::cout << "begin1: " << (*itStellarMatches).begin1 << std::endl;

                    // Alignment Score using _analyzeAlignment(row0, row1, alignLen, matchNumber) from stellar_output.h
                    // matchNumber and alignLen reset within _analyzeAlignment
                    _analyzeAlignment((*itStellarMatches).row1, (*itStellarMatches).row2, alignLen, matchNumber);
                    alignDistance = alignLen - matchNumber;
                    file2 << queryIDs[i] << "\t"
                          << "0\t"
                    // << databaseIDs[0] << "\t" //<< matches[i].matches[stellarmatches].id << "\t"
                          << (*itStellarMatches).id << "\t"
                    // 1-based position in DB/ref!
                          << "begin2: " << (*itStellarMatches).begin2 + 1 << "\t"
                          << "end2: " << (*itStellarMatches).end2 + 1 << "\t"
                          << "begin1: " << (*itStellarMatches).begin1 + 1 << "\t"
                          << "end1: " << (*itStellarMatches).end1 + 1 << "\t"
                          << ((*itStellarMatches).orientation ? '+' : '-') << "\t"
                          << "0\t"
                          << cigar.str() << "\t"
                          << "distance: " << alignDistance << "\t"
                          << "# matches: " << matchNumber << "\t"
                          << "alignLen: " << alignLen << "\t"
                    // << "*\t0\t0\t"
                    // << "MatchSeq" << "\t"
                    << "*\n";
                }
                ++itStellarMatches;                // better goNext()?
            }

            reverseComplement(databases);
            itStellarMatches = begin(stellarMatches[i].matches);

            while (itStellarMatches < itEndStellarMatches)
            {
                // for(unsigned stellarmatches = 0; stellarmatches < length(matches[i].matches); ++stellarmatches){

                if (!(*itStellarMatches).orientation)
                {

                    std::stringstream cigar, mutations;
                    _getCigarLine((*itStellarMatches).row1, (*itStellarMatches).row2, cigar, mutations);
                    //                                std::cout << "cigar string: " << cigar.str() << std::endl;
//                                    std::cout << "begin1: " << (*itStellarMatches).begin1 << std::endl;

                    // Alignment Score using _analyzeAlignment(row0, row1, alignLen, matchNumber) from stellar_output.h
                    // matchNumber and alignLen reset within _analyzeAlignment
                    _analyzeAlignment((*itStellarMatches).row1, (*itStellarMatches).row2, alignLen, matchNumber);
                    alignDistance = alignLen - matchNumber;
                    file2 << queryIDs[i] << "\t"
                          << "0\t"
                          << databaseIDs[0] << "\t"                   // << matches[i].matches[stellarmatches].id << "\t"
                          << (*itStellarMatches).id << "\t"
                          << "begin2: " << (*itStellarMatches).begin2 + 1 << "\t"                   // << matches[i].matches[stellarmatches].begin2 << "\t"
                          << "end2: " << (*itStellarMatches).end2 + 1 << "\t"                   // << matches[i].matches[stellarmatches].begin2 << "\t"
                          << "begin1: " << (*itStellarMatches).begin1 + 1 << "\t"                   // << matches[i].matches[stellarmatches].begin2 << "\t"
                          << "end1: " << (*itStellarMatches).end1 + 1 << "\t"
                          << ((*itStellarMatches).orientation ? '+' : '-') << "\t"
                          << "0\t"
                          << cigar.str() << "\t"
                          << "distance: " << alignDistance << "\t"
                          << "# matches: " << matchNumber << "\t"
                          << "alignLen: " << alignLen << "\t"
                          << "*\n";
                }
                ++itStellarMatches;                // better goNext()?
            }
            reverseComplement(databases);
            file2 << "\n";
        }
        file2.close();
    }
}

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_OUT_H_
