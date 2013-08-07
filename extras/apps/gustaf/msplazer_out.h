// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011-2013, Kathrin Trappe, FU Berlin
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

#ifndef SEQAN_EXTRAS_APPS_GUSTAF_MSPLAZER_OUT_H_
#define SEQAN_EXTRAS_APPS_GUSTAF_MSPLAZER_OUT_H_

#include <iostream>
#include <fstream>
#include <seqan/gff_io.h>
#include <seqan/file.h>

#include "../../../core/apps/stellar/stellar.h"
#include "msplazer.h"

using namespace seqan;

/**
.Tag.DotDrawingMSplazer
..cat:Input/Output
..summary:Switch to trigger drawing in dot format.
..value.DotDrawingMSplazer:Graphs in dot format.
..include:seqan/msplazer.h
*/

struct DotDrawingMSplazer_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazer;

/**
.Tag.DotDrawingMSplazerBestChain
..cat:Input/Output
..summary:Switch to trigger drawing in dot format.
..value.DotDrawingMSplazerBestChain:Best chain graphs in dot format.
..include:seqan/msplazer.h
*/

struct DotDrawingMSplazerBestChain_;
typedef Tag<DotDrawingMSplazer_> DotDrawingMSplazerBestchain;


/**
.Function.write:
..signature:write(file, msplazerchain, stellarmatches, tag)
..param.file:The file to write to.
...type:Class.Graph
..param.tag:A tag to select the output format.
...type:Tag.DotDrawingMSplazer
..include:seqan/msplazer.h
 */
template <typename TGraph, typename TVertexDescriptor, typename TScoreAlloc, typename TMatch, // typename TFile, 
          typename TBreakpoint, typename TPos, typename TMatchAlloc>
// typename TBreakpointAlloc, typename TMatchAlloc> // Requires Value<SparsePropertyMap> specialisation in msplazer.h
void
write(std::ostream & out, // TFile & file,  // std::ostream & out,  // std::fstream f; f.open(..); if (!f.good()) ... ; write(f, ...);  // write(std::cerr/std::cout, ...
      // MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, TBreakpointAlloc, // Requires Value<SparsePropertyMap> specialisation in msplazer.h
      MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, SparsePropertyMap<TBreakpoint, TPos>,
                    TMatchAlloc> const & msplazerchain,
      TMatch const & queryMatches,
      unsigned const & queryLength,
      DotDrawingMSplazer const &)
{
    // IOREV _doc_ _batchreading_
    SEQAN_CHECKPOINT
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
            out << queryMatches[i].end1 + 1;
            out << "  ";
            out << (queryMatches[i].orientation ? '+' : '-');
            out << "\\n read: ";
            out << queryMatches[i].begin2 + 1;
            out << "...";
            out << queryMatches[i].end2 + 1;
            out << "\\n";
            out << i;
            out << "\"];\n";
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
            out << queryLength + 1;
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
        << "}\n";
}

// DotWriting call for read graphs
template <typename TSequence, typename TId, typename TMSplazerChain>
void _writeDotfiles(StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & stellarMatches,
                    StringSet<TSequence> const & queries,
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
            write(f, queryChains[i], stellarMatches[i].matches, length(queries[i]), DotDrawingMSplazer());
            f.close();
        }
    }
}

template <typename TBreakpoint>
void _fillGffRecord(GffRecord & record, TBreakpoint & bp, unsigned id)
{
    typedef typename TBreakpoint::TId TId;
    TId sId;
    typedef typename TBreakpoint::TPos TPos;
    _getShortId(sId, bp.startSeqId);
    record.ref = sId;
    record.source = "GUSTAF";
    record.type = bp.svtype;
    record.beginPos = bp.startSeqPos;
    if (bp.svtype == "deletion" || bp.svtype == "inversion")
        record.endPos = bp.endSeqPos;
    else
        record.endPos = bp.startSeqPos + 1;
    if (bp.startSeqStrand)
        record.strand = '+';
    else
        record.strand = '-';

    appendValue(record.tagName, "ID");
    appendValue(record.tagValue, toString(id));
    if (bp.svtype == "insertion")
    {
        appendValue(record.tagName, "size");
        appendValue(record.tagValue, '-' + toString(length(bp.insertionSeq)));
        appendValue(record.tagName, "seq");
        appendValue(record.tagValue, bp.insertionSeq);
    }
    else if (bp.svtype == "deletion")
    {
        appendValue(record.tagName, "size");
        appendValue(record.tagValue, toString(static_cast<TPos>(bp.endSeqPos - bp.startSeqPos)));
    }
    else
    {
        appendValue(record.tagName, "endChr");
        _getShortId(sId, bp.endSeqId);
        appendValue(record.tagValue, sId);
        appendValue(record.tagName, "endPos");
        appendValue(record.tagValue, toString(bp.endSeqPos));
        appendValue(record.tagName, "endStrand");
        if (bp.endSeqStrand)
            appendValue(record.tagValue, '+');
        else
            appendValue(record.tagValue, '-');
    }
    appendValue(record.tagName, "support");
    appendValue(record.tagValue, toString(bp.support));
    appendValue(record.tagName, "supportIds");

    std::stringstream s;
    for (unsigned i = 0; i < length(bp.supportIds); ++i)
    {
        s << bp.supportIds[i];
        s << ',';
    }

    appendValue(record.tagValue, s.str());
    appendValue(record.tagName, "breakpoint");
    appendValue(record.tagValue, toString(bp.readStartPos + 1));
}

// Breakpoint writing call
template <typename TBreakpoint>
bool _writeGlobalBreakpoints(String<TBreakpoint> const & globalBreakpoints,
                             MSplazerOptions const & msplazerOptions)
{
    // Creating GFF record and writing to gff file
    std::string fn_gff = toString(msplazerOptions.breakpointOutFile);
    GffStream gffOut(fn_gff.c_str(), GffStream::WRITE);
    if (!isGood(gffOut))
        std::cerr << "Error while opening breakpoint file!" << std::endl;

    GffRecord gff_record;
    for (unsigned i = 0; i < length(globalBreakpoints); ++i)
    {
        if (globalBreakpoints[i].svtype != "none" && globalBreakpoints[i].support >= msplazerOptions.support)
        {
            // Fill record
            _fillGffRecord(gff_record, globalBreakpoints[i], i);
            // Write record
            if (writeRecord(gffOut, gff_record) != 0)
                std::cerr << "Error while writing breakpoint gff record!" << std::endl;
            clear(gff_record);
        }
    }

    std::cout << " completed writing " << msplazerOptions.breakpointOutFile << std::endl;
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
    std::cout << "  inital gap threshold (ith)       : " << options.initGapThresh << std::endl;
    std::cout << "Penalties:" << std::endl;
    std::cout << "  translocation penalty (tp)       : " << options.diffDBPen << std::endl;
    std::cout << "  inversion penalty (ip)           : " << options.diffStrandPen << std::endl;
    std::cout << "  order penalty (op)               : " << options.diffOrderPen << std::endl;

    std::cout << "  required read support (st)       : " << options.support << std::endl;
}

/**TODO(ktrappe): Add write functionality for SAM/VCF format
.Function.write:
..signature:write(file, msplazerchain, stellarmatches, nodeMap, edgeMap, tag)
..param.graph:The graph to write out.
...type:Class.Graph
..param.nodeMap:A mapping from vertex descriptor to vertex label.
..param.edgeMap:A mapping from edge descriptor to edge label.
..param.tag:A tag to select the output format.
...type:Tag.SamMSplazer
..include:seqan/msplazer.h
 */

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
