// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011, Knut Reinert, FU Berlin
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
#include <seqan/file.h>

#include "../../../core/apps/stellar/stellar.h"
#include "msplazer.h"

using namespace seqan;

/**
.Function.write:
..signature:write(file, msplazerchain, stellarmatches, tag)
..param.file:The file to write to.
...type:Class.Graph
..param.tag:A tag to select the output format.
...type:Tag.DotDrawingMSplazer
..include:seqan/msplazer.h
 */
template <typename TFile, typename TGraph, typename TVertexDescriptor, typename TScoreAlloc, typename TMatch,
          typename TBreakpoint, typename TPos, typename TMatchAlloc>
// typename TBreakpointAlloc, typename TMatchAlloc> // Requires Value<SparsePropertyMap> specialisation in msplazer.h
void
write(TFile & file,
      // MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, TBreakpointAlloc, // Requires Value<SparsePropertyMap> specialisation in msplazer.h
      MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, SparsePropertyMap<TBreakpoint, TPos>,
                    TMatchAlloc> const & msplazerchain,
      TMatch const & queryMatches,
      unsigned const & queryLength,
      DotDrawingMSplazer const &)
{
    // IOREV _doc_ _batchreading_
    SEQAN_CHECKPOINT
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    // typedef typename Value<TBreakpointAlloc>::Type TBreakpoint; // Requires Value<SparsePropertyMap> specialisation in msplazer.h
    typedef typename TBreakpoint::TId TId;

    // _writeGraphType(file,g,DotDrawing());
    _streamWrite(file, "digraph G {\n");
    _streamPut(file, '\n');
    _streamWrite(file, "/* Graph Attributes */\n");
    _streamWrite(file, "graph [rankdir = LR, clusterrank = local];\n");
    _streamPut(file, '\n');
    _streamWrite(file, "/* Node Attributes */\n");
    _streamWrite(file,
                 "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n");
    _streamPut(file, '\n');
    _streamWrite(file, "/* Edge Attributes */\n");
    _streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
    _streamPut(file, '\n');

    _streamWrite(file, "/* Nodes */\n");
    typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
    TConstIter it(msplazerchain.graph);
    unsigned i = 0;
    bool atEndV = false;
    for (; !atEnd(it); ++it)
    {
        TId sId;
        if (i < length(queryMatches))
        {
            _streamPutInt(file, *it);
            _streamWrite(file, " [label = \"");
            _streamWrite(file, "chr: ");
            _getShortId(sId, queryMatches[i].id);
            _streamWrite(file, sId);
            _streamWrite(file, "\\n db: ");
            _streamPutInt(file, queryMatches[i].begin1 + 1);
            _streamWrite(file, "...");
            _streamPutInt(file, queryMatches[i].end1 + 1);
            _streamWrite(file, "  ");
            _streamWrite(file, (queryMatches[i].orientation ? '+' : '-'));
            _streamWrite(file, "\\n read: ");
            _streamPutInt(file, queryMatches[i].begin2 + 1);
            _streamWrite(file, "...");
            _streamPutInt(file, queryMatches[i].end2 + 1);
            _streamWrite(file, "\\n");
            _streamPutInt(file, i);
            _streamWrite(file, "\"];\n");
        }
        else if (!atEndV)
        {
            _streamPutInt(file, *it);
            _streamWrite(file, " [label = \"start");
            _streamWrite(file, "\"];\n");
            atEndV = true;
        }
        else if (atEndV)
        {
            _streamPutInt(file, *it);
            _streamWrite(file, " [label = \"end");
            _streamWrite(file, "\\n");
            _streamPutInt(file, queryLength + 1);
            _streamWrite(file, "\"];\n");
        }
        else
            std::cerr << "in writing dot: default vertex??" << std::endl;
        ++i;
    }

    _streamPut(file, '\n');
    _streamWrite(file, "/* Edges */\n");
    typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
    TConstEdIter itEd(msplazerchain.graph);
    // bpC = 0;
    for (; !atEnd(itEd); ++itEd)
    {
        TVertexDescriptor sc = sourceVertex(itEd);
        TVertexDescriptor tr = targetVertex(itEd);
        _streamPutInt(file, sc);
        _writeEdgeType(file, msplazerchain.graph, DotDrawing());
        _streamPutInt(file, tr);
        _streamWrite(file, " [label = \"");
        _streamPutInt(file, getCargo(*itEd));
        TBreakpoint bp;
        bool foundBP = getProperty(msplazerchain.breakpoints, value(itEd), bp);
        if (foundBP)
        {
            _streamWrite(file, "*");
            // _streamPutInt(file, bpC);
            // ++bpC;
        }
        _streamWrite(file, "\"];\n");
    }
    _streamPut(file, '\n');
    _streamWrite(file, "}\n");
}

/**TODO(ktrappe)
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
            // if(length(stellarMatches[i].matches) > 45){
            // String< char > fn = toString(msplazerOptions.outDir) + "read" + toString(i+1) + '_'
            // + toString(msplazerOptions.queryFile) + '_' + toString(msplazerOptions.databaseFile) + ".dot";
            String<char> fn = "read" + toString(i + 1) + '_' + toString(msplazerOptions.jobName) + ".dot";
            // std::cerr  << fn << std::endl;
            FILE * strmWrite = fopen(toCString(fn), "w");
            // write(strmWrite, queryChains[i].graph, DotDrawing());
            write(strmWrite, queryChains[i], stellarMatches[i].matches, length(queries[i]), DotDrawingMSplazer());

            fclose(strmWrite);
            // std::cerr << " completed dot write in: " << i << std::endl;
        }
    }
}

// Breakpoint writing call
template <typename TBreakpoint>
bool _writeGlobalBreakpoints(String<TBreakpoint> const & globalBreakpoints,
                             MSplazerOptions const & msplazerOptions,
                             unsigned const & support)
{
    // String<char> fn = toString(msplazerOptions.outDir) + toString(msplazerOptions.jobName) + toString(fileName);
    String<char> fn = toString(msplazerOptions.breakpointOutFile);
    // std::cerr  << fn << std::endl;
    FILE * strmWrite = fopen(toCString(fn), "w");

    // std::ofstream file;
    /*file.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::app);
    if (!strmWrite.is_open()) {
        std::cerr << "Could not open output file." << std::endl;
        return 1;
    }

    if (!open(f, filename, "rb"))
    {
        std::cerr << "ERROR: GZip file has the wrong format!" << std::endl;
        return 1;
    }
    */
    /*
    _streamWrite(strmWrite, "Global breakpoints found on best Chains sorted according to genome Ids\n");
    _streamWrite(strmWrite, "Database file: ");
    _streamWrite(strmWrite, msplazerOptions.databaseFile);
    _streamPut(strmWrite, '\n');
    */
    // _streamPut(strmWrite, BreakpointFileHeader());
    // print bps
    for (unsigned i = 0; i < length(globalBreakpoints); ++i)
    {
        // if(globalBreakpoints[i].svtype != "none" && globalBreakpoints[i].support >= msplazerOptions.support)
        if (globalBreakpoints[i].svtype != "none" && globalBreakpoints[i].support >= support)
            _streamWrite(strmWrite, globalBreakpoints[i], i);
        // _streamPut(strmWrite, '\n');
    }
    fclose(strmWrite);
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

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_OUT_H_