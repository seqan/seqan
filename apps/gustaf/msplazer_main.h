// ==========================================================================
//                                    Gustaf
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

#ifndef SEQAN_APPS_GUSTAF_MSPLAZER_MAIN_H_
#define SEQAN_APPS_GUSTAF_MSPLAZER_MAIN_H_


#include "../stellar/stellar.h"
#include "../stellar/stellar_output.h"

#include "msplazer.h"
#include "msplazer_out.h"
#include "msplazer_algorithms.h"
#include "gustaf_matepairs.h"
#include "stellar_routines.h"
#include "create_stellarmatches_from_file.h"
using namespace seqan;

// /////////////////////////////////////////////////////////////////////////////
// MSplazer Wrapper
int msplazer(StellarOptions & stellarOptions, MSplazerOptions & msplazerOptions)
{
    // Finder
    typedef String<Dna5> TSequence;
    //  typedef FragmentStore<void>::TReadSeqStore TSequence;

    // Database and query ID

    typedef CharString TId;

    // import query sequences using _importSequences from Stellar
    StringSet<TSequence> queries;
    StringSet<TId> queryIDs;
    String<unsigned> readJoinPositions;
    // TODO (ktrappe) distinguish between paired and single and, call appropriate
    // importSeq function (and preprocess query files)
    if (msplazerOptions.pairedEndMode)
    {
        std::cout << "Loading paired-end read sequences... ";
        if (!_importSequences(msplazerOptions.queryFile[0], msplazerOptions.queryFile[1], msplazerOptions.revCompl, queries, queryIDs, readJoinPositions))
            return 1;
    }else
    {
        std::cout << "Loading query sequences... ";
        if (!_importSequences(stellarOptions.queryFile, "query", queries, queryIDs))
            return 1;
    }
    StringSet<TId> shortQueryIDs = queryIDs;

    /*
    unsigned readLength = 0;
    for(unsigned i = 0; i < length(queries); ++i)
        readLength += length(queries[i]);
        */
    // std::cerr << "Loaded read seq: " << queries[i] << std::endl;
    std::cout << "done" << std::endl;

    // import database sequence using _importSequences from Stellar
    StringSet<TSequence> databases;
    StringSet<TId> databaseIDs;
    std::cout << "Loading reference sequences... ";
    if (!_importSequences(stellarOptions.databaseFile, "database", databases, databaseIDs))
        return 1;

    StringSet<TId> shortDatabaseIDs = databaseIDs;

    for (unsigned i = 0; i < length(databases); ++i)
    {
        std::cout << "Loaded db seq with length: " << length(databases[i]) << std::endl;
        std::cout << "Loaded db ID: " << databaseIDs[i] << std::endl;
    }

    std::cout << "done" << std::endl;


    // /////////////////////////////////////////////////////////////////////////
    // Compute Stellar Matches and their score
    // Note: Matches will be sorted when calling score function _getMatchDistanceScore
    // Note: Matches of the reverse strand will be modified in the sense that the match positions refer to the forward
    // strand (and not the reverse strand) get distance scores for all matches

    // Stellar Match Space Allocator
    typedef StringSet<QueryMatches<StellarMatch<TSequence, TId> > > TMatches;
    // FragmentStore
    /*
       TMatches matches;
       resize(matches, length(fragments.readSeqStore));
       */

    TMatches stellarMatches;
    TMatches testStMatches;
    resize(stellarMatches, length(queries));

    // Score space allocator
    typedef StellarMatch<TSequence, TId> TMatch;
    typedef String<int> TScoreAlloc;

    String<TScoreAlloc> distanceScores;
    resize(distanceScores, length(stellarMatches));

    // get Stellar matches
    bool doStellar = true;
    if (msplazerOptions.stellarInputFile != "")
    {
        doStellar = false;
        std::cout << " Getting STELLAR matches from file, not calling STELLAR" << std::endl;
    }
    if (doStellar)
    {
        std::cout << "Calling STELLAR..." << std::endl;
        std::cout << "Stellar options:" << std::endl;
        _writeFileNames(stellarOptions);
        _writeSpecifiedParams(stellarOptions);
        _writeCalculatedParams(stellarOptions);
        _getStellarMatches(queries, databases, databaseIDs, stellarOptions, stellarMatches);
        std::cout << "done" << std::endl;
    }
    else
    {
        std::cout << "Importing STELLAR matches from file " << msplazerOptions.stellarInputFile << std::endl;
        double startST = sysTime();
        // TODO (ktrappe) distinguish call with queryIDs and shortQueryIDs in case of mate pairs? stellar writes out short
        // query IDs anyway...
        if (!_getStellarMatchesFromFile(queries, shortQueryIDs, databases, databaseIDs, msplazerOptions.stellarInputFile,
                                        stellarMatches, msplazerOptions.numThreads))
            return 1;

        std::cout << "done" << std::endl;
        std::cout << "TIME importing stellar matches " <<   (sysTime() - startST) << "s" << std::endl;
    }
    /*
    for(unsigned i = 0; i < length(queries); ++i){
        for(unsigned j = 0; j < length(stellarMatches[i].matches); ++j)
            std::cout << stellarMatches[i].matches[j] << std::endl;
    }
    */
    double startDist = sysTime();
    std::cout << "Getting match distance..." << std::endl;
    _getMatchDistanceScore(stellarMatches, distanceScores, msplazerOptions.numThreads);
    std::cout << "TIME getting match distance " <<   (sysTime() - startDist) << "s" << std::endl;


    // Graph statistics
    /*
    unsigned maxSize = 0;
    for(unsigned i = 0; i < length(queries); ++i){
        maxSize = std::max(static_cast<unsigned>(length(stellarMatches[i].matches)), maxSize);
        for(unsigned j = 0; j < length(stellarMatches[i].matches); ++j)
            std::cout << stellarMatches[i].matches[j] << std::endl;
    }
    String<unsigned> sizeCount;
    String<unsigned> chainCount;
    String<unsigned> partialChainCount;
    resize(sizeCount, maxSize+1);
    resize(chainCount, maxSize+1);
    resize(partialChainCount, maxSize+1);
    for(unsigned i = 0; i < length(sizeCount); ++i){
        sizeCount[i] = 0;
        chainCount[i] = 0;
        partialChainCount[i] = 0;
    }
    for(unsigned i = 0; i < length(queries); ++i){
        ++sizeCount[length(stellarMatches[i].matches)];
    }
    for(unsigned i = 0; i < length(sizeCount); ++ i){
        std::cout << "Number of occurrences for matches " << i << " : " << sizeCount[i] << std::endl;
    }
    for(unsigned i = 0; i < length(distanceScores); ++i){
        TScoreAlloc scores = distanceScores[i];
        std::cerr << "read score: " << i << std::endl;
        for(unsigned j = 0; j < length(scores); ++j){
            std::cerr << "score: " << scores[j] << std::endl;
        }
        std::cerr << std::endl;
    }
    */

    // ///////////////////////////////////////////////////////////////////////
    // Create chains/graphs for all queries using Stellar matches


    // Graph structure for stellar matches
    typedef Graph<Directed<int> > TGraph; // TRowSize> > TGraph;
    // static_cast<Nothing>(TGraph());
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    // Breakpoint property map
    typedef SparsePropertyMap<Breakpoint<TSequence, TId>, unsigned> TSparsePropertyMap;
    typedef String<TMatch> TMatchAlloc;

    // Container for graph/chain, scores and start and end vertices
    typedef MSplazerChain<TGraph, TVertexDescriptor, TScoreAlloc, TSparsePropertyMap, TMatchAlloc> TMSplazerChain;
    String<TMSplazerChain> queryChains;

    std::cout << "Constructing graphs... ";
    // TODO distinguish call with queryIDs and shortQueryIDs for mate pairs?
    double startGraphs = sysTime();
    std::cout << "Building graphs... ";
    _chainQueryMatches(stellarMatches, distanceScores, queryChains, queryIDs, queries, readJoinPositions, msplazerOptions);
    std::cout << "...done... " << (sysTime() - startGraphs) << "s"  << std::endl;

    // ///////////////////////////////////////////////////////////////////////
    // Analyze chains

    double startAnalysis = sysTime();
    std::cout << "Analyzing graphs... ";
    _analyzeChains(queryChains);
    std::cout << "...done... " << (sysTime() - startAnalysis) << "s" << std::endl;

    // ///////////////////////////////////////////////////////////////////////
    // Breakpoints

    typedef Breakpoint<TSequence, TId> TBreakpoint;
    String<TBreakpoint> globalBreakpoints;
    // String<TBreakpoint> globalStellarIndels;
    double startBP = sysTime();
    std::cout << "Extracting breakpoints... ";
    _findAllBestChains(queryChains, stellarMatches, globalBreakpoints, msplazerOptions);
    std::cout << "...done... " << (sysTime() - startBP) << "s" << std::endl;
    // _findAllBestChains(queryChains, stellarMatches, queries, queryIDs, globalBreakpoints, globalStellarIndels, msplazerOptions);
    // _findAllChains(queryChains, stellarMatches, queries, queryIDs, globalBreakpoints, globalStellarIndels, msplazerOptions);
    // _findAllChains(queryChains);

    // Graph statistics output
    /*
    for(unsigned i = 0; i < length(queryChains); ++i)
        if(!queryChains[i].isEmpty)
            _findPartialChains(queryChains[i]);

    for(unsigned i = 0; i < length(queryChains); ++i){
        if(!queryChains[i].isEmpty && !queryChains[i].isPartial)
            ++chainCount[length(stellarMatches[i].matches)];
        else if(queryChains[i].isPartial)
            ++partialChainCount[length(stellarMatches[i].matches)];
    }

    for(unsigned i = 0; i < length(chainCount); ++ i){
        std::cout << "Number of complete chains for matches " << i << " : " << chainCount[i] << std::endl;
    }
    for(unsigned i = 0; i < length(partialChainCount); ++ i){
        std::cout << "Number of partial chains for matches " << i << " : " << partialChainCount[i] << std::endl;
    }
    */


    // std sort in ascending order
    double startWriting = sysTime();
    std::cout << "Sorting and writing breakpoints... ";
    std::sort(begin(globalBreakpoints), end(globalBreakpoints));
    // std::sort(begin(globalStellarIndels), end(globalStellarIndels));
    _writeGlobalBreakpoints(globalBreakpoints, msplazerOptions, Gff());
    _writeGlobalBreakpoints(globalBreakpoints, databases, databaseIDs, msplazerOptions, Vcf());
    std::cout << "...done " << (sysTime() - startWriting) << "s" << std::endl;
    // _writeGlobalBreakpoints(globalStellarIndels, msplazerOptions, msplazerOptions.support);

    // ///////////////////////////////////////////////////////////////////////
    // SAM Output

    // SAM out of Stellar matches
    // String <char> samOutfile = toString(msplazerOptions.outDir) + toString(msplazerOptions.samOutFile);//"example/chr19.sam";
    // _writeSamFile(toCString(samOutfile), stellarMatches, databases, databaseIDs, queryIDs);
    // std::cerr << " completed sam " << std::endl;



    // ///////////////////////////////////////////////////////////////////////
    // Write dot files
    //
    // if(length(queries) < 200)
    if (msplazerOptions.dotOut)
        _writeDotfiles(stellarMatches, queries, queryIDs, queryChains, msplazerOptions);
    return 0;
}

#endif  // #ifndef SEQAN_APPS_GUSTAF_MSPLAZER_MAIN_H_
