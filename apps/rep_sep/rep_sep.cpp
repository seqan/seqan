/*=========================================================================
  Copyright (C) 2009 by Stephan Aiche

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================
  $Id$
 ==========================================================================*/

#define SEQAN_PROFILE

#include <iostream>
#include <sstream>
#include <fstream>

// load command line parser
#include <seqan/arg_parse.h>

#include <seqan/store.h>

#include "utils.h"
#include "assembly_parser.h"
#include "column_scanner.h"
#include "rgraph.h"

using namespace seqan;
using namespace std;

struct RepSepOptions
{
    CharString assembly;
    int contig;
    int copyNumber;
    bool noClean;
    bool createDotFile;
    CharString outputPrefix;
    //bool useDNP;
    //bool useSimpleColumns;
    double error;
    bool hmce;
    bool hsce;

    RepSepOptions() :
        assembly(""),
        contig(0),
        copyNumber(2),
        noClean(false),
        createDotFile(false),
        outputPrefix(""),
        //useDNP(true),
        //useSimpleColumns(false),
        error(0.0),
        hmce(true),
        hsce(false)
    {}
};

seqan::ArgumentParser::ParseResult
parseCommandLine(RepSepOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser;

    addUsageLine(parser, "[OPTION]... --assembly <input file> --output-prefix <prefix>");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setShortDescription(parser, "Repeat Separation Tool -- Copyright (c) 2009, Stephan Aiche");

    // needed input file
    addOption(parser, ArgParseOption("a", "assembly", "Input assembly filename.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "assembly", "afg");
    setRequired(parser, "assembly");

    addOption(parser, ArgParseOption("c", "contig", "Index of the contig in the assembly that should be analyzed (NOTE: the index is 0 based).", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "contig", 0);
    setMinValue(parser, "contig", "0");

    addOption(parser, ArgParseOption("n", "copy-number", "Number of compressed repeat copies in the given contig.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "copy-number", 2);
    setMinValue(parser, "copy-number", "2");

    addOption(parser, ArgParseOption("", "no-clean", "Disable automatic graph cleaning."));
    addOption(parser, ArgParseOption("", "dotfile", "Write constructed graph as dotfile to visualize in Graphviz."));

    addOption(parser, ArgParseOption("p", "output-prefix", "Filename prefix for the result files. Files for the ILP and the Result will be named PREFIX (lp|rs|dot|heu).", ArgParseArgument::STRING, "PREFIX"));
    setRequired(parser, "output-prefix");

    addSection(parser, "Column Detection Strategy");
    //addOption(parser, ArgParseOption("d", "dnp", "Use DNP strategy for column detection [DEFAULT]"));
    //addOption(parser, ArgParseOption("s", "simple-columns", "Use the simple (normal distributed) column detection strategy"));
    addOption(parser, ArgParseOption("e", "error", "Expected sequencing error.", ArgParseArgument::DOUBLE, "ERROR"));

    addSection(parser, "Heuristic selection");
    addOption(parser, ArgParseOption("", "hmce", "Solve the problem with the multi-component-expansion (mce) heuristic [DEFAULT]."));
    addOption(parser, ArgParseOption("", "hsce", "Solve the problem with the single-component-expansion (sce) heuristic."));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // extract option values
    getOptionValue(options.assembly, parser, "assembly");

    getOptionValue(options.contig, parser, "contig");
    getOptionValue(options.copyNumber, parser, "copy-number");

    options.noClean = isSet(parser, "no-clean");
    options.createDotFile = isSet(parser, "dotfile");

    getOptionValue(options.outputPrefix, parser, "output-prefix");

    //options.useDNP = isSet(parser, "dnp");
    //options.useSimpleColumns = isSet(parser, "simple-columns") && !options.useDNP; // we prefer useDNP
    getOptionValue(options.error, parser, "error");

    options.hmce = isSet(parser, "hmce");
    options.hsce = isSet(parser, "hsce") && !options.hmce; // we prefer hmce

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, const char * argv[])
{
    RepSepOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // run starts
    SEQAN_PROTIMESTART(profileTime);

    typedef FragmentStore<> TFragmentStore;
    typedef Size<TFragmentStore>::Type TSize;
    TFragmentStore fragStore;

    cout << "loading data from " << options.assembly << " .. " << flush;

    TSize numberOfContigs = 0;
    // TSize numberOfReads = 0;

    ifstream strmReads;
    if (!open(strmReads, toCString(options.assembly)))
    {
        cout << "Could not open " << options.assembly << endl;
        return 1;
    }
    read(fragStore, strmReads, Amos());
    close(strmReads);
    numberOfContigs = length(fragStore.contigStore);
    // numberOfReads = length(fragStore.readStore);

    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;

    /*
    cout << endl << "------------------------------- " << endl;
    cout << "stats: " << endl;
    cout << "number of contigs: " << numberOfContigs << endl;
    cout << "number of reads:   " << numberOfReads << endl;
    cout << "------------------------------- " << endl;
    */
    if (!(static_cast<TSize>(options.contig) < numberOfContigs))
    {
        cout << "You have selected an invalid contig! Only " << numberOfContigs << " different contig(s) were found" << endl;
        cout << "in the currently assembly" << endl;
        return 1;
    }

    cout << "#INFO: you have selected contig nr. " << options.contig << endl;

    // construct matrix for parsing
    typedef char TAlpahbet;
    typedef Id<TFragmentStore::TAlignedReadStore>::Type TId;
    typedef Size<TFragmentStore::TReadSeq>::Type TReadPos;

    typedef Triple<TAlpahbet, TId, TReadPos> TMatrixValue;
    typedef String<TMatrixValue> TMatrix;

    TMatrix matrix;
    typedef Triple<char, TId, TReadPos> TColumnAlphabet;
    typedef String<TColumnAlphabet> TCandidateColumn;
    typedef Pair<TReadPos, TCandidateColumn> TAnnotatedCandidateColumn;
    String<TAnnotatedCandidateColumn> candidates;

    SEQAN_PROTIMEUPDATE(profileTime);
    cout << "parsing contig for candidate columns .. " << flush;

    // TODO: add code for selection of scanner
    SimpleColumn algoSpec;
    algoSpec.parameters.error_probability = options.error;
    parseContig(fragStore, options.contig, candidates, algoSpec);

    for (TSize x = 0; x < length(candidates); ++x)
    {
        cout << candidates[x].i1 << endl;
    }

    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;
    cout << "#INFO: <" << length(candidates) << "> possible candidate columns were identified" << endl;


    ReadGraph<TColumnAlphabet, Value<TFragmentStore::TAlignedReadStore>::Type, TReadPos> rgraph;

    cout << "adding candidates to graph .. " << flush;
    construct(rgraph, candidates, fragStore, options.contig);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;

    cout << "adding mate pairs to graph .. " << flush;
    add_mate_pairs(rgraph, fragStore, options.contig);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;

    GraphScoring scoring_scheme;

    cout << "compute scores for graph edges .. " << flush;
    scoreGraph_(rgraph, fragStore, static_cast<TSize>(options.contig), scoring_scheme);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;

    if (!options.noClean)
    {
        cout << "cleaning graph .. " << flush;
        // TODO: cleaning needs to be implemented
        cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;
    }
    else
    {
        cout << "#INFO: skipped graph cleaning" << endl;
    }

    cout << "analyzing constructed graph .. " << flush;
    bool hasMultiComp = hasMultipleComponents(rgraph);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;

    if (hasMultiComp)
    {
        cout << "This graph has multiple components!!" << endl;
    }

    // TODO: add some statistics for the graph
    cout << endl << "################ graph info ################" << endl;

    cout << " |edges| " << numEdges(rgraph) << endl;
    cout << " |vertices| " << numVertices(rgraph) << endl;

    cout << "############################################" << endl << endl << flush;

    // TODO: implement heurisitics
    typedef SelectGraph_<ReadGraph<TColumnAlphabet, Value<TFragmentStore::TAlignedReadStore>::Type, TReadPos> >::Type TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    typedef String<TVertexDescriptor> TComponent;
    typedef String<TComponent> TComponentList;

    TComponentList components;

    if (options.hmce)
    {
        cout << "separate the graph using (GuidedParallelComponentMerge) .. " << flush;
        solve(rgraph, components, options.copyNumber, GuidedParallelComponentMerge());
        cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;
    }
    else
    {
        cout << "separate the graph using (SingleComponentExpansion) .. " << flush;
        solve(rgraph, components, options.copyNumber, SingleComponentExpansion());
        cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;
    }

    for (TSize c = 0; c < length(components); ++c)
    {
        cout << "Component <" << c << "> contains " << length(components[c]) << " reads " << endl;
        for (TSize r = 0; r < length(components[c]); ++r)
        {
            cout << value(rgraph.vertexCargo, components[c][r]).alignedRead.readId << " (" << fragStore.readNameStore[value(rgraph.vertexCargo, value(components[c], r)).alignedRead.readId] << ")" << endl;
        }
    }
    // TODO: implement IO stuff
}
