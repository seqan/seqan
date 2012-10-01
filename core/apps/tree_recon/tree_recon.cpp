/*==========================================================================
               SeqAn - The Library for Sequence Analysis
                         http://www.seqan.de 
============================================================================
Copyright (C) 2007

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
==========================================================================*/

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <fstream>


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TMat, typename TNames>
inline void 
_readPhylipMatrix(TFile& file,
				  TMat& matrix,
				  TNames& names)
{
//IOREV can probably stay here since its too unimportant for the rest of seqan
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TMat>::Type TDistance;
	typedef typename Size<TMat>::Type TSize;
	typedef typename Value<TNames>::Type TName;
	typedef typename Iterator<TMat, Standard>::Type TMatIter;

	// Parse the file and convert the internal ids
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;
		_parseSkipWhitespace(file, c);
		TSize nseq = _parseReadNumber(file, c);
		_parseSkipWhitespace(file, c);
		resize(matrix, nseq * nseq);
		resize(names, nseq);
		TMatIter it = begin(matrix, Standard());
		for(TSize row = 0; row<nseq; ++row) {
			_parseReadIdentifier(file, names[row], c);
			_parseSkipWhitespace(file, c);
			for(TSize col = 0; col<nseq; ++col, ++it) {
				*it = _parseReadDouble(file, c);
				_parseSkipWhitespace(file, c);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[])
{
    // -----------------------------------------------------------------------
    // Setup argument parser
    // -----------------------------------------------------------------------
    seqan::ArgumentParser parser("tree_recon");

    // Set short description, version, date.
    setShortDescription(parser, "Tree reconstruction");
    setVersion(parser, "1.02");
    setDate(parser, "July 17, 2012");

    // Usage line and description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-m\\fP \\fIIN.DIST\\fP");
    addDescription(parser, "Reconstruct phylogenetic tree from Phylip matrix \\fIIN.DIST\\fP.");

	addSection(parser, "Input / Output");
    addOption(parser, seqan::ArgParseOption("m", "matrix", "Name Phylip distance matrix file.  Must contain at least three species.", seqan::ArgParseArgument::INPUTFILE, "FILE"));
    setRequired(parser, "matrix");
	addOption(parser, seqan::ArgParseOption("o", "out-file", "Path to write output to.", seqan::ArgParseArgument::OUTPUTFILE, "FILE"));
    setDefaultValue(parser, "out-file", "tree.dot");
	addOption(parser, seqan::ArgParseOption("f", "format", "The output format.", seqan::ArgParseArgument::STRING, "FORMAT"));
    setValidValues(parser, "format", "dot newick");
    setDefaultValue(parser, "format", "dot");

    addSection(parser, "Algorithm Options");
    addOption(parser, seqan::ArgParseOption("b", "build", "Tree building method. \\fInj\\fP: neighbour-joining, \\fImin\\fP: UPGMA single linkage, \\fImax\\fP: UPGMA complete linkage, \\fIavg\\fP: UPGMA average linkage, \\fIwavg\\fP: UPGMA weighted average linkage.  Neighbour-joining creates an unrooted tree.  We root that tree at the least joined pair.", seqan::ArgParseArgument::STRING, "METHOD"));
    setValidValues(parser, "build", "nj min max avg wavg");
    setDefaultValue(parser, "build", "nj");

    addTextSection(parser, "Contact and References");
    addListItem(parser, "For questions or comments, contact:", "Tobias Rausch <rausch@embl.de>");
    addListItem(parser, "SeqAn Homepage:", "http://www.seqan.de");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

	// Tree reconstruction
	typedef double TDistanceValue;
	typedef String<char> TName;
	typedef Size<TName>::Type TSize;

	// Read the options	
	::std::string infile;
    getOptionValue(infile, parser, "matrix");
	::std::string outfile;
    getOptionValue(outfile, parser, "out-file");
	TSize build = 0;
	String<char> meth;
    getOptionValue(meth, parser, "build");
	if (meth == "nj") build = 0;
	else if (meth == "min") build = 1;
	else if (meth == "max") build = 2;
	else if (meth == "avg") build = 3;
	else if (meth == "wavg") build = 4;
	String<char> format;
	getOptionValue(format, parser, "format");

	// Read the distance matrix
	String<TName> names;
	String<TDistanceValue> matrix;
	FILE* strmMat = fopen(infile.c_str(), "rb");
    if (strmMat == 0)
    {
        std::cerr << "Could not open file " << infile.c_str() << std::endl;
        return 1;
    }
	_readPhylipMatrix(strmMat, matrix, names);	
	fclose(strmMat);

	// Create the tree
	Graph<Tree<TDistanceValue> > tree;
	if (build == 0) njTree(matrix, tree);
	else if (build == 1) upgmaTree(matrix, tree, UpgmaMin());
	else if (build == 2) upgmaTree(matrix, tree, UpgmaMax());
	else if (build == 3) upgmaTree(matrix, tree, UpgmaAvg());
	else if (build == 4) upgmaTree(matrix, tree, UpgmaWeightAvg());
	
	if (format == "dot") {
		TSize nameLen = length(names);
		resize(names, numVertices(tree));
		// Add the label prefix for leaves
		for(TSize i = 0;i < nameLen; ++i) {
			TName tmpName = "label = \"";
			append(tmpName, names[i], Generous());
			append(tmpName, '"');
			names[i] = tmpName;
		}
		// Append emty names for internal vertices
		for(;nameLen < length(names); ++nameLen) {
			names[nameLen] = "label = \"\"";
		}

		// Write the result
		FILE* strmDot;
		strmDot = fopen(outfile.c_str(), "w");
		write(strmDot, tree, names, DotDrawing());
		fclose(strmDot);
	} else if (format == "newick") {
		FILE* strmDot;
		strmDot = fopen(outfile.c_str(), "w");
		// If nj tree collapse the root
		if (build == 0) write(strmDot, tree, names, true, NewickFormat());
		else write(strmDot, tree, names, false, NewickFormat());
		fclose(strmDot);
	}
	return 0;
}
