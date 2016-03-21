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
#include <seqan/stream.h>

#include <iostream>


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

template <typename TMat, typename TNames, typename TForwardIter>
void readPhylipMatrix(TMat & matrix,
                      TNames & names,
                      TForwardIter & iter)
{
    //typedef typename Value<TFile>::Type TValue;
    typedef typename Value<TMat>::Type TFloat;
    typedef typename Size<TMat>::Type TSize;
    //typedef typename Value<TNames>::Type TName;
    typedef typename Iterator<TMat, Standard>::Type TMatIter;

    if (atEnd(iter))
        throw UnexpectedEnd();

    // Parse the file and convert the internal ids.
    CharString buffer;
    while (!atEnd(iter))
    {
        clear(buffer);

        // Read "<blanks><num>\n"
        skipUntil(iter, NotFunctor<IsBlank>());
        readUntil(buffer, iter, OrFunctor<NotFunctor<IsGraph>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Nothing> >());
        TSize nseq = lexicalCast<TSize>(buffer);
        skipLine(iter);

        // Allocate memory in matrix and names.
        resize(matrix, nseq * nseq);
        resize(names, nseq);

        TMatIter it = begin(matrix, Standard());
        for (TSize row = 0; row < nseq; ++row)
        {
            readUntil(names[row], iter, NotFunctor<IsGraph>());
            skipUntil(iter, OrFunctor<NotFunctor<IsSpace>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Nothing> >());

            for (TSize col = 0; col < nseq; ++col, ++it)
            {
                clear(buffer);
                readUntil(buffer, iter, NotFunctor<IsGraph>());
                *it = lexicalCast<TFloat>(buffer);

                if (col + 1 < nseq)
                    skipUntil(iter, OrFunctor<NotFunctor<IsSpace>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Nothing> >());
                else
                    skipLine(iter);
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
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "Phylogeneny");

    // Usage line and description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-m\\fP \\fIIN.DIST\\fP");
    addDescription(parser, "Reconstruct phylogenetic tree from Phylip matrix \\fIIN.DIST\\fP.");

    addSection(parser, "Input / Output");
    addOption(parser, seqan::ArgParseOption("m", "matrix", "Name Phylip distance matrix file.  Must contain at least three species.", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
    setRequired(parser, "matrix");
    setValidValues(parser, "matrix", "dist");
    addOption(parser, seqan::ArgParseOption("o", "out-file", "Path to write output to.", seqan::ArgParseArgument::OUTPUT_FILE, "FILE"));
    setDefaultValue(parser, "out-file", "tree.dot");
    setValidValues(parser, "out-file", "dot newick");

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
        return res == seqan::ArgumentParser::PARSE_ERROR;

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
    String<char> tmp = outfile;
    toLower(tmp);
    if (endsWith(tmp, ".dot"))
        format = "dot";
    else
        format = "newick";

    // Read the distance matrix
    String<TName> names;
    String<TDistanceValue> matrix;
    VirtualStream<char, Input> inPhylip(infile.c_str());
    DirectionIterator<VirtualStream<char, Input>, Input>::Type iter(directionIterator(inPhylip, Input()));
    readPhylipMatrix(matrix, names, iter);

    // Create the tree
    Graph<Tree<TDistanceValue> > tree;
    switch (build)
    {
        case 0:
            njTree(matrix, tree);
            break;
        case 1:
            upgmaTree(matrix, tree, UpgmaMin());
            break;
        case 2:
            upgmaTree(matrix, tree, UpgmaMax());
            break;
        case 3:
            upgmaTree(matrix, tree, UpgmaAvg());
            break;
        case 4:
            upgmaTree(matrix, tree, UpgmaWeightAvg());
            break;
        default:
            SEQAN_FAIL("unknown build method.");
    }

    VirtualStream<char, Output> oStream(outfile.c_str());

    if (format == "dot")
    {
        TSize nameLen = length(names);
        resize(names, numVertices(tree));
        // Add the label prefix for leaves
        for (TSize i = 0;i < nameLen; ++i)
        {
            TName tmpName = "label = \"";
            append(tmpName, names[i], Generous());
            appendValue(tmpName, '"');
            names[i] = tmpName;
        }

        // Append emty names for internal vertices
        for (; nameLen < length(names); ++nameLen)
            names[nameLen] = "label = \"\"";

        // Write the result
        writeRecords(oStream, tree, names, DotDrawing());
    }
    else if (format == "newick")
    {
        // If nj tree collapse the root
        if (build == 0)
            writeRecords(oStream, tree, names, true, NewickFormat());
        else
            writeRecords(oStream, tree, names, false, NewickFormat());
    }

    return 0;
}
