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
  ==========================================================================
  SeqCons -- Read alignment via realignment or MSA.
  ==========================================================================
  Author: Tobias Rausch <rausch@embl.de>
  ==========================================================================
*/

#include <seqan/basic.h>
#include <seqan/consensus.h>
#include <seqan/modifier.h>

#include <iostream>
#include <fstream>

#include <seqan/arg_parse.h>

using namespace seqan;

void
_setVersion(ArgumentParser & parser)
{
	::std::string rev = "$Revision: 4663 $";
	setVersion(parser, "0.23 Revision: " + rev.substr(11, 4) + "");
	setDate(parser, "Nov 21, 2012");
}

int parseCommandLine(ConsensusOptions & consOpt, int argc, const char * argv[])
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("seqcons");
    _setVersion(parser);
    setShortDescription(parser, "Multi-read alignment.");
    addDescription(parser, "(c) Copyright 2009 by Tobias Rausch");

	addUsageLine(parser, "-r <FASTA file with reads> [Options]");
	addUsageLine(parser, "-a <AMOS message file> [Options]");
	addUsageLine(parser, "-s <Sam file> [-c <FASTA contigs file>] [Options]");


	addSection(parser, "Main Options:");
	addOption(parser, ArgParseOption("r", "reads", "file with reads", ArgParseArgument::STRING, "<FASTA reads file>"));
	addOption(parser, ArgParseOption("a", "afg", "message file", ArgParseArgument::STRING, "<AMOS afg file>"));
	addOption(parser, ArgParseOption("s", "sam", "Sam file", ArgParseArgument::STRING, "<Sam file>"));
	addOption(parser, ArgParseOption("c", "contigs", "FASTA file with contigs, ignored if not Sam input", ArgParseArgument::STRING, "<FASTA contigs file>"));
	addOption(parser, ArgParseOption("o", "outfile", "output filename", ArgParseArgument::STRING, "<Filename>"));
	setValidValues(parser, "outfile", "afg seqan afg am");
	setDefaultValue(parser, "outfile", "align.txt");

	addOption(parser, ArgParseOption("m", "method", "alignment method", ArgParseArgument::STRING, "realign", "[realign | msa]"));
	setDefaultValue(parser, "method", "realign");
	addOption(parser, ArgParseOption("b", "bandwidth", "bandwidth", ArgParseArgument::INTEGER, "<Int>"));
	setDefaultValue(parser, "bandwidth", "8");
    addOption(parser, ArgParseOption("n", "noalign", "no align, only convert input", ArgParseArgument::INTEGER, "<Bool>"));
    setDefaultValue(parser, "noalign", "0");

	addOption(parser, ArgParseOption("ma", "matchlength", "min. overlap length", ArgParseArgument::INTEGER, "<Int>"));
	setDefaultValue(parser, "matchlength", "15");
	addOption(parser, ArgParseOption("qu", "quality", "min. overlap precent identity", ArgParseArgument::INTEGER, "<Int>"));
	setDefaultValue(parser, "quality", "80");
	addOption(parser, ArgParseOption("ov", "overlaps", "min. number of overlaps per read", ArgParseArgument::INTEGER, "<Int>"));
	setDefaultValue(parser, "overlaps", "3");
	addOption(parser, ArgParseOption("wi", "window", "The window size.  If this parameter is greater than 0 then all overlaps within a given window are computed.", ArgParseArgument::INTEGER,"<Int>"));
	setDefaultValue(parser, "window", "0");
 
	addSection(parser, "ReAlign Method Options:");
	addOption(parser, ArgParseOption("in", "include", "include contig sequence", ArgParseArgument::INTEGER, "<Bool>"));
	setDefaultValue(parser, "include", "0");
	addOption(parser, ArgParseOption("rm", "rmethod", "realign method", ArgParseArgument::STRING, "[nw | gotoh]"));
	setDefaultValue(parser, "rmethod", "gotoh");

	if (argc == 1)
	{
		printShortHelp(parser, std::cerr);	// print short help and exit
		return 1;
	}

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        return res;// == seqan::ArgumentParser::PARSE_ERROR;
    }

	// Main options
	getOptionValue(consOpt.readsfile, parser, "reads");
	getOptionValue(consOpt.afgfile, parser, "afg");
	getOptionValue(consOpt.samfile, parser, "sam");
	getOptionValue(consOpt.contigsfile, parser, "contigs");
	getOptionValue(consOpt.outfile, parser, "outfile");

    if (empty(consOpt.samfile) && !empty(consOpt.contigsfile))
        std::cerr << "WARNING: Contigs specified by input is not FASTA, ignoring --contigs parameters." << std::endl;

    // Get lower case of the output file name.  File endings are accepted in both upper and lower case.
    CharString tmp = consOpt.outfile;
    toLower(tmp);

    if (endsWith(tmp, ".seqan"))
        consOpt.output = 0;
    else if (endsWith(tmp, ".afg"))
        consOpt.output = 1;
    else if (endsWith(tmp, ".frg"))
        consOpt.output = 2;
    else if (endsWith(tmp, ".cgb"))
        consOpt.output = 3;
    else if (endsWith(tmp, ".sam"))
        consOpt.output = 4;

    String<char> optionVal;
	getOptionValue(optionVal, parser, "method");
	if (optionVal == "realign")
        consOpt.method = 0;
	else if (optionVal == "msa")
        consOpt.method = 1;

	getOptionValue(consOpt.bandwidth, parser, "bandwidth");
#ifdef CELERA_OFFSET
	if (!isSetLong(parser, "bandwidth")) consOpt.bandwidth = 15;
#endif
    getOptionValue(consOpt.noalign, parser, "noalign");

	// MSA options
	getOptionValue(consOpt.matchlength, parser, "matchlength");
	getOptionValue(consOpt.quality, parser, "quality");
	getOptionValue(consOpt.overlaps, parser, "overlaps");
#ifdef CELERA_OFFSET
	if (!isSetLong(parser, "overlaps")) consOpt.overlaps = 5;
#endif
	getOptionValue(consOpt.window, parser, "window");

	// ReAlign options
	getOptionValue(consOpt.include, parser, "include");
	getOptionValue(optionVal, parser, "rmethod");
	if (optionVal == "nw")
        consOpt.rmethod = 0;
	else if (optionVal == "gotoh")
        consOpt.rmethod = 1;
    return 0;

}

// Load the reads and layout positions
template <typename TFragmentStore, typename TSize>
int loadFiles(TFragmentStore & fragStore, TSize & numberOfContigs, ConsensusOptions const & consOpt)
{
//IOREV
    std::cerr << "Reading input..." << std::endl;
	if (!empty(consOpt.readsfile)) {
		// Load simple read file
        std::fstream strmReads(consOpt.readsfile.c_str(), std::fstream::in | std::fstream::binary);
		bool moveToFront = false;
		if (consOpt.noalign) moveToFront = true;
		if (_convertSimpleReadFile(strmReads, fragStore, consOpt.readsfile, moveToFront) != 0)
			return 1;
		numberOfContigs = 1;
	} else if (!empty(consOpt.afgfile)) {
		// Load Amos message file
        std::fstream strmReads(consOpt.afgfile.c_str(), std::fstream::in | std::fstream::binary);
		read(strmReads, fragStore, Amos());	
		numberOfContigs = length(fragStore.contigStore);
	} else if (!empty(consOpt.samfile)) {
        // Possibly load contigs into fragment store.
        if (!empty(consOpt.contigsfile)) {
            if (!loadContigs(fragStore, consOpt.contigsfile.c_str())) {
                std::cerr << "Could not load contigs file " << consOpt.contigsfile.c_str() << std::endl;
                return 1;
            }
        }
		// Load Sam message file
        std::fstream strmReads(consOpt.samfile.c_str(), std::fstream::in | std::fstream::binary);
		read(strmReads, fragStore, Sam());
		numberOfContigs = length(fragStore.contigStore);
	} else {
		return 1;
	}

    return 0;
}

// Write resulting alignment.
template <typename TFragmentStore>
int writeOutput(TFragmentStore /*const*/ & fragStore, ConsensusOptions const & consOpt)
{
//IOREV
    std::cerr << "Writing output..." << std::endl;
	if (consOpt.output == 0) {
		// Write old SeqAn multi-read alignment format
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		write(strmWrite, fragStore, FastaReadFormat());	
		fclose(strmWrite);
	} else if (consOpt.output == 1) {
		// Write Amos
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		write(strmWrite, fragStore, Amos());	
		fclose(strmWrite);
	} else if (consOpt.output == 2) {
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		_writeCeleraFrg(strmWrite, fragStore);	
		fclose(strmWrite);
	} else if (consOpt.output == 3) {
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		_writeCeleraCgb(strmWrite, fragStore);	
		fclose(strmWrite);
	} else if (consOpt.output == 4) {
		// Write out resulting MSA in a Sam file.
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		write(strmWrite, fragStore, Sam());
		fclose(strmWrite);
		// Write out resulting consensus sequence.
		char buffer[10*1024];
		buffer[0] = '\0';
		strcat(buffer, consOpt.outfile.c_str());
		strcat(buffer, ".consensus.fasta");
		strmWrite = fopen(buffer, "w");
		writeContigs(strmWrite, fragStore, Fasta());
		fclose(strmWrite);
	}
    return 0;
}

int main(int argc, const char *argv[]) {
	// Command line parsing
	ArgumentParser parser;
	ConsensusOptions consOpt;
    int ret = parseCommandLine(consOpt, argc, argv);
    if (ret != 0)
        return ret;

	// Create a new fragment store
	typedef FragmentStore<> TFragmentStore;
	typedef Size<TFragmentStore>::Type TSize;
	TFragmentStore fragStore;

	// Load the reads and layout positions
	TSize numberOfContigs = 0;
    ret = loadFiles(fragStore, numberOfContigs, consOpt);

	// Multi-realignment desired or just conversion of the input
	if (!consOpt.noalign) {
		// Iterate over all contigs
        if (consOpt.method == 0)
            std::cerr << "Performing realignment..." << std::endl;
        else
            std::cerr << "Performing consensus alignment..." << std::endl;
		for (TSize currentContig = 0; currentContig < numberOfContigs; ++currentContig) {
            std::cerr << "contig " << (currentContig + 1) << "/" << numberOfContigs << std::endl;
			if (consOpt.method == 0) {
				Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > combinedScore;
				reAlign(fragStore, combinedScore, currentContig, consOpt.rmethod, consOpt.bandwidth, consOpt.include);
				if (consOpt.include) reAlign(fragStore, combinedScore, currentContig, consOpt.rmethod, consOpt.bandwidth, false);
			} else {
                std::cerr << "Performing consensus alignment..." << std::endl;
				// Import all reads of the given contig
				typedef TFragmentStore::TReadSeq TReadSeq;
				StringSet<TReadSeq, Owner<> > readSet;
				String<Pair<TSize, TSize> > begEndPos;
				_loadContigReads(readSet, begEndPos, fragStore, currentContig);
				if (!length(readSet)) continue;

				// Align the reads
				typedef StringSet<TReadSeq, Dependent<> > TStringSet;
				typedef Graph<Alignment<TStringSet, void, WithoutEdgeId> > TAlignGraph;
				TAlignGraph gOut(readSet);
				consensusAlignment(gOut, begEndPos, consOpt);
			
				// Update the contig in the fragment store
				if (!empty(gOut)) updateContig(fragStore, gOut, currentContig);
				clear(gOut);

				//// Debug code for CA
				//mtRandInit();
				//String<char> fileTmp1 = "tmp1";
				//String<char> fileTmp2 = "tmp2";
				//for(int i = 0; i<10; ++i) {
				//	int file = (mtRand() % 20) + 65;
				//	appendValue(fileTmp1, char(file));
				//	appendValue(fileTmp2, char(file));
				//}
				//std::fstream strm3;
				//strm3.open(toCString(fileTmp2), std::ios_base::out | std::ios_base::trunc);
				//for(int i = 0;i<(int) length(origStrSet); ++i) {
				//	std::stringstream name;
				//	name << value(begEndPos, i).i1 << "," << value(begEndPos, i).i2;
				//	String<char> myTitle = name.str();
				//	write(strm3, origStrSet[i], myTitle, Fasta());			
				//	if (value(begEndPos, i).i1 > value(begEndPos, i).i2) reverseComplement(origStrSet[i]);
				//}
				//strm3.close();
			}
		} // end loop over all contigs
	}

    // Write result.
    ret = writeOutput(fragStore, consOpt);
    return ret;
}
