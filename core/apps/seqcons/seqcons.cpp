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
#include <seqan/misc/misc_cmdparser.h>

#include <iostream>
#include <fstream>


using namespace seqan;

void
_addVersion(CommandLineParser & parser)
{
	::std::string rev = "$Revision: 4663 $";
	addVersionLine(parser, "Version 0.22 (06. August 2009) Revision: " + rev.substr(11, 4) + "");
}

void setUpCommandLineParser(CommandLineParser & parser)
{
	_addVersion(parser);

	addTitleLine(parser, "***************************************");
	addTitleLine(parser, "* Multi-read alignment - SeqCons      *");
	addTitleLine(parser, "* (c) Copyright 2009 by Tobias Rausch *");
	addTitleLine(parser, "***************************************");

	addUsageLine(parser, "-r <FASTA file with reads> [Options]");
	addUsageLine(parser, "-a <AMOS message file> [Options]");
	addUsageLine(parser, "-s <Sam file> [-c <FASTA contigs file>] [Options]");

	addSection(parser, "Main Options:");
	addOption(parser, addArgumentText(CommandLineOption("r", "reads", "file with reads", OptionType::String), "<FASTA reads file>"));
	addOption(parser, addArgumentText(CommandLineOption("a", "afg", "message file", OptionType::String), "<AMOS afg file>"));
	addOption(parser, addArgumentText(CommandLineOption("s", "sam", "Sam file", OptionType::String), "<Sam file>"));
	addOption(parser, addArgumentText(CommandLineOption("c", "contigs", "FASTA file with contigs, ignored if not Sam input", OptionType::String), "<FASTA contigs file>"));
	addOption(parser, addArgumentText(CommandLineOption("o", "outfile", "output filename", (int)OptionType::String, "align.txt"), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("f", "format", "output format", (int)OptionType::String, "afg"), "[seqan | afg | sam]"));
	addOption(parser, addArgumentText(CommandLineOption("m", "method", "alignment method", (int)OptionType::String, "realign"), "[realign | msa]"));
	addOption(parser, addArgumentText(CommandLineOption("b", "bandwidth", "bandwidth", (int)OptionType::Int, 8), "<Int>"));
	addOption(parser, CommandLineOption("n", "noalign", "no align, only convert input", OptionType::Boolean));

	addSection(parser, "MSA Method Options:");
	addOption(parser, addArgumentText(CommandLineOption("ma", "matchlength", "min. overlap length", (int)OptionType::Int, 15), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("qu", "quality", "min. overlap precent identity", (int)OptionType::Int, 80), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("ov", "overlaps", "min. number of overlaps per read", (int)OptionType::Int, 3), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("wi", "window", "window size", (int)OptionType::Int, 0), "<Int>"));
	addHelpLine(parser, "/*If this parameter is > 0 then all");
	addHelpLine(parser, "  overlaps within a given window");
	addHelpLine(parser, "  are computed.*/");

	addSection(parser, "ReAlign Method Options:");
	addOption(parser, CommandLineOption("in", "include", "include contig sequence", OptionType::Boolean));
	addOption(parser, addArgumentText(CommandLineOption("rm", "rmethod", "realign method", (int)OptionType::String, "gotoh"), "[nw | gotoh]"));
}

int parseCommandLine(ConsensusOptions & consOpt, CommandLineParser & parser, int argc, const char * argv[])
{
    if (argc == 1)
	{
		shortHelp(parser, std::cerr);	// print short help and exit
		return 1;
	}

	if (!parse(parser, argc, argv, ::std::cerr))
        return 1;
	if (isSetLong(parser, "help") || isSetLong(parser, "version"))
        return 0;	// print help or version and exit

	// Main options
	getOptionValueLong(parser, "reads", consOpt.readsfile);
	getOptionValueLong(parser, "afg", consOpt.afgfile);
	getOptionValueLong(parser, "sam", consOpt.samfile);
	getOptionValueLong(parser, "contigs", consOpt.contigsfile);
	getOptionValueLong(parser, "outfile", consOpt.outfile);

    if (empty(consOpt.samfile) && !empty(consOpt.contigsfile))
        std::cerr << "WARNING: Contigs specified by input is not FASTA, ignoring --contigs parameters." << std::endl;

	String<char> optionVal;
	getOptionValueLong(parser, "format", optionVal);
	if (optionVal == "seqan")
        consOpt.output = 0;
	else if (optionVal == "afg")
        consOpt.output = 1;
	else if (optionVal == "frg")
        consOpt.output = 2;
	else if (optionVal == "cgb")
        consOpt.output = 3;
	else if (optionVal == "sam")
        consOpt.output = 4;

	getOptionValueLong(parser, "method", optionVal);
	if (optionVal == "realign")
        consOpt.method = 0;
	else if (optionVal == "msa")
        consOpt.method = 1;

	getOptionValueLong(parser, "bandwidth", consOpt.bandwidth);
#ifdef CELERA_OFFSET
	if (!isSetLong(parser, "bandwidth")) consOpt.bandwidth = 15;	
#endif
	getOptionValueLong(parser, "noalign", consOpt.noalign);

	// MSA options
	getOptionValueLong(parser, "matchlength", consOpt.matchlength);
	getOptionValueLong(parser, "quality", consOpt.quality);
	getOptionValueLong(parser, "overlaps", consOpt.overlaps);
#ifdef CELERA_OFFSET
	if (!isSetLong(parser, "overlaps")) consOpt.overlaps = 5;	
#endif
	getOptionValueLong(parser, "window", consOpt.window);
	
	// ReAlign options
	getOptionValueLong(parser, "include", consOpt.include);
	getOptionValueLong(parser, "rmethod", optionVal);
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
		bool success = _convertSimpleReadFile(strmReads, fragStore, consOpt.readsfile, moveToFront);
		if (!success)
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
	CommandLineParser parser;
    setUpCommandLineParser(parser);

	// Get all command line options
	ConsensusOptions consOpt;
    int ret = parseCommandLine(consOpt, parser, argc, argv);
    if (ret != 0)
        return ret;

	// Create a new fragment store
	typedef FragmentStore<> TFragmentStore;
	typedef Size<TFragmentStore>::Type TSize;
	TFragmentStore fragStore;

	// Load the reads and layout positions
	TSize numberOfContigs = 0;
    ret = loadFiles(fragStore, numberOfContigs, consOpt);
    if (ret != 0) {
        shortHelp(parser, std::cerr);
        return ret;
    }

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
