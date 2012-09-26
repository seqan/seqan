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
#include <seqan/align.h>
#include <seqan/modifier.h>
#include <seqan/misc/misc_cmdparser.h>

#include <iostream>
#include <fstream>


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 4566 $";
	addVersionLine(parser, "Version 1.0 (15. July 2009) Revision: " + rev.substr(11, 4) + "");
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TSeqSet, typename TNameSet>
bool _loadSequences(TSeqSet& sequences, 
					TNameSet& fastaIDs,
					const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);
	unsigned seqCount = length(multiFasta);
	resize(sequences, seqCount, Exact());
	resize(fastaIDs, seqCount, Exact());
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		assignSeqId(fastaIDs[i], multiFasta[i], format);
		assignSeq(sequences[i], multiFasta[i], format);
	}
	return (seqCount > 0);
}

// TODO(holtgrew): Make publically available.
template<typename TStringSet, typename TCargo, typename TSpec>
inline int
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Lcs)
{
	return globalAlignment(g, stringSet(g), Lcs());
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TAlignConfig, typename TScore, typename TSeqFile, typename TMethod, typename TDiag, typename TOutputFormat, typename TOutfile>
inline void
pairwise_align(TScore const& sc,
			   TSeqFile& seqfile,
			   TMethod method,
			   TDiag low,
			   TDiag high,
			   bool banded,
			   TOutputFormat outputFormat,
			   TOutfile& outfile) 
{
	// Load the 2 sequences
	typedef String<TAlphabet> TSequence;
	StringSet<TSequence, Owner<> > sequenceSet;
	StringSet<String<char> > sequenceNames;
	_loadSequences(sequenceSet, sequenceNames, seqfile.c_str());

	// Fix low and high diagonal.
	low = _max(low, -1 * (int) length(sequenceSet[1]));
	high = _min(high, (int) length(sequenceSet[0]));

	// Align the sequences
	Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > gAlign(sequenceSet);
	
	int aliScore = 0;
	// Banded alignment?
	if (!banded) {
		if (method == 0) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), NeedlemanWunsch());
		else if (method == 1) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), Gotoh());
		else if (method == 2) aliScore = localAlignment(gAlign, sc);
		else if (method == 3) aliScore = globalAlignment(gAlign, Lcs());
	} else {
		if (method == 0) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), low, high, NeedlemanWunsch());
		else if (method == 1) aliScore = globalAlignment(gAlign, sc, TAlignConfig(), low, high, Gotoh());
	}
	
	// Alignment output
	std::cout << "Alignment score: " << aliScore << std::endl;
	if (outputFormat == 0) {
		FILE* strmWrite = fopen(outfile.c_str(), "w");
		write(strmWrite, gAlign, sequenceNames, FastaFormat());
		fclose(strmWrite);
	} else if (outputFormat == 1) {
		FILE* strmWrite = fopen(outfile.c_str(), "w");
		write(strmWrite, gAlign, sequenceNames, MsfFormat());
		fclose(strmWrite);
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_setMatchScore(TScore&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_setMismatchScore(TScore&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSc>
inline void
_setMatchScore(Score<int, Simple>& sc, TSc msc) {
	sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TSc>
inline void
_setMismatchScore(Score<int, Simple>& sc, TSc mmsc) {
	sc.data_mismatch = mmsc;
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
_initAlignParams(CommandLineParser& parser, TScore& sc) {
	// Set options
	getOptionValueLong(parser, "gop", sc.data_gap_open);
	getOptionValueLong(parser, "gex", sc.data_gap_extend);
	int msc = 0;
	getOptionValueLong(parser, "msc", msc);
	_setMatchScore(sc, msc);
	int mmsc = 0;
	getOptionValueLong(parser, "mmsc", mmsc);
	_setMismatchScore(sc, mmsc);
	::std::string seqfile;
	getOptionValueLong(parser, "seq", seqfile);
	::std::string outfile = "out.fasta";
	getOptionValueLong(parser, "outfile", outfile);
	unsigned int method = 0;
	String<char> meth;
	getOptionValueLong(parser, "method", meth);
	if (meth == "nw") method = 0;
	else if (meth == "gotoh") method = 1;
	else if (meth == "sw") method = 2;
	else if (meth == "lcs") method = 3;
	unsigned int outputFormat = 0;
	String<char> format;
	getOptionValueLong(parser, "format", format);
	if (format == "fasta") outputFormat = 0;
	else if (format == "msf") outputFormat = 1;
	int low = 0;
	int high = 0;
	bool banded = false;
	if (isSetLong(parser, "low")) {
		getOptionValueLong(parser, "low", low);
		banded = true;
	}
	if (isSetLong(parser, "high")) {
		getOptionValueLong(parser, "high", high);
		banded = true;
	}

	// Check options
	if (!isSetLong(parser, "seq")) { help(parser); exit(1); }
	if (low > high) banded = false;
	
	// Do pairwise alignment
	if (isSetLong(parser, "config")) {
		String<char> config;
		getOptionValueLong(parser, "config", config);
		if (config == "tttt") pairwise_align<TAlphabet, AlignConfig<true, true, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "tttf") pairwise_align<TAlphabet, AlignConfig<true, true, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "ttft") pairwise_align<TAlphabet, AlignConfig<true, true, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "ttff") pairwise_align<TAlphabet, AlignConfig<true, true, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "tftt") pairwise_align<TAlphabet, AlignConfig<true, false, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "tftf") pairwise_align<TAlphabet, AlignConfig<true, false, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "tfft") pairwise_align<TAlphabet, AlignConfig<true, false, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "tfff") pairwise_align<TAlphabet, AlignConfig<true, false, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "fttt") pairwise_align<TAlphabet, AlignConfig<false, true, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "fttf") pairwise_align<TAlphabet, AlignConfig<false, true, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "ftft") pairwise_align<TAlphabet, AlignConfig<false, true, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "ftff") pairwise_align<TAlphabet, AlignConfig<false, true, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "fftt") pairwise_align<TAlphabet, AlignConfig<false, false, true, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "fftf") pairwise_align<TAlphabet, AlignConfig<false, false, true, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "ffft") pairwise_align<TAlphabet, AlignConfig<false, false, false, true> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
		else if (config == "ffff") pairwise_align<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
	} else pairwise_align<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, outputFormat, outfile);
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(CommandLineParser& parser, Dna5 const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initAlignParams<Dna5>(parser, sc);
	} else {
		Score<int> sc;
		_initAlignParams<Dna5>(parser, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(CommandLineParser& parser, char const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initAlignParams<char>(parser, sc);
	} else {
		Score<int> sc;
		_initAlignParams<char>(parser, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(CommandLineParser& parser, Rna5 const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initAlignParams<Rna5>(parser, sc);
	} else {
		Score<int> sc;
		_initAlignParams<Rna5>(parser, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(CommandLineParser& parser, AminoAcid const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initAlignParams<AminoAcid>(parser, sc);
	} else {
		Blosum62 sc;
		_initAlignParams<AminoAcid>(parser, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {

	// Command line parsing
	CommandLineParser parser;
	_addVersion(parser);
	
	addTitleLine(parser, "***************************************");
	addTitleLine(parser, "* Pairwise alignment - PairAlign      *");
	addTitleLine(parser, "* (c) Copyright 2009 by Tobias Rausch *");
	addTitleLine(parser, "***************************************");

	addUsageLine(parser, "-s <FASTA sequence file> [Options]");

	addSection(parser, "Main Options:");
	addOption(parser, addArgumentText(CommandLineOption("s", "seq", "file with 2 sequences", OptionType::String), "<FASTA Sequence File>"));
	addOption(parser, addArgumentText(CommandLineOption("a", "alphabet", "sequence alphabet", (int)OptionType::String, "protein"), "[protein | dna | rna | text]"));
	addOption(parser, addArgumentText(CommandLineOption("m", "method", "alignment method", (int)OptionType::String, "gotoh"), "[nw, gotoh, sw, lcs]"));
	addHelpLine(parser, "nw = Needleman-Wunsch");
	addHelpLine(parser, "gotoh = Gotoh");
	addHelpLine(parser, "sw = Smith-Waterman");
	addHelpLine(parser, "lcs = Longest common subsequence");
	addOption(parser, addArgumentText(CommandLineOption("o", "outfile", "output filename", (int)OptionType::String, "out.fasta"), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("f", "format", "output format", (int)OptionType::String, "fasta"), "[fasta | msf]"));
	
	addSection(parser, "Scoring Options:");
	addOption(parser, addArgumentText(CommandLineOption("g", "gop", "gap open penalty", (int)OptionType::Int, -11), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("e", "gex", "gap extension penalty", (int)OptionType::Int, -1), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("ma", "matrix", "score matrix", (int)OptionType::String, "Blosum62"), "<Matrix file>"));
	addOption(parser, addArgumentText(CommandLineOption("ms", "msc", "match score", (int)OptionType::Int, 5), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("mm", "mmsc", "mismatch penalty", (int)OptionType::Int, -4), "<Int>"));
	
	addSection(parser, "Banded Alignment Options:");
	addOption(parser, addArgumentText(CommandLineOption("lo", "low", "lower diagonal", OptionType::Int), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("hi", "high", "upper diagonal", OptionType::Int), "<Int>"));
			
	addSection(parser, "DP Matrix Configuration Options:");
	addOption(parser, addArgumentText(CommandLineOption("c", "config", "alignment configuration", (int)OptionType::String, "ffff"), "[ffff | ... | tttt]"));
	addHelpLine(parser, "tfff = First row with 0's");
	addHelpLine(parser, "ftff = First column with 0's");
	addHelpLine(parser, "fftf = Search last column for max");
	addHelpLine(parser, "ffft = Search last row for max");
	addHelpLine(parser, "All combinations are allowed.");

	if (argc == 1)
	{
		shortHelp(parser, std::cerr);	// print short help and exit
		return 0;
	}

	if (!parse(parser, argc, argv, ::std::cerr)) return 1;
	if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit

	// Basic command line options
	String<char> alphabet;
	getOptionValueLong(parser, "alphabet", alphabet);
	
	// Initialize scoring matrices
	if (alphabet == "dna") _initScoreMatrix(parser, Dna5());
	else if (alphabet == "rna") _initScoreMatrix(parser, Rna5());
	else if (alphabet == "protein") _initScoreMatrix(parser, AminoAcid());
	else _initScoreMatrix(parser, char());

	return 0;
}
