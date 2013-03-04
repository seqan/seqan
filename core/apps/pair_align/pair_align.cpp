/*==========================================================================
               SeqAn - The Library for Sequence Analysis
                         http://www.seqan.de 
============================================================================
Copyright (C) 2007-2012

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
#include <seqan/arg_parse.h>

#include <iostream>
#include <fstream>


using namespace seqan;

// --------------------------------------------------------------------------
// Class Options
// --------------------------------------------------------------------------

struct Options
{
    static int const INVALID_DIAGONAL;

    seqan::CharString inputFile;
    seqan::CharString outputFile;
    seqan::CharString alphabet;
    seqan::CharString method;
    int outputFormat;
    int gop;
    int gex;
    seqan::CharString matrix;
    int msc;
    int mmsc;
    int low;
    int high;
    seqan::CharString config;

    Options() : gop(0), gex(0), msc(0), mmsc(0), low(INVALID_DIAGONAL), high(INVALID_DIAGONAL)
    {}
};

int const Options::INVALID_DIAGONAL = seqan::MaxValue<int>::VALUE;

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
_initAlignParams(Options const & options, TScore& sc) {
    // Set options
    sc.data_gap_open = options.gop;
    sc.data_gap_extend = options.gex;
    int msc = options.msc;
    _setMatchScore(sc, msc);
    int mmsc = options.mmsc;
    _setMismatchScore(sc, mmsc);
    ::std::string seqfile = toCString(options.inputFile);
    ::std::string outfile = toCString(options.outputFile);
    unsigned int method = 0;
    String<char> meth = options.method;
    if (meth == "nw") method = 0;
    else if (meth == "gotoh") method = 1;
    else if (meth == "sw") method = 2;
    else if (meth == "lcs") method = 3;
    int low = 0;
    int high = 0;
    bool banded = false;
    if (options.low != Options::INVALID_DIAGONAL)
    {
        low = options.low;
        banded = true;
    }
    if (options.high != Options::INVALID_DIAGONAL)
    {
        high = options.high;
        banded = true;
    }

    // Check options
    if (low > high) banded = false;
    
    // Do pairwise alignment
    String<char> config = options.config;
    if (!empty(config))
    {
        if (config == "tttt")
            pairwise_align<TAlphabet, AlignConfig<true, true, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tttf")
            pairwise_align<TAlphabet, AlignConfig<true, true, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ttft")
            pairwise_align<TAlphabet, AlignConfig<true, true, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ttff")
            pairwise_align<TAlphabet, AlignConfig<true, true, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tftt")
            pairwise_align<TAlphabet, AlignConfig<true, false, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tftf")
            pairwise_align<TAlphabet, AlignConfig<true, false, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tfft")
            pairwise_align<TAlphabet, AlignConfig<true, false, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "tfff")
            pairwise_align<TAlphabet, AlignConfig<true, false, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fttt")
            pairwise_align<TAlphabet, AlignConfig<false, true, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fttf")
            pairwise_align<TAlphabet, AlignConfig<false, true, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ftft")
            pairwise_align<TAlphabet, AlignConfig<false, true, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ftff")
            pairwise_align<TAlphabet, AlignConfig<false, true, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fftt")
            pairwise_align<TAlphabet, AlignConfig<false, false, true, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "fftf")
            pairwise_align<TAlphabet, AlignConfig<false, false, true, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ffft")
            pairwise_align<TAlphabet, AlignConfig<false, false, false, true> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
        else if (config == "ffff")
            pairwise_align<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
    }
    else
    {
        pairwise_align<TAlphabet, AlignConfig<false, false, false, false> >(sc, seqfile, method, low, high, banded, options.outputFormat, outfile);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, Dna5 const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<Dna5>(options, sc);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<Dna5>(options, sc);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, char const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<char>(options, sc);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<char>(options, sc);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, Rna5 const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {    
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<Rna5>(options, sc);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<Rna5>(options, sc);
    }
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(Options const & options, AminoAcid const) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, options.matrix);
        _initAlignParams<AminoAcid>(options, sc);
    }
    else
    {
        Blosum62 sc;
        _initAlignParams<AminoAcid>(options, sc);
    }
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("pair_align");
    // Set short description, version, and date.
    setShortDescription(parser, "Pairwise alignment");
    setVersion(parser, "1.1");
    setDate(parser, "November 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-s\\fP \\fIIN.fa\\fP");
	setCategory(parser, "Sequence Alignment");
    addDescription(parser,
                   "The program allows to align two sequences using dyamic programming alignment algorithms while "
                   "tweaking various parameters.");

    addSection(parser, "Main Options");
    addOption(parser, seqan::ArgParseOption("s", "seq", "FASTA file with two sequences.", seqan::ArgParseOption::INPUTFILE, "IN.fa"));
    setRequired(parser, "seq");
    setValidValues(parser, "seq", "fasta fa");
    addOption(parser, seqan::ArgParseOption("a", "alphabet", "Sequence alphabet.", seqan::ArgParseOption::STRING, "ALPHABET"));
    setValidValues(parser, "alphabet", "protein dna rna text");
    setDefaultValue(parser, "alphabet", "protein");
    addOption(parser, seqan::ArgParseOption("m", "method",
                                            "DP alignment method: Needleman-Wunsch, Gotoh, Smith-Waterman, "
                                            "Longest Common Subsequence",
                                            seqan::ArgParseOption::STRING, "METHOD"));
    setValidValues(parser, "method", "nw gotoh sw lcs");
    setDefaultValue(parser, "method", "gotoh");
    addOption(parser, seqan::ArgParseOption("o", "outfile", "Output filename.", seqan::ArgParseOption::OUTPUTFILE, "OUT"));
    setDefaultValue(parser, "outfile", "out.fasta");
	setValidValues(parser, "outfile", "fa fasta msf");
	//TODO(rmaerker): We removed this option. The file format is derived from the outfile format.
    //addOption(parser, seqan::ArgParseOption("f", "format", "Output format.", seqan::ArgParseOption::STRING));
    //setValidValues(parser, "format", "fa fasta msf");
    //setDefaultValue(parser, "format", "fasta");

    addSection(parser, "Scoring Options");
    addOption(parser, seqan::ArgParseOption("g", "gop", "Gap open penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gop", "-11");
    addOption(parser, seqan::ArgParseOption("e", "gex", "Gap extension penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gex", "-1");
    addOption(parser, seqan::ArgParseOption("ma", "matrix", "Score matrix.", seqan::ArgParseOption::STRING, "MATRIX_FILE"));
    addOption(parser, seqan::ArgParseOption("ms", "msc", "Match score.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msc", "5");
    addOption(parser, seqan::ArgParseOption("mm", "mmsc", "Mismatch penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "mmsc", "-4");

    addSection(parser, "Banded Alignment Options");
    addOption(parser, seqan::ArgParseOption("lo", "low", "Lower diagonal.", seqan::ArgParseOption::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("hi", "high", "Upper diagonal.", seqan::ArgParseOption::INTEGER, "INT"));

    addSection(parser, "DP Matrix Configuration Options");
    addOption(parser, seqan::ArgParseOption("c", "config", "Alignment configuration.", seqan::ArgParseOption::STRING, "CONF"));
    setValidValues(parser, "config", "ffff ffft fftf fftt ftff ftft fttf fttt tfff tfft tftf tftt ttff ttft tttf tttt");

    addTextSection(parser, "Alignment configuration");
    addText(parser,
            "The alignment configuration is a string of four characters, each being either t or f. All "
            "combinations are allowed. The meaning is as follows.");
    addListItem(parser, "tfff", "First row initialized with 0s.");
    addListItem(parser, "ftff", "First column initialized with 0s.");
    addListItem(parser, "fftf", "Search last column for maximum.");
    addListItem(parser, "ffft", "Search last row for maximum.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.inputFile, parser, "seq");
    getOptionValue(options.outputFile, parser, "outfile");
	// Guess file format based on extension of file.
	CharString tmp = options.outputFile;
	if (endsWith(tmp, ".fa") || endsWith(tmp, "fasta"))
	    options.outputFormat = 0;
	else if (endsWith(tmp, ".msf"))
	    options.outputFormat = 1;

    getOptionValue(options.alphabet, parser, "alphabet");
    getOptionValue(options.method, parser, "method");

    //getOptionValue(options.format, parser, "format");
    getOptionValue(options.gop, parser, "gop");
    getOptionValue(options.gex, parser, "gex");
    getOptionValue(options.matrix, parser, "matrix");
    getOptionValue(options.msc, parser, "msc");
    getOptionValue(options.mmsc, parser, "mmsc");
    getOptionValue(options.low, parser, "low");
    getOptionValue(options.high, parser, "high");
    getOptionValue(options.config, parser, "config");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, const char *argv[])
{
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    
    // Initialize scoring matrices
    if (options.alphabet == "dna") _initScoreMatrix(options, Dna5());
    else if (options.alphabet == "rna") _initScoreMatrix(options, Rna5());
    else if (options.alphabet == "protein") _initScoreMatrix(options, AminoAcid());
    else _initScoreMatrix(options, char());

    return 0;
}
