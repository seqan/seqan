// ==========================================================================
//                            pair_align_lib.h
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: René Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implementation of the pairwise align tool.
// ==========================================================================

#ifndef APPS_PAIR_ALIGN_PAIR_ALIGN_LIB_H_
#define APPS_PAIR_ALIGN_PAIR_ALIGN_LIB_H_

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

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

//////////////////////////////////////////////////////////////////////////////////

template <typename TSeqSet, typename TNameSet>
bool _loadSequences(TSeqSet& sequences,
                    TNameSet& fastaIDs,
                    const char *fileName)
{
    SeqFileIn inFile;
    if (!open(inFile, fileName))
    {
        std::cerr << "Could not open file " << fileName << " for reading!" << std::endl;
        return false;
    }

    try
    {
        readRecords(fastaIDs, sequences, inFile);
    }
    catch (seqan::ParseError const & e)
    {
        std::cerr << "ERROR: Problem parsing input file: " << e.what() << "\n";
        return false;
    }

    if (length(fastaIDs) > 2u)  // Limit number of sequences to 2.
    {
        resize(fastaIDs, 2, Exact());
        resize(sequences, 2, Exact());
        return true;
    }
    return (length(fastaIDs) == 2u);
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

// Helper interfaces to select correct public interface given the method specialization.
template <typename TAlign, typename TScore, typename TConfig, typename TBandSize, typename TIsBand, typename TMethod>
inline typename Value<TScore>::Type
_alignHelper(TAlign & gAlign,
             TScore const & sc,
             TConfig const & config,
             TBandSize low,
             TBandSize high,
             TIsBand banded,
             TMethod const & /*globalAlign*/)
{
    if (banded)
        return globalAlignment(gAlign, sc, config, low, high, TMethod());
    return globalAlignment(gAlign, sc, config, TMethod());
}

template <typename TAlign, typename TScore, typename TConfig, typename TBandSize, typename TIsBand>
inline typename Value<TScore>::Type
_alignHelper(TAlign & gAlign,
             TScore const & sc,
             TConfig const & /*config*/,
             TBandSize low,
             TBandSize high,
             TIsBand banded,
             SmithWaterman const & /*tag*/)
{
    if (banded)
        return localAlignment(gAlign, sc, low, high);
    return localAlignment(gAlign, sc);
}

template <typename TAlign, typename TScore, typename TConfig, typename TBandSize, typename TIsBand>
inline typename Value<TScore>::Type
_alignHelper(TAlign & gAlign,
             TScore const &,
             TConfig,
             TBandSize,
             TBandSize,
             TIsBand,
             Lcs const &)
{
    return globalAlignment(gAlign, Lcs());
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TAlignConfig, typename TScore, typename TSeqFile, typename TMethod, typename TDiag, typename TOutfile>
inline void
pairwise_align(TScore const& sc,
               TSeqFile& seqfile,
               TMethod method,
               TDiag low,
               TDiag high,
               bool banded,
               TOutfile& outfile)
{
    typedef VirtualStream<char, Output> TOutputStream;

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

    int aliScore = _alignHelper(gAlign, sc, TAlignConfig(), low, high, banded, method);

    // Alignment output
    std::cout << "Alignment score: " << aliScore << std::endl;
    TOutputStream outputFile;
    if (!open(outputFile, toCString(outfile)))
    {
        std::cerr << "Could not open " << outfile << " for writing!" << std::endl;
        return;
    }

    if (guessFormatFromFilename(outfile, Fasta()))
        write(outputFile, gAlign, sequenceNames, FastaFormat());
    else
        write(outputFile, gAlign, sequenceNames, MsfFormat());
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

template<typename TAlphabet, typename TScore, typename TAlgorithm, typename TAlignSpec>
inline void
_initAlignParams(Options const & options, TScore& sc, TAlgorithm const & alg, TAlignSpec) {
    // Set options
    sc.data_gap_open = options.gop;
    sc.data_gap_extend = options.gex;
    int msc = options.msc;
    _setMatchScore(sc, msc);
    int mmsc = options.mmsc;
    _setMismatchScore(sc, mmsc);
    ::std::string seqfile = toCString(options.inputFile);
    ::std::string outfile = toCString(options.outputFile);
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
    pairwise_align<TAlphabet, TAlignSpec>(sc, seqfile, alg, low, high, banded, outfile);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithm, typename TAlignSpec>
inline void
_initScoreMatrix(Options const & options, Dna5 const, TAlgorithm const & alg, TAlignSpec const & spec) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, toCString(options.matrix));
        _initAlignParams<Dna5>(options, sc, alg, spec);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<Dna5>(options, sc, alg, spec);
    }
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithm, typename TAlignSpec>
inline void
_initScoreMatrix(Options const & options, char const, TAlgorithm const & alg, TAlignSpec const & spec) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, toCString(options.matrix));
        _initAlignParams<char>(options, sc, alg, spec);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<char>(options, sc, alg, spec);
    }
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithm, typename TAlignSpec>
inline void
_initScoreMatrix(Options const & options, Rna5 const, TAlgorithm const & alg, TAlignSpec const & spec) {
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, toCString(options.matrix));
        _initAlignParams<Rna5>(options, sc, alg, spec);
    }
    else
    {
        Score<int> sc;
        _initAlignParams<Rna5>(options, sc, alg, spec);
    }
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithm, typename TAlignSpec>
inline void
_initScoreMatrix(Options const & options, AminoAcid const, TAlgorithm const & alg, TAlignSpec const & spec)
{
    String<char> matrix = options.matrix;
    if (!empty(options.matrix))
    {
        Score<int, ScoreMatrix<> > sc;
        loadScoreMatrix(sc, toCString(options.matrix));
        _initAlignParams<AminoAcid>(options, sc, alg, spec);
    }
    else
    {
        Blosum62 sc;
        _initAlignParams<AminoAcid>(options, sc, alg, spec);
    }
}

// --------------------------------------------------------------------------
// Function pairAlignConfig()
// --------------------------------------------------------------------------

template <typename TAlgorithm, typename TSpec>
inline void pairAlignConfig(Options const & options, TAlgorithm const & alg, TSpec const & spec)
{
    // Initialize scoring matrices
    if (options.alphabet == "dna") _initScoreMatrix(options, Dna5(), alg, spec);
    else if (options.alphabet == "rna") _initScoreMatrix(options, Rna5(), alg, spec);
    else if (options.alphabet == "protein") _initScoreMatrix(options, AminoAcid(), alg, spec);
    else _initScoreMatrix(options, char(), alg, spec);
}

// ==========================================================================
// Declaration of external interfaces.
// ==========================================================================

void pairAlignLocal(Options const &);
void pairAlignLcs(Options const &);
void pairAlignGlobal_1111(Options const &);
void pairAlignGlobal_1110(Options const &);
void pairAlignGlobal_1101(Options const &);
void pairAlignGlobal_1100(Options const &);
void pairAlignGlobal_1011(Options const &);
void pairAlignGlobal_1010(Options const &);
void pairAlignGlobal_1001(Options const &);
void pairAlignGlobal_1000(Options const &);
void pairAlignGlobal_0111(Options const &);
void pairAlignGlobal_0110(Options const &);
void pairAlignGlobal_0101(Options const &);
void pairAlignGlobal_0100(Options const &);
void pairAlignGlobal_0011(Options const &);
void pairAlignGlobal_0010(Options const &);
void pairAlignGlobal_0001(Options const &);
void pairAlignGlobal_0000(Options const &);

inline void pairAlignMain(Options const & options)
{

    if (options.method == "lcs")
        pairAlignLcs(options);
    else if (options.method == "sw")
        pairAlignLocal(options);
    else
    {
        if (options.config == "tttt")
            pairAlignGlobal_1111(options);
        else if (options.config == "tttf")
            pairAlignGlobal_1110(options);
        else if (options.config == "ttft")
            pairAlignGlobal_1101(options);
        else if (options.config == "ttff")
            pairAlignGlobal_1100(options);
        else if (options.config == "tftt")
            pairAlignGlobal_1011(options);
        else if (options.config == "tftf")
            pairAlignGlobal_1010(options);
        else if (options.config == "tfft")
            pairAlignGlobal_1001(options);
        else if (options.config == "tfff")
            pairAlignGlobal_1000(options);
        else if (options.config == "fttt")
        	pairAlignGlobal_0111(options);
        else if (options.config == "fttf")
        	pairAlignGlobal_0110(options);
        else if (options.config == "ftft")
        	pairAlignGlobal_0101(options);
        else if (options.config == "ftff")
        	pairAlignGlobal_0100(options);
        else if (options.config == "fftt")
        	pairAlignGlobal_0011(options);
        else if (options.config == "fftf")
        	pairAlignGlobal_0010(options);
        else if (options.config == "ffft")
        	pairAlignGlobal_0001(options);
        else
            pairAlignGlobal_0000(options);
    }
}

#endif // APPS_PAIR_ALIGN_PAIR_ALIGN_LIB_H_
