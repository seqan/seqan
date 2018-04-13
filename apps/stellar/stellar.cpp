// ==========================================================================
//                    STELLAR - SwifT Exact LocaL AligneR
//                   http://www.seqan.de/projects/stellar/
// ==========================================================================
// Copyright (C) 2010-2012 by Birte Kehr
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

#include "stellar.h"
#include "stellar_output.h"

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
// Initializes a Finder object for a database sequence,
//  calls stellar, and writes matches to file
template <typename TSequence, typename TId, typename TPattern, typename TMatches>
inline bool
_stellarOnOne(TSequence & database,
              TId & databaseID,
              TPattern & swiftPattern,
              bool databaseStrand,
              TMatches & matches,
              StellarOptions & options)
{
    std::cout << "  " << databaseID;
    if (!databaseStrand)
        std::cout << ", complement";
    std::cout << std::flush;

    // finder
    typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
    TFinder swiftFinder(database, options.minRepeatLength, options.maxRepeatPeriod);

    // stellar
    if (options.fastOption == CharString("exact"))
        stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop,
                options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
                databaseID, databaseStrand, matches, AllLocal());
    else if (options.fastOption == "bestLocal")
        stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop,
                options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
                databaseID, databaseStrand, matches, BestLocal());
    else if (options.fastOption == "bandedGlobal")
        stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop,
                options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
                databaseID, databaseStrand, matches, BandedGlobal());
    else if (options.fastOption == "bandedGlobalExtend")
        stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop,
                options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
                databaseID, databaseStrand, matches, BandedGlobalExtend());
    else
    {
        std::cerr << "\nUnknown verification strategy: " << options.fastOption << std::endl;
        return false;
    }

    std::cout << std::endl;
    return true;
}

//////////////////////////////////////////////////////////////////////////////
namespace seqan {

template <typename TStringSet, typename TShape, typename TSpec>
struct Cargo<Index<TStringSet, IndexQGram<TShape, TSpec> > >
{
    typedef struct
    {
        double      abundanceCut;
    } Type;
};

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TStringSet, typename TShape, typename TSpec>
inline bool _qgramDisableBuckets(Index<TStringSet, IndexQGram<TShape, TSpec> > & index)
{
    typedef typename Fibre<TStringSet, QGramDir>::Type      TDir;
    typedef typename Iterator<TDir, Standard>::Type         TDirIterator;
    typedef typename Value<TDir>::Type                      TSize;

    TDir & dir    = indexDir(index);
    bool result  = false;
    unsigned counter = 0;
    TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
    if (thresh < 100)
        thresh = 100;

    TDirIterator it = begin(dir, Standard());
    TDirIterator itEnd = end(dir, Standard());
    for (; it != itEnd; ++it)
        if (*it > thresh)
        {
            *it = (TSize) - 1;
            result = true;
            ++counter;
        }

    if (counter > 0)
        std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

    return result;
}

template <>
struct FunctorComplement<AminoAcid>:
    public::std::unary_function<AminoAcid, AminoAcid>
{
    inline AminoAcid operator()(AminoAcid x) const
    {
        return x;
    }

};

}

///////////////////////////////////////////////////////////////////////////////
// Initializes a Pattern object with the query sequences,
//  and calls _stellarOnOne for each database sequence
template <typename TSequence, typename TId>
inline bool
_stellarOnAll(StringSet<TSequence> & databases,
              StringSet<TId> & databaseIDs,
              StringSet<TSequence> & queries,
              StringSet<TId> & queryIDs,
              StellarOptions & options)
{
    // pattern
    typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    TQGramIndex qgramIndex(queries);
    resize(indexShape(qgramIndex), options.qGram);
    cargo(qgramIndex).abundanceCut = options.qgramAbundanceCut;
    Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

    if (options.verbose)
        swiftPattern.params.printDots = true;

    // Construct index
    std::cout << "Constructing index..." << std::endl;
    indexRequire(qgramIndex, QGramSADir());
    std::cout << std::endl;

    // container for eps-matches
    StringSet<QueryMatches<StellarMatch<TSequence, TId> > > matches;
    resize(matches, length(queries));

    std::cout << "Aligning all query sequences to database sequence..." << std::endl;
    for (unsigned i = 0; i < length(databases); ++i)
    {
        // positive database strand
        if (options.forward)
        {
            if (!_stellarOnOne(databases[i], databaseIDs[i], swiftPattern, true, matches, options))
                return 1;
        }
        // negative (reverse complemented) database strand
        if (options.reverse && options.alphabet != "protein" && options.alphabet != "char")
        {
            reverseComplement(databases[i]);
            if (!_stellarOnOne(databases[i], databaseIDs[i], swiftPattern, false, matches, options))
                return 1;

            reverseComplement(databases[i]);
        }
    }
    std::cout << std::endl;

    // file output
    if (options.disableThresh != std::numeric_limits<unsigned>::max())
    {
        if (!_outputMatches(matches, queries, queryIDs, databases, options.verbose,
                            options.outputFile, options.outputFormat, options.disabledQueriesFile))
            return 1;
    }
    else
    {
        if (!_outputMatches(matches, queryIDs, databases, options.verbose,
                            options.outputFile, options.outputFormat))
            return 1;
    }

    return 0;
}

template <typename TId>
inline bool
_checkUniqueId(std::set<TId> & uniqueIds, TId & id)
{
    TId shortId;
    typedef typename Iterator<TId>::Type TIterator;

    TIterator it = begin(id);
    TIterator itEnd = end(id);
    while (it != itEnd && *it > 32)
    {
        appendValue(shortId, *it);
        ++it;
    }

    if (uniqueIds.count(shortId) == 0)
    {
        uniqueIds.insert(shortId);
        return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Imports sequences from a file,
//  stores them in the StringSet seqs and their identifiers in the StringSet ids
template <typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileName,
                 CharString const & name,
                 StringSet<TSequence> & seqs,
                 StringSet<TId> & ids)
{
    SeqFileIn inSeqs;
    if (!open(inSeqs, (toCString(fileName))))
    {
        std::cerr << "Failed to open " << name << " file." << std::endl;
        return false;
    }

    std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
    bool idsUnique = true;

    TSequence seq;
    TId id;
    unsigned seqCount = 0;
    for (; !atEnd(inSeqs); ++seqCount)
    {
        readRecord(id, seq, inSeqs);

        idsUnique &= _checkUniqueId(uniqueIds, id);

        appendValue(seqs, seq, Generous());
        appendValue(ids, id, Generous());
    }

    std::cout << "Loaded " << seqCount << " " << name << " sequence" << ((seqCount > 1) ? "s." : ".") << std::endl;
    if (!idsUnique)
        std::cerr << "WARNING: Non-unique " << name << " ids. Output can be ambiguous.\n";
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and from sequences and writes them to std::cout
template <typename TStringSet>
void _writeMoreCalculatedParams(StellarOptions & options, TStringSet & databases, TStringSet & queries)
{
//IOREV _notio_
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;

    if (options.qgramAbundanceCut != 1)
    {
        std::cout << "Calculated parameters:" << std::endl;
    }

    TSize queryLength = length(concat(queries));
    if (options.qgramAbundanceCut != 1)
    {
        std::cout << "  q-gram expected abundance : ";
        std::cout << queryLength / (double)((long)1 << (options.qGram << 1)) << std::endl;
        std::cout << "  q-gram abundance threshold: ";
        std::cout << _max(100, (int)(queryLength * options.qgramAbundanceCut)) << std::endl;
        std::cout << std::endl;
    }

    if (IsSameType<TAlphabet, Dna5>::VALUE || IsSameType<TAlphabet, Rna5>::VALUE)
    {
        // Computation of maximal E-value for this search

        TSize maxLengthQueries = 0;
        TSize maxLengthDatabases = 0;

        typename Iterator<TStringSet>::Type dbIt = begin(databases);
        typename Iterator<TStringSet>::Type dbEnd = end(databases);
        while (dbIt != dbEnd)
        {
            if (length(*dbIt) > maxLengthDatabases)
            {
                maxLengthDatabases = length(*dbIt);
            }
            ++dbIt;
        }

        typename Iterator<TStringSet>::Type queriesIt = begin(queries);
        typename Iterator<TStringSet>::Type queriesEnd = end(queries);
        while (queriesIt != queriesEnd)
        {
            if (length(*queriesIt) > maxLengthQueries)
            {
                maxLengthQueries = length(*queriesIt);
            }
            ++queriesIt;
        }

        TSize errors = static_cast<TSize>(options.minLength * options.epsilon);
        TSize minScore = options.minLength - 3 * errors; // #matches - 2*#errors // #matches = minLenght - errors,

        std::cout << "All matches matches resulting from your search have an E-value of: " << std::endl;
        std::cout << "        " << _computeEValue(minScore, maxLengthQueries, maxLengthDatabases) << " or smaller";
        std::cout << "  (match score = 1, error penalty = -2)" << std::endl;

        std::cout << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and writes them to std::cout
void _writeCalculatedParams(StellarOptions & options)
{
//IOREV _notio_
    int errMinLen = (int) floor(options.epsilon * options.minLength);
    int n = (int) ceil((errMinLen + 1) / options.epsilon);
    int errN = (int) floor(options.epsilon * n);
    unsigned smin = (unsigned) _min(ceil((double)(options.minLength - errMinLen) / (errMinLen + 1)),
                                    ceil((double)(n - errN) / (errN + 1)));

    std::cout << "Calculated parameters:" << std::endl;
    if (options.qGram == (unsigned)-1)
    {
        options.qGram = (unsigned)_min(smin, 32u);
        std::cout << "  k-mer length : " << options.qGram << std::endl;
    }

    int threshold = (int) _max(1, (int) _min((n + 1) - options.qGram * (errN + 1),
                                             (options.minLength + 1) - options.qGram * (errMinLen + 1)));
    int overlap = (int) floor((2 * threshold + options.qGram - 3) / (1 / options.epsilon - options.qGram));
    int distanceCut = (threshold - 1) + options.qGram * overlap + options.qGram;
    int logDelta = _max(4, (int) ceil(log((double)overlap + 1) / log(2.0)));
    int delta = 1 << logDelta;

    std::cout << "  s^min        : " << smin << std::endl;
    std::cout << "  threshold    : " << threshold << std::endl;
    std::cout << "  distance cut : " << distanceCut << std::endl;
    std::cout << "  delta        : " << delta << std::endl;
    std::cout << "  overlap      : " << overlap << std::endl;
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes user specified parameters from options object to std::cout
template <typename TOptions>
void
_writeSpecifiedParams(TOptions & options)
{
//IOREV _notio_
    // Output user specified parameters
    std::cout << "User specified parameters:" << std::endl;
    std::cout << "  minimal match length             : " << options.minLength << std::endl;
    std::cout << "  maximal error rate (epsilon)     : " << options.epsilon << std::endl;
    std::cout << "  maximal x-drop                   : " << options.xDrop << std::endl;
    if (options.qGram != (unsigned)-1)
        std::cout << "  k-mer (q-gram) length            : " << options.qGram << std::endl;
    std::cout << "  search forward strand            : " << ((options.forward) ? "yes" : "no") << std::endl;
    std::cout << "  search reverse complement        : " << ((options.reverse) ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    std::cout << "  verification strategy            : " << options.fastOption << std::endl;
    if (options.disableThresh != (unsigned)-1)
    {
        std::cout << "  disable queries with more than   : " << options.disableThresh << " matches" << std::endl;
    }
    std::cout << "  maximal number of matches        : " << options.numMatches << std::endl;
    std::cout << "  duplicate removal every          : " << options.compactThresh << std::endl;
    if (options.maxRepeatPeriod != 1 || options.minRepeatLength != 1000)
    {
        std::cout << "  max low complexity repeat period : " << options.maxRepeatPeriod << std::endl;
        std::cout << "  min low complexity repeat length : " << options.minRepeatLength << std::endl;
    }
    if (options.qgramAbundanceCut != 1)
    {
        std::cout << "  q-gram abundance cut ratio       : " << options.qgramAbundanceCut << std::endl;
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes file name from options object to std::cout
template <typename TOptions>
void
_writeFileNames(TOptions & options)
{
//IOREV _notio_
    std::cout << "I/O options:" << std::endl;
    std::cout << "  database file   : " << options.databaseFile << std::endl;
    std::cout << "  query file      : " << options.queryFile << std::endl;
    std::cout << "  alphabet        : " << options.alphabet << std::endl;
    std::cout << "  output file     : " << options.outputFile << std::endl;
    std::cout << "  output format   : " << options.outputFormat << std::endl;
    if (options.disableThresh != (unsigned)-1)
    {
        std::cout << "  disabled queries: " << options.disabledQueriesFile << std::endl;
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
template <typename TOptions>
ArgumentParser::ParseResult
_parseOptions(ArgumentParser & parser, TOptions & options)
{
    getArgumentValue(options.databaseFile, parser, 0);
    getArgumentValue(options.queryFile, parser, 1);

    // output options
    getOptionValue(options.outputFile, parser, "out");
    getOptionValue(options.disabledQueriesFile, parser, "outDisabled");
    getOptionValue(options.noRT, parser, "no-rt");

    CharString tmp = options.outputFile;
    toLower(tmp);

    if (endsWith(tmp, ".gff"))
        options.outputFormat = "gff";
    else if (endsWith(tmp, ".txt"))
        options.outputFormat = "txt";

    // main options
    getOptionValue(options.qGram, parser, "kmer");
    getOptionValue(options.minLength, parser, "minLength");
    getOptionValue(options.epsilon, parser, "epsilon");
    getOptionValue(options.xDrop, parser, "xDrop");
    getOptionValue(options.alphabet, parser, "alphabet");

    if (isSet(parser, "forward") && !isSet(parser, "reverse"))
        options.reverse = false;
    if (!isSet(parser, "forward") && isSet(parser, "reverse"))
        options.forward = false;

    getOptionValue(options.fastOption, parser, "verification");
    getOptionValue(options.disableThresh, parser, "disableThresh");
    getOptionValue(options.numMatches, parser, "numMatches");
    getOptionValue(options.compactThresh, parser, "sortThresh");
    getOptionValue(options.maxRepeatPeriod, parser, "repeatPeriod");
    getOptionValue(options.minRepeatLength, parser, "repeatLength");
    getOptionValue(options.qgramAbundanceCut, parser, "abundanceCut");

    getOptionValue(options.verbose, parser, "verbose");

    if (isSet(parser, "kmer") && options.qGram >= 1 / options.epsilon)
    {
        std::cerr << "Invalid parameter value: Please choose q-gram length lower than 1/epsilon." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    if (options.numMatches > options.compactThresh)
    {
        std::cerr << "Invalid parameter values: Please choose numMatches <= sortThresh." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

///////////////////////////////////////////////////////////////////////////////
// Set-Up of Argument Parser
void _setParser(ArgumentParser & parser)
{
    setShortDescription(parser, "the SwifT Exact LocaL AligneR");
    setDate(parser, SEQAN_DATE);
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setCategory(parser, "Local Alignment");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIFASTA FILE 1\\fP> <\\fIFASTA FILE 2\\fP>");

    addDescription(parser,
                   "STELLAR implements the SWIFT filter algorithm (Rasmussen et al., 2006) "
                   "and a verification step for the SWIFT hits that applies local alignment, "
                   "gapped X-drop extension, and extraction of the longest epsilon-match.");
    addDescription(parser,
                   "Input to STELLAR are two files, each containing one or more sequences "
                   "in FASTA format. Each sequence from file 1 will be compared to each "
                   "sequence in file 2. The sequences from file 1 are used as database, the "
                   "sequences from file 2 as queries.");
    addDescription(parser, "(c) 2010-2012 by Birte Kehr");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE 1"));
    setValidValues(parser, 0, "fa fasta");  // allow only fasta files as input
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE 2"));
    setValidValues(parser, 1, "fa fasta");  // allow only fasta files as input

    addSection(parser, "Main Options");

    addOption(parser, ArgParseOption("e", "epsilon", "Maximal error rate (max 0.25).", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e", "0.05");
    setMinValue(parser, "e", "0.0000001");
    setMaxValue(parser, "e", "0.25");
    addOption(parser, ArgParseOption("l", "minLength", "Minimal length of epsilon-matches.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "l", "100");
    setMinValue(parser, "l", "0");
    addOption(parser, ArgParseOption("f", "forward", "Search only in forward strand of database."));
    addOption(parser, ArgParseOption("r", "reverse", "Search only in reverse complement of database."));
    addOption(parser, ArgParseOption("a", "alphabet",
                                     "Alphabet type of input sequences (dna, rna, dna5, rna5, protein, char).",
                                     ArgParseArgument::STRING));
    setValidValues(parser, "a", "dna dna5 rna rna5 protein char");
    addOption(parser, ArgParseOption("v", "verbose", "Set verbosity mode."));

    addSection(parser, "Filtering Options");

    addOption(parser, ArgParseOption("k", "kmer", "Length of the q-grams (max 32).", ArgParseArgument::INTEGER));
    setMinValue(parser, "k", "1");
    setMaxValue(parser, "k", "32");
    addOption(parser, ArgParseOption("rp", "repeatPeriod",
                                     "Maximal period of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "rp", "1");
    addOption(parser, ArgParseOption("rl", "repeatLength",
                                     "Minimal length of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "rl", "1000");
    addOption(parser, ArgParseOption("c", "abundanceCut", "k-mer overabundance cut ratio.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "c", "1");
    setMinValue(parser, "c", "0");
    setMaxValue(parser, "c", "1");

    addSection(parser, "Verification Options");

    addOption(parser, ArgParseOption("x", "xDrop", "Maximal x-drop for extension.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "x", "5");
    addOption(parser, ArgParseOption("vs", "verification", "Verification strategy: exact or bestLocal or bandedGlobal",
                                     ArgParseArgument::STRING));
    //addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
    //addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
    //addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");
    setDefaultValue(parser, "vs", "exact");
    setValidValues(parser, "vs", "exact bestLocal bandedGlobal");
    addOption(parser, ArgParseOption("dt", "disableThresh",
                                     "Maximal number of verified matches before disabling verification for one query "
                                     "sequence (default infinity).", ArgParseArgument::INTEGER));
    setMinValue(parser, "dt", "0");
    addOption(parser, ArgParseOption("n", "numMatches",
                                     "Maximal number of kept matches per query and database. If STELLAR finds more matches, "
                                     "only the longest ones are kept.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "n", "50");
    addOption(parser, ArgParseOption("s", "sortThresh",
                                     "Number of matches triggering removal of duplicates. Choose a smaller value for saving "
                                     "space.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "s", "500");

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "out", "Name of output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "o", "gff txt");
    setDefaultValue(parser, "o", "stellar.gff");
    addOption(parser, ArgParseOption("od", "outDisabled",
                                     "Name of output file for disabled query sequences.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "outDisabled", seqan::SeqFileOut::getFileExtensions());
    setDefaultValue(parser, "od", "stellar.disabled.fasta");
    addOption(parser, ArgParseOption("t", "no-rt", "Suppress printing running time."));
    hideOption(parser, "t");

    addTextSection(parser, "References");
    addText(parser, "Kehr, B., Weese, D., Reinert, K.: STELLAR: fast and exact local alignments. BMC Bioinformatics, "
                    "12(Suppl 9):S15, 2011.");
}

// TODO(holtgrew): Move this into a SeqAn misc module.

// not supported anymore in vc2015
// https://msdn.microsoft.com/en-us/library/bb531344.aspx
class ScientificNotationExponentOutputNormalizer
{
public:
    unsigned _oldExponentFormat;

    ScientificNotationExponentOutputNormalizer() :
        _oldExponentFormat(0)
    {
    }

    ~ScientificNotationExponentOutputNormalizer()
    {
    }
};

///////////////////////////////////////////////////////////////////////////////
// Parses and outputs parameters, calls _stellarOnAll().
template <typename TOptions, typename TAlphabet>
int mainWithOptions(TOptions & options, String<TAlphabet>)
{
    typedef String<TAlphabet> TSequence;

    // output file names
    _writeFileNames(options);

    // output parameters
    _writeSpecifiedParams(options);
    _writeCalculatedParams(options);

    // import query sequences
    StringSet<TSequence> queries;
    StringSet<CharString> queryIDs;
    if (!_importSequences(options.queryFile, "query", queries, queryIDs))
        return 1;

    // import database sequence
    StringSet<TSequence> databases;
    StringSet<CharString> databaseIDs;
    if (!_importSequences(options.databaseFile, "database", databases, databaseIDs))
        return 1;

    std::cout << std::endl;
    _writeMoreCalculatedParams(options, databases, queries);

    // open output files
    std::ofstream file;
    file.open(toCString(options.outputFile));
    if (!file.is_open())
    {
        std::cerr << "Could not open output file." << std::endl;
        return 1;
    }
    file.close();

    if (options.disableThresh != std::numeric_limits<unsigned>::max())
    {
        std::ofstream daFile;
        daFile.open(toCString(options.disabledQueriesFile));
        if (!daFile.is_open())
        {
            std::cerr << "Could not open file for disabled queries." << std::endl;
            return 1;
        }
        daFile.close();
    }

    // stellar on all databases and queries writing results to file

    double startTime = sysTime();
    if (!_stellarOnAll(databases, databaseIDs, queries, queryIDs, options))
        return 1;

    if (options.verbose && options.noRT == false)
        std::cout << "Running time: " << sysTime() - startTime << "s" << std::endl;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Program entry point.
int main(int argc, const char * argv[])
{
    // Makes sure that printing doubles in scientific notation is normalized on all platforms.
    ScientificNotationExponentOutputNormalizer scientificNotationNormalizer;

    // command line parsing
    ArgumentParser parser("stellar");

    StellarOptions options = StellarOptions();
    _setParser(parser);
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == ArgumentParser::PARSE_OK)
        res = _parseOptions(parser, options);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    if (options.alphabet == "dna")
        mainWithOptions(options, String<Dna>());
    else if (options.alphabet == "dna5")
        mainWithOptions(options, String<Dna5>());
    else if (options.alphabet == "rna")
        mainWithOptions(options, String<Rna>());
    else if (options.alphabet == "rna5")
        mainWithOptions(options, String<Rna5>());
    else if (options.alphabet == "protein")
        mainWithOptions(options, String<AminoAcid>());
    else if (options.alphabet == "char")
        mainWithOptions(options, String<char>());
}
