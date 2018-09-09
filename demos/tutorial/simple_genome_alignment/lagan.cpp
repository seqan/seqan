#include <cstdlib>

//![include_seq_io]
#include <seqan/seq_io.h>
#include <seqan/stream.h>
//![include_seq_io]

//![include_align]
#include <seqan/align.h>
//![include_align]
//![include_arg_parse]
#include <seqan/arg_parse.h>
//![include_arg_parse]
//![include_index]
#include <seqan/index.h>
//![include_index]
//![include_seeds]
#include <seqan/seeds.h>
//![include_seeds]

using namespace seqan;

std::stringstream quiet_output;

//![lagan_option]
struct LaganOption
{
    std::string seq1Filename;
    std::string seq2Filename;
    std::string outFilename;
    unsigned qGramSize;
    unsigned distanceCriteria;
    unsigned gapCriteria;

};
//![lagan_option]
/*
//![parse_arguments]
auto parseCommandLine(LaganOption & options, int const argc, char * argv[])
{
    ArgumentParser parser("lagan");

    // Set short description, version, and date.
    setShortDescription(parser, "Glocal alignment computation.");
    setVersion(parser, "0.1");
    setDate(parser, "September 2017");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fISEQUENCE_1\\fP\" \"\\fISEQUENCE_2\\fP\"");
    addDescription(parser,
                  "lagan is a large-scale sequence alignment tool "
                  "using a glocal alignment approach. It first filters for good seeding matches which are "
                  "chained together using CHAOS chaining and finally a global chain is computed and "
                  "the alignment around this chain.");

    // We require two arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "SEQUENCE_1"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "SEQUENCE_2"));

    addOption(parser, seqan::ArgParseOption("o", "output", "Output file to write the alignment to.",
                                            ArgParseArgument::OUTPUT_FILE, "FILE"));
    setValidValues(parser, "output", ".align");
    setDefaultValue(parser, "output", "out.align");

    addOption(parser, seqan::ArgParseOption("q", "qgram", "Size of the qGram.", ArgParseArgument::INTEGER, "INT"));
    addDefaultValue(parser, "qgram", "12");

    addOption(parser, seqan::ArgParseOption("d", "distance", "Distance criteria for CHAOS chaining.",
                                            ArgParseArgument::INTEGER, "INT"));
    addDefaultValue(parser, "distance", "10");

    addOption(parser, seqan::ArgParseOption("g", "gap", "Gap criteria for CHAOS chaining.",
                                            ArgParseArgument::INTEGER, "INT"));
    addDefaultValue(parser, "gap", "15");

    auto parseResult = parse(parser, argc, argv);

    if (parseResult == ArgumentParser::PARSE_OK)
    {
        getArgumentValue(options.seq1Filename, parser, 0);
        getArgumentValue(options.seq2Filename, parser, 1);
        getOptionValue(options.outFilename, parser, "output");
        getOptionValue(options.qGramSize, parser, "qgram");
        getOptionValue(options.distanceCriteria, parser, "distance");
        getOptionValue(options.gapCriteria, parser, "gap");
    }
    return parseResult;
}
//![parse_arguments]
*/

auto parseCommandLine(LaganOption & options, int const argc, char * argv[])
{
    ArgumentParser parser("lagan");

    // Set short description, version, and date.
    setShortDescription(parser, "Glocal alignment computation.");
    setVersion(parser, "0.1");
    setDate(parser, "September 2017");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fISEQUENCE_1\\fP\" \"\\fISEQUENCE_2\\fP\"");
    addDescription(parser,
                  "lagan is a large-scale sequence alignment tool "
                  "using a glocal alignment approach. It first filters for good seeding matches which are "
                  "chained together using CHAOS chaining and finally a global chain is computed and "
                  "the alignment around this chain.");

    // We require two arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "SEQUENCE_1"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "SEQUENCE_2"));

    addOption(parser, seqan::ArgParseOption("o", "output", "Output file to write the alignment to.",
                                            ArgParseArgument::OUTPUT_FILE, "FILE"));
    setValidValues(parser, "output", ".align");
    setDefaultValue(parser, "output", "out.align");

    addOption(parser, seqan::ArgParseOption("q", "qgram", "Size of the qGram.", ArgParseArgument::INTEGER, "INT"));
    addDefaultValue(parser, "qgram", "12");

    addOption(parser, seqan::ArgParseOption("d", "distance", "Distance criteria for CHAOS chaining.",
                                            ArgParseArgument::INTEGER, "INT"));
    addDefaultValue(parser, "distance", "10");

    addOption(parser, seqan::ArgParseOption("g", "gap", "Gap criteria for CHAOS chaining.",
                                            ArgParseArgument::INTEGER, "INT"));
    addDefaultValue(parser, "gap", "15");

    auto parseResult = parse(parser, argc, argv, quiet_output, quiet_output);

    if (parseResult == ArgumentParser::PARSE_OK)
    {
        getArgumentValue(options.seq1Filename, parser, 0);
        getArgumentValue(options.seq2Filename, parser, 1);
        getOptionValue(options.outFilename, parser, "output");
        getOptionValue(options.qGramSize, parser, "qgram");
        getOptionValue(options.distanceCriteria, parser, "distance");
        getOptionValue(options.gapCriteria, parser, "gap");
    }
    return parseResult;
}

/* load the sequences */
/*
//![load_sequence_template]
inline std::pair<Dna5String, CharString>
loadSequence(std::string const & fileName)
{
    // Write your code here.
}
//![load_sequence_template]
*/

//![load_sequence_solution]
inline std::pair<CharString, Dna5String>
loadSequence(std::string const & fileName)
{
    if (empty(fileName))
    {
        std::cerr << "Invalid file name!\n";
        std::exit(EXIT_FAILURE);
    }
    try {
        SeqFileIn seqFile(fileName.c_str());

        Dna5String seq;
        CharString seqId;
        readRecord(seqId, seq, seqFile);
        return {std::move(seqId), std::move(seq)};
    }
    catch (ParseError & error)
    {
        std::cerr << "Error while parsing " << fileName << "!\n";
        std::cerr << error.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
    catch (IOError & error)
    {
        std::cerr << "Could not read " << fileName << "!\n";
        std::cerr << error.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
}
//![load_sequence_solution]

/*
//![initial_main]
int main(int const argc, char * argv[])
{
    LaganOption options;
    if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
    {
        std::exit(EXIT_FAILURE);  // terminate the program.
    }
}
//![initial_main]
*/
int main(int const argc, char * argv[])
{
    LaganOption options;
    if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
    {
        std::exit(EXIT_SUCCESS);  // terminate the program.
    }

//![load_seq_1]
    /* load the sequences */
    CharString seq1Id;
    Dna5String seq1;

    std::tie(seq1Id, seq1) = loadSequence(options.seq1Filename);
    std::cout << "Loaded " << seq1Id << " (" << length(seq1) << "bp)\n";
//![load_seq_1]

//![load_seq_2]
    CharString seq2Id;
    Dna5String seq2;

    std::tie(seq2Id, seq2) = loadSequence(options.seq2Filename);
    std::cout << "Loaded " << seq2Id << " (" << length(seq2) << "bp)\n";
//![load_seq_2]

//![declare_index]
    using TIndex    = Index<Dna5String, IndexQGram<Shape<Dna5, SimpleShape>, OpenAddressing>>;
//![declare_index]

//![declare_seed_set]
    using TSeed     = Seed<ChainedSeed>;
    using TSeedSize = typename Size<TSeed>::Type;
    using TSeedSet  = SeedSet<TSeed>;
//![declare_seed_set]

//![init_seed_set]
    TSeedSet seedSet;
    Score<int, Simple> scoreScheme{2, -1, -2};
//![init_seed_set]

//![init_index]
    TIndex qGramIndex{seq1};

    // Initialize the shape.
    resize(indexShape(qGramIndex), options.qGramSize);
    hashInit(indexShape(qGramIndex), begin(seq2, Standard()));
//![init_index]

/*
//![solution_assignment2]
for (auto seqIter = begin(seq2, Standard());
     seqIter != end(seq2, Standard()) - length(indexShape(qGramIndex)) + 1;
     ++seqIter)
{
    hashNext(indexShape(qGramIndex), seqIter);
    std::cout << getOccurrences(qGramIndex, indexShape(qGramIndex));
}
//![solution_assignment2]
*/

//![solution_seeding]
    for (auto seqIter = begin(seq2, Standard());
         seqIter != end(seq2, Standard()) - length(indexShape(qGramIndex)) + 1;
         ++seqIter)
    {
        hashNext(indexShape(qGramIndex), seqIter);
        for (TSeedSize seq1Pos : getOccurrences(qGramIndex, indexShape(qGramIndex)))
        {
            /* filter the seeds */
            if (!addSeed(seedSet,
                         TSeed{seq1Pos, static_cast<TSeedSize>(seqIter - begin(seq2, Standard())), options.qGramSize},
                         options.distanceCriteria, options.gapCriteria, scoreScheme, seq1, seq2, Chaos{}))
            {
                addSeed(seedSet,
                        TSeed{seq1Pos, static_cast<TSeedSize>(seqIter - begin(seq2, Standard())), options.qGramSize},
                        Single{});
            }
        }
    }
    std::cout << "Finished seeding: " << length(seedSet) << " seeds!\n";
//![solution_seeding]

//![chain_seeds]
    /* chain seeds */
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    std::cout << "Finished chaining: " << length(seedChain) << " seeds!\n";
//![chain_seeds]

//![build_alignment]
    /* build alignment */
    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seq1);
    assignSource(row(alignment, 1), seq2);

    int score = bandedChainAlignment(alignment, seedChain, Score<int, Simple>{4, -5, -1, -11}, 15);
    std::cout << "Finished alignment: Score = " << score << "\n";
//![build_alignment]

//![output_alignment]
    /* output alignment */
    std::ofstream outStream;
    outStream.open(options.outFilename.c_str());
    if (!outStream.good())
    {
        std::cerr << "Could not open " << options.outFilename << "!\n";
        std::exit(EXIT_FAILURE);
    }

    outStream << "#seq1: " << seq1Id << "\n";
    outStream << "#seq2: " << seq2Id << "\n";
    outStream << "#score: " << score << "\n";
    outStream << "#alignment\n";
    outStream << alignment;

    std::exit(EXIT_SUCCESS);
//![output_alignment]
}
