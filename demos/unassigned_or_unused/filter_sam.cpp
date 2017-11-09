// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Filter the alignments in a SAM file.
//
// Exmple: Sort alignments by distance and limit the number of
// alignments to output per read to 100.
//
// filter_sam INPUT.sam -r REFERENCE.fasta --sort-distance --limit 100
//   -o OUT.sam
// ==========================================================================

#include <iostream>
#include <fstream>
#include <random>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/store.h>
#include <seqan/arg_parse.h>

using namespace seqan;

struct Options
{
    CharString samFilename;
    CharString referenceFilename;
    CharString outputFilename;
    bool randomTieBreak;
    bool sortDistance;
    unsigned limit;
    int seed;
};

struct SortAlignmentDistance_;
typedef Tag<SortAlignmentDistance_> SortAlignmentDistance;

template <typename TAlignedRead, typename TDistanceString>
struct LessAlignedReadDistance
{
    TDistanceString const & distances;

    LessAlignedReadDistance(TDistanceString const & _distances) :
        distances(_distances)
    {}

    inline bool
    operator()(TAlignedRead const & lhs, TAlignedRead const & rhs) const
    {
        return distances[lhs.id] < distances[rhs.id];
    }

};

template <typename TAlignedReadStore, typename TDistanceString>
void sortAlignedReads(TAlignedReadStore & alignedReadStore,
                      TDistanceString const & distances,
                      SortAlignmentDistance const &)
{
    std::stable_sort(
        begin(alignedReadStore, Standard()),
        end(alignedReadStore, Standard()),
        LessAlignedReadDistance<typename Value<TAlignedReadStore>::Type, TDistanceString>(distances));
}

void performWork(Options const & options)
{
    typedef FragmentStore<>::TContigSeq TContigSeq;
    FragmentStore<> fragmentStore;

    if (!empty(options.referenceFilename))
    {
        std::cerr << "Loading contigs..." << std::endl;
        loadContigs(fragmentStore, options.referenceFilename);
        std::cerr << "  loaded " << length(fragmentStore.contigStore) << " contigs" << std::endl;
    }
    std::cerr << "Loading alignments..." << std::endl;
    BamFileIn ins(toCString(options.samFilename));
    readRecords(fragmentStore, ins);
    std::cerr << "  loaded " << length(fragmentStore.alignedReadStore) << " alignments" << std::endl;

    // Compute distances with reference.
    std::cerr << "Computing distances..." << std::endl;
    String<int> distances;
    resize(distances, length(fragmentStore.alignedReadStore), Exact());
    typedef Iterator<FragmentStore<>::TAlignedReadStore, Standard>::Type TIterator;
    for (TIterator it = begin(fragmentStore.alignedReadStore), itEnd = end(fragmentStore.alignedReadStore); it != itEnd; ++it)
    {
        Align<TContigSeq> align;
        resize(rows(align), 2);
        size_t beginPos = it->beginPos;
        size_t endPos = it->endPos;
        if (beginPos > endPos)
            std::swap(beginPos, endPos);
        typedef Infix<TContigSeq>::Type TInfix;
        TInfix infixCopy(fragmentStore.contigStore[it->contigId].seq, beginPos, endPos);
        assignSource(row(align, 0), infixCopy);
        assignSource(row(align, 1), fragmentStore.readSeqStore[it->readId]);
        int s = globalAlignment(align, Score<int, EditDistance>(), NeedlemanWunsch());
        distances[it->id] = -s;
    }

    // Sort reads, ties are broken randomly.
    std::cerr << "Sorting alignments..." << std::endl;
    std::mt19937 rng;
    if (options.randomTieBreak)
        shuffle(fragmentStore.alignedReadStore, rng);
    if (options.sortDistance)
        sortAlignedReads(fragmentStore.alignedReadStore, distances, SortAlignmentDistance());
    sortAlignedReads(fragmentStore.alignedReadStore, SortReadId());

    // Copy over at most options.limit alignments per read.
    std::cerr << "Filtering reads..." << std::endl;
    FragmentStore<>::TAlignedReadStore rsCopy;
    size_t readId = std::numeric_limits<size_t>::max();
    size_t alignmentCount = 0;
    for (TIterator it = begin(fragmentStore.alignedReadStore), itEnd = end(fragmentStore.alignedReadStore); it != itEnd; ++it)
    {
        if (readId != it->readId)
        {
            alignmentCount = 1;
            readId = it->readId;
        }
        else
        {
            alignmentCount += 1;
        }
        if (alignmentCount <= options.limit)
            appendValue(rsCopy, *it);
    }

    // Copy the filtered reads back into fragment store.
    std::swap(fragmentStore.alignedReadStore, rsCopy);

    // Write out the fragment store.
    std::cerr << "Writing alignments..." << std::endl;
    BamFileOut outFile;
    if (empty(options.outputFilename))
        open(outFile, std::cout);
    else
        open(outFile, toCString(options.outputFilename));
    writeRecords(outFile, fragmentStore);
}

int main(int argc, char const ** argv)
{
    Options options;
    options.randomTieBreak = true;
    options.sortDistance = false;
    options.limit = 100;
    options.seed = 42;

    // Setup command line parser.
    ArgumentParser parser;
    addUsageLine(parser, "filter_sam [OPTIONS] INPUT.sam");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "INPUT"));
    addOption(parser, ArgParseOption("r", "reference", "Reference sequence in FASTA file.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("o", "output-filename", "Filename of the output.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("sd", "sort-distance", "Sort alignments by distance."));
    addOption(parser, ArgParseOption("l", "limit", "Number of alignments to output per read."));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    getOptionValue(options.referenceFilename, parser, "reference");
    getOptionValue(options.outputFilename, parser, "output-filename");
    getOptionValue(options.limit, parser, "limit");
    getOptionValue(options.sortDistance, parser, "sort-distance");

    getArgumentValue(options.samFilename, parser, 0);

    if (options.sortDistance && empty(options.referenceFilename))
    {
        std::cerr << "Reference has to be given if --sort-distance- is specified." << std::endl;
        return 1;
    }

    performWork(options);

    return 0;
}
