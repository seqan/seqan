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

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/random.h>
#include <seqan/misc/misc_cmdparser.h>

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

    LessAlignedReadDistance(TDistanceString const & _distances)
            : distances(_distances)
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

    if (!empty(options.referenceFilename)) {
        std::cerr << "Loading contigs..." << std::endl;
        loadContigs(fragmentStore, options.referenceFilename);
        std::cerr << "  loaded " << length(fragmentStore.contigStore) << " contigs" << std::endl;
    }
    std::cerr << "Loading alignments..." << std::endl;
    std::ifstream ins(toCString(options.samFilename));
    read(ins, fragmentStore, Sam());
    std::cerr << "  loaded " << length(fragmentStore.alignedReadStore) << " alignments" << std::endl;

    // Compute distances with reference.
    std::cerr << "Computing distances..." << std::endl;
    String<int> distances;
    resize(distances, length(fragmentStore.alignedReadStore), Exact());
    typedef Iterator<FragmentStore<>::TAlignedReadStore, Standard>::Type TIterator;
    for (TIterator it = begin(fragmentStore.alignedReadStore), itEnd = end(fragmentStore.alignedReadStore); it != itEnd; ++it) {
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
    Rng<> rng;
    if (options.randomTieBreak)
        shuffle(fragmentStore.alignedReadStore, rng);
    if (options.sortDistance)
        sortAlignedReads(fragmentStore.alignedReadStore, distances, SortAlignmentDistance());
    sortAlignedReads(fragmentStore.alignedReadStore, SortReadId());

    // Copy over at most options.limit alignments per read.
    std::cerr << "Filtering reads..." << std::endl;
    FragmentStore<>::TAlignedReadStore rsCopy;
    size_t readId = MaxValue<size_t>::VALUE;
    size_t alignmentCount = 0;
    for (TIterator it = begin(fragmentStore.alignedReadStore), itEnd = end(fragmentStore.alignedReadStore); it != itEnd; ++it) {
        if (readId != it->readId) {
            alignmentCount = 1;
            readId = it->readId;
        } else {
            alignmentCount += 1;
        }
        if (alignmentCount <= options.limit)
            appendValue(rsCopy, *it);
    }

    // Copy the filtered reads back into fragment store.
    std::swap(fragmentStore.alignedReadStore, rsCopy);

    // Write out the fragment store.
    std::cerr << "Writing alignments..." << std::endl;
    std::ostream * outStream;
    if (empty(options.outputFilename))
        outStream = &std::cout;
    else
        outStream = new std::ofstream(toCString(options.outputFilename));
    write(*outStream, fragmentStore, Sam());
    if (empty(options.outputFilename))
        delete outStream;
}

int main(int argc, char const ** argv)
{
    Options options;
    options.randomTieBreak = true;
    options.sortDistance = false;
    options.limit = 100;
    options.seed = 42;

    // Setup command line parser.
    CommandLineParser parser;
    addUsageLine(parser, "filter_sam [OPTIONS] INPUT.sam");
    addOption(parser, CommandLineOption("r", "reference", "Reference sequence in FASTA file.", OptionType::String, options.referenceFilename));
    addOption(parser, CommandLineOption("o", "output-filename", "Filename of the output.", OptionType::String, options.outputFilename));
    addOption(parser, CommandLineOption("sd", "sort-distance", "Sort alignments by distance.", OptionType::Bool, options.sortDistance));
    addOption(parser, CommandLineOption("l", "limit", "Number of alignments to output per read.", OptionType::Int, options.limit));
    requiredArguments(parser, 1);

    // Parse command line.
    bool stop = !parse(parser, argc, argv, std::cerr);
    if (stop)
        return !isSetShort(parser, "h");
    getOptionValueLong(parser, "reference", options.referenceFilename);
    getOptionValueLong(parser, "output-filename", options.outputFilename);
    getOptionValueLong(parser, "limit", options.limit);
    if (isSetLong(parser, "sort-distance"))
        options.sortDistance = true;
    options.samFilename = getArgumentValue(parser, 0);

    int res = 0;
    if (options.sortDistance && empty(options.referenceFilename)) {
        std::cerr << "Reference has to be given if --sort-distance- is specified." << std::endl;
        stop = true;
        res = 1;
    }

    if (!stop)
        performWork(options);

    return res;
}
