#include <iostream>
#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/interval_tree.h>
#include <seqan/parallel.h>

using namespace seqan;


// define used types
typedef FragmentStore<>                         TStore;
typedef Value<TStore::TAnnotationStore>::Type   TAnnotation;
typedef TAnnotation::TId                        TId;
typedef TAnnotation::TPos                       TPos;
typedef IntervalAndCargo<TPos, TId>             TInterval;
typedef IntervalTree<TPos, TId>                 TIntervalTree;
typedef Value<TStore::TAlignedReadStore>::Type  TAlignedRead;

//
// 1. Load annotations and alignments from files
//
bool loadFiles(TStore & store, std::string const & annotationFileName,  std::string const & alignmentFileName)
{
    BamFileIn alignmentFile;
    if (!open(alignmentFile, alignmentFileName.c_str()))
    {
        std::cerr << "Couldn't open alignment file " << alignmentFileName << std::endl;
        return false;
    }
    std::cout << "Loading read alignments ..... " << std::flush;
    readRecords(store, alignmentFile);
    std::cout << "[" << length(store.alignedReadStore) << "]" << std::endl;

    // load annotations
    GffFileIn annotationFile;
    if (!open(annotationFile, toCString(annotationFileName)))
    {
        std::cerr << "Couldn't open annotation file" << annotationFileName << std::endl;
        return false;
    }
    std::cout << "Loading genome annotation ... " << std::flush;
    readRecords(store, annotationFile);
    std::cout << "[" << length(store.annotationStore) << "]" << std::endl;

    return true;
}

//
// 2. Extract intervals from gene annotations (grouped by contigId)
//
void extractGeneIntervals(String<String<TInterval> > & intervals, TStore const & store)
{
    // extract intervals from gene annotations (grouped by contigId)
    resize(intervals, length(store.contigStore));

    Iterator<TStore const, AnnotationTree<> >::Type it = begin(store, AnnotationTree<>());

    if (!goDown(it))
        return;

    do
    {
        SEQAN_ASSERT_EQ(getType(it), "gene");
        TPos beginPos = getAnnotation(it).beginPos;
        TPos endPos = getAnnotation(it).endPos;
        TId contigId = getAnnotation(it).contigId;

        if (beginPos > endPos)
            std::swap(beginPos, endPos);

        // insert forward-strand interval of the gene and its annotation id
        appendValue(intervals[contigId], TInterval(beginPos, endPos, value(it)));
    }
    while (goRight(it));
}

//
// 3. Construct interval trees
//
void constructIntervalTrees(String<TIntervalTree> & intervalTrees,
                            String<String<TInterval> > & intervals)
{
    int numContigs = length(intervals);
    resize(intervalTrees, numContigs);

    SEQAN_OMP_PRAGMA(parallel for)
    for (int i = 0; i < numContigs; ++i)
        createIntervalTree(intervalTrees[i], intervals[i]);
}

//
// 4. Count reads per gene
//
void countReadsPerGene(String<unsigned> & readsPerGene, String<TIntervalTree> const & intervalTrees, TStore const & store)
{
    resize(readsPerGene, length(store.annotationStore), 0);
    String<TId> result;
    int numAlignments = length(store.alignedReadStore);

    // iterate aligned reads and get search their begin and end positions
    SEQAN_OMP_PRAGMA(parallel for private (result))
    for (int i = 0; i < numAlignments; ++i)
    {
        TAlignedRead const & ar = store.alignedReadStore[i];
        TPos queryBegin = _min(ar.beginPos, ar.endPos);
        TPos queryEnd = _max(ar.beginPos, ar.endPos);

        // search read-overlapping genes
        findIntervals(result, intervalTrees[ar.contigId], queryBegin, queryEnd);

        // increase read counter for each overlapping annotation given the id in the interval tree
        for (unsigned j = 0; j < length(result); ++j)
        {
            SEQAN_OMP_PRAGMA(atomic)
            readsPerGene[result[j]] += 1;
        }
    }
}

//![solution]
//
// 5. Output RPKM values
//
void outputGeneCoverage(String<unsigned> const & readsPerGene, TStore const & store)
{
    // output abundances for covered genes
    Iterator<TStore const, AnnotationTree<> >::Type transIt = begin(store, AnnotationTree<>());
    Iterator<TStore const, AnnotationTree<> >::Type exonIt;
    double millionMappedReads = length(store.alignedReadStore) / 1000000.0;

    std::cout << "#gene name\tRPKM value" << std::endl;
    for (unsigned j = 0; j < length(readsPerGene); ++j)
    {
        if (readsPerGene[j] == 0)
            continue;

        unsigned mRNALengthMax = 0;
        goTo(transIt, j);

        // determine maximal mRNA length (which we use as gene length)
        SEQAN_ASSERT_NOT(isLeaf(transIt));
        goDown(transIt);

        do
        {
            exonIt = nodeDown(transIt);
            unsigned mRNALength = 0;

            // determine mRNA length, sum up the lengths of its exons
            do
            {
                if (getAnnotation(exonIt).typeId == store.ANNO_EXON)
                    mRNALength += abs((int)getAnnotation(exonIt).beginPos - (int)getAnnotation(exonIt).endPos);
            }
            while (goRight(exonIt));

            if (mRNALengthMax < mRNALength)
                mRNALengthMax = mRNALength;
        }
        while (goRight(transIt));

        // RPKM is number of reads mapped to a gene divided by its gene length in kbps
        // and divided by millions of total mapped reads
        std::cout << store.annotationNameStore[j] << '\t';
        std::cout << readsPerGene[j] / (mRNALengthMax / 1000.0) / millionMappedReads << std::endl;
    }
}
//![solution]

//![main]
int main(int argc, char const * argv[])
{
    TStore store;
    String<String<TInterval> > intervals;
    String<TIntervalTree> intervalTrees;
    String<unsigned> readsPerGene;

    std::string annotationFileName = getAbsolutePath("demos/tutorial/simple_rna_seq/example.gtf");
    std::string alignmentFileName = getAbsolutePath("demos/tutorial/simple_rna_seq/example.sam");

    if (!loadFiles(store, annotationFileName, alignmentFileName))
        return 1;

    extractGeneIntervals(intervals, store);
    constructIntervalTrees(intervalTrees, intervals);
    countReadsPerGene(readsPerGene, intervalTrees, store);
    outputGeneCoverage(readsPerGene, store);

    return 0;
}
//![main]
