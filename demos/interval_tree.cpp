///An example for using interval trees.
#include <iostream>
#include <seqan/graph_align.h>

using namespace seqan;

int main()
{

    typedef CharString TCargo;  // id type
    typedef int TValue;         // position type

    typedef IntervalAndCargo<TValue, TCargo> TInterval;
    typedef IntervalTree<TValue, TCargo> TIntervalTree;

    String<TInterval> intervals;
    resize(intervals, 5);

///Store gene annotation in intervals.
    intervals[0].i1 = 5;   intervals[0].i2 = 1000;
    intervals[0].cargo = "gene";

    intervals[1].i1 = 50;  intervals[1].i2 = 200;
    intervals[1].cargo = "exon";

    intervals[2].i1 = 600; intervals[2].i2 = 800;
    intervals[2].cargo = "exon";

    intervals[3].i1 = 100; intervals[3].i2 = 200;
    intervals[3].cargo = "coding";

    intervals[4].i1 = 600; intervals[4].i2 = 700;
    intervals[4].cargo = "coding";

    TIntervalTree tree(intervals);

///Add another interval.
    TInterval interval;
    interval.i1 = 200; interval.i2 = 600;
    interval.cargo = "intron";

    addInterval(tree, interval);

///Query a genomic region.
    TValue delBegin = 300;
    TValue delEnd   = 500;
    String<TCargo> results;

    findIntervals(results, tree, delBegin, delEnd);

    std::cout << "Deletion " << delBegin << ".." << delEnd << " overlaps with ";
    for (unsigned i = 0; i < length(results); ++i)
        std::cout << results[i] << ",";
    std::cout << std::endl;

///Query a single position.
    TValue snpPos = 150;
    findIntervals(results, tree, snpPos);

    std::cout << "SNP " << snpPos << " overlaps with ";
    for (unsigned i = 0; i < length(results); ++i)
        std::cout << results[i] << ",";
    std::cout << std::endl;

    CharString iCargo("exon");
    bool res = removeInterval(tree, 50, 200, iCargo);
    if (res)
        std::cout << "Removed exon interval 50..200.\n";

///Now, redo the query. This time one interval less should be returned.
    String<TCargo> results2;
    findIntervals(results2, tree, snpPos);

    std::cout << "SNP " << snpPos << " overlaps with ";
    for (unsigned i = 0; i < length(results2); ++i)
        std::cout << results2[i] << ",";
    std::cout << std::endl;

    return 0;
}
