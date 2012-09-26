#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

#include "curve_smoothing.h"
#include "witio.h"
#include "find_myers_ukkonen_reads.h"

using namespace seqan;

int main(int argc, char **argv) {
    if (argc != 3) return 1;

    Dna5String contig = argv[1];
    Dna5String read = argv[2];

    Finder<Dna5String> finder(contig);
    Pattern<Dna5String, MyersUkkonenReads> pattern(read, -length(read));

    EditDistanceScore scoring;

    while (find(finder, pattern)) {
        if (endPosition(finder) < length(read))
            continue;
        if (endPosition(finder) > length(contig))
            continue;
        while (findBegin(finder, pattern, getScore(pattern)))
            continue;

        std::cout << "end = " << endPosition(finder) << ", begin = " << beginPosition(finder) << ", last = " << endPosition(finder) - 1 << ", score = " << getScore(pattern) << std::endl;
        Align<Segment<Dna5String, InfixSegment> > ali;
        appendValue(rows(ali), infix(finder));
        appendValue(rows(ali), infix(read, 0, length(read)));
        int scoreValue = globalAlignment(ali, scoring, NeedlemanWunsch());
        std::cout << "Global alignment with score " << scoreValue << std::endl;
        std::cout << ali << std::endl;

    }

    return 0;
}
