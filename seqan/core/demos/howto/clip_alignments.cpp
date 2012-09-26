// FRAGMENT(headers)
#include <seqan/file.h>
#include <seqan/align.h>

// TODO(holtgrew): Remove, not required for new Gaps objects.

using namespace seqan;

// FRAGMENT(print-function)
/*
template<typename TRow, typename TString>
void printDataArray(TRow & row, TString & name) {
    typedef typename Size<typename TRow::TArr>::Type TSize;
    typename TRow::TArr dataArray = row.data_arr;

    std::cout << "data array " << name << ": [";
    if (length(dataArray) > 0) std::cout << value(dataArray, 0);
    for(TSize i = 1; i < length(dataArray); ++i) {
        std::cout << "," << value(dataArray, i);
    }
    std::cout << "]" << std::endl;
}
*/

// FRAGMENT(main)
int main(int, const char *[]) {
/*

// FRAGMENT(unclipped)
    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "acgtttacgaat");
    assignSource(row(align, 1), "agtttatcggt");

    globalAlignment(align, Score<int>(1, -1, -1));
    std::cout << align;

    printDataArray(row(align, 0), "row 0");
    printDataArray(row(align, 1), "row 1");
    std::cout << std::endl;

    std::cout << "toViewPosition(row0, 5): " << toViewPosition(row(align, 0), 5) << std::endl;
    std::cout << "toViewPosition(row1, 5): " << toViewPosition(row(align, 1), 5) << std::endl;
    std::cout << std::endl;

    std::cout << "toSourcePosition(row0, 4): " << toSourcePosition(row(align, 0), 4) << std::endl;
    std::cout << "toSourcePosition(row1, 4): " << toSourcePosition(row(align, 1), 4) << std::endl;
    std::cout << std::endl;

// FRAGMENT(clipping)

    setClippedBeginPosition(row(align, 0), toSourcePosition(row(align, 0), 2));
    setClippedBeginPosition(row(align, 1), toSourcePosition(row(align, 1), 2));
    setClippedEndPosition(row(align, 0), toSourcePosition(row(align, 0), 10));
    setClippedEndPosition(row(align, 1), toSourcePosition(row(align, 1), 10));

    std::cout << align;

    printDataArray(row(align, 0), "row 0");
    printDataArray(row(align, 1), "row 1");
    std::cout << std::endl;

    std::cout << "clippedBeginPosition(row0): " << clippedBeginPosition(row(align, 0)) << std::endl;
    std::cout << "clippedBeginPosition(row1): " << clippedBeginPosition(row(align, 1)) << std::endl;
    std::cout << std::endl;

// FRAGMENT(clipped)

    std::cout << "toViewPosition(row0, 5): " << toViewPosition(row(align, 0), 5) << std::endl;
    std::cout << "toViewPosition(row1, 5): " << toViewPosition(row(align, 1), 5) << std::endl;
    std::cout << std::endl;

    std::cout << "toSourcePosition(row0, 4): " << toSourcePosition(row(align, 0), 4) << std::endl;
    std::cout << "toSourcePosition(row1, 4): " << toSourcePosition(row(align, 1), 4) << std::endl;
    std::cout << std::endl;

// FRAGMENT(tasks)

    std::cout << "TASK 1 (clipped view pos cvp of clipped source pos csp in row 0): " << std::endl;
    std::cout << "  csp = 2 -> cvp = ";
    std::cout << toViewPosition(row(align, 0), 2 + clippedBeginPosition(row(align, 0))) - toViewPosition(row(align, 0), clippedBeginPosition(row(align, 0))) << std::endl;
    std::cout << "  csp = 6 -> cvp = ";
    std::cout << toViewPosition(row(align, 0), 6 + clippedBeginPosition(row(align, 0))) - toViewPosition(row(align, 0), clippedBeginPosition(row(align, 0))) << std::endl;
    std::cout << std::endl;

    std::cout << "TASK 2 (source pos sp of clipped view pos cvp in row 0): " << std::endl;
    std::cout << "  cvp = 2 -> sp = ";
    std::cout << toSourcePosition(row(align, 0), 2 + toViewPosition(row(align, 0), clippedBeginPosition(row(align, 0)))) << std::endl;
    std::cout << "  cvp = 6 -> sp = ";
    std::cout << toSourcePosition(row(align, 0), 6 + toViewPosition(row(align, 0), clippedBeginPosition(row(align, 0)))) << std::endl;
    std::cout << std::endl;

    std::cout << "TASK 3 (clipped source pos csp of source pos sp in row 0): " << std::endl;
    std::cout << "  sp = 4 -> csp = ";
    std::cout << 4 - clippedBeginPosition(row(align, 0)) << std::endl;
    std::cout << "  sp = 8 -> csp = ";
    std::cout << 8 - clippedBeginPosition(row(align, 0)) << std::endl;
    std::cout << std::endl;
*/
}

