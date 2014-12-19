//![solution]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    // Defining all types that are needed.
    typedef String<char> TSequence;
    typedef Align<TSequence, ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    typedef Iterator<TRow>::Type TRowIterator;

    TSequence seq1 = "ACGTCACCTC";
    TSequence seq2 = "ACGGGCCTATC";

    // Initializing the align object.
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    // Use references to the rows of align.
    TRow & row1 = row(align, 0);
    TRow & row2 = row(align, 1);

    // Insert gaps.
    insertGaps(row1, 2, 2);
    insertGap(row1, 7);  // We need to pass the view position which is changed due to the previous insertion.
    insertGaps(row2, 9, 2);

    // Initialize the row iterators.
    TRowIterator itRow1 = begin(row1);
    TRowIterator itEndRow1 = end(row1);
    TRowIterator itRow2 = begin(row2);

    // Iterate over both rows simultaneously.
    int gapCount = 0;
    for (; itRow1 != itEndRow1; ++itRow1, ++itRow2)
    {
        if (isGap(itRow1))
        {
            gapCount += countGaps(itRow1);  // Count the number of consecutive gaps from the current position in row1.
            itRow1 += countGaps(itRow1);    // Jump to next position to check for gaps.
            itRow2 += countGaps(itRow1);    // Jump to next position to check for gaps.
        }
        if (isGap(itRow2))
        {
            gapCount += countGaps(itRow2);  // Count the number of consecutive gaps from the current position in row2.
            itRow1 += countGaps(itRow2);    // Jump to next position to check for gaps.
            itRow2 += countGaps(itRow2);    // Jump to next position to check for gaps.
        }
    }
    // Print the result.
    std::cout << "Number of gaps: " << gapCount << std::endl;
}
//![solution]
