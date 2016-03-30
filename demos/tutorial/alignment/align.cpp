//![main]
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
//![main]
//![typedefs]
    typedef char TChar;                             // character type
    typedef String<TChar> TSequence;                // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;     // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
//![typedefs]

std::cout << "//![output_init]\n";
//![init]
    TSequence seq1 = "CDFGDC";
    TSequence seq2 = "CDEFGAHGC";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
    std::cout << align;
//![init]
std::cout << "//![output_init]\n";

std::cout << "//![output_manipulation]\n";
//![manipulation]
    TRow & row1 = row(align, 0);
    TRow & row2 = row(align, 1);
    insertGap(row1, 2);
    std::cout << align;
    insertGaps(row1, 5, 2);
    std::cout << align;
//![manipulation]
std::cout << "//![output_manipulation]\n";

std::cout << "//![output_source_positions]\n";
//![printingSourcePos]
    std::cout << std::endl << "ViewToSource1: " << std::endl;
    for (auto c: row1)
        std::cout << c << " ";
    std::cout << std::endl;

    for (unsigned i = 0; i < length(row1); ++i)
        std::cout << toSourcePosition(row1, i) << " ";
    std::cout << std::endl;


    std::cout << std::endl << "ViewToSource2: " << std::endl;
    for (auto c: row2)
        std::cout << c << " ";
    std::cout << std::endl;

    for (unsigned i = 0; i < length(row2); ++i)
        std::cout << toSourcePosition(row2, i) << " ";
    std::cout << std::endl;
//![printingSourcePos]
std::cout << "//![output_source_positions]\n";

std::cout << "//![output_view_positions]\n";
//![printingViewPos]
    std::cout << std::endl << "SourceToView1: " << std::endl;
    for (auto c: source(row1))
        std::cout << c << " ";
    std::cout << std::endl;

    for (unsigned i = 0; i < length(source(row1)); ++i)
        std::cout << toViewPosition(row1, i) << " ";
    std::cout << std::endl;


    std::cout << std::endl << "SourceToView2: " << std::endl;
    for (auto c: source(row2))
        std::cout << c << " ";
    std::cout << std::endl;

    for (unsigned i = 0; i < length(source(row2)); ++i)
        std::cout << toViewPosition(row2, i) << " ";
    std::cout << std::endl;
//![printingViewPos]
std::cout << "//![output_view_positions]\n";

std::cout << "//![output_clipping]\n";
//![clipping]
    std::cout << std::endl << "Before clipping:\n" << align;
    setClippedBeginPosition(row1, 1);
    setClippedEndPosition(row1, 7);
    setClippedBeginPosition(row2, 1);
    setClippedEndPosition(row2, 7);
    std::cout << std::endl << "After clipping:\n" << align;
//![clipping]
std::cout << "//![output_clipping]\n";

std::cout << "//![output_clipping_positions]\n";
//![clipping_positions]
    std::cout << std::endl << "ViewToSource1: ";
    for (unsigned i = 0; i < length(row1); ++i)
        std::cout << toSourcePosition(row1, i) << " ";

    std::cout << std::endl << "ViewToSource2: ";
    for (unsigned i = 0; i < length(row2); ++i)
        std::cout << toSourcePosition(row2, i) << " ";
    std::cout << std::endl;

    std::cout << std::endl << "SourceToView1: ";
    for (unsigned i = 0; i < length(source(row1)); ++i)
        std::cout << toViewPosition(row1, i) << " ";

    std::cout << std::endl << "SourceToView2: ";
    for (unsigned i = 0; i < length(source(row2)); ++i)
        std::cout << toViewPosition(row2, i) << " ";
    std::cout << std::endl;
//![clipping_positions]
std::cout << "//![output_clipping_positions]\n";

std::cout << "//![output_iteratingRowClipped]\n";
//![iteratingRowClipped]
    typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row1),
                 itEnd = end(row1);
    for (; it != itEnd; ++it)
    {
        TChar c = isGap(it) ? gapValue<TChar>() : *it;
        std::cout << c << " ";
    }
    std::cout << std::endl;
//![iteratingRowClipped]
std::cout << "//![output_iteratingRowClipped]\n";

std::cout << "//![output_iteratingRowClipped2]\n";
//![iteratingRowClipped2]
    clearClipping(row1);

    it = begin(row1);
    itEnd = end(row1);
    for (; it != itEnd; ++it)
    {
        TChar c = isGap(it) ? gapValue<TChar>() : *it;
        std::cout << c << " ";
    }
    std::cout << std::endl;
//![iteratingRowClipped2]
std::cout << "//![output_iteratingRowClipped2]\n";
//![return]
    return 0;
}
//![return]
