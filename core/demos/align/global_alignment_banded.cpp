#include <iostream>

#include <seqan/base.h>
#include <seqan/align.h>
#include <seqan/file.h>  // for printint strings

int main()
{
    seqan::Dna5String seqH = "CGATT";
    seqan::Dna5String seqV = "CGAAATT";

    seqan::Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seqH);
    assignSource(row(align, 0), seqV);

    seqan::Score<int, seqan::Simple> scoringScheme(2, -1, -2);
    seqan::AlignConfig<> alignConfig;

    int lDiag = -2;
    int uDiag = -2;

    int result = globalAlignment(align, scoringScheme, alignConfig, lDiag, uDiag);

    std::cerr << "The resulting alignment is\n"
              << align;

    return 0;
}
