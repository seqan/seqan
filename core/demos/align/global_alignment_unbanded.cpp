#include <iostream>

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/file.h>  // for printint strings

using namespace seqan;

int main()
{
    Dna5String seqH = "CGATT";
    Dna5String seqV = "CGAAATT";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seqH);
    assignSource(row(align, 1), seqV);

    Score<int, Simple> scoringScheme(2, -1, -2, -1);
    AlignConfig<> alignConfig;

    int result = globalAlignment(align, scoringScheme, alignConfig);

    std::cout << "Score: " << result << std::endl;
    std::cout << "The resulting alignment is\n"
              << align << std::endl;

    return 0;
}
