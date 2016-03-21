#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/align.h>


using namespace seqan;

int main()
{
    typedef String<Dna> TSequence;                   // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;      // align type

//![lookup_rule]
    seqan::String<char> s = "example";
    unsigned i = length(s);
//![lookup_rule]
    ignoreUnusedVariableWarning(i);

    TAlign align;
//![seqan_tags]
    globalAlignment(align, seqan::MyersHirschberg());
//![seqan_tags]
    return 0;
}
