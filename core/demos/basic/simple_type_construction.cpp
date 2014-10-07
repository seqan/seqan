#include <iostream>
#include <seqan/basic.h>

using namespace seqan;

int main()
{
//![simple type construction and assignment]
    seqan::Dna a('C');  // => a.value == 1
    seqan::Dna b(3);    // => b.value == 3
    seqan::Dna c('T');  // => c.value == 3
    seqan::Dna d('t');  // => d.value == 3

    seqan::Dna e;  // => e.value == 0
    e = 'N';       // => e.value == 0
    e = 'c';       // => e.value == c
//![simple type construction and assignment]
    
    return 0;
}
