//![swap-headers]
#include <iostream>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace seqan;
//![swap-headers]

//![swap-declaration]
template <typename T>
void swap(T & container, int i, int j, int k)
{
//![swap-declaration]

//![swap-metafunction]
    // define helper variable
    T help;
    resize(help, k);

    for (int x = 0; x < k; ++x)
        value(help, x) = container[i + x];
//![swap-metafunction]

//![swap-work]
    for (int x = 0; x < k; ++x)
        value(container, i + x) = value(container, j + x);
    for (int x = 0; x < k; ++x)
        value(container, j + x) = help[x];

    return;
}
//![swap-work]

//![swap-main]
int main()
{
    typedef String<Dna> TDnaString;
    TDnaString dna = "AAAATTTT";

    typedef String<int> TIntString;
    typedef Iterator<String<int>, Rooted>::Type TIntIterator;

    TIntString numbers;
    appendValue(numbers, 1);   appendValue(numbers, 1);   appendValue(numbers, 1);
    appendValue(numbers, 1);   appendValue(numbers, 1);   appendValue(numbers, 1);
    appendValue(numbers, 3);   appendValue(numbers, 3);   appendValue(numbers, 3);
    appendValue(numbers, 3);   appendValue(numbers, 3);   appendValue(numbers, 3);
//![swap-main]

//![swap-apply]
    swap(dna, 1, 4, 2);
    std::cout << dna << "\n";

    swap(numbers, 1, 7, 2);
    for (TIntIterator it = begin(numbers); !atEnd(it); goNext(it))
        std::cout << *it;
    std::cout << "\n";

    return 0;
}
//![swap-apply]
