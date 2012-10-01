// FRAGMENT(header)
#include <seqan/seeds.h>

using namespace seqan;

// FRAGMENT(writeSeed)
template<typename TSeed, typename TStringSet>
void
writeSeed(TSeed & seed, TStringSet const & seqs) {
    for (unsigned dim = 0; dim < dimension(seed); ++dim) {
        std::cout << "Seed from position " << leftPosition(seed, dim);
        std::cout << " to " << rightPosition(seed, dim);
        std::cout << " in Sequence " << dim << " = " << value(seqs, dim) << ": ";
        std::cout << infix(value(seqs, dim), leftPosition(seed, dim), rightPosition(seed, dim)+1);
        std::cout << std::endl;
    }
}

//FRAGMENT(main)
int main() {
    StringSet<DnaString> seqs;
    appendValue(seqs, "ATCCGCAT");
    appendValue(seqs, "TATGCAGCC");
    appendValue(seqs, "CACTCCGCT");

    typedef Seed<int, MultiSeed> TSeed;
    TSeed seed(3);

    setLeftPosition(seed, 0, 2);
    setRightPosition(seed, 0, 5);
    setLeftPosition(seed, 1, 4);
    setRightPosition(seed, 1, 7);
    setLeftPosition(seed, 2, 4);
    setRightPosition(seed, 2, 7);

    writeSeed(seed, seqs);

    return 0;
}
