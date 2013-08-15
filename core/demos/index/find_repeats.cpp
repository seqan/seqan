#include <fstream>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

int main()
{
    using namespace seqan;

    // Get path to file to search for repeats in.
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/core/demos/index/ref.fa");

    // Load first sequence from file.
    Dna5String seq;
    if (readFasta(seq, toCString(path)) != 0)
    {
        std::cerr << "Could not open " << path << "\n";
        return 1;
    }

    // Find repeats and print them.
    String<Repeat<unsigned, unsigned> > repeats;
    findRepeats(repeats, seq, 3);

    std::cerr << "# of repeats: " << length(repeats) << "\n";
    for (unsigned i = 0; i < length(repeats); ++i)
        std::cerr << "i == " << i << ", beginPosition = " << repeats[i].beginPosition
                  << ", endPosition = " << repeats[i].endPosition
                  << ", period = " << repeats[i].period << "\n";

    return 0;
}
