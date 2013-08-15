#include <iostream>
#include <fstream>

#include <seqan/index.h>
#include <seqan/statistics.h>
#include <seqan/seq_io.h>

using namespace std;
using namespace seqan;

int main()
{
    // Build path to background FASTA file.
    CharString bgPath = SEQAN_PATH_TO_ROOT();
    append(bgPath, "/extras/demos/statistics/background.fa");

    // Read the background from a file into X.
    StringSet<DnaString> X;
    SequenceStream seqStream(toCString(bgPath));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open " << bgPath << "\n";
        return 1;
    }
    StringSet<CharString> ids;  // will be ignored
    if (readAll(ids, X, seqStream) != 0)
    {
        std::cerr << "ERROR: Problem reading from " << bgPath << "\n";
        return 1;
    }

    // Create MarkovModel of order 3 from the background.
    MarkovModel<Dna> mm(3);
    buildMarkovModel(mm, X);

    // Build set of words that we want to compute the zscore of.
    StringSet<DnaString> W;
    appendValue(W, "CCCAAAGC");
    appendValue(W, "CCCAAAGTAAATT");

    // Compute and print zscore.
    std::cout << "zscore: " << zscore(W, X, mm, AhoCorasick()) << "\n";

// //TODO his path has to be set explicitely when calling the demo
//  FILE *fd = fopen("projects/library/demos/zscore_human_mm.3","r");
//  read(fd, mm);
//  fclose(fd);

    //std::cout << zscore(W, X, mm, WuManber()) << std::endl;

    return 0;
}
