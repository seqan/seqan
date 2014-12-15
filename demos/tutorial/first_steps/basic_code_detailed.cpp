// FRAGMENT(includes)
#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;

// FRAGMENT(main)
int main()
{
    // FRAGMENT(sequences)
    // Initialization
    String<char> text = "This is an awesome tutorial to get to know SeqAn!";
    String<char> pattern = "tutorial";

    // FRAGMENT(score)
    String<int> score;
    // FRAGMENT(resize)
    resize(score, length(text) - length(pattern) + 1);

    // FRAGMENT(similarity)
    // Computation of the similarities
    // Iteration over the text (outer loop)
    for (unsigned i = 0; i < length(text) - length(pattern) + 1; ++i)
    {
        int localScore = 0;
        // Iteration over the pattern for character comparison
        for (unsigned j = 0; j < length(pattern); ++j)
        {
            if (text[i + j] == pattern[j])
                ++localScore;
        }
        score[i] = localScore;
    }

    // FRAGMENT(print)
    // Printing the result
    for (unsigned i = 0; i < length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;

    // FRAGMENT(end)
    return 0;
}
