//![includes]
#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;
//![includes]

//![main]
int main()
{
//![main]
//![sequences]
    // Initialization
    String<char> text = "This is an awesome tutorial to get to know SeqAn!";
    String<char> pattern = "tutorial";
//![sequences]

//![score]
    String<int> score;
//![score]
//![resize]
    resize(score, length(text) - length(pattern) + 1);
//![resize]

//![similarity]
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
//![similarity]

//![print]
    // Printing the result
    for (unsigned i = 0; i < length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;
//![print]

//![end]
    return 0;
}
//![end]
