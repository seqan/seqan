// FRAGMENT(all)
#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

int computeLocalScore(seqan::String<char> subText, seqan::String<char> pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}

seqan::String<int> computeScore(seqan::String<char> text, seqan::String<char> pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text) - seqan::length(pattern) + 1, 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    
    return score;
}

int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to know SeqAn!";
    seqan::String<char> pattern = "tutorial";
    seqan::String<int> score = computeScore(text, pattern);

    for (unsigned i = 0; i < seqan::length(score); ++i)
        std::cout << score[i] << " ";
    std::cout << std::endl;

    return 0;
}
