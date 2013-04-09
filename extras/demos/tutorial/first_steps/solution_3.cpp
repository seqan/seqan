// Adjust your current code to be more memory and time efficient by using references in the function header.

#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

int computeLocalScore(seqan::String<char> const & subText, seqan::String<char> const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}

seqan::String<int> computeScore(seqan::String<char> const & text, seqan::String<char> const & pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text) - seqan::length(pattern) + 1, 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    
    return score;
}

void print(seqan::String<int> const & text)
{
    for (unsigned i = 0; i < seqan::length(text); ++i)
        std::cout << text[i] << " ";
    std::cout << std::endl;
}

int main()
{
    seqan::String<char> text = "This is an awesome tutorial to get to now SeqAn!";
    seqan::String<char> pattern = "tutorial";
    seqan::String<int> score = computeScore(text, pattern);

    print(score);

    return 0;
}
