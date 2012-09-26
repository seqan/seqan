// Provide a generic print function which is used when the input type is not String<int>. 

#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

template <typename TText>
int computeLocalScore(TText const & subText, seqan::String<seqan::AminoAcid> const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        localScore += seqan::score(seqan::Blosum62(), subText[i], pattern[i]);
    
    return localScore;
}

template <typename TText, typename TPattern>
int computeLocalScore(TText const & subText, TPattern const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < seqan::length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;
    
    return localScore;
}

template <typename TText, typename TPattern>
seqan::String<int> computeScore(TText const & text, TPattern const & pattern)
{
    seqan::String<int> score;
    seqan::resize(score, seqan::length(text), 0);

    for (unsigned i = 0; i < seqan::length(text) - seqan::length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + seqan::length(pattern)), pattern);
    
    return score;
}

template <typename TText>
void print(TText const & text)
{
    std::cout << text << std::endl;
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

    print(text);
    // > This is an awesome tutorial to get to now SeqAn!
    print(pattern);
    // > tutorial
    print(score);
    // > 1 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 8 0 1 0 0 0 0 2 0 1 0 0 1 0 3 0 1 1 0 0 0 0 0 0 0 0 0 0 0

    return 0;
}
