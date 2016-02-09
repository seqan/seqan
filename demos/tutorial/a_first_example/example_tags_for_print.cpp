#include <iostream>
#include <seqan/sequence.h>
#include <seqan/score.h>

template <typename TText, typename TSpec>
void print(TText const & text, TSpec const & /*tag*/)
{
    for (unsigned i = 0; i < seqan::length(text); ++i)
        std::cout << text[i] << " ";
    std::cout << std::endl;
}

struct MaxOnly {};

template <typename TText>
void print(TText const & score, MaxOnly const & /*tag*/)
{
    int maxScore = score[0];
    seqan::String<int> output;
    appendValue(output, 0);
    for (unsigned i = 1; i < seqan::length(score); ++i)
    {
        if (score[i] > maxScore)
        {
            maxScore = score[i];
            clear(output);
            resize(output, 1, i);
        }
        else if (score[i] == maxScore)
        {
            appendValue(output, i);
        }
    }

    for (unsigned i = 0; i < seqan::length(output); ++i)
        std::cout << output[i] << " ";
    std::cout << std::endl;
}

int main()
{
    return 0;
}
