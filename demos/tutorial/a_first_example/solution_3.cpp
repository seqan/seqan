//![all]
#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;

//![head_local]
int computeLocalScore(String<char> const & subText, String<char> const & pattern)
//![head_local]
{
    int localScore = 0;
    for (unsigned i = 0; i < length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;

    return localScore;
}

//![head]
String<int> computeScore(String<char> const & text, String<char> const & pattern)
//![head]
{
    String<int> score;
    resize(score, length(text) - length(pattern) + 1, 0);

    for (unsigned i = 0; i < length(text) - length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + length(pattern)), pattern);

    return score;
}

void print(String<int> const & text)
{
    for (unsigned i = 0; i < length(text); ++i)
        std::cout << text[i] << " ";
    std::cout << std::endl;
}

int main()
{
    String<char> text = "This is an awesome tutorial to get to know SeqAn!";
    String<char> pattern = "tutorial";
    String<int> score = computeScore(text, pattern);

    print(score);

    return 0;
}
//![all]
