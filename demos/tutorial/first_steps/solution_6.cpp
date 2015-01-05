// Provide a print function that prints pairs of positions and their score if the score is greater than 0.

#include <iostream>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

using namespace seqan;

template <typename TText>
int computeLocalScore(TText const & subText, String<AminoAcid> const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < length(pattern); ++i)
        localScore += score(Blosum62(), subText[i], pattern[i]);

    return localScore;
}

template <typename TText, typename TPattern>
int computeLocalScore(TText const & subText, TPattern const & pattern)
{
    int localScore = 0;
    for (unsigned i = 0; i < length(pattern); ++i)
        if (subText[i] == pattern[i])
            ++localScore;

    return localScore;
}

template <typename TText, typename TPattern>
String<int> computeScore(TText const & text, TPattern const & pattern)
{
    String<int> score;
    resize(score, length(text) - length(pattern) + 1, 0);

    for (unsigned i = 0; i < length(text) - length(pattern) + 1; ++i)
        score[i] = computeLocalScore(infix(text, i, i + length(pattern)), pattern);

    return score;
}

template <typename TText>
void print(TText const & text)
{
    std::cout << text << std::endl;
}

void print(String<int> const & text)
{
    for (unsigned i = 0; i < length(text); ++i)
        std::cout << text[i] << " ";
    std::cout << std::endl;
}

template <typename TText, typename TSpec>
void print(TText const & text, TSpec const & /*tag*/)
{
    print(text);
}

struct MaxOnly {};

template <typename TText>
void print(TText const & score, MaxOnly const & /*tag*/)
{
    int maxScore = score[0];
    String<int> output;
    appendValue(output, 0);
    for (unsigned i = 1; i < length(score); ++i)
    {
        if (score[i] > maxScore)
        {
            maxScore = score[i];
            clear(output);
            resize(output, 1, i);
        }
        else if (score[i] == maxScore)
            appendValue(output, i);
    }

    print(output);
}

struct GreaterZero {};

template <typename TText>
void print(TText const & score, GreaterZero const & /*tag*/)
{
    String<Pair<int> > output;
    for (unsigned i = 1; i < length(score); ++i)
        if (score[i] > 0)
            appendValue(output, Pair<int>(i, score[i]));

    for (unsigned i = 0; i < length(output); ++i)
        std::cout << "(" << output[i].i1 << "; " << output[i].i2 << ") ";
    std::cout << std::endl;
}

int main()
{
    String<char> text = "This is an awesome tutorial to get to now SeqAn!";
    String<char> pattern = "tutorial";
    String<int> score = computeScore(text, pattern);

    print(text);
    // > This is an awesome tutorial to get to now SeqAn!
    print(pattern);
    // > tutorial
    print(score);
    // > 1 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 8 0 1 0 0 0 0 2 0 1 0 0 1 0 3 0 1 1 0 0 0 0
    print(score, MaxOnly());
    // > 19
    print(score, GreaterZero());
    // > (2; 1) (5; 1) (12; 1) (17; 1) (19; 8) (21; 1) (26; 2) (28; 1) (31; 1) (33; 3) (35; 1) (36; 1)

    // And now for a protein pattern
    String<AminoAcid> protein = "tutorial";
    String<int> proteinScore = computeScore(text, protein);

    print(text);
    // > This is an awesome tutorial to get to now SeqAn!
    print(protein);
    // > TXTXRIAL
    print(proteinScore);
    // > 6 -9 -3 -6 -6 0 -9 -8 -7 -3 -9 -5 -8 -4 -5 -6 -6 1 -6 25 -7 2 -6 -6 -9 -6 -5 -7 1 -7 -5 -4 -6 2 -6 -3 -8 -9 -10 -4 -6 0 0 0 0 0 0 0
    print(proteinScore, MaxOnly());
    // > 19
    print(proteinScore, GreaterZero());
    // > (17; 1) (19; 25) (21; 2) (28; 1) (33; 2)

    return 0;
}
