///A tutorial about global alignments.
//#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/align.h>

int main()
{
	using namespace seqan;
	typedef Value<Gaps<Dna5String, ArrayGaps> >::Type TValue;
	/*
	using namespace seqan;

///Two DNA sequences that shall be aligned.
    typedef String<Dna> TSequence;
    TSequence seq1 = "atcgaatgcgga";
    TSequence seq2 = "actcgttgca";
///Scoring objects are used to define a scoring scheme.
///In this case, affine gap costs with match = 0, mismatch = -1, gapextend = -1 and gapopen = -2.
    Score<int> scoringScheme(0, -1, -1, -2);
///Example 1: We use @Class.Align@ to align the two sequences.   
///Since we do not specify an @Tag.Pairwise Global Alignment Algorithms|algorithm tag@ when we call @Function.globalAlignment@, 
///a suitable algorithm (@Tag.Pairwise Global Alignment Algorithms|Gotoh@) is automatically choosen.
    Align<TSequence, ArrayGaps> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    int score = globalAlignment(align, scoringScheme);
    std::cout << "Score = " << score << std::endl;
    std::cout << align << std::endl;
///Example 2: We now choose explicitely the algorithm @Tag.Pairwise Global Alignment Algorithms|MyersHirschberg@.
///Since this algorithm always works on Levenshtein distance, we do not need score.
    score = globalAlignment(align, MyersHirschberg());
    std::cout << "Score = " << score << std::endl;
    std::cout << align << std::endl;
///Example 3: We now do the same as in case 1, but now we use an @Spec.Alignment Graph@ for storing the alignment.
///Here we use @Tag.Pairwise Global Alignment Algorithms|Gotoh's algorithm@.
    typedef StringSet<TSequence, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    TStringSet string_set;
    appendValue(string_set, seq1);
    appendValue(string_set, seq2);
    TAlignmentGraph alignment_graph(string_set);

    score = globalAlignment(alignment_graph, scoringScheme, Gotoh());
    std::cout << "Score = " << score << std::endl;
    std::cout << alignment_graph << std::endl;
	*/
    return 0;
}
