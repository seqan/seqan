#include <seqan/align.h>
#include <seqan/seeds.h>

using namespace seqan;

int main() {
	typedef String<Dna>					TSequence;
	typedef Infix<TSequence>::Type		TInfix;
	typedef Position<TSequence>::Type	TPos;
	typedef int							TScoreValue;
	
	//// read sequences
	//TSequence seq_0 = String<Dna, FileReader<Fasta> >("../../../demos/sequence_1a.fa");
	//TSequence seq_1 = String<Dna, FileReader<Fasta> >("../../../demos/sequence_2a.fa");

	//// set input infixes
	//TInfix infix_0 = infix(seq_0, 681, 878);
	//TInfix infix_1 = infix(seq_1, 949, 1166);
	
	// read sequences
	TSequence seq_0 = String<Dna, FileReader<Fasta> >("../../demos/sequence_1.fa");
	TSequence seq_1 = String<Dna, FileReader<Fasta> >("../../demos/sequence_2.fa");

	// set input infixes
	TInfix infix_0 = infix(seq_0, 1033, 1147);
	TInfix infix_1 = infix(seq_1, 1029, 1163);

	// set the scoring scheme
	TScoreValue matchScore = 1;
	TScoreValue penalty = -19;
	Score<TScoreValue> scoreMatrix(matchScore, penalty, penalty);

	// align object for the complete sequences
	Align<TSequence> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), seq_0);
	assignSource(row(align, 1), seq_1);

	// ---------- banded local alignment on infixes ----------

	// align object for the infixes
	Align<TInfix> localAlign;
    resize(rows(localAlign), 2);
    assignSource(row(localAlign, 0), infix_0);
    assignSource(row(localAlign, 1), infix_1);

	// banded local alignment and integration of gaps in large alignment
	int highDiag = 0;
	int lowDiag = length(infix_0) - length(infix_1);
	localAlignment(localAlign, scoreMatrix, lowDiag, highDiag);
	integrateAlign(align, localAlign);

	// ---------- seed extension with local alignment as seed ----------

	// begin and end position of local alignment in sequences
	TPos locAliBegin_0 = clippedBeginPosition(row(localAlign, 0)) + beginPosition(infix_0);
	TPos locAliBegin_1 = clippedBeginPosition(row(localAlign, 1)) + beginPosition(infix_1);
	TPos locAliEnd_0 = clippedEndPosition(row(localAlign, 0)) + beginPosition(infix_0);
	TPos locAliEnd_1 = clippedEndPosition(row(localAlign, 1)) + beginPosition(infix_1);

	// seed extension (gapped X-drop)
	Seed<int, SimpleSeed> seed(locAliBegin_0, locAliBegin_1, locAliEnd_0 - 1, locAliEnd_1 - 1);
	extendSeed(seed, 5*(-penalty), scoreMatrix, seq_0, seq_1, 2, GappedXDrop());

	// alignment on left extension
	if (leftPosition(seed, 0) < (int)locAliBegin_0 || leftPosition(seed, 1) < (int)locAliBegin_1) {
		StringSet<TInfix> leftExtensions;
		appendValue(leftExtensions, infix(seq_0, leftPosition(seed, 0), locAliBegin_0));
		appendValue(leftExtensions, infix(seq_1, leftPosition(seed, 1), locAliBegin_1));

		Align<TInfix> leftAlign;
		resize(rows(leftAlign), 2);
		assignSource(row(leftAlign, 0), leftExtensions[0]);
		assignSource(row(leftAlign, 1), leftExtensions[1]);

		lowDiag = startDiagonal(seed) - leftDiagonal(seed);
		highDiag = startDiagonal(seed) - rightDiagonal(seed);

		globalAlignment(leftAlign, scoreMatrix, lowDiag, highDiag, NeedlemanWunsch());
		integrateAlign(align, leftAlign);
	}

	// alignment on rigth extension
	if (rightPosition(seed, 0) + 1 > (int)locAliEnd_0 || rightPosition(seed, 1) + 1 > (int)locAliEnd_1) {
		StringSet<TInfix> rightExtensions;
		appendValue(rightExtensions, infix(seq_0, locAliEnd_0, rightPosition(seed, 0) + 1));
		appendValue(rightExtensions, infix(seq_1, locAliEnd_1, rightPosition(seed, 1) + 1));

		Align<TInfix> rightAlign;
		resize(rows(rightAlign), 2);
		assignSource(row(rightAlign, 0), rightExtensions[0]);
		assignSource(row(rightAlign, 1), rightExtensions[1]);

		int startDiag = locAliEnd_1 - locAliEnd_0;
		lowDiag = startDiag - leftDiagonal(seed);
		highDiag = startDiag - rightDiagonal(seed);

		globalAlignment(rightAlign, scoreMatrix, lowDiag, highDiag, NeedlemanWunsch());
		integrateAlign(align, rightAlign);
	}

	// set extended begin and end positions of align
	setClippedBeginPosition(row(align, 0), leftPosition(seed, 0));
	setClippedBeginPosition(row(align, 1), leftPosition(seed, 1));
	setBeginPosition(row(align, 0), 0);
	setBeginPosition(row(align, 1), 0);
	setClippedEndPosition(row(align, 0), rightPosition(seed, 0) + 1);
	setClippedEndPosition(row(align, 1), rightPosition(seed, 1) + 1);

	std::cout << align;
	return 0;
}
