// FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
// FRAGMENT(typedefs)
	typedef String<char> TSequence;					// sequence type
	typedef Align<TSequence,ArrayGaps> TAlign;		// align type
	typedef Row<TAlign>::Type TRow;					// gapped sequence type

// FRAGMENT(init)
	TSequence seq1 = "CDFGHC";
	TSequence seq2 = "CDEFGAHC";

	TAlign align;
	resize(rows(align), 2); 
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

// FRAGMENT(manipulation)
	std::cout << align;
	TRow &row1 = row(align,0);
	TRow &row2 = row(align,1);
	insertGap(row1,2);
	insertGap(row1,5);
	std::cout << align;
	
// FRAGMENT(printingViewPos)
	std::cout << std::endl << "ViewToSource1: ";
	for(unsigned i = 0; i < length(row1); ++i)
		std::cout << toSourcePosition(row1, i) << ",";
	
	std::cout << std::endl << "ViewToSource2: ";
	for(unsigned i = 0; i < length(row2); ++i)
		std::cout << toSourcePosition(row2, i) << ",";
	std::cout << std::endl;

// FRAGMENT(printingSourcePos)
	std::cout << std::endl << "SourceToView1: ";
	for(unsigned i = 0; i < length(source(row1)); ++i)
		std::cout << toViewPosition(row1, i) << ",";

	std::cout << std::endl << "SourceToView2: ";
	for(unsigned i = 0; i < length(source(row2)); ++i)
		std::cout << toViewPosition(row2, i) << ",";
	std::cout << std::endl;



	return 0;
}
