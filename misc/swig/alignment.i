%define DOCSTRING
"Demo"
%enddef

%module(docstring=DOCSTRING) alignment


%{
#define SWIG_FILE_WITH_INIT
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/alignment.h>
#include "alignment.h"
%}
%include "alignment.h"

%feature("autodoc", """alignDna(\"seq1\",\"seq2\",[scores],\"alignment algorithms\",binary value (0-15),lowDiag=0,highDiag=0)->(AlignmentObject,AlignmentScore)\n Aligns two dna sequences, returns an AlignmentObject and the corresponding alignment score.\n The Parameters are:\n seq1= \t\t\t\t string sequence \n seq1= \t\t\t\t string sequence \n [scores] = \t\t\t an integer list containing either 4 values describing scores for a match, mismatch, opening gaps and extending gaps or a list with size of the alphabet^2 which describes the score for each substitution. \n alignment algorithms= \t\t MyersHirschberg \n \t\t\t\t\t Hirschberg \n \t\t\t\t\t NeedlemanWunsch \n \t\t\t\t\t Gotoh \n \t\t\t\t\t SmithWaterman \n \t\t\t\t\t BandedNeedlemanWunsch \n \t\t\t\t\t BandedGotoh \n binary value= \t\t\t\t Is a binary representation of an integer from an intervall 0-15 which indicates how the DP matrix is initialised and what ends gaps are free. The original SeqAn object (AlignConfig) use 4 bool values to change this, bool TTop, bool TLeft, bool TRight and bool TBottom. If TTop is true the first row of the DP Matrix is initialised with 0's. If TLeft is true the first column is initialized with 0's. If TRight is true, the maximum is search in the last column. If TBottom is true, the maximum is searched in the last row. All options can be combined in all possible ways, which results in 16 different combinations. Each combination is mapped uniquely to the interval of 0-15 (false,false,false,false)->0, (false,false,false,true)->1, (false,false,true,false)->2 and so on. This feature is not yet supported for all alignment algorithms.\n Example: AlignmentObject,AlignmentScore=alignment.alignDna(\"ACGTC\",\"AACGTCCC\",[3,1,-1,0],\"NeedlemanWunsch\",bin(1)) """) seqan::alignSequence<seqan::Dna5>;

%template(alignDna)seqan::alignSequence<seqan::Dna5>;
%template(alignChar)seqan::alignSequence<char>;
%template(alignAmino)seqan::alignSequence<seqan::AminoAcid>;
%template(printDnaAlignment)seqan::printAlignment<seqan::Dna5>;
%template(printCharAlignment)seqan::printAlignment<char>;
%template(printAminoAcidAlignment)seqan::printAlignment<seqan::AminoAcid>;

%newobject seqan::alignSequence();
