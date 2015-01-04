/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  A simple tool for pairwise local and alignment given at the command line.
  Currently works for type DnaString only.
  ===========================================================================*/

#include <iostream>

#include <seqan/align.h>

// Banner to print at top of program output.
const char * kBanner = (
    "********************************\n"
    " arg_align  -- Simple Alignment \n"
    "********************************\n"
    "\n"
    "We use edit distance scoring\n");

// Usage help message.
const char * kUsageStr = (
    "Align two strings, given on the command line.\n"
    "\n"
    "Usage: local_align SEQ1 SEQ2");

// Program return codes.
const int kRetOk = 0;      // No error.
const int kRetArgErr = 1;  // Argument error.

using namespace seqan;


typedef Dna5String TSequenceType;  // The sequence type used below.


// Parse parameters.
//
// Return true iff the parameters were valid.
bool parseParameters(int argc, char ** argv, TSequenceType & seq1, TSequenceType & seq2)
{
    if (argc != 3)
        return false;

    // Simply assign parameters.
    // TODO(holtgrew): Checking for valid characters.
    seq1 = argv[1];
    seq2 = argv[2];
    return true;
}

int main(int argc, char ** argv)
{
    std::cout << kBanner << std::endl;

    // Get parameter from args.
    TSequenceType seq1, seq2;
    if (!parseParameters(argc, argv, seq1, seq2))
    {
        std::cerr << kUsageStr << std::endl;
        return kRetArgErr;
    }

    // Perform local alignments.
    //
    // First, setup the book-keeping objects.
    Align<TSequenceType> ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), seq1);
    assignSource(row(ali, 1), seq2);
    Score<int> scoring; //(0, -1, -1);

    // Then, find the best local alignment.
    // TODO(holtgrew): Score always seems to be 0.
    {
        int scoreValue = localAlignment(ali, scoring);
        std::cout << "Local alignment with score " << scoreValue << std::endl;
        std::cout << ali << std::endl;
    }

    // Then, find the best local alignment.
    {
        int scoreValue = globalAlignment(ali, scoring, NeedlemanWunsch());
        std::cout << "Global alignment with score " << scoreValue << std::endl;
        std::cout << ali << std::endl;
    }

    return kRetOk;
}
