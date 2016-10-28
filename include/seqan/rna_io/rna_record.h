// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Authors: Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <winkler@molgen.mpg.de>
// ==========================================================================


#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_RECORD_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_RECORD_H_

namespace seqan {

// ============================================================================
// Typedefs
// ============================================================================

typedef Graph<Undirected<double> > TRnaRecordGraph;

typedef typename Iterator<TRnaRecordGraph, AdjacencyIterator>::Type TAdjacencyIterator;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RnaRecord
// ----------------------------------------------------------------------------

class RnaRecord
{
public:
    static const int INVALID_POS = -1;

    // Amount of records.
    unsigned amount;

    //beginning and ending positions of the sequence
    int begPos;
    int endPos;

    // Energy
    float energy;

    // Record's name.
    CharString name;

    // Sequence identifier
    StringSet<CharString> seq_id;
    
    //string of base at each position in Rna strand
    StringSet<Rna5String, Owner<JournaledSet> > sequence;

    // Undirected graph for connected bases and probabilities
    TRnaRecordGraph graph;

    bool isInjective; // each base has at most 1 connection

    ////////RDAT FILES
    CharString qual; //I think?

    int offset;

    String<CharString> seqpos;

    String<CharString> annotation;

    CharString comment;

    //Annotation data 1
    //annotation data 2

    String<float> reactivity;

    String<float> reactivity_error;

    String<float> xsel;

    String<float> xsel_refine;

    //mutpos
  
    // Default constructor.
    RnaRecord() : amount(0), begPos(INVALID_POS), endPos(INVALID_POS), energy(0), name(" "), isInjective(true),
        offset(0), comment("")
    {}                                                                                      

};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------


inline void clear(RnaRecord & record)
{
    clear(record.name);
    clear(record.sequence);
    clear(record.graph);
    clear(record.qual);

    clear(record.seqpos);
    clear(record.annotation);
    clear(record.comment);
    clear(record.reactivity);
    clear(record.reactivity_error);
    clear(record.xsel);    
    clear(record.xsel_refine);
}

}  // namespace seqan

#endif  //SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_RECORD_H_
