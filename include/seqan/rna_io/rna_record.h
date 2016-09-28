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
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================


#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_RECORD_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_RECORD_H_

#include <climits>

namespace seqan {

// ============================================================================
// Typedefs
// ============================================================================

typedef Graph<Undirected<double> > TRnaRecordGraph;

typedef typename Iterator<TRnaRecordGraph, AdjacencyIterator>::Type TRnaAdjacencyIterator;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RnaInterGraph
// ----------------------------------------------------------------------------

class RnaInterGraph {
public:
    // Specs of the Method used to compute the bpp matrix or the structure
    CharString specs;
    // graph storing the weight of the RNA structure interactions
    TRnaRecordGraph inter;
    // Default constructor.
    RnaInterGraph() : specs("") {}
    // Constructor with given TRnaRecordGraph
    RnaInterGraph(TRnaRecordGraph const & _inter) : specs(""), inter(_inter) {}
};

// ----------------------------------------------------------------------------
// Class RnaRecord
// ----------------------------------------------------------------------------

class RnaRecord
{
public:
    unsigned const UNDEF = UINT_MAX;

    // identification of record
    unsigned recordID;

    // Amount of records.
    unsigned seqLen;

    // Start position of the sequence
    unsigned offset;

    // Energy
    float energy;

    // Record's name.
    CharString name;

    // string of base at each position in Rna strand, ONLY SINGLE-SEQUENCE RECORDS
    Rna5String sequence;

    // sequence identifier for aligned sequences, ONLY ALIGNMENT RECORDS
    StringSet<CharString> seqID;

    // alignment of several sequences (gaps allowed), ONLY ALIGNMENT RECORDS
    Align<Rna5String, ArrayGaps> align;

    // Undirected graph for base pairings
    // vertices: sequence/alignment column index, edges: base pair with assigned probability
    //String<TRnaRecordGraph> graph;

    // Vector of base pair probability graphs extracted from the input files
    String<RnaInterGraph> bppMatrGraphs;
    // Vector of fixed structure graphs extracted from the input files
    String<RnaInterGraph> fixedGraphs;

    CharString quality;

    // indices for type attribute (see header)
    String<unsigned> typeID;

    ////////RDAT FILES
    //CharString qual; //I think?

    //String<CharString> seqpos;

    //String<CharString> annotation;

    CharString comment;

    //Annotation data 1
    //annotation data 2

    StringSet<String<float> > reactivity;

    StringSet<String<float> > reactError;

    //String<float> xsel;

    //String<float> xsel_refine;

    //mutpos
  
    // Default constructor.
    RnaRecord() : recordID(UNDEF), seqLen(0), offset(1), energy(0.0f), name(""), sequence(""), quality(""), comment("")
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
    record.recordID = 0;
    record.seqLen = 0;
    record.offset = 1;
    record.energy = 0.0f;
    clear(record.name);
    clear(record.sequence);
    clear(record.seqID);
    clearGaps(record.align);
    clear(record.fixedGraphs);
    clear(record.bppMatrGraphs);
    clear(record.quality);
    clear(record.typeID);
    clear(record.comment);
    clear(record.reactivity);
    clear(record.reactError);
}

}  // namespace seqan

#endif  //SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_RECORD_H_
