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

#include <cstdint>

namespace seqan {

// ============================================================================
// Typedefs
// ============================================================================

/*!
 * @typedef RnaAdjacencyIterator
 * @headerfile <seqan/rna_io.h>
 * @brief Iterator for adjacent vertices in a @link RnaStructureGraph @endlink.
 * @signature typedef typename Iterator<Graph<Undirected<double> >, AdjacencyIterator>::Type RnaAdjacencyIterator;
 */
typedef typename Iterator<Graph<Undirected<double> >, AdjacencyIterator>::Type RnaAdjacencyIterator;

/*!
 * @typedef TSizeRna5String
 * @headerfile <seqan/rna_io.h>
 * @brief Type for the length of @link Rna5String @endlink.
 * @signature typedef typename Size<Rna5String>::Type TSizeRna5String;
 */
typedef typename Size<Rna5String>::Type TSizeRna5String;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RnaStructureGraph
// ----------------------------------------------------------------------------

/*!
 * @class RnaStructureGraph
 * @headerfile <seqan/rna_io.h>
 * @brief An undirected graph representing an RNA structure.
 *
 * @signature class RnaStructureGraph;
 *
 * Storage of a graph that contains base pairs and their probabilities of an RNA structure.
 * The associated <i>specs</i> string describes the method that was used to obtain the RNA structure.
 */
class RnaStructureGraph {
public:
    /*!
     * @var Graph<Undirected<double>> RnaStructureGraph::inter
     * @brief Undirected graph for base pairings of the RNA structure.
     *
     * The vertices denote the sequence/alignment column index.
     * Edges are drawn among all pairing nucleotides.
     * An edge's cargo represents the probability for the respective pairing.
     */
    Graph<Undirected<double> > inter;

    /*!
     * @var CharString RnaStructureGraph::specs
     * @brief Specs of the Method used to compute the bpp matrix or the structure.
     */
    CharString specs{};

    /*!
     * @var float RnaStructureGraph::energy
     * @brief Energy of the RNA structure.
     */
    float energy{};
};

// ----------------------------------------------------------------------------
// Class RnaRecord
// ----------------------------------------------------------------------------

/*!
 * @class RnaRecord
 * @headerfile <seqan/rna_io.h>
 * @brief A container for RNA structure data.
 *
 * @signature class RnaRecord;
 *
 * The container stores all kinds of data that can be obtained by reading RNA structure file records.
 */
class RnaRecord
{
private:
    // Constant for an undefined ID.
    static std::uint32_t const undef = UINT32_MAX;

public:
    /*!
     * @var std::uint32_t RnaRecord::recordID
     * @brief Identification of the record.
     *
     * In an RNA structure file the first record gets ID 0, the following ID 1 and so on.
     */
    std::uint32_t recordID{undef};

    /*!
     * @var TSizeRna5String RnaRecord::seqLen
     * @brief Length of the sequence or alignment stored in this record.
     */
    TSizeRna5String seqLen{};

    /*!
     * @var unsigned RnaRecord::offset
     * @brief Start index of the sequence.
     */
    unsigned offset{1u};

    /*!
     * @var CharString RnaRecord::name
     * @brief Sequence name.
     */
    CharString name{};

    /*!
     * @var CharString RnaRecord::quality
     * @brief Quality values for the sequence.
     */
    CharString quality{};

    /*!
     * @var String<RnaStructureGraph> RnaRecord::bppMatrGraphs
     * @brief Vector of base pair probability graphs extracted from the input files.
     */
    String<RnaStructureGraph> bppMatrGraphs;

    /*!
     * @var String<RnaStructureGraph> RnaRecord::fixedGraphs
     * @brief Vector of fixed structure graphs extracted from the input files.
     */
    String<RnaStructureGraph> fixedGraphs;

    /*!
     * @var CharString RnaRecord::comment
     * @brief Comment to be stored together with the record.
     */
    CharString comment{};

    /*!
     * @var StringSet<String<float>> RnaRecord::reactivity
     * @brief The area peak/likelihood estimate that represents the <b>reactivity</b> at each position.
     *
     * This member variable is used only if biological validated data is found (T.. fields in EBPSEQ are set).
     */
    StringSet<String<float> > reactivity;

    /*!
     * @var StringSet<String<float>> RnaRecord::reactError
     * @brief Error of the @link RnaRecord::reactivity @endlink.
     *
     * This member variable is used only if biological validated data is found (T.. fields in EBPSEQ are set).
     *
     * If @link RnaRecord::reactivity @endlink was derived as a consensus of different replicates,
     * this indicates the standard error between samples used.
     * If only one experiment was done, this may be some measure of variation between the data in that experiment,
     * e.g. 0.2 of the standard deviation.
     */
    StringSet<String<float> > reactError;

    /*!
     * @var String<std::uint32_t> RnaRecord::typeID
     * @brief Indices of the assigned type (T..) attributes records from EBPSEQ files.
     *
     * This member variable is used only if biological validated data is found (T.. fields in EBPSEQ are set).
     */
    String<std::uint32_t> typeID;

    /*!
     * @var Rna5String RnaRecord::sequence
     * @brief String of bases in RNA strand.
     *
     * This member variable is only used in sequence-based records (from CT, DBN, DBV, BPSEQ, EBPSEQ files).
     */
    Rna5String sequence{};

    /*!
     * @var StringSet<CharString> RnaRecord::seqID
     * @brief Sequence identifier for aligned sequences.
     *
     * This member variable is only used in alignment-based records (from STH files).
     */
    StringSet<CharString> seqID;

    /*!
     * @var Align<Rna5String,ArrayGaps> RnaRecord::align
     * @brief Alignment of several sequences (including gaps).
     *
     * This member variable is only used in alignment-based records (from STH files).
     */
    Align<Rna5String, ArrayGaps> align;

    /*!
     * @fn RnaRecord::hasUndefinedID()
     * @brief Test for an undefined @link RnaRecord::recordID @endlink value.
     * @return bool True if @link RnaRecord::recordID @endlink is not set.
     * @signature bool RnaRecord::hasUndefinedID()
     */
    bool hasUndefinedID() const
    {
        return recordID == undef;
    }

    /*!
     * @fn RnaRecord::clearID()
     * @brief Clear value of @link RnaRecord::recordID @endlink and set to undefined.
     * @signature void RnaRecord::clearID()
     */
    void clearID()
    {
        recordID = undef;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

inline void clear(RnaRecord & record)
{
    record.clearID();
    record.seqLen = 0;
    record.offset = 1;
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

#endif  // SEQAN_INCLUDE_SEQAN_RNA_IO_RNA_RECORD_H_
