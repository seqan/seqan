// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Module for handling NCBI Blast I/O and E-Value computation
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_BLAST_BASE_H_
#define SEQAN_EXTRAS_BLAST_BLAST_BASE_H_


#include <cstdlib>
#include <seqan/version.h>
#include <seqan/align.h>

namespace seqan {


// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @enum BlastFormatFile
 * @brief Enum with BLAST file format specs
 * @signature enum class BlastFormatFile : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * SeqAn can currently write Pairwise, Tabular and TabularWithHeader. It can
 * read Tabular and TabularWithHeader.
 *
 * @val BlastFormatFile BlastFormatFile::Pairwise;
 * @brief default blast output (blastall -m 0 &amp; blast* -outfmt 0)
 *
 * @val BlastFormatFile BlastFormatFile::MasterSlaveIdent;
 * @brief master-slave showing identities (blastall -m 1 &amp; blast* -outfmt 1)
 *
 * @val BlastFormatFile BlastFormatFile::MasterSlaveNoIdent;
 * @brief master-slave without identities (blastall -m 2 &amp; blast* -outfmt 2)
 *
 * @val BlastFormatFile BlastFormatFile::FlatMasterSlaveIdent;
 * @brief flat master-slave showing identities (blastall -m 3 &amp; blast* -outfmt 3)
 *
 * @val BlastFormatFile BlastFormatFile::FlatMasterSlaveNoIdent;
 * @brief master-slave without identities, with blunt ends (blastall -m 5, not available with Blast+)
 *
 * @val BlastFormatFile BlastFormatFile::MasterSlaveBluntEnds;
 * @brief flat master-slave without identities, with blunt ends (blastall -m 6, not available with Blast+)
 *
 * @val BlastFormatFile BlastFormatFile::XML;
 * @brief XML (blastall -m 7 &amp; blast* -outfmt 8)
 *
 * @val BlastFormatFile BlastFormatFile::Tabular;
 * @brief tab-seperated (blastall -m 8 &amp; blast* -outfmt 6)
 *
 * @val BlastFormatFile BlastFormatFile::TabularWithHeader;
 * @brief tab-seperated with Header / comments (blastall -m 9 &amp; blast* -outfmt 7)
 *
 * @val BlastFormatFile BlastFormatFile::TextASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 8, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::BinASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 9, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::CSV;
 * @brief comma-seperated values (blast* -outfmt 10, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::BlastArASN1;
 * @brief Blast Archive Format, Abstract Syntax Notation One (blast* -outfmt 11, not available in traditional BLAST)
 */

enum class BlastFormatFile : uint8_t
{
    Pairwise = 0,
    MasterSlaveIdent = 1,
    MasterSlaveNoIdent = 2,
    FlatMasterSlaveIdent = 3,
    FlatMasterSlaveNoIdent = 4,
    MasterSlaveBluntEnds = 5,       // only available in Generation==Blast
    FlatMasterSlaveBluntEnds = 6,   // only available in Generation==Blast
    XML = 7,
    Tabular = 8,
    TabularWithHeader = 9,
    TextASN1 = 10,                 // only available in Generation==Blast+
    BinASN1 = 11,                  // only available in Generation==Blast+
    CSV = 12,                      // only available in Generation==Blast+
    BlastArASN1 = 13,              // only available in Generation==Blast+
    INVALID_File=255
};

/*!
 * @enum BlastFormatProgram
 * @brief Enum with BLAST program spec
 * @signature enum class BlastFormatProgram : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * @val BlastFormatProgram BlastFormatProgram::BlastN
 * @brief Nucleotide Query VS Nucleotide Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::BlastP
 * @brief Protein Query VS Protein Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::BlastX
 * @brief translated Nucleotide Query VS Protein Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::TBlastN
 * @brief Protein Query VS translated Nucleotide Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::TBlastX
 * @brief translated Nucleotide Query VS translated Nucleotide Subject
 *
 */
enum class BlastFormatProgram : uint8_t
{
    BlastN,         //              nucl VS             nucl
    BlastP,         //              prot VS             prot
    BlastX,         // translated   nucl VS             prot
    TBlastN,        //              prot VS translated  nucl
    TBlastX,        // translated   nucl VS translated  nucl
    INVALID_Program=255
};

/*!
 * @enum BlastFormatGeneration
 * @brief Enum with BLAST program generation
 * @signature enum class BlastFormatGeneration : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * @val BlastFormatGeneration BlastFormatGeneration::Blast
 * @brief traditional NCBI Blast, written in C ("blastall")
 *
 * @val BlastFormatGeneration BlastFormatGeneration::BlastPlus
 * @brief NCBI Blast+, written in C++
 *
 */

enum class BlastFormatGeneration : uint8_t
{
    Blast,
    BlastPlus,
    INVALID_Generation=255
};


/*!
 * @class BlastFormat
 *
 * @brief Blast Format specifier
 *
 * @signature template <BlastFormatFile _f, BlastFormatProgram _p, BlastFormatGeneration _g>
 *            struct BlastFormat;
 *
 * @headerfile seqan/blast.h
 * @tparam _f    File Type Format
 * @see BlastFormatFile
 * @tparam _p    Program Type Format
 * @see BlastFormatProgram
 * @tparam _g    Program Generation
 * @see BlastFormatGeneration
 */


template <BlastFormatFile _f, BlastFormatProgram _p, BlastFormatGeneration _g>
struct BlastFormat
{
    // have static members for easy and run-time access to "type"
    static constexpr BlastFormatFile        f = _f;
    static constexpr BlastFormatProgram     p = _p;
    static constexpr BlastFormatGeneration  g = _g;

};

// template <BlastFormatFile _f, BlastFormatProgram _p, BlastFormatGeneration _g>
// using BlastFormat = Tag<BlastFormat_<_f, _p, _g>>;



// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// getBlastProgramType()
// ----------------------------------------------------------------------------

/*!
 * @mfn getBlastProgramType
 * @headerfile seqan/blast.h
 * @brief for given query and subject alphabets, return the @link BlastFormatProgram @endlink
 * @signature  getBlastProgramType(TQueryAlph const &, TSubjAlph const &)
 *
 * @tparam          TQueryAlph  Alphabet of of the query sequences
 * @tparam          TSubjAlph   Alphabet of of the subejct sequences
 * @return          Value of @link BlastFormatProgram @endlink . For nucl VS nucl BlastFormatProgram::BlastN is returned, although TBlastX would also be legal.
 *
 * This is a convenience function for determining the Blast-Program type
 * corresponding to input type combinations.
 *
 * This is a constexpr function.
 */

template< typename TQueryAlph, typename TSubjAlph>
constexpr
BlastFormatProgram
getBlastProgramType(TQueryAlph const &, TSubjAlph const &)
{
    return BlastFormatProgram::INVALID_Program;
}

template<typename TQueryAlph, typename TSubjAlph,
         typename TSpec, typename TSpec2>
constexpr
BlastFormatProgram
getBlastProgramType(String<TQueryAlph, TSpec> const &,
                    String<TSubjAlph, TSpec2> const &)
{
    // needs constexpr constructors of Alphabet types
    return getBlastProgramType(TQueryAlph(), TSubjAlph());
}

// --- DNA vs DNA ---
// NOTE that Dna VS Dna could also be TBlastX, but BlastN is more common
constexpr
BlastFormatProgram
getBlastProgramType(Dna const &, Dna const &)
{
    return BlastFormatProgram::BlastN;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna const &, Dna5 const &)
{
    return BlastFormatProgram::BlastN;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna5 const &, Dna const &)
{
    return BlastFormatProgram::BlastN;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna5 const &, Dna5 const &)
{
    return BlastFormatProgram::BlastN;
}

// --- Protein vs Protein ---
constexpr
BlastFormatProgram
getBlastProgramType(AminoAcid const &, AminoAcid const &)
{
    return BlastFormatProgram::BlastP;
}

// --- Dna vs Protein ---
constexpr
BlastFormatProgram
getBlastProgramType(Dna const &, AminoAcid const &)
{
    return BlastFormatProgram::BlastX;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna5 const &, AminoAcid const &)
{
    return BlastFormatProgram::BlastX;
}

// --- Protein vs Dna ---
constexpr
BlastFormatProgram
getBlastProgramType(AminoAcid const &, Dna const &)
{
    return BlastFormatProgram::TBlastX;
}

constexpr
BlastFormatProgram
getBlastProgramType(AminoAcid const &, Dna5 const &)
{
    return BlastFormatProgram::TBlastX;
}

// ============================================================================
// Functions
// ============================================================================

/*!
 * @defgroup BlastScoringScheme
 * @brief functions for converting to and from Blast's scoring-scheme behaviour
 *
 * Blast (and many other tools) compute scores of a stretch of gaps as
 * <tt>s = gO + n * gE</tt>
 * where gO is the gapOpen score, gE is the gap extend score and n ist the
 * total number of gaps.
 *
 * SeqAn, however, computes them as as <tt>s = gO + (n-1) * gE</tt>. The
 * functions below convert between the behaviours by adjusting the
 * gapOpen score.
 *
 * For more information, see <a href="https://trac.seqan.de/ticket/1091">https://trac.seqan.de/ticket/1091</a>.
 *
 * Please note that independent of this issue, SeqAn always works with
 * scores, never with penalties, i.e. a penalty is represented by a negative
 * score.
 *
 * @fn BlastScoringScheme#seqanScoringScheme2blastScoringScheme
 * @signature void seqanScoringScheme2blastScoringScheme(scoringScheme);
 * @brief convert to Blast's behaviour
 * @param[in,out]      scoringScheme      The @link Score @endlink object to modify.
 *
 * @fn BlastScoringScheme#blastScoringScheme2seqanScoringScheme
 * @signature void blastScoringScheme2seqanScoringScheme(scoringScheme);
 * @brief convert from Blast's behaviour
 * @param[in,out]      scoringScheme      The @link Score @endlink object to modify.
 *
 */

template <typename TValue, typename TSpec>
inline void
seqanScoringScheme2blastScoringScheme(Score<TValue, TSpec> & scheme)
{
    setScoreGapOpen(scheme, scoreGapOpen(scheme) - scoreGapExtend(scheme));
}

template <typename TValue, typename TSpec>
inline void
blastScoringScheme2seqanScoringScheme(Score<TValue, TSpec> & scheme)
{
    setScoreGapOpen(scheme, scoreGapOpen(scheme) + scoreGapExtend(scheme));
}

// ----------------------------------------------------------------------------
// _programTagToString()
// ----------------------------------------------------------------------------

template <BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             p,
                                             g> const & /*tag*/)
{
    return "UNKOWN BLAST PROGRAM";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::BlastN,
                                             g> const & /*tag*/)
{
    return "BlastN";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::BlastP,
                                             g> const & /*tag*/)
{
    return "BLASTP";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::BlastX,
                                             g> const & /*tag*/)
{
    return "BLASTX";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::TBlastN,
                                             g> const & /*tag*/)
{
    return "TBLASTN";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::TBlastX,
                                             g> const & /*tag*/)
{
    return "TBLASTX";
}

// ----------------------------------------------------------------------------
// _defaultFields()
// ----------------------------------------------------------------------------

template <BlastFormatGeneration g>
constexpr
const char * _defaultFields()
{
    return "ERROR Fields not specializied for this type";
}

template <>
constexpr
const char * _defaultFields<BlastFormatGeneration::Blast>()
{
    return "Fields: Query id, Subject id, % identity, alignment length," \
           " mismatches, gap openings, q. start, q. end, s. start, s." \
           " end, e-value, bit score";
}

template <>
constexpr
const char * _defaultFields<BlastFormatGeneration::BlastPlus>()
{
    return "Fields: query id, subject id, % identity, alignment " \
           "length, mismatches, gap opens, q. start, q. end, s. " \
           "start, s. end, evalue, bit score";
}

template <BlastFormatProgram p, BlastFormatGeneration g>
constexpr
const char * _defaultFields(BlastFormat<BlastFormatFile::TabularWithHeader,
                                        p,
                                        g> const &)
{
    return _defaultFields<g>();
}


} // namespace seqan

#endif // header guard
