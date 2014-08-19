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
 * SeqAn can currently write PAIRWISE, TABULAR and TABULAR_WITH_HEADER. It can
 * read TABULAR and TABULAR_WITH_HEADER.
 *
 * @val BlastFormatFile BlastFormatFile::PAIRWISE;
 * @brief default blast output (blastall -m 0 &amp; blast* -outfmt 0)
 *
 * @val BlastFormatFile BlastFormatFile::MASTER_SLAVE_IDENT;
 * @brief master-slave showing identities (blastall -m 1 &amp; blast* -outfmt 1)
 *
 * @val BlastFormatFile BlastFormatFile::MASTER_SLAVE_NO_IDENT;
 * @brief master-slave without identities (blastall -m 2 &amp; blast* -outfmt 2)
 *
 * @val BlastFormatFile BlastFormatFile::FLAT_MASTER_SLAVE_IDENT;
 * @brief flat master-slave showing identities (blastall -m 3 &amp; blast* -outfmt 3)
 *
 * @val BlastFormatFile BlastFormatFile::FLAT_MASTER_SLAVE_NO_IDENT;
 * @brief master-slave without identities, with blunt ends (blastall -m 5, not available with Blast+)
 *
 * @val BlastFormatFile BlastFormatFile::MASTER_SLAVE_BLUNT_ENDS;
 * @brief flat master-slave without identities, with blunt ends (blastall -m 6, not available with Blast+)
 *
 * @val BlastFormatFile BlastFormatFile::XML;
 * @brief XML (blastall -m 7 &amp; blast* -outfmt 8)
 *
 * @val BlastFormatFile BlastFormatFile::TABULAR;
 * @brief tab-seperated (blastall -m 8 &amp; blast* -outfmt 6)
 *
 * @val BlastFormatFile BlastFormatFile::TABULAR_WITH_HEADER;
 * @brief tab-seperated with Header / comments (blastall -m 9 &amp; blast* -outfmt 7)
 *
 * @val BlastFormatFile BlastFormatFile::TEXT_ASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 8, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::BIN_ASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 9, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::CSV;
 * @brief comma-seperated values (blast* -outfmt 10, not available in traditional BLAST)
 *
 * @val BlastFormatFile BlastFormatFile::BLAST_AR_ASN1;
 * @brief Blast Archive Format, Abstract Syntax Notation One (blast* -outfmt 11, not available in traditional BLAST)
 */

enum class BlastFormatFile : uint8_t
{
    PAIRWISE = 0,
    MASTER_SLAVE_IDENT = 1,
    MASTER_SLAVE_NO_IDENT = 2,
    FLAT_MASTER_SLAVE_IDENT = 3,
    FLAT_MASTER_SLAVE_NO_IDENT = 4,
    MASTER_SLAVE_BLUNT_ENDS = 5,       // only available in Generation==Blast
    FLAT_MASTER_SLAVE_BLUNT_ENDS = 6,   // only available in Generation==Blast
    XML = 7,
    TABULAR = 8,
    TABULAR_WITH_HEADER = 9,
    TEXT_ASN1 = 10,                 // only available in Generation==Blast+
    BIN_ASN1 = 11,                  // only available in Generation==Blast+
    CSV = 12,                      // only available in Generation==Blast+
    BLAST_AR_ASN1 = 13,              // only available in Generation==Blast+
    INVALID_File=255
};

/*!
 * @enum BlastFormatProgram
 * @brief Enum with BLAST program spec
 * @signature enum class BlastFormatProgram : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * @val BlastFormatProgram BlastFormatProgram::BLASTN
 * @brief Nucleotide Query VS Nucleotide Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::BLASTP
 * @brief Protein Query VS Protein Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::BLASTX
 * @brief translated Nucleotide Query VS Protein Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::TBLASTN
 * @brief Protein Query VS translated Nucleotide Subject
 *
 * @val BlastFormatProgram BlastFormatProgram::TBLASTX
 * @brief translated Nucleotide Query VS translated Nucleotide Subject
 *
 */
enum class BlastFormatProgram : uint8_t
{
    BLASTN,         //              NUCL VS             NUCL
    BLASTP,         //              PROT VS             PROT
    BLASTX,         // TRANSLATED   NUCL VS             PROT
    TBLASTN,        //              PROT VS TRANSLATED  NUCL
    TBLASTX,        // TRANSLATED   NUCL VS TRANSLATED  NUCL
    INVALID_Program=255
};

/*!
 * @enum BlastFormatGeneration
 * @brief Enum with BLAST program generation
 * @signature enum class BlastFormatGeneration : uint8_t { ... };
 *
 * @headerfile seqan/blast.h
 *
 * @val BlastFormatGeneration BlastFormatGeneration::BLAST
 * @brief traditional Blast, written in C ("blastall" binary); all behaviour related
 * to this Format is based on NCBI BLAST-2.2.26
 *
 * @val BlastFormatGeneration BlastFormatGeneration::BLAST_PLUS
 * @brief Blast+, written in C++; all behaviour related
 * to this Format is based on NCBI BLAST-2.2.27+
 *
 */

enum class BlastFormatGeneration : uint8_t
{
    BLAST,
    BLAST_PLUS,
    INVALID_Generation=255
};

/*!
 * @class BlastFormat
 * @headerfile seqan/blast.h
 * @brief Blast Format specifier
 *
 * @signature template <BlastFormatFile _f, BlastFormatProgram _p, BlastFormatGeneration _g>
 *            struct BlastFormat;
 *
 * @section Data structures
 *
 * The easiest way to use the BLAST module is to organize your data in the
 * structures provided, i.e. create records (@link BlastRecord @endlink)
 * that contain all matches (@link BlastMatch @endlink) of one query
 * sequence. Store some general information of your database in a
 * @link BlastDbSpecs @endlink object.
 *
 *
 * @section E-Value statistics
 *
 * This module provides
 * 
 * @section Output
 *
 * When writing Blast(-like) results to a file the easiest way to proceed
 * is to organize your results into
 *
 * When writing this data, call @link BlastFormat#writeTop @endlink
 * first, than iterate over the records, calling
 * @link BlastRecord#writeRecord @endlink on each, and finishing with
 * @link BlastFormat#writeBottom @endlink .
 * Certain BlastFormats have additional more fine-grained functions, e.g.
 * @link BlastRecord#writeHeader @endlink and
 * @link BlastMatch#writeMatch @endlink .
 *
 * @tparam _f    File Type Format
 * @see BlastFormatFile
 * @tparam _p    Program Type Format
 * @see BlastFormatProgram
 * @tparam _g    Program Generation
 * @see BlastFormatGeneration
 */

template <BlastFormatFile _f, BlastFormatProgram _p, BlastFormatGeneration _g>
struct BlastFormat_
{
    // have static members for easy and run-time access to "type"
    static constexpr BlastFormatFile        f = _f;
    static constexpr BlastFormatProgram     p = _p;
    static constexpr BlastFormatGeneration  g = _g;

};

template <BlastFormatFile _f, BlastFormatProgram _p, BlastFormatGeneration _g>
using BlastFormat = Tag<BlastFormat_<_f, _p, _g>>;

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
 * @return          Value of @link BlastFormatProgram @endlink . For nucl VS nucl BlastFormatProgram::BLASTN is returned, although TBlastX would also be legal.
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
    return BlastFormatProgram::BLASTN;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna const &, Dna5 const &)
{
    return BlastFormatProgram::BLASTN;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna5 const &, Dna const &)
{
    return BlastFormatProgram::BLASTN;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna5 const &, Dna5 const &)
{
    return BlastFormatProgram::BLASTN;
}

// --- Protein vs Protein ---
constexpr
BlastFormatProgram
getBlastProgramType(AminoAcid const &, AminoAcid const &)
{
    return BlastFormatProgram::BLASTP;
}

// --- Dna vs Protein ---
constexpr
BlastFormatProgram
getBlastProgramType(Dna const &, AminoAcid const &)
{
    return BlastFormatProgram::BLASTX;
}

constexpr
BlastFormatProgram
getBlastProgramType(Dna5 const &, AminoAcid const &)
{
    return BlastFormatProgram::BLASTX;
}

// --- Protein vs Dna ---
constexpr
BlastFormatProgram
getBlastProgramType(AminoAcid const &, Dna const &)
{
    return BlastFormatProgram::TBLASTX;
}

constexpr
BlastFormatProgram
getBlastProgramType(AminoAcid const &, Dna5 const &)
{
    return BlastFormatProgram::TBLASTX;
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
                                             BlastFormatProgram::BLASTN,
                                             g> const & /*tag*/)
{
    return "BlastN";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::BLASTP,
                                             g> const & /*tag*/)
{
    return "BLASTP";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::BLASTX,
                                             g> const & /*tag*/)
{
    return "BLASTX";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::TBLASTN,
                                             g> const & /*tag*/)
{
    return "TBLASTN";
}

template <BlastFormatFile f, BlastFormatGeneration g>
constexpr
const char * _programTagToString(BlastFormat<f,
                                             BlastFormatProgram::TBLASTX,
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
const char * _defaultFields<BlastFormatGeneration::BLAST>()
{
    return "Query id, Subject id, % identity, alignment length," \
           " mismatches, gap openings, q. start, q. end, s. start, s." \
           " end, e-value, bit score";
}

template <>
constexpr
const char * _defaultFields<BlastFormatGeneration::BLAST_PLUS>()
{
    return "query id, subject id, % identity, alignment " \
           "length, mismatches, gap opens, q. start, q. end, s. " \
           "start, s. end, evalue, bit score";
}

template <BlastFormatProgram p, BlastFormatGeneration g>
constexpr
const char * _defaultFields(BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        p,
                                        g> const &)
{
    return _defaultFields<g>();
}

template <BlastFormatProgram p, BlastFormatGeneration g>
constexpr
const char * _defaultFields(BlastFormat<BlastFormatFile::TABULAR,
                                        p,
                                        g> const &)
{
    return "";
}


// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#writeTop
 * @headerfile seqan/blast.h
 * @brief write the top-most section of a BLAST output file (NO-OP for tabular formats)
 * @signature int writeTop(stream, blastDbSpecs, blastFormatTag)
 *
 * @param stream            The file to write to (FILE, fstream, @link Stream @endlink ...)
 * @param blastDbSpecs      The @link BlastDbSpecs @endlink of your database-
 * @param blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastMatch
 */

template <typename TStream,
          typename TDbSpecs,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
writeTop(TStream                    & /**/,
         TDbSpecs             const & /**/,
         BlastFormat<f, p, g> const & /*tag*/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastRecord#writeRecord
 * @headerfile seqan/blast.h
 * @brief write a @link BlastRecord @endlink inluding it's @link BlastMatch @endlink es to a file.
 * @signature int writeRecord(stream, blastRecord, blastDbSpecs. blastFormatTag)
 *
 * @param stream        The file to write to (FILE, fstream, @link Stream @endlink ...)
 * @param blastRecord   The @link BlastRecord @endlink you wish to print.
 * @param blastDbSpecs  The @link BlastDbSpecs @endlink .
 * @param blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 */


// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#writeBottom
 * @headerfile seqan/blast.h
 * @brief write the top-most section of a BLAST output file (NO-OP for tabular formats)
 * @signature int writeBottom(stream, blastDbSpecs, scoringAdapter, blastFormatTag)
 *
 * @param stream            The file to write to (FILE, fstream, @link Stream @endlink ...)
 * @param scoringAdapter    A @link BlastScoringAdapter @endlink with relevant information.
 * @param blastDbSpecs      The @link BlastDbSpecs @endlink of your database-
 * @param blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastMatch
 */

template <typename TStream,
          typename TDbSpecs,
          typename TBlastScoringAdapater,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
writeBottom(TStream                           & /**/,
            TDbSpecs                    const & /**/,
            TBlastScoringAdapater       const & /**/,
            BlastFormat<f, p,g>         const & /*tag*/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function _untranslatePositions() -- retransform positions
// ----------------------------------------------------------------------------

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & /**/,
                      signed char const /**/,
                      False const & /*hasReverseComplement*/,
                      False const & /*hasFrames*/)
{
    // BLAST is 1-indexed, but end positions are "on" instead of behind
    // so only the begin positions need adapting
    ++effectiveStart;
}

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & effectiveEnd,
                      signed char const frameShift,
                      True const & /*hasReverseComplement*/,
                      False const & /*hasFrames*/)
{
    // BLAST is 1-indexed, but end positions are "on" instead of behind
    // so only the begin positions need adapting
    ++effectiveStart;
    // reverse strand symoblized by swapped begin and end
    if (frameShift < 0)
        std::swap(effectiveStart, effectiveEnd);
}

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & effectiveEnd,
                      signed char const frameShift,
                      True const & /*hasReverseComplement*/,
                      True const & /*hasFrames*/)
{
    // correct for codon translation and frameshift
    // subtract 1 because frameshift is 1-indexed
    effectiveStart = effectiveStart * 3 + std::abs(frameShift) - 1;
    effectiveEnd = effectiveEnd * 3 + std::abs(frameShift) - 1;

    _untranslatePositions(effectiveStart, effectiveEnd, frameShift, True(), False());
}

// ----------------------------------------------------------------------------
// Function _nextPos() -- increment/decrement position
// ----------------------------------------------------------------------------

constexpr int8_t
_step(signed char const /**/,
      False const & /*hasReverseComplement*/,
      False const & /*hasFrames*/)
{
    return 1;
}

constexpr int8_t
_step(signed char const frameShift,
      True const & /*hasReverseComplement*/,
      False const & /*hasFrames*/)
{
    // iterate backwards for reverse frames
    return (frameShift < 0) ? -1 : 1;
}

constexpr int8_t
_step(signed char const frameShift,
      True const & /*hasReverseComplement*/,
      True const & /*hasFrames*/)
{
    // iterate three nucleotides per amino acid
    return (frameShift < 0) ? -3 : 3;
}

} // namespace seqan

#endif // header guard
