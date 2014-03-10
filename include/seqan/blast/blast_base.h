// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
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

#ifndef SEQAN_CXX11_STANDARD //C++98
#define constexpr inline
#endif

namespace seqan {


// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BlastFormatOptions
{

/*!
 * @enum BlastFormatOptions::M
 * @brief Enum with BLAST file format specs
 *
 * @headerfile seqan/blast.h
 *
 * @var BlastFormatOptions::M Pairwise;
 * @brief default blast output (blastall -m 0 & blast* -outfmt 0)
 *
 * @var BlastFormatOptions::M MasterSlaveIdent;
 * @brief master-slave showing identities (blastall -m 1 & blast* -outfmt 1)
 *
 * @var BlastFormatOptions::M MasterSlaveNoIdent;
 * @brief master-slave without identities (blastall -m 2 & blast* -outfmt 2)
 *
 * @var BlastFormatOptions::M FlatMasterSlaveIdent;
 * @brief flat master-slave showing identities (blastall -m 3 & blast* -outfmt 3)
 *
 * @var BlastFormatOptions::M FlatMasterSlaveNoIdent;
 * @brief master-slave without identities, with blunt ends (blastall -m 5, not available with Blast+)
 *
 * @var BlastFormatOptions::M MasterSlaveBluntEnds;
 * @brief flat master-slave without identities, with blunt ends (blastall -m 6, not available with Blast+)
 *
 * @var BlastFormatOptions::M XML;
 * @brief XML (blastall -m 7 & blast* -outfmt 8)
 *
 * @var BlastFormatOptions::M Tabular;
 * @brief tab-seperated (blastall -m 8 & blast* -outfmt 6)
 *
 * @var BlastFormatOptions::M TabularWithHeader;
 * @brief tab-seperated with Header / comments (blastall -m 9 & blast* -outfmt 7)
 *
 * @var BlastFormatOptions::M TextASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 8, not available in traditional BLAST)
 *
 * @var BlastFormatOptions::M BinASN1;
 * @brief Abstract Syntax Notation One (blast* -outfmt 9, not available in traditional BLAST)
 *
 * @var BlastFormatOptions::M CSV;
 * @brief comma-seperated values (blast* -outfmt 10, not available in traditional BLAST)
 *
 * @var BlastFormatOptions::M BlastArASN1;
 * @brief Blast Archive Format, Abstract Syntax Notation One (blast* -outfmt 11, not available in traditional BLAST)
 *
 * @section Remarks
 *
 * SeqAn currently implements Pairwise, Tabular and TabularWithHeader
 */
    enum M
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
        INVALID_M=1023
    };

/*!
 * @enum BlastFormatOptions::Program
 * @brief Enum with BLAST program spec
 *
 * @headerfile seqan/blast.h
 *
 * @var BlastFormatOptions::Program BlastN
 * @brief Nucleotide Query VS Nucleotide Subject
 *
 * @var BlastFormatOptions::Program BlastP
 * @brief Protein Query VS Protein Subject
 *
 * @var BlastFormatOptions::Program BlastX
 * @brief translated Nucleotide Query VS Protein Subject
 *
 * @var BlastFormatOptions::Program TBlastN
 * @brief Protein Query VS translated Nucleotide Subject
 *
 * @var BlastFormatOptions::Program TBlastX
 * @brief translated Nucleotide Query VS translated Nucleotide Subject
 *
 */
    enum Program
    {
        BlastN,         //              nucl VS             nucl
        BlastP,         //              prot VS             prot
        BlastX,         // translated   nucl VS             prot
        TBlastN,        //              prot VS translated  nucl
        TBlastX,        // translated   nucl VS translated  nucl
        INVALID_Program=1023
    };

/*!
 * @enum BlastFormatOptions::Generation
 * @brief Enum with BLAST program generation
 *
 * @headerfile seqan/blast.h
 *
 * @var BlastFormatOptions::Generation Blast
 * @brief traditional NCBI Blast, written in C ("blastall")
 *
 * @var BlastFormatOptions::Generation BlastPlus
 * @brief NCBI Blast+, written in C++
 *
 */
    enum Generation
    {
        Blast,
        BlastPlus,
        INVALID_Generation=1023
    };
};

/*!
 * @class BlastFormat
 *
 * @brief Blast Format specifier
 *
 * @signature template <BlastFormatOptions::M _m,BlastFormatOptions::Program _p, BlastFormatOptions::Generation _g>
 *            struct BlastFormat;
 *
 * @headerfile seqan/blast.h
 * @tparam _m    File Type Format
 * @see BlastFormatOptions::M
 * @tparam _p    Program Type Format
 * @see BlastFormatOptions::Program
 * @tparam _g    Program Generation
 * @see BlastFormatOptions::Generation
 */

/*TODO(C++11): change struct BlastFormat to struct BlastFormat_ and wrap
a type-dependent typedef Tag around it */
template <BlastFormatOptions::M            _m,
          BlastFormatOptions::Program      _p,
          BlastFormatOptions::Generation   _g>
struct BlastFormat
{
    // have static members for run-time access to "type"
    static const BlastFormatOptions::M          m = _m;
    static const BlastFormatOptions::Program    p = _p;
    static const BlastFormatOptions::Generation g = _g;

};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// getBlastProgramType()
// ----------------------------------------------------------------------------

//TODO(h4nn3s): test, document
template< typename TQueryAlph, typename TSubjAlph>
constexpr
BlastFormatOptions::Program
getBlastProgramType(TQueryAlph const &, TSubjAlph const &)
{
    return BlastFormatOptions::INVALID_Program;
}

template<typename TQueryAlph, typename TSubjAlph,
         typename TSpec, typename TSpec2>
constexpr
BlastFormatOptions::Program
getBlastProgramType(String<TQueryAlph, TSpec> const &,
                    String<TSubjAlph, TSpec2> const &)
{
    // needs constexpr constructors of Alphabet types
    return getBlastProgramType(TQueryAlph(), TSubjAlph());
}

// --- DNA vs DNA ---
// NOTE that Dna VS Dna could also be TBlastX, but BlastN is more common
constexpr
BlastFormatOptions::Program
getBlastProgramType(Dna const &, Dna const &)
{
    return BlastFormatOptions::BlastN;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(Dna const &, Dna5 const &)
{
    return BlastFormatOptions::BlastN;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(Dna5 const &, Dna const &)
{
    return BlastFormatOptions::BlastN;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(Dna5 const &, Dna5 const &)
{
    return BlastFormatOptions::BlastN;
}

// --- Protein vs Protein ---
constexpr
BlastFormatOptions::Program
getBlastProgramType(AminoAcid const &, AminoAcid const &)
{
    return BlastFormatOptions::BlastP;
}

// --- Dna vs Protein ---
constexpr
BlastFormatOptions::Program
getBlastProgramType(Dna const &, AminoAcid const &)
{
    return BlastFormatOptions::BlastX;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(Dna5 const &, AminoAcid const &)
{
    return BlastFormatOptions::BlastX;
}

// --- Protein vs Dna ---
constexpr
BlastFormatOptions::Program
getBlastProgramType(AminoAcid const &, Dna const &)
{
    return BlastFormatOptions::TBlastX;
}

constexpr
BlastFormatOptions::Program
getBlastProgramType(AminoAcid const &, Dna5 const &)
{
    return BlastFormatOptions::TBlastX;
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// _programTagToString()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             p,
                                             g> const & /*tag*/)
{
    return "UNKOWN BLAST PROGRAM";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::BlastN,
                                             g> const & /*tag*/)
{
    return "BlastN";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::BlastP,
                                             g> const & /*tag*/)
{
    return "BLASTP";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::BlastX,
                                             g> const & /*tag*/)
{
    return "BLASTX";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::TBlastN,
                                             g> const & /*tag*/)
{
    return "TBLASTN";
}

template <BlastFormatOptions::M m, BlastFormatOptions::Generation g>
constexpr
const char * _programTagToString(BlastFormat<m,
                                             BlastFormatOptions::TBlastX,
                                             g> const & /*tag*/)
{
    return "TBLASTX";
}

// ----------------------------------------------------------------------------
// _defaultFields()
// ----------------------------------------------------------------------------

template <BlastFormatOptions::Generation g>
constexpr
const char * _defaultFields()
{
    return "ERROR Fields not specializied for this type";
}

template <>
constexpr
const char * _defaultFields<BlastFormatOptions::Blast>()
{
    return "Fields: Query id, Subject id, % identity, alignment length," \
           " mismatches, gap openings, q. start, q. end, s. start, s." \
           " end, e-value, bit score";
}

template <>
constexpr
const char * _defaultFields<BlastFormatOptions::BlastPlus>()
{
    return "Fields: query id, subject id, % identity, alignment " \
           "length, mismatches, gap opens, q. start, q. end, s. " \
           "start, s. end, evalue, bit score";
}

template <BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
constexpr
const char * _defaultFields(BlastFormat<BlastFormatOptions::TabularWithHeader,
                                        p,
                                        g> const &)
{
    return _defaultFields<g>();
}


} // namespace seqan

#ifndef SEQAN_CXX11_STANDARD //C++98
#undef constexpr
#endif

#endif // header guard
