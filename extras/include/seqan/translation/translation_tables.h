// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Code for Dna(5) to AminoAcid Translation - Conversion Tables
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_TABLES_H_
#define EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_TABLES_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @enum GeneticCodeSpec GeneticCode Specs
 * @brief Specialization values for @link GeneticCode @endlink
 * @headerfile seqan/translation.h
 *
 * @signature enum class  GeneticCodeSpec : uint8_t { ...};
 *
 * @see translate
 * @see GeneticCode
 *
 * The numeric values of the enums correspond to the genbank transl_table values
 * (see <a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi">http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi</a>).
 * Some genetic codes have been deprecated, so not all numeric values are available.
 *
 * Please not that this is part of the translation module which requires C++11.
 *
 * @val GeneticCodeSpec GeneticCodeSpec::CANONICAL = 1
 * @brief The Standard Genetic Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::VERT_MITOCHONDRIAL = 2
 * @brief The Vertebrate Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::YEAST_MITOCHONDRIAL = 3
 * @brief The Yeast Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::MOLD_MITOCHONDRIAL = 4
 * @brief The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::INVERT_MITOCHONDRIAL = 5
 * @brief The Invertebrate Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::CILIATE = 6
 * @brief The CILIATE, Dasycladacean and Hexamita Nuclear Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::FLATWORM_MITOCHONDRIAL = 9
 * @brief The Echinoderm and Flatworm Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::EUPLOTID = 10
 * @brief The EUPLOTID Nuclear Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::PROKARYOTE = 11
 * @brief The Bacterial, Archaeal and Plant Plastid Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::ALT_YEAST = 12
 * @brief The Alternative Yeast Nuclear Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::ASCIDIAN_MITOCHONDRIAL = 13
 * @brief The Ascidian Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::ALT_FLATWORM_MITOCHONDRIAL = 14
 * @brief The Alternative Flatworm Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::BLEPHARISMA = 15
 * @brief BLEPHARISMA Nuclear Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::CHLOROPHYCEAN_MITOCHONDRIAL = 16
 * @brief Chlorophycean Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::TREMATODE_MITOCHONDRIAL = 21
 * @brief Trematode Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::SCENEDESMUS_MITOCHONDRIAL = 22
 * @brief Scenedesmus obliquus mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::THRAUSTOCHYTRIUM_MITOCHONDRIAL = 23
 * @brief Thraustochytrium Mitochondrial Code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::PTEROBRANCHIA_MITOCHONDRIAL = 24
 * @brief Pterobranchia mitochondrial code
 *
 * @val GeneticCodeSpec GeneticCodeSpec::GRACILIBACTERIA = 25
 * @brief Candidate Division SR1 and GRACILIBACTERIA Code
 */

enum class GeneticCodeSpec : uint8_t
{
    CANONICAL=1,
    VERT_MITOCHONDRIAL,
    YEAST_MITOCHONDRIAL,
    MOLD_MITOCHONDRIAL,
    INVERT_MITOCHONDRIAL,
    CILIATE,
    FLATWORM_MITOCHONDRIAL = 9,
    EUPLOTID,
    PROKARYOTE,
    ALT_YEAST,
    ASCIDIAN_MITOCHONDRIAL,
    ALT_FLATWORM_MITOCHONDRIAL,
    BLEPHARISMA,
    CHLOROPHYCEAN_MITOCHONDRIAL,
    TREMATODE_MITOCHONDRIAL = 21,
    SCENEDESMUS_MITOCHONDRIAL,
    THRAUSTOCHYTRIUM_MITOCHONDRIAL,
    PTEROBRANCHIA_MITOCHONDRIAL,
    GRACILIBACTERIA
};

// -----------------------------------------------------------------------
// Tag GeneticCode
// -----------------------------------------------------------------------

/*!
 * @tag GeneticCode
 * @brief DNA/RNA to AminoAcid translation code, needs to be spec'ed by one of @link GeneticCodeSpec @endlink.
 * @headerfile seqan/translation.h
 * @signature GeneticCode<GeneticCodeSpec::value>
 *
 * @see translate
 * @see GeneticCodeSpec
 *
 */

template <GeneticCodeSpec = GeneticCodeSpec::CANONICAL>
struct GeneticCode
{};


// -----------------------------------------------------------------------
// struct TranslateTableDnaToAminoAcid_
// -----------------------------------------------------------------------


template <typename TCodeSpec, typename TVoidSpec = void>
struct TranslateTableDnaToAminoAcid_;

/* Tables according to:
 * http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 */

//  DEFINITIONS

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::CANONICAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::VERT_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE [4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { '*', 'S', '*', 'S' }, // ag?
            { 'M', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::YEAST_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'M', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'T', 'T', 'T', 'T' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::MOLD_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::INVERT_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'S', 'S', 'S', 'S' }, // ag?
            { 'M', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::CILIATE>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { 'Q', 'Y', 'Q', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::FLATWORM_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'N', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'S', 'S', 'S', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::EUPLOTID>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'C', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::PROKARYOTE>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::ALT_YEAST>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'S', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::ASCIDIAN_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'G', 'S', 'G', 'S' }, // ag?
            { 'M', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::ALT_FLATWORM_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'N', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'S', 'S', 'S', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { 'Y', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::BLEPHARISMA>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', 'Q', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::CHLOROPHYCEAN_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', 'L', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::TREMATODE_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'N', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'S', 'S', 'S', 'S' }, // ag?
            { 'M', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::SCENEDESMUS_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', 'L', 'Y' }, // ua?
            { '*', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::THRAUSTOCHYTRIUM_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { '*', 'C', 'W', 'C' }, // ug?
            { '*', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::PTEROBRANCHIA_MITOCHONDRIAL>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'S', 'S', 'K', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'W', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

template <typename TVoidSpec>
struct TranslateTableDnaToAminoAcid_<GeneticCode<
        GeneticCodeSpec::GRACILIBACTERIA>, TVoidSpec>
{
    static constexpr char VALUE[4][4][4]
    {
        { // a??
            { 'K', 'N', 'K', 'N' }, // aa?
            { 'T', 'T', 'T', 'T' }, // ac?
            { 'R', 'S', 'R', 'S' }, // ag?
            { 'I', 'I', 'M', 'I' }  // au?
        }, { // c??
            { 'Q', 'H', 'Q', 'H' }, // ca?
            { 'P', 'P', 'P', 'P' }, // cc?
            { 'R', 'R', 'R', 'R' }, // cg?
            { 'L', 'L', 'L', 'L' }  // cu?
        }, { // g??
            { 'E', 'D', 'E', 'D' }, // ga?
            { 'A', 'A', 'A', 'A' }, // gc?
            { 'G', 'G', 'G', 'G' }, // gg?
            { 'V', 'V', 'V', 'V' }  // gu?
        }, { // u??
            { '*', 'Y', '*', 'Y' }, // ua?
            { 'S', 'S', 'S', 'S' }, // uc?
            { 'G', 'C', 'W', 'C' }, // ug?
            { 'L', 'F', 'L', 'F' }  // uu?
        }
        //     a    c    g    u
    };
};

//  DECLARATIONS

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::CANONICAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::VERT_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::YEAST_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::MOLD_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::INVERT_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::CILIATE>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::FLATWORM_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::EUPLOTID>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::PROKARYOTE>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::ALT_YEAST>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::ASCIDIAN_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::ALT_FLATWORM_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::BLEPHARISMA>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::CHLOROPHYCEAN_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::TREMATODE_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::SCENEDESMUS_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::THRAUSTOCHYTRIUM_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::PTEROBRANCHIA_MITOCHONDRIAL>,
TVoidSpec>::VALUE[4][4][4];

template <typename TVoidSpec>
constexpr char
TranslateTableDnaToAminoAcid_<GeneticCode<GeneticCodeSpec::GRACILIBACTERIA>,
TVoidSpec>::VALUE[4][4][4];

}

#endif // EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_TABLES_H_
