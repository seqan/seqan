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
 * @see translate
 * @see GeneticCode
 *
 * Some genetic codes have been deprecated, so not all transl_table-numbers are present.
 *
 * This is a scoped enum, only available with C++11.
 *
 * @val GeneticCodeSpec GeneticCodeSpec::Canonical = 1
 * @brief The Standard Genetic Code (genbank transl_table=1)
 *
 * @val GeneticCodeSpec GeneticCodeSpec::VertMitochondrial = 2
 * @brief The Vertebrate Mitochondrial Code (genbank transl_table=2)
 *
 * @val GeneticCodeSpec GeneticCodeSpec::YeastMitochondrial = 3
 * @brief The Yeast Mitochondrial Code (genbank transl_table=3)
 *
 * @val GeneticCodeSpec GeneticCodeSpec::MoldMitochondrial = 4
 * @brief The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
 *        (genbank transl_table=4).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::InvertMitochondrial = 5
 * @brief The Invertebrate Mitochondrial Code (genbank transl_table=5)
 *
 * @val GeneticCodeSpec GeneticCodeSpec::Ciliate = 6
 * @brief The Ciliate, Dasycladacean and Hexamita Nuclear Code (genbank transl_table=6).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::FlatwormMitochondrial = 9
 * @brief The Echinoderm and Flatworm Mitochondrial Code (genbank transl_table=9).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::Euplotid = 10
 * @brief The Euplotid Nuclear Code (genbank transl_table=10).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::Prokaryote = 11
 * @brief The Bacterial, Archaeal and Plant Plastid Code (genbank transl_table=11).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::AltYeast = 12
 * @brief The Alternative Yeast Nuclear Code (genbank transl_table=12).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::AscidianMitochondrial = 13
 * @brief The Ascidian Mitochondrial Code (genbank transl_table=13).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::AltFlatwormMitochondrial = 14
 * @brief The Alternative Flatworm Mitochondrial Code (genbank transl_table=14).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::Blepharisma = 15
 * @brief Blepharisma Nuclear Code (genbank transl_table=15).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::ChlorophyceanMitochondrial = 16
 * @brief Chlorophycean Mitochondrial Code (genbank transl_table=16).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::TrematodeMitochondrial = 21
 * @brief Trematode Mitochondrial Code (genbank transl_table=21).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::ScenedesmusMitochondrial = 22
 * @brief Scenedesmus obliquus mitochondrial Code (genbank transl_table=22).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::ThraustochytriumMitochondrial = 23
 * @brief Thraustochytrium Mitochondrial Code (genbank transl_table=23).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::PterobranchiaMitochondrial = 24
 * @brief Pterobranchia mitochondrial code (genbank transl_table=24).
 *
 * @val GeneticCodeSpec GeneticCodeSpec::Gracilibacteria = 25
 * @brief Candidate Division SR1 and Gracilibacteria Code (genbank transl_table=25).
 */

enum class GeneticCodeSpec : uint8_t
{
    Canonical=1,
    VertMitochondrial,
    YeastMitochondrial,
    MoldMitochondrial,
    InvertMitochondrial,
    Ciliate,
    FlatwormMitochondrial = 9,
    Euplotid,
    Prokaryote,
    AltYeast,
    AscidianMitochondrial,
    AltFlatwormMitochondrial,
    Blepharisma,
    ChlorophyceanMitochondrial,
    TrematodeMitochondrial = 21,
    ScenedesmusMitochondrial,
    ThraustochytriumMitochondrial,
    PterobranchiaMitochondrial,
    Gracilibacteria
};

// -----------------------------------------------------------------------
// Tag GeneticCode
// -----------------------------------------------------------------------

/*!
 * @tag GeneticCode
 * @brief Dna to AminoAcid translation code, needs to be spec'ed by one of @link GeneticCodeSpec @endlink.
 * @signature GeneticCode<GeneticCodeSpec::value>
 *
 * @see translate
 * @see GeneticCodeSpec
 *
 */

template <GeneticCodeSpec = GeneticCodeSpec::Canonical>
struct GeneticCode
{};

// -----------------------------------------------------------------------
// struct TranslateTableDnaToAminoAcid_
// -----------------------------------------------------------------------


template <typename T = void>
struct TranslateTableDnaToAminoAcid_
{
    static char const VALUE[4][4][4];
};

/* Tables according to:
 * http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 */

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::Canonical> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::VertMitochondrial> >::VALUE [4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::YeastMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::MoldMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::InvertMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::Ciliate> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::FlatwormMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::Euplotid> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::Prokaryote> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::AltYeast> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::AscidianMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::AltFlatwormMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::Blepharisma> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::ChlorophyceanMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::TrematodeMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::ScenedesmusMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::ThraustochytriumMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::PterobranchiaMitochondrial> >::VALUE[4][4][4] =
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

template <>
char const TranslateTableDnaToAminoAcid_<
    GeneticCode<
        GeneticCodeSpec::Gracilibacteria> >::VALUE[4][4][4] =
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

}

#endif // EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_TABLES_H_
