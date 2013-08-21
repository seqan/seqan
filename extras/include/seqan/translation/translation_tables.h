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
 * @defgroup GeneticCodeSpecs GeneticCode Specs
 * @brief Specialization tags for @link GeneticCode @endlink
 *
 * @tag GeneticCodeSpecs#Canonical
 * @brief The Standard Genetic Code (genbank transl_table=1)
 *
 * @tag GeneticCodeSpecs#VertMitochondrial
 * @brief The Vertebrate Mitochondrial Code (genbank transl_table=2)
 *
 * @tag GeneticCodeSpecs#YeastMitochondrial
 * @brief The Yeast Mitochondrial Code (genbank transl_table=3)
 *
 * @tag GeneticCodeSpecs#MoldMitochondrial
 * @brief The Mold, Protozoan, and Coelenterate Mitochondrial Code and the
 * Mycoplasma/Spiroplasma Code (genbank transl_table=4)
 *
 * @tag GeneticCodeSpecs#InvertMitchondrial
 * @brief The Invertebrate Mitochondrial Code (genbank transl_table=5)
 *
 * @tag GeneticCodeSpecs#Ciliate
 * @brief The Ciliate, Dasycladacean and Hexamita Nuclear Code
 * (genbank transl_table=6)
 *
 * @tag GeneticCodeSpecs#FlatwormMitochondrial
 * @brief The Echinoderm and Flatworm Mitochondrial Code (genbank transl_table=9)
 *
 * @tag GeneticCodeSpecs#Euplotid
 * @brief The Euplotid Nuclear Code (genbank transl_table=10)
 *
 * @tag GeneticCodeSpecs#Prokaryote
 * @brief The Bacterial, Archaeal and Plant Plastid Code (genbank
 * transl_table=11)
 *
 * @tag GeneticCodeSpecs#AltYeast
 * @brief The Alternative Yeast Nuclear Code (genbank transl_table=12)
 *
 * @tag GeneticCodeSpecs#AscidianMitochondrial
 * @brief The Ascidian Mitochondrial Code (genbank transl_table=13)
 *
 * @tag GeneticCodeSpecs#AltFlatwormMitochondrial
 * @brief The Alternative Flatworm Mitochondrial Code (genbank transl_table=14)
 *
 * @tag GeneticCodeSpecs#Blepharisma
 * @brief Blepharisma Nuclear Code (genbank genbank transl_table=15)
 *
 * @tag GeneticCodeSpecs#ChlorophyceanMitochondrial
 * @brief Chlorophycean Mitochondrial Code (genbank genbank transl_table=16)
 *
 * @tag GeneticCodeSpecs#TrematodeMitochondrial
 * @brief Trematode Mitochondrial Code (genbank transl_table=21)
 *
 * @tag GeneticCodeSpecs#ScenedesmusMitochondrial
 * @brief Scenedesmus obliquus mitochondrial Code (genbank transl_table=22)
 *
 * @tag GeneticCodeSpecs#ThraustochytriumMitochondrial
 * @brief Thraustochytrium Mitochondrial Code (genbank transl_table=23)
 *
 * @tag GeneticCodeSpecs#PterobranchiaMitochondrial
 * @brief Pterobranchia mitochondrial code (genbank transl_table=24)
 *
 * @tag GeneticCodeSpecs#Gracilibacteria
 * @brief Candidate Division SR1 and Gracilibacteria Code (genbank
 * transl_table=25)
 *
 * @headerfile seqan/translation.h
 *
 * @see translate
 * @see GeneticCode
 *
 * @section Remarks
 *
 * Some genetic codes have been deprecated, so not all transl_table-numbers
 * are present.
 */

// -----------------------------------------------------------------------
// Tag Canonical
// -----------------------------------------------------------------------

struct Canonical_
{};

typedef Tag<Canonical_> Canonical;

// -----------------------------------------------------------------------
// Tag VertMitochondrial
// -----------------------------------------------------------------------

struct VertMitochondrial_
{};

typedef Tag<VertMitochondrial_> VertMitochondrial;

// -----------------------------------------------------------------------
// Tag YeastMitochondrial
// -----------------------------------------------------------------------

struct YeastMitochondrial_
{};

typedef Tag<YeastMitochondrial_> YeastMitochondrial;

// -----------------------------------------------------------------------
// Tag MoldMitochondrial
// -----------------------------------------------------------------------

struct MoldMitochondrial_
{};

typedef Tag<MoldMitochondrial_> MoldMitochondrial;

// -----------------------------------------------------------------------
// Tag InvertMitochondrial
// -----------------------------------------------------------------------

struct InvertMitochondrial_
{};

typedef Tag<InvertMitochondrial_> InvertMitochondrial;

// -----------------------------------------------------------------------
// Tag Ciliate
// -----------------------------------------------------------------------

struct Ciliate_
{};

typedef Tag<Ciliate_> Ciliate;

// -----------------------------------------------------------------------
// Tag FlatwormMitochondrial
// -----------------------------------------------------------------------

struct FlatwormMitochondrial_
{};

typedef Tag<FlatwormMitochondrial_> FlatwormMitochondrial;

// -----------------------------------------------------------------------
// Tag Euplotid
// -----------------------------------------------------------------------

struct Euplotid_
{};

typedef Tag<Euplotid_> Euplotid;

// -----------------------------------------------------------------------
// Tag Prokaryote
// -----------------------------------------------------------------------

struct Prokaryote_
{};

typedef Tag<Prokaryote_> Prokaryote;

// -----------------------------------------------------------------------
// Tag AltYeast
// -----------------------------------------------------------------------

struct AltYeast_
{};

typedef Tag<AltYeast_> AltYeast;

// -----------------------------------------------------------------------
// Tag AscidianMitochondrial
// -----------------------------------------------------------------------

struct AscidianMitochondrial_
{};

typedef Tag<AscidianMitochondrial_> AscidianMitochondrial;

// -----------------------------------------------------------------------
// Tag AltFlatwormMitochondrial
// -----------------------------------------------------------------------

struct AltFlatwormMitochondrial_
{};

typedef Tag<AltFlatwormMitochondrial_> AltFlatwormMitochondrial;

// -----------------------------------------------------------------------
// Tag Blepherisma
// -----------------------------------------------------------------------

struct Blepherisma_
{};

typedef Tag<Blepherisma_> Blepherisma;

// -----------------------------------------------------------------------
// Tag ChlorophyceanMitochondrial
// -----------------------------------------------------------------------

struct ChlorophyceanMitochondrial_
{};

typedef Tag<ChlorophyceanMitochondrial_> ChlorophyceanMitochondrial;

// -----------------------------------------------------------------------
// Tag TrematodeMitochondrial
// -----------------------------------------------------------------------

struct TrematodeMitochondrial_
{};

typedef Tag<TrematodeMitochondrial_> TrematodeMitochondrial;

// -----------------------------------------------------------------------
// Tag ScenedesmusMitochondrial
// -----------------------------------------------------------------------

struct ScenedesmusMitochondrial_
{};

typedef Tag<ScenedesmusMitochondrial_> ScenedesmusMitochondrial;

// -----------------------------------------------------------------------
// Tag ThraustochytriumMitochondrial
// -----------------------------------------------------------------------

struct ThraustochytriumMitochondrial_
{};

typedef Tag<ThraustochytriumMitochondrial_> ThraustochytriumMitochondrial;

// -----------------------------------------------------------------------
// Tag PterobranchiaMitochondrial
// -----------------------------------------------------------------------

struct PterobranchiaMitochondrial_
{};

typedef Tag<PterobranchiaMitochondrial_> PterobranchiaMitochondrial;

// -----------------------------------------------------------------------
// Tag Gracilibacteria
// -----------------------------------------------------------------------

struct Gracilibacteria_
{};

typedef Tag<Gracilibacteria_> Gracilibacteria;


// -----------------------------------------------------------------------
// Tag GeneticCode
// -----------------------------------------------------------------------

/*!
 * @tag GeneticCode
 * @brief Dna to AminoAcid translation code, needs to be spec'ed by
 * one of @link GeneticCodeSpecs @endlink
 * @signature GeneticCode<GeneticCodeSpec>
 *
 * @see translate
 * @see GeneticCodeSpecs
 *
 * @headerfile seqan/translation.h
 */

template <typename T = Canonical>
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
        Canonical> >::VALUE[4][4][4] =
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
        VertMitochondrial> >::VALUE [4][4][4] =
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
        YeastMitochondrial> >::VALUE[4][4][4] =
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
        MoldMitochondrial> >::VALUE[4][4][4] =
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
        InvertMitochondrial> >::VALUE[4][4][4] =
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
        Ciliate> >::VALUE[4][4][4] =
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
        FlatwormMitochondrial> >::VALUE[4][4][4] =
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
        Euplotid> >::VALUE[4][4][4] =
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
        Prokaryote> >::VALUE[4][4][4] =
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
        AltYeast> >::VALUE[4][4][4] =
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
        AscidianMitochondrial> >::VALUE[4][4][4] =
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
        AltFlatwormMitochondrial> >::VALUE[4][4][4] =
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
        Blepherisma> >::VALUE[4][4][4] =
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
        ChlorophyceanMitochondrial> >::VALUE[4][4][4] =
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
        TrematodeMitochondrial> >::VALUE[4][4][4] =
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
        ScenedesmusMitochondrial> >::VALUE[4][4][4] =
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
        ThraustochytriumMitochondrial> >::VALUE[4][4][4] =
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
        PterobranchiaMitochondrial> >::VALUE[4][4][4] =
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
        Gracilibacteria> >::VALUE[4][4][4] =
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
