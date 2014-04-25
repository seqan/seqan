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
// Code for Dna(5) to AminoAcid Translation
// ==========================================================================


#ifndef EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_
#define EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Enum TranslationFrames
// -----------------------------------------------------------------------

/*!
 * @enum TranslationFrames
 * @headerfile seqan/translation.h
 * @brief Class Enum with frames for @link translate @endlink()
 *
 * @signature enum class TranslationFrames : uint8_t { ... };
 *
 * @val TranslationFrames TranslationFrames::SingleFrame = 0;
 * @brief Translate the sequence(s) "as is", n input sequences result in n output sequences.
 *
 * @val TranslationFrames TranslationFrames::WithReverseComplement = 1;
 * @brief Translate the sequence(s) as well as their reverse complements (n -> * 2n).
 *
 * @val TranslationFrames TranslationFrames::WithFrameShifts = 2;
 * @brief Translate the sequence(s) as well as their shifted frames (n -> 3n).
 *
 * @val TranslationFrames TranslationFrames::SixFrame = 3;
 * @brief Equals (TranslationFrames::WithReverseComplement | TranslationFrames::WithFrameShifts); shifted frames of original and reverse complement are
 *        translated (n -> 6n).
 */

enum class TranslationFrames : uint8_t
{
    SingleFrame             = 0,
    WithReverseComplement   = 1,
    WithFrameShifts         = 2,
    SixFrame                = 3
};

// -----------------------------------------------------------------------
// Tag Frames_ (internal)
// -----------------------------------------------------------------------

template <uint8_t num>
struct Frames_
{};

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction ReverseComplement_
// -----------------------------------------------------------------------

// type of the reverse complement of a string, also works with infixes and
// ModifiedStrings

template <typename T>
struct ReverseComplement_
{
    typedef ModifiedString<
        ModifiedString<T, ModView<
                           FunctorComplement<
                               typename Value<T>::Type > > >, ModReverse> Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function _ord()
// --------------------------------------------------------------------------

// returns ordValue of a DNA(5) or RNA(5) character
// for everything else (e.g. char) the character is converted to Dna5 first
// TODO(C++11): when ordValue is constexpr, this should be, too.

template <typename T>
inline typename ValueSize<T>::Type
_ord(T const & c)
{
    return ordValue(Dna5(c));
}

inline ValueSize<Dna>::Type
_ord(Dna const & c)
{
    return ordValue(c);
}

inline ValueSize<Dna5>::Type
_ord(Dna5 const & c)
{
    return ordValue(c);
}

inline ValueSize<Rna>::Type
_ord(Rna const & c)
{
    return ordValue(c);
}

inline ValueSize<Rna5>::Type
_ord(Rna5 const & c)
{
    return ordValue(c);
}

// --------------------------------------------------------------------------
// Function _translateTriplet()
// --------------------------------------------------------------------------


template <typename TOrd, GeneticCodeSpec codeSpec>
constexpr AminoAcid
_translateTriplet(TOrd const & c1,
                  TOrd const & c2,
                  TOrd const & c3,
                  GeneticCode<codeSpec> const & /**/)
{
    return (( c1 > 3 ) || ( c2 > 3 ) || ( c3 > 3 ) )
            ? 'X'
            : TranslateTableDnaToAminoAcid_<
                GeneticCode<codeSpec> >::VALUE[c1][c2][c3];
}

// --------------------------------------------------------------------------
// Function _translateString()
// --------------------------------------------------------------------------

template <typename TOutString, typename TInString, GeneticCodeSpec codeSpec>
inline void
_translateString(TOutString & target,
                 TInString const & source,
                 GeneticCode<codeSpec> const & /**/)
{
    SEQAN_ASSERT_EQ(length(source)/3, length(target));
    typedef typename Position<TInString>::Type TPos;

    for (TPos i = 0; i+2 < length(source); i+=3)
    {
        target[i/3] = _translateTriplet(_ord(value(source, i  )),
                                        _ord(value(source, i+1)),
                                        _ord(value(source, i+2)),
                                        GeneticCode<codeSpec>());
    }
}

template <typename TOutString, typename TSpec, typename TInString, GeneticCodeSpec codeSpec>
inline void
_translateString(Segment<TOutString, TSpec> && target,
                 TInString const & source,
                 GeneticCode<codeSpec> const & /**/)
{
    SEQAN_ASSERT_EQ(length(source)/3, length(target));
    typedef typename Position<TInString>::Type TPos;

    for (TPos i = 0; i+2 < length(source); i+=3)
    {
        target[i/3] = _translateTriplet(_ord(value(source, i  )),
                                        _ord(value(source, i+1)),
                                        _ord(value(source, i+2)),
                                        GeneticCode<codeSpec>());
    }
}

// --------------------------------------------------------------------------
// Function _translateImplLoop()
// --------------------------------------------------------------------------

// single frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec codeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<codeSpec> const & /**/,
                   Frames_<1u> const & /**/)
{
    typedef GeneticCode<codeSpec> TCode;
    _translateString(target[i], source[i], TCode());
}

// with reverse complement
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec codeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<codeSpec> const & /**/,
                   Frames_<2u> const & /**/)
{
    typedef typename Value<StringSet<TInString, TSpec3> const>::Type TVal;
    typedef typename ReverseComplement_<TVal>::Type TRevComp;
    typedef GeneticCode<codeSpec> TCode;

    if (i % 2)
    {
        TRevComp revComp(value(source, i/2));
        _translateString(target[i], revComp, TCode());
    }
    else
    {
        _translateString(target[i], source[i/2], TCode());
    }
}

// three frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec codeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<codeSpec> const & /**/,
                   Frames_<3u> const & /**/)
{
    typedef GeneticCode<codeSpec> TCode;
    _translateString(target[i], suffix(source[i/3], i % 3), TCode());
}

// six frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          GeneticCodeSpec codeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   unsigned const i,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<codeSpec> const & /**/,
                   Frames_<6u> const & /**/)
{
    typedef typename Prefix<
        typename Value<
            StringSet<TInString, TSpec3> const>::Type>::Type TVal;
    typedef typename ReverseComplement_<TVal>::Type TRevComp;
    typedef GeneticCode<codeSpec> TCode;

    if (i % 6 > 2)
    {
        TRevComp revComp(prefix(value(source, i/6),
                                length(value(source,i/6)) - (i % 3)));
        _translateString(target[i], revComp, TCode());
    }
    else
    {
        _translateString(target[i], suffix(source[i/6], i % 3), TCode());
    }

}

// --------------------------------------------------------------------------
// Function _translateImplLoopOMPWrapper()
// --------------------------------------------------------------------------


template <typename TSource, typename TTarget, uint8_t frames,
          GeneticCodeSpec codeSpec>
inline void
_translateImplLoopOMPWrapper(TTarget & target,
                             TSource const & source,
                             GeneticCode<codeSpec> const & /**/,
                             Frames_<frames> const & /**/,
                             Parallel const & /**/)
{
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    for (unsigned i = 0; i < length(target); ++i)
        _translateImplLoop(target, i, source, GeneticCode<codeSpec>(),
                           Frames_<frames>());
}

template <typename TSource, typename TTarget, uint8_t frames,
          GeneticCodeSpec codeSpec>
inline void
_translateImplLoopOMPWrapper(TTarget & target,
                             TSource const & source,
                             GeneticCode<codeSpec> const & /**/,
                             Frames_<frames> const & /**/,
                             Serial const & /**/)
{
    for (unsigned i = 0; i < length(target); ++i)
        _translateImplLoop(target, i, source, GeneticCode<codeSpec>(),
                           Frames_<frames>());
}

// --------------------------------------------------------------------------
// Function _translateImpl()
// --------------------------------------------------------------------------

// general case
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TParallelism, GeneticCodeSpec codeSpec, unsigned char n>
inline void
_translateImpl(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
               StringSet<TInString, TSpec3> const & source,
               GeneticCode<codeSpec> const & /**/,
               Frames_<n> const & /**/,
               TParallelism const & /**/)
{
    typedef typename Position<StringSet<TInString, TSpec3> >::Type TPos;

    resize(target, length(source) * n, Exact());
    for (TPos i = 0; i < length(target); ++i)
    {
        // current dnastring's length / 3 (3DNA -> 1 AA)
        TPos len = length(source[i/n]) / 3;
        // shorten for shifted frames
        if (( n > 2 ) && ( length(source[i/n]) % 3 ) < ( i%3 ))
            --len;
        resize(target[i], len, Exact());
    }

    _translateImplLoopOMPWrapper(target, source, GeneticCode<codeSpec>(),
                                 Frames_<n>(),
                                 TParallelism());
}

// ConcatDirect specialization
template <typename TSpec1, typename TSpec3, typename TInString,
          typename TParallelism, GeneticCodeSpec codeSpec, unsigned char n>
inline void
_translateImpl(StringSet<String<AminoAcid,
                                TSpec1>, Owner<ConcatDirect<> > > & target,
               StringSet<TInString, TSpec3> const & source,
               GeneticCode<codeSpec> const & /**/,
               Frames_<n> const & /**/,
               TParallelism const & /**/)
{
    typedef typename Position<StringSet<TInString, TSpec3> >::Type TPos;

    resize(target.limits, length(source) * n + 1, Exact());
    target.limits[0] = 0;
    for (TPos i = 0; i+1 < length(target.limits); ++i)
    {
        // current dnastring's length / 3 (3DNA -> 1 AA)
        TPos len = length(source[i/n]) / 3;
        // shorten for shifted frames
        if (( n > 2 ) && ( length(source[i/n]) % 3 ) < ( i%3 ))
            --len;
        target.limits[i+1] = target.limits[i] + len;
    }

    resize(target.concat, back(target.limits), Exact());

    _translateImplLoopOMPWrapper(target, source, GeneticCode<codeSpec>(),
                                 Frames_<n>(),
                                 TParallelism());
}

// --------------------------------------------------------------------------
// Function _translateInputWrap()
// --------------------------------------------------------------------------

// stringset to stringset
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TParallelism, GeneticCodeSpec codeSpec, unsigned char n>
inline int
_translateInputWrap(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                    StringSet<TInString, TSpec3> const & source,
                    GeneticCode<codeSpec> const & /**/,
                    Frames_<n> const & /**/,
                    TParallelism const & /**/)
{
    _translateImpl(target, source, GeneticCode<codeSpec>(), Frames_<n>(),
                   TParallelism());
    return 0;
}

// single string to stringset conversion
template <typename TSpec1, typename TSpec2, typename TInString,
          typename TParallelism, GeneticCodeSpec codeSpec, unsigned char n>
inline int
_translateInputWrap(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                    TInString const & source,
                    GeneticCode<codeSpec> const & /**/,
                    Frames_<n> const & /**/,
                    TParallelism const & /**/)
{
    StringSet<TInString, Dependent<> > set;
    appendValue(set, source);
    _translateImpl(target, set, GeneticCode<codeSpec>(), Frames_<n>(),
                   TParallelism());

    return 0;
}


// bail out because multiple frames don't fit in one string
template <typename TSpec1, typename TInString, typename TParallelism,
          GeneticCodeSpec codeSpec, unsigned char n>
inline int
_translateInputWrap(String<AminoAcid, TSpec1> & /**/,
                    TInString const & /**/,
                    GeneticCode<codeSpec> const & /**/,
                    Frames_<n> const & /**/,
                    TParallelism const & /**/)
{
    return -1;
}
// single string to single string conversion
template <typename TSpec1, typename TInString, typename TParallelism,
          GeneticCodeSpec codeSpec>
inline int
_translateInputWrap(String<AminoAcid, TSpec1> & target,
                    TInString const & source,
                    GeneticCode<codeSpec> const & /**/,
                    Frames_<1> const & /**/,
                    TParallelism const & /**/)
{
    resize(target, length(source)/3, Exact());
    _translateString(target, source, GeneticCode<codeSpec>());

    return 0;
}

// --------------------------------------------------------------------------
// Function translate()
// --------------------------------------------------------------------------

/*!
 * @fn translate
 * @headerfile seqan/translation.h
 * @brief translate sequences of Dna or Rna into amino acid alphabet, optionally with frames
 * @signature int translate(target, source[, frames][, geneticCode][, TParallelism])
 * @signature int translate(target, source[, frames][, geneticCodeSpec][, TParallelism])
 *
 * @param[out]      target      The amino acid sequence(s).  @link StringSet @endlink of @link AminoAcid @endlink
 *                              or @link String @endlink of @link AminoAcid @endlink if source is a single string
 *                              and frames is <tt>SingleFrame</tt>.
 * @param[in]       source      Source sequences @link String @endlink or @link StringSet @endlink.
 *                              If the value type is not Dna, Dna5, Rna, Rna5 then it is converted
 *                              to Dna5.
 * @param[in]       frame       The @link TranslationFrames @endlink, defaults to SingleFrame.
 * @param[in]       geneticCode The @link GeneticCode @endlink to use, defaults to GeneticCode<GeneticCodeSpec::Canonical> (this is compile-time constant)
 * @param[in]       geneticCodeSpec The @link GeneticCodeSpec @endlink to use, this is a run-time parameter ("run-time argument")
 * @param[in]       TParallelism Whether to use SMP or not, see @link ParallelismTags @endlink .
 *
 * @return int 0 on success, and -1 on incompatible parameters (e.g. multiple frames but target type not StringSet).
 *
 * If OpenMP is supported by platform and TParallelism is not specified as
 * "Serial", translation will be parallelized. The only exception is when doing
 * single-frame translation of a single string -- which is never parallelized.
 *
 * The translation process is fastest when using ConcatDirect-StringSets for
 * both input and output StringSets and when not having to convert the alphabet
 * of the source (see below).
 *
 * @section Example
 *
 * @code{.cpp}
 * StringSet<Dna5String> dnaSeqs;
 *
 * // do something that fills up dnaSeqs, e.g. read from file or assign
 *
 * StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > aaSeqs;
 *
 * translate(aaSeqs, dnaSeqs, TranslationFrames::SixFrame);
 *
 * // do something with the aaSeqs
 * @endcode
 *
 * @see TranslationFrames
 * @see GeneticCode
 */

template <typename TTarget, typename TSource, typename TParallelism,
          GeneticCodeSpec codeSpec>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCode<codeSpec> const & /**/,
          TParallelism const & /**/)
{
    static_assert(std::is_same<TParallelism, Parallel>::value or
                  std::is_same<TParallelism, Serial>::value,
                  "TParallelism must either be Parallel() or Serial().");

    typedef GeneticCode<codeSpec> TCode;
    switch (frames)
    {
    case TranslationFrames::SingleFrame:
        return _translateInputWrap(target, source, TCode(), Frames_<1>(),
                                   TParallelism());
    case TranslationFrames::WithReverseComplement:
        return _translateInputWrap(target, source, TCode(), Frames_<2>(),
                                   TParallelism());
    case TranslationFrames::WithFrameShifts:
        return _translateInputWrap(target, source, TCode(), Frames_<3>(),
                                   TParallelism());
    case TranslationFrames::SixFrame:
        return _translateInputWrap(target, source, TCode(), Frames_<6>(),
                                   TParallelism());
    default:
        return -1;
    }
    return 0;
}

template <typename TTarget, typename TSource, GeneticCodeSpec codeSpec>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCode<codeSpec> const & /**/)
{
    return translate(target, source, frames, GeneticCode<codeSpec>(),
                     Parallel());
}

template <typename TTarget, typename TSource, GeneticCodeSpec codeSpec>
inline int
translate(TTarget & target,
          TSource const & source,
          GeneticCode<codeSpec> const & /**/)
{
    return translate(target, source, TranslationFrames::SingleFrame,
                     GeneticCode<codeSpec>(), Parallel());
}

template <typename TTarget, typename TSource>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames)
{
    return translate(target, source, frames, GeneticCode<>(), Parallel());
}

template <typename TTarget, typename TSource>
inline int
translate(TTarget & target,
          TSource const & source)
{
    return translate(target, source, TranslationFrames::SingleFrame,
                     GeneticCode<>(), Parallel());
}

template <typename TTarget, typename TSource, typename TParallelism>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCodeSpec const & geneticCode,
          TParallelism const & /**/)
{
    static_assert(std::is_same<TParallelism, Parallel>::value or
                  std::is_same<TParallelism, Serial>::value,
                  "TParallelism must either be Parallel() or Serial().");
    switch(geneticCode)
    {
        case GeneticCodeSpec::Canonical:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::Canonical>(),
                             TParallelism());
        case GeneticCodeSpec::VertMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::VertMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::YeastMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::YeastMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::MoldMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::MoldMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::InvertMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::InvertMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::Ciliate:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::Ciliate>(),
                             TParallelism());
        case GeneticCodeSpec::FlatwormMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::FlatwormMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::Euplotid:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::Euplotid>(),
                             TParallelism());
        case GeneticCodeSpec::Prokaryote:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::Prokaryote>(),
                             TParallelism());
        case GeneticCodeSpec::AltYeast:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::AltYeast>(),
                             TParallelism());
        case GeneticCodeSpec::AscidianMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::AscidianMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::AltFlatwormMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::AltFlatwormMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::Blepharisma:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::Blepharisma>(),
                             TParallelism());
        case GeneticCodeSpec::ChlorophyceanMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::ChlorophyceanMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::TrematodeMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::TrematodeMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::ScenedesmusMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::ScenedesmusMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::ThraustochytriumMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::ThraustochytriumMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::PterobranchiaMitochondrial:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::PterobranchiaMitochondrial>(),
                             TParallelism());
        case GeneticCodeSpec::Gracilibacteria:
            return translate(target, source, frames,
                             GeneticCode<GeneticCodeSpec::Gracilibacteria>(),
                             TParallelism());
    }

    std::cerr << "Invalid genetic code translation table selected."
              << std::endl;
    return -1;
}

template <typename TTarget, typename TSource>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          GeneticCodeSpec const & geneticCode)
{
    return translate(target, source, frames, geneticCode, Parallel());
}

template <typename TTarget, typename TSource, typename TParallelism>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationFrames const frames,
          TParallelism const & /**/)
{
    return translate(target, source, frames, GeneticCode<>(), TParallelism());
}

}

#endif // EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_
