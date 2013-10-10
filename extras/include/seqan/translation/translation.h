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

#ifdef SEQAN_CXX11_STANDARD
#define RVREF   &&
#else
#define RVREF
#endif

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Enum TranslationOptions
// -----------------------------------------------------------------------

/*!
 * @enum TranslationOptions
 * @brief Enum with options for @link translate @endlink()
 *
 * @signature enum TranslationOptions { ... };
 *
 * @headerfile seqan/translation.h
 *
 * @var TranslationOptions SINGLE_FRAME
 * @brief translate the sequence(s) "as is", n input sequences result in n
 * output sequences
 *
 * @var TranslationOptions WITH_REV_COMP
 * @brief translate the sequence(s) as well as their reverse complements (n ->
 * 2n)
 *
 * @var TranslationOptions WITH_FRAME_SHIFT
 * @brief translate the sequence(s) as well as their shifted frames (n -> 3n)
 *
 * @var TranslationOptions SIX_FRAME
 * @brief equals (WITH_REV_COMP | WITH_FRAME_SHIFT); shifted frames of original
 * and reverse complement are translated (n -> 6n)
 */

enum TranslationOptions
{
    SINGLE_FRAME     = 0,
    WITH_REV_COMP    = 1,
    WITH_FRAME_SHIFT = 2,
    SIX_FRAME        = 3
};

// -----------------------------------------------------------------------
// Tag Frames_ (internal)
// -----------------------------------------------------------------------

template <unsigned char num>
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


template <typename TOrd, typename TCodeSpec>
inline AminoAcid
_translateTriplet(TOrd const & c1,
                  TOrd const & c2,
                  TOrd const & c3,
                  GeneticCode<TCodeSpec> const & /**/)
{
    if (( c1 > 3 ) || ( c2 > 3 ) || ( c3 > 3 ) )
        return 'X';
    return TranslateTableDnaToAminoAcid_<
               GeneticCode<
                   TCodeSpec> >::VALUE[c1][c2][c3];
}

// --------------------------------------------------------------------------
// Function _translateString()
// --------------------------------------------------------------------------

template <typename TOutString, typename TInString, typename TCodeSpec>
inline void
_translateString(TOutString & target,
                 TInString const & source,
                 GeneticCode<TCodeSpec> const & /**/)
{
    SEQAN_ASSERT_EQ(length(source)/3, length(target));
    typedef typename Position<TInString>::Type TPos;

    for (TPos i = 0; i+2 < length(source); i+=3)
    {
        target[i/3] = _translateTriplet(_ord(value(source, i  )),
                                        _ord(value(source, i+1)),
                                        _ord(value(source, i+2)),
                                        GeneticCode<TCodeSpec>());
    }
}

// the above doesn't work for elements of ConcatDirect, because they are infixes
// in the æther thereby RVALUES. Here we use a Hack to get around it:
// * if SeqAn is compiled with CPP11 Support, we use an R-Value-Reference
//   in the definition, which works, because the implementation only operates
//   on the host of the infix and doesn't change the properties of the infix
//   itself (these changes would be lost)
// * if there is no CPP11 Support then we simply copy the infix, which is
//   suboptimal, but not too expensive
// * NOTE that if the infix had a move constructor, we could transparently
//   just use the above function, but with an && and remove the specialization
//   below. This is due to reference collapsing in C++11 ( 'A& &&' is translated
//   to 'A&', NOT 'A&&'),
template <typename TOutString, typename TSpec, typename TInString, typename TCodeSpec>
inline void
_translateString(Segment<TOutString, TSpec> RVREF target,
                 TInString const & source,
                 GeneticCode<TCodeSpec> const & /**/)
{
    SEQAN_ASSERT_EQ(length(source)/3, length(target));
    typedef typename Position<TInString>::Type TPos;

    for (TPos i = 0; i+2 < length(source); i+=3)
    {
        target[i/3] = _translateTriplet(_ord(value(source, i  )),
                                        _ord(value(source, i+1)),
                                        _ord(value(source, i+2)),
                                        GeneticCode<TCodeSpec>());
    }
}

// --------------------------------------------------------------------------
// Function _translateImplLoop()
// --------------------------------------------------------------------------

// single frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TCodeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<TCodeSpec> const & /**/,
                   Frames_<1u> const & /**/)
{
    typedef GeneticCode<TCodeSpec> TCode;
    #ifndef SEQAN_TRANSLATION_NO_PARALLEL
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    #endif
    for (unsigned i = 0; i < length(target); ++i)
        _translateString(target[i], source[i], TCode());
}

// with reverse complement
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TCodeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<TCodeSpec> const & /**/,
                   Frames_<2u> const & /**/)
{
    typedef typename Value<StringSet<TInString, TSpec3> const>::Type TVal;
    typedef typename ReverseComplement_<TVal>::Type TRevComp;
    typedef GeneticCode<TCodeSpec> TCode;

    #ifndef SEQAN_TRANSLATION_NO_PARALLEL
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    #endif
    for (unsigned i = 0; i < length(target); ++i)
    {
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
}

// three frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TCodeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<TCodeSpec> const & /**/,
                   Frames_<3u> const & /**/)
{
    typedef GeneticCode<TCodeSpec> TCode;
    #ifndef SEQAN_TRANSLATION_NO_PARALLEL
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    #endif
    for (unsigned i = 0; i < length(target); ++i)
        _translateString(target[i], suffix(source[i/3], i % 3), TCode());
}

// six frame
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TCodeSpec>
inline void
_translateImplLoop(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                   StringSet<TInString, TSpec3> const & source,
                   GeneticCode<TCodeSpec> const & /**/,
                   Frames_<6u> const & /**/)
{
    typedef typename Prefix<
        typename Value<
            StringSet<TInString, TSpec3> const>::Type>::Type TVal;
    typedef typename ReverseComplement_<TVal>::Type TRevComp;
    typedef GeneticCode<TCodeSpec> TCode;
    #ifndef SEQAN_TRANSLATION_NO_PARALLEL
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic))
    #endif
    for (unsigned i = 0; i < length(target); ++i)
    {
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
}

// --------------------------------------------------------------------------
// Function _translateImpl()
// --------------------------------------------------------------------------

// general case
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TCodeSpec, unsigned char n>
inline void
_translateImpl(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
               StringSet<TInString, TSpec3> const & source,
               GeneticCode<TCodeSpec> const & /**/,
               Frames_<n> const & /**/)
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

    _translateImplLoop(target, source, GeneticCode<TCodeSpec>(), Frames_<n>());
}

// ConcatDirect specialization
template <typename TSpec1, typename TSpec3, typename TInString,
          typename TCodeSpec, unsigned char n>
inline void
_translateImpl(StringSet<String<AminoAcid,
                                TSpec1>, Owner<ConcatDirect<> > > & target,
               StringSet<TInString, TSpec3> const & source,
               GeneticCode<TCodeSpec> const & /**/,
               Frames_<n> const & /**/)
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

    _translateImplLoop(target, source, GeneticCode<TCodeSpec>(), Frames_<n>());
}

// --------------------------------------------------------------------------
// Function _translateInputWrap()
// --------------------------------------------------------------------------

// stringset to stringset
template <typename TSpec1, typename TSpec2, typename TSpec3, typename TInString,
          typename TCodeSpec, unsigned char n>
inline int
_translateInputWrap(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                    StringSet<TInString, TSpec3> const & source,
                    GeneticCode<TCodeSpec> const & /**/,
                    Frames_<n> const & /**/)
{
    _translateImpl(target, source, GeneticCode<TCodeSpec>(), Frames_<n>());
    return 0;
}

// single string to stringset conversion
template <typename TSpec1, typename TSpec2, typename TInString,
          typename TCodeSpec, unsigned char n>
inline int
_translateInputWrap(StringSet<String<AminoAcid, TSpec1>, TSpec2> & target,
                    TInString const & source,
                    GeneticCode<TCodeSpec> const & /**/,
                    Frames_<n> const & /**/)
{
    StringSet<TInString, Dependent<> > set;
    appendValue(set, source);
    _translateImpl(target, set, GeneticCode<TCodeSpec>(), Frames_<n>());

    return 0;
}


// bail out because multiple frames don't fit in one string
template <typename TSpec1, typename TInString, typename TCodeSpec,
          unsigned char n>
inline int
_translateInputWrap(String<AminoAcid, TSpec1> & /**/,
                    TInString const & /**/,
                    GeneticCode<TCodeSpec> const & /**/,
                    Frames_<n> const & /**/)
{
    return -1;
}
// single string to single string conversion
template <typename TSpec1, typename TInString, typename TCodeSpec>
inline int
_translateInputWrap(String<AminoAcid, TSpec1> & target,
                    TInString const & source,
                    GeneticCode<TCodeSpec> const & /**/,
                    Frames_<1> const & /**/)
{
    resize(target, length(source)/3, Exact());
    _translateString(target, source, GeneticCode<TCodeSpec>());

    return 0;
}

// --------------------------------------------------------------------------
// Function translate()
// --------------------------------------------------------------------------

/*!
 * @fn translate
 * @brief translate sequences of Dna or Rna into amino acid alphabet, optionally with frames
 * @signature int translate(target, source[, options][, codeSpec])
 *
 * @param[out]      target      the amino acid sequence(s)
 * [StringSet&lt;String&lt;AminoAcid&gt;, *&gt; or String&lt;AminoAcid*gt; iff source is a single
 * string and options == SINGLE_FRAME]
 * @param[in]       source      source sequences [String<Dna|Dna5|Rna|Rna5> or
 * StringSet thereof, other input will be converted to Dna5]
 * @param[in]       options     the @link TranslationOptions @endlink, defaults to SINGLE_FRAME
 * @param[in]       codeSpec    the @link GeneticCode @endlink to use, defaults to canonical
 *
 * @return 0 on success, and -1 on incompatible parameters (e.g. multiple frames
 * but target type not StringSet)
 *
 * @headerfile seqan/translation.h
 * @section Remarks
 *
 * This call uses OpenMP internally, if supported by platform. If you want to
 * disable this, #define SEQAN_TRANSLATION_NO_PARALLEL.
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
 * translate(aaSeqs, dnaSeqs, SIX_FRAME);
 *
 * // do something with the aaSeqs
 * @endcode
 *
 * @see TranslationOptions
 * @see GeneticCode
 */

template <typename TTarget, typename TSource, typename TCodeSpec>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationOptions const options,
          GeneticCode<TCodeSpec> const & /**/)
{
    typedef GeneticCode<TCodeSpec> TCode;
    switch (options)
    {
    case SINGLE_FRAME:
        return _translateInputWrap(target, source, TCode(), Frames_<1>());
    case WITH_REV_COMP:
        return _translateInputWrap(target, source, TCode(), Frames_<2>());
    case WITH_FRAME_SHIFT:
        return _translateInputWrap(target, source, TCode(), Frames_<3>());
    case SIX_FRAME:
        return _translateInputWrap(target, source, TCode(), Frames_<6>());
    default:
        return -1;
    }
    return 0;
}


template <typename TTarget, typename TSource>
inline int
translate(TTarget & target,
          TSource const & source,
          TranslationOptions const options)
{
    return translate(target, source, options, GeneticCode<>());
}

template <typename TTarget, typename TSource, typename TCodeSpec>
inline int
translate(TTarget & target,
          TSource const & source,
          GeneticCode<TCodeSpec> const & /**/)
{
    return translate(target, source, SINGLE_FRAME, GeneticCode<TCodeSpec>());
}

template <typename TTarget, typename TSource>
inline int
translate(TTarget & target,
          TSource const & source)
{
    return translate(target, source, SINGLE_FRAME, GeneticCode<>());
}



}

#endif // EXTRAS_INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_
