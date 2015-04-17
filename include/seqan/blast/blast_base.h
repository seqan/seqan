// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Module for handling NCBI Blast I/O and E-Value computation
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_BLAST_BASE_H_
#define SEQAN_EXTRAS_BLAST_BLAST_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Enum BlastProgram
// ----------------------------------------------------------------------------

/*!
 * @enum BlastProgram
 * @brief Enum with BLAST program spec
 * @signature enum class BlastProgram : uint8_t { ... };
 *
 * @headerfile <seqan/blast.h>
 * @see BlastProgramTag
 *
 * @val BlastProgram BlastProgram::BLASTN
 * @brief Nucleotide Query VS Nucleotide Subject
 *
 * @val BlastProgram BlastProgram::BLASTP
 * @brief Protein Query VS Protein Subject
 *
 * @val BlastProgram BlastProgram::BLASTX
 * @brief translated Nucleotide Query VS Protein Subject
 *
 * @val BlastProgram BlastProgram::TBLASTN
 * @brief Protein Query VS translated Nucleotide Subject
 *
 * @val BlastProgram BlastProgram::TBLASTX
 * @brief translated Nucleotide Query VS translated Nucleotide Subject
 *
 * @val BlastProgram BlastProgram::UNKNOWN
 * @brief Unkown type. Used internally and to signify that the type could not be inferred from the file.
 */
enum class BlastProgram : uint8_t
{
    BLASTN,         //              NUCL VS             NUCL
    BLASTP,         //              PROT VS             PROT
    BLASTX,         // TRANSLATED   NUCL VS             PROT
    TBLASTN,        //              PROT VS TRANSLATED  NUCL
    TBLASTX,        // TRANSLATED   NUCL VS TRANSLATED  NUCL
    UNKNOWN=255
};

/*!
 * @defgroup BlastProgramTag BlastProgram interface
 * @brief Integral constants of @link BlastProgram @endlink
 *
 * @tag BlastProgramTag#BlastProgramTagBlastN
 * @brief Tag for @link BlastProgram::BLASTN @endlink
 * @headerfile <seqan/blast.h>
 *
 * @tag BlastProgramTag#BlastProgramTagBlastP
 * @brief Tag for @link BlastProgram::BLASTP @endlink
 * @headerfile <seqan/blast.h>
 *
 * @tag BlastProgramTag#BlastProgramTagBlastX
 * @brief Tag for @link BlastProgram::BLASTX @endlink
 * @headerfile <seqan/blast.h>
 *
 * @tag BlastProgramTag#BlastProgramTagTBlastN
 * @brief Tag for @link BlastProgram::TBLASTN @endlink
 * @headerfile <seqan/blast.h>
 *
 * @tag BlastProgramTag#BlastProgramTagTBlastX
 * @brief Tag for @link BlastProgram::TBLASTX @endlink
 * @headerfile <seqan/blast.h>
 *
 * @tag BlastProgramTag#BlastProgramTagUnknown
 * @brief Tag for @link BlastProgram::UNKNOWN @endlink
 * @headerfile <seqan/blast.h>
 */

template <BlastProgram p>
using BlastProgramTag = std::integral_constant<BlastProgram, p>;

typedef BlastProgramTag<BlastProgram::BLASTN>  BlastProgramTagBlastN;
typedef BlastProgramTag<BlastProgram::BLASTP>  BlastProgramTagBlastP;
typedef BlastProgramTag<BlastProgram::BLASTX>  BlastProgramTagBlastX;
typedef BlastProgramTag<BlastProgram::TBLASTN> BlastProgramTagTBlastN;
typedef BlastProgramTag<BlastProgram::TBLASTX> BlastProgramTagTBlastX;
typedef BlastProgramTag<BlastProgram::UNKNOWN> BlastProgramTagUnknown;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// getBlastProgram()
// ----------------------------------------------------------------------------

/*
 * Given both an enum value and a tag, return the tag's value unless it is UNKNOWN
 * in which case the enum's value is returned. Can be used in implementations that
 * want to cover compile-time constants and run-time settings both.
 */

// this is always evaluated at compile-time
template <BlastProgram p>
constexpr BlastProgram
getBlastProgram(BlastProgram const, BlastProgramTag<p> const &)
{
    return p;
}

// this will likely be evaluated at run-time
constexpr BlastProgram
getBlastProgram(BlastProgram const p, BlastProgramTag<BlastProgram::UNKNOWN> const &)
{
    return p;
}

// ----------------------------------------------------------------------------
// qHasRevComp()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramTag#qHasRevComp
 * @brief Whether the reverse complement of the <b>query</b> sequence(s) is searched with the given BlastProgram
 * @signature constexpr bool qHasRevComp(BlastProgram const p);
 * template <BlastProgram p>
 * constexpr bool qHasRevComp(BlastProgramTag<p> const &);
 * template <BlastProgram p>
 * constexpr bool qHasRevComp(BlastProgram const _p, BlastProgramTag<p> const &);
 *
 * @section Remarks
 *
 * This function can be used to check whether the defined BlastProgam implies the existence/creation/... of a
 * reverse complement sequence of the query sequence.
 *
 * If both a tag and an enum are specified, the tag is given preference, unless it is set to
 * @link BlastProgram::UNKNOWN @endlink.
 *
 * @return false for @link BlastProgram::BLASTP @endlink and @link BlastProgram::TBLASTN @endlink
 * @return true otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
qHasRevComp(BlastProgram const p)
{
    return ((p==BlastProgram::BLASTP) || (p==BlastProgram::TBLASTN))
            ? false
            : true;
}

template <BlastProgram p>
constexpr bool
qHasRevComp(BlastProgramTag<p> const &)
{
    return qHasRevComp(p);
}

template <BlastProgram p>
constexpr bool
qHasRevComp(BlastProgram const _p, BlastProgramTag<p> const &)
{
    return qHasRevComp(getBlastProgram(_p, BlastProgramTag<p>()));
}

// ----------------------------------------------------------------------------
// qIsTranslated()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramTag#qIsTranslated
 * @brief Whether the <b>query</b> sequence is translated in the given BlastProgram
 * @signature constexpr bool qIsTranslated(BlastProgram const p);
 * template <BlastProgram p>
 * constexpr bool qIsTranslated(BlastProgramTag<p> const &);
 * template <BlastProgram p>
 * constexpr bool qIsTranslated(BlastProgram const _p, BlastProgramTag<p> const &);
 *
 * @section Remarks
 *
 * This function can be used to check whether the defined BlastProgam implies the
 * nucleotide to aminoacid translation of the query sequence.
 *
 * If both a tag and an enum are specified, the tag is given preference, unless it is set to
 * @link BlastProgram::UNKNOWN @endlink.
 *
 * @return true for @link BlastProgram::BLASTX @endlink and @link BlastProgram::TBLASTX @endlink
 * @return false otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
qIsTranslated(BlastProgram const p)
{
    return ((p==BlastProgram::BLASTX) || (p==BlastProgram::TBLASTX))
            ? true
            : false;
}

template <BlastProgram p>
constexpr bool
qIsTranslated(BlastProgramTag<p> const &)
{
    return qIsTranslated(p);
}

template <BlastProgram p>
constexpr bool
qIsTranslated(BlastProgram const _p, BlastProgramTag<p> const &)
{
    return qIsTranslated(getBlastProgram(_p, BlastProgramTag<p>()));
}

// ----------------------------------------------------------------------------
// qNumFrames()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramTag#qNumFrames
 * @brief The number of frames per <b>query</b> sequence implied by the given BlastProgram
 * @signature constexpr uint8_t qNumFrames(BlastProgram const p);
 * template <BlastProgram p>
 * constexpr uint8_t qNumFrames(BlastProgramTag<p> const &);
 * template <BlastProgram p>
 * constexpr uint8_t qNumFrames(BlastProgram const _p, BlastProgramTag<p> const &);
 *
 * @section Remarks
 *
 * If both a tag and an enum are specified, the tag is given preference, unless it is set to
 * @link BlastProgram::UNKNOWN @endlink.
 *
 * @return 6 for @link BlastProgram::BLASTX @endlink and @link BlastProgram::TBLASTX @endlink
 * @return 2 for @link BlastProgram::BLASTN @endlink
 * @return 1 for @link BlastProgram::BLASTP @endlink and @link BlastProgram::TBLASTN @endlink
 * @headerfile <seqan/blast.h>
 */

constexpr uint8_t
qNumFrames(BlastProgram p)
{
    return (qIsTranslated(p)
            ? 6
            : (qHasRevComp(p)
                ? 2
                : 1));
}

template <BlastProgram p>
constexpr bool
qNumFrames(BlastProgramTag<p> const &)
{
    return qNumFrames(p);
}

template <BlastProgram p>
constexpr bool
qNumFrames(BlastProgram const _p, BlastProgramTag<p> const &)
{
    return qNumFrames(getBlastProgram(_p, BlastProgramTag<p>()));
}

// ----------------------------------------------------------------------------
// sHasRevComp()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramTag#sHasRevComp
 * @brief Whether the reverse complement of the <b>subject</b> sequence(s) is searched with the given BlastProgram
 * @signature constexpr bool sHasRevComp(BlastProgram const p);
 * template <BlastProgram p>
 * constexpr bool sHasRevComp(BlastProgramTag<p> const &);
 * template <BlastProgram p>
 * constexpr bool sHasRevComp(BlastProgram const _p, BlastProgramTag<p> const &);
 *
 * @section Remarks
 *
 * This function can be used to check whether the defined BlastProgam implies the existence/creation/... of a
 * reverse complement sequence of the subject sequence.
 *
 * If both a tag and an enum are specified, the tag is given preference, unless it is set to
 * @link BlastProgram::UNKNOWN @endlink.
 *
 * @return true for @link BlastProgram::TBLASTX @endlink and @link BlastProgram::TBLASTN @endlink
 * @return false otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
sHasRevComp(BlastProgram const p)
{
    return ((p==BlastProgram::TBLASTX) ||
            (p==BlastProgram::TBLASTN))
            ? true
            : false;
}

template <BlastProgram p>
constexpr bool
sHasRevComp(BlastProgramTag<p> const &)
{
    return sHasRevComp(p);
}

template <BlastProgram p>
constexpr bool
sHasRevComp(BlastProgram const _p, BlastProgramTag<p> const &)
{
    return sHasRevComp(getBlastProgram(_p, BlastProgramTag<p>()));
}

// ----------------------------------------------------------------------------
// sIsTranslated()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramTag#sIsTranslated
 * @brief Whether the <b>subject</b> sequence is translated in the given BlastProgram
 * @signature constexpr bool sIsTranslated(BlastProgram const p);
 * template <BlastProgram p>
 * constexpr bool sIsTranslated(BlastProgramTag<p> const &);
 * template <BlastProgram p>
 * constexpr bool sIsTranslated(BlastProgram const _p, BlastProgramTag<p> const &);
 *
 * @section Remarks
 *
 * This function can be used to check whether the defined BlastProgam implies the
 * nucleotide to aminoacid translation of the subject sequence.
 *
 * If both a tag and an enum are specified, the tag is given preference, unless it is set to
 * @link BlastProgram::UNKNOWN @endlink.
 *
 * @return true for @link BlastProgram::TBLASTX @endlink and @link BlastProgram::TBLASTN @endlink
 * @return false otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
sIsTranslated(BlastProgram const p)
{
    return ((p==BlastProgram::TBLASTX) ||
            (p==BlastProgram::TBLASTN))
            ? true
            : false;
}

template <BlastProgram p>
constexpr bool
sIsTranslated(BlastProgramTag<p> const &)
{
    return sIsTranslated(p);
}

template <BlastProgram p>
constexpr bool
sIsTranslated(BlastProgram const _p, BlastProgramTag<p> const &)
{
    return sIsTranslated(getBlastProgram(_p, BlastProgramTag<p>()));
}

// ----------------------------------------------------------------------------
// sNumFrames()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramTag#sNumFrames
 * @brief The number of frames per <b>query</b> sequence implied by the given BlastProgram
 * @signature constexpr uint8_t sNumFrames(BlastProgram const p);
 * template <BlastProgram p>
 * constexpr uint8_t sNumFrames(BlastProgramTag<p> const &);
 * template <BlastProgram p>
 * constexpr uint8_t sNumFrames(BlastProgram const _p, BlastProgramTag<p> const &);
 *
 * @section Remarks
 *
 * If both a tag and an enum are specified, the tag is given preference, unless it is set to
 * @link BlastProgram::UNKNOWN @endlink.
 *
 * @return 6 for @link BlastProgram::TBLASTX @endlink and @link BlastProgram::TBLASTN @endlink
 * @return 1 otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr uint8_t
sNumFrames(BlastProgram p)
{
    return (sIsTranslated(p)
            ? 6
            : (sHasRevComp(p)
                ? 2
                : 1));
}

template <BlastProgram p>
constexpr bool
sNumFrames(BlastProgramTag<p> const &)
{
    return sNumFrames(p);
}

template <BlastProgram p>
constexpr bool
sNumFrames(BlastProgram const _p, BlastProgramTag<p> const &)
{
    return sNumFrames(getBlastProgram(_p, BlastProgramTag<p>()));
}

// ----------------------------------------------------------------------------
// _programTagToString()
// ----------------------------------------------------------------------------

template <BlastProgram p>
constexpr const char *
_programTagToString()
{
    return "UNKOWN BLAST PROGRAM";
}

template <>
constexpr const char *
_programTagToString<BlastProgram::BLASTN>()
{
    return "BLASTN";
}

template <>
constexpr const char *
_programTagToString<BlastProgram::BLASTP>()
{
    return "BLASTP";
}

template <>
constexpr const char *
_programTagToString<BlastProgram::BLASTX>()
{
    return "BLASTX";
}

template <>
constexpr const char *
_programTagToString<BlastProgram::TBLASTN>()
{
    return "TBLASTN";
}

template <>
constexpr const char *
_programTagToString<BlastProgram::TBLASTX>()
{
    return "TBLASTX";
}

// if not compile time fixed, than we can operate with const char * and need std::string instead
inline std::string const
_programTagToString(BlastProgram const _p)
{
    switch(_p)
    {
        case BlastProgram::BLASTN:
            return std::string(_programTagToString<BlastProgram::BLASTN>());
        case BlastProgram::BLASTP:
            return std::string(_programTagToString<BlastProgram::BLASTP>());
        case BlastProgram::BLASTX:
            return std::string(_programTagToString<BlastProgram::BLASTX>());
        case BlastProgram::TBLASTN:
            return std::string(_programTagToString<BlastProgram::TBLASTN>());
        case BlastProgram::TBLASTX:
            return std::string(_programTagToString<BlastProgram::TBLASTX>());
        case BlastProgram::UNKNOWN:
            break;
    }
    return std::string(_programTagToString<BlastProgram::UNKNOWN>());
}

// if known at compile-time, deduce at compile-time
template <BlastProgram p>
constexpr const char *
_programTagToString(BlastProgram const, BlastProgramTag<p> const &)
{
    return _programTagToString<p>();
}

// if SPEC == unkown, deduce from run-time parameter
inline std::string const
_programTagToString(BlastProgram const _p, BlastProgramTag<BlastProgram::UNKNOWN> const &)
{
    return _programTagToString(_p);
}

// ----------------------------------------------------------------------------
// _programStringToTag()
// ----------------------------------------------------------------------------

template <typename TString>
inline BlastProgram
_programStringToTag(TString const & str)
{
    if (str == _programTagToString<BlastProgram::BLASTN>())
        return BlastProgram::BLASTN;
    else if (str == _programTagToString<BlastProgram::BLASTP>())
        return BlastProgram::BLASTP;
    else if (str == _programTagToString<BlastProgram::BLASTX>())
        return BlastProgram::BLASTX;
    else if (str == _programTagToString<BlastProgram::TBLASTN>())
        return BlastProgram::TBLASTN;
    else if (str == _programTagToString<BlastProgram::TBLASTX>())
        return BlastProgram::TBLASTX;

    return BlastProgram::UNKNOWN;
}

// ----------------------------------------------------------------------------
// Function _untranslatePositions() -- retransform positions
// ----------------------------------------------------------------------------

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & /**/,
                      signed char const /**/,
                      TPos const & /**/,
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
                      TPos const & length,
                      True const & /*hasReverseComplement*/,
                      False const & /*hasFrames*/)
{
    if (frameShift > 0)
    {
        // BLAST is 1-indexed, but end positions are "on" instead of behind
        // so only the begin positions need adapting
        ++effectiveStart;
    } else
    {
        // reverse strand coordinates have to be transformed
        effectiveStart = length - effectiveStart;
        effectiveEnd = length - effectiveEnd + 1;
        // end is incremented instead of start
    }
}

template <typename TPos>
inline void
_untranslatePositions(TPos & effectiveStart,
                      TPos & effectiveEnd,
                      signed char const frameShift,
                      TPos const & length,
                      True const & /*hasReverseComplement*/,
                      True const & /*hasFrames*/)
{
    // correct for codon translation and frameshift
    // subtract 1 because frameshift is 1-indexed
    effectiveStart = effectiveStart * 3 + std::abs(frameShift) - 1;
    effectiveEnd = effectiveEnd * 3 + std::abs(frameShift) - 1;

    _untranslatePositions(effectiveStart, effectiveEnd, frameShift,
                          length, True(), False());
}

template <typename TPos, BlastProgram p>
inline void
_untranslateQPositions(TPos & effectiveStart,
                      TPos & effectiveEnd,
                      int8_t const frameShift,
                      TPos const & length,
                      BlastProgram const _p,
                      BlastProgramTag<p> const &)
{
     if (qIsTranslated(_p, BlastProgramTag<p>()))
         _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, True(), True());
     else if (qHasRevComp(_p, BlastProgramTag<p>()))
         _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, True(), False());
     else
        _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, False(), False());
}

template <typename TPos, BlastProgram p>
inline void
_untranslateSPositions(TPos & effectiveStart,
                      TPos & effectiveEnd,
                      int8_t const frameShift,
                      TPos const & length,
                      BlastProgram const _p,
                      BlastProgramTag<p> const &)
{
     if (sIsTranslated(_p, BlastProgramTag<p>()))
         _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, True(), True());
     else if (sHasRevComp(_p, BlastProgramTag<p>()))
         _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, True(), False());
     else
        _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, False(), False());
}

} // namespace seqan

#endif // header guard
