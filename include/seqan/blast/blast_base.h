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
 * @see BlastProgramSelector
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
 * @brief Unkown type. Used to signify that the type could not be inferred from the file.
 *
 * @val BlastProgram BlastProgram::DYNAMIC
 * @brief This can only be used when defining a @link BlastProgramSelector @endlink
 *
 */
enum class BlastProgram : uint8_t
{
    BLASTN,         //              NUCL VS             NUCL
    BLASTP,         //              PROT VS             PROT
    BLASTX,         // TRANSLATED   NUCL VS             PROT
    TBLASTN,        //              PROT VS TRANSLATED  NUCL
    TBLASTX,        // TRANSLATED   NUCL VS TRANSLATED  NUCL
    UNKNOWN=254,
    DYNAMIC=255
};

/*!
 * @class BlastProgramSelector
 * @brief A datatype that can act as either a @link BlastProgram @endlink or as an constexpr integral constant
 * thereof.
 *
 * @signature template <BlastProgram p> struct BlastProgramSelector { ... };
 * @headerfile <seqan/blast.h>
 *
 * This is a proxy datatype that enables compile-time optimizations through constexpressions iff the value
 * is known at compile time. You will rarely need to instantiate objects of this type yourself, but they
 * are used in the @link BlastIOContext @endlink.
 *
 * The interface functions work on regular enum values of @link BlastProgram @endlink, as well, but are
 * gathered here for convenience.
 *
 * Please note that the default value for <tt>BlastProgramSelector<BlastProgram::DYNAMIC></tt> is
 * <tt>BlastProgram::UNKNOWN</tt>, not <tt>BlastProgram::DYNAMIC</tt>.
 *
 * @subsection Example
 *
 * mutable variable:
 * @code{.cpp}
 * BlastProgramSelector<BlastProgram::DYNAMIC> myProgram = BlastProgram::BLASTN;
 * // same as:
 * // BlastProgram myProgram = BlastProgram::BLASTN;
 *
 * SEQAN_ASSERT(myProgram == BlastProgram::BLASTN); // assertion is checked at run-time
 * myProgram = BlastProgram::BLASTP; // works without problems
 * @endcode
 *
 * compile time integral constant:
 * @code{.cpp}
 * BlastProgramSelector<BlastProgram::BLASTN> myProgram;
 * static_assert(myProgram == BlastProgram::BLASTN, ""); // assertion is checked at compile time
 * myProgram = BlastProgram::BLASTP; // would fail, because value is fixed (and different)
 * @endcode
 */
template <BlastProgram _p>
struct BlastProgramSelector
{
    constexpr operator BlastProgram() const
    {
        return _p;
    }

    BlastProgramSelector operator=(BlastProgram const p)
    {
        if (p != _p)
            SEQAN_FAIL("ERROR: Tried to set blastProgram on context, but was already defined at compile time (and set to a "
                       "different value)!");
        return *this;
    }
};

template <>
struct BlastProgramSelector<BlastProgram::DYNAMIC>
{
    BlastProgram _runtimeValue = BlastProgram::UNKNOWN;

    operator BlastProgram() const
    {
        return _runtimeValue;
    }

    BlastProgramSelector operator=(BlastProgram const p)
    {
        _runtimeValue = p;
        return *this;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// qHasRevComp()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramSelector#qHasRevComp
 * @brief Whether the reverse complement of the <b>query</b> sequence(s) is searched with the given BlastProgram
 * @signature constexpr bool qHasRevComp(BlastProgram const p);
 *
 * @return false for @link BlastProgram::BLASTP @endlink and @link BlastProgram::TBLASTN @endlink
 * @return true otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
qHasRevComp(BlastProgram const p)
{
    return ((p!=BlastProgram::BLASTP) && (p!=BlastProgram::TBLASTN));
}

// ----------------------------------------------------------------------------
// qIsTranslated()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramSelector#qIsTranslated
 * @brief Whether the <b>query</b> sequence is translated in the given BlastProgram
 * @signature constexpr bool qIsTranslated(BlastProgram const p);
 *
 * @return true for @link BlastProgram::BLASTX @endlink and @link BlastProgram::TBLASTX @endlink
 * @return false otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
qIsTranslated(BlastProgram const p)
{
    return ((p==BlastProgram::BLASTX) || (p==BlastProgram::TBLASTX));
}

// ----------------------------------------------------------------------------
// qNumFrames()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramSelector#qNumFrames
 * @brief The number of frames per <b>query</b> sequence implied by the given BlastProgram
 * @signature constexpr uint8_t qNumFrames(BlastProgram const p);
 *
 * @return 6 for @link BlastProgram::BLASTX @endlink and @link BlastProgram::TBLASTX @endlink
 * @return 2 for @link BlastProgram::BLASTN @endlink
 * @return 1 for @link BlastProgram::BLASTP @endlink and @link BlastProgram::TBLASTN @endlink
 * @headerfile <seqan/blast.h>
 */

constexpr uint8_t
qNumFrames(BlastProgram p)
{
    return (qIsTranslated(p) ? 6 : (qHasRevComp(p) ? 2 : 1));
}

// ----------------------------------------------------------------------------
// sHasRevComp()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramSelector#sHasRevComp
 * @brief Whether the reverse complement of the <b>subject</b> sequence(s) is searched with the given BlastProgram
 * @signature constexpr bool sHasRevComp(BlastProgram const p);
 *
 * @return true for @link BlastProgram::TBLASTX @endlink and @link BlastProgram::TBLASTN @endlink
 * @return false otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
sHasRevComp(BlastProgram const p)
{
    return ((p==BlastProgram::TBLASTX) || (p==BlastProgram::TBLASTN));
}

// ----------------------------------------------------------------------------
// sIsTranslated()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramSelector#sIsTranslated
 * @brief Whether the <b>subject</b> sequence is translated in the given BlastProgram
 * @signature constexpr bool sIsTranslated(BlastProgram const p);
 *
 * @return true for @link BlastProgram::TBLASTX @endlink and @link BlastProgram::TBLASTN @endlink
 * @return false otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr bool
sIsTranslated(BlastProgram const p)
{
    return ((p==BlastProgram::TBLASTX) || (p==BlastProgram::TBLASTN));
}

// ----------------------------------------------------------------------------
// sNumFrames()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastProgramSelector#sNumFrames
 * @brief The number of frames per <b>query</b> sequence implied by the given BlastProgram
 * @signature constexpr uint8_t sNumFrames(BlastProgram const p);
 *
 * @return 6 for @link BlastProgram::TBLASTX @endlink and @link BlastProgram::TBLASTN @endlink
 * @return 1 otherwise
 * @headerfile <seqan/blast.h>
 */

constexpr uint8_t
sNumFrames(BlastProgram p)
{
    return (sIsTranslated(p) ? 6 : (sHasRevComp(p) ? 2 : 1));
}

// ----------------------------------------------------------------------------
// _programTagToString()
// ----------------------------------------------------------------------------

template <typename TVoidSpec = void>
struct BlastProgramStrings_
{
    static constexpr char const * const VALUE [6] =
    {
        "BLASTN",
        "BLASTP",
        "BLASTX",
        "TBLASTN",
        "TBLASTX",
        "UNKOWN BLAST PROGRAM"
    };
};

template <typename TVoidSpec>
constexpr char const * const BlastProgramStrings_<TVoidSpec>::VALUE[6];

constexpr const char *
_programTagToString(BlastProgram const _p)
{
    return (uint8_t(_p) < 5) ? BlastProgramStrings_<>::VALUE[uint8_t(_p)] : BlastProgramStrings_<>::VALUE[5];
}

// ----------------------------------------------------------------------------
// _programStringToTag()
// ----------------------------------------------------------------------------

template <typename TString>
inline BlastProgram
_programStringToTag(TString const & str)
{
    for (uint8_t i = 0; i < 5; ++i)
        if (str == BlastProgramStrings_<>::VALUE[i])
            return BlastProgram(i);

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
                      BlastProgramSelector<p> const & selector)
{
     if (qIsTranslated(selector))
         _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, True(), True());
     else if (qHasRevComp(selector))
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
                      BlastProgramSelector<p> const & selector)
{
     if (sIsTranslated(selector))
         _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, True(), True());
     else if (sHasRevComp(selector))
         _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, True(), False());
     else
        _untranslatePositions(effectiveStart, effectiveEnd, frameShift, length, False(), False());
}

} // namespace seqan

#endif // header guard
