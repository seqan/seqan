// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de
// ==========================================================================
// Code for dynamically selectable protein matrixes
// ==========================================================================

#ifndef SEQAN_SCORE_MATRIX_DYN_H_
#define SEQAN_SCORE_MATRIX_DYN_H_

namespace seqan
{

// ============================================================================
// Classes, Tags, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// helpers
// ----------------------------------------------------------------------------

struct ScoreSpecSelectable {};

namespace impl
{
namespace score
{
// ATTENTION: always keep in sync with the enum below!
using MatrixTags = TagList<ScoreSpecBlosum30,
                   TagList<ScoreSpecBlosum45,
                   TagList<ScoreSpecBlosum62,
                   TagList<ScoreSpecBlosum80,
                   TagList<ScoreSpecPam40,
                   TagList<ScoreSpecPam120,
                   TagList<ScoreSpecPam200,
                   TagList<ScoreSpecPam250,
                   TagList<ScoreSpecVtml200
                   > > > > > > > > >;
} // namespace impl::score
} // namespace impl

// ----------------------------------------------------------------------------
// Enum AminoAcidScoreMatrixID
// ----------------------------------------------------------------------------

/*!
 * @enum AminoAcidScoreMatrixID
 * @brief Enum with aliases for the different Score specializations
 * @signature enum class AminoAcidScoreMatrixID : uint8_t { ... };
 *
 * @headerfile <seqan/score.h>
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::BLOSUM30
 * @brief Blosum30, see also @link Blosum30 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::BLOSUM45
 * @brief Blosum45, see also @link Blosum45 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::BLOSUM62
 * @brief Blosum62, see also @link Blosum62 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::BLOSUM80
 * @brief Blosum80, see also @link Blosum80 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::PAM40
 * @brief Pam40, see also @link Pam40 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::PAM120
 * @brief Pam120, see also @link Pam120 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::PAM200
 * @brief Pam200, see also @link Pam200 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::PAM250
 * @brief Pam250, see also @link Pam250 @endlink
 *
 * @val AminoAcidScoreMatrixID AminoAcidScoreMatrixID::VTML200
 * @brief Vtml200, see also @link Vtml200 @endlink
 */

enum class AminoAcidScoreMatrixID : std::underlying_type_t<decltype(Find<impl::score::MatrixTags, ScoreSpecBlosum30>::VALUE)>
{
    BLOSUM30 = Find<impl::score::MatrixTags, ScoreSpecBlosum30>::VALUE,
    BLOSUM45 = Find<impl::score::MatrixTags, ScoreSpecBlosum45>::VALUE,
    BLOSUM62 = Find<impl::score::MatrixTags, ScoreSpecBlosum62>::VALUE,
    BLOSUM80 = Find<impl::score::MatrixTags, ScoreSpecBlosum80>::VALUE,
    PAM40    = Find<impl::score::MatrixTags, ScoreSpecPam40>::VALUE,
    PAM120   = Find<impl::score::MatrixTags, ScoreSpecPam120>::VALUE,
    PAM200   = Find<impl::score::MatrixTags, ScoreSpecPam200>::VALUE,
    PAM250   = Find<impl::score::MatrixTags, ScoreSpecPam250>::VALUE,
    VTML200  = Find<impl::score::MatrixTags, ScoreSpecVtml200>::VALUE
};

// ----------------------------------------------------------------------------
// Class SelectableAminoAcidMatrix
// ----------------------------------------------------------------------------

/*!
 * @typedef SelectableAminoAcidMatrix
 * @headerfile <seqan/score.h>
 * @brief An AminoAcid score matrix that can be "specialized" at run-time
 *
 * @signature using SelectableAminoAcidMatrix =  Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> >;
 */

template <typename TValue>
class Score<TValue, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > :
    public Score<TValue, ScoreMatrix<AminoAcid, Default> >
{
public:

    using TBaseScore = Score<TValue, ScoreMatrix<AminoAcid, Default>>;

    AminoAcidScoreMatrixID _ident = static_cast<AminoAcidScoreMatrixID>(0u);

    explicit Score(TValue _gap_extend = -1) :
        TBaseScore(_gap_extend)
    {}

    Score(TValue _gap_extend, TValue _gap_open) :
        TBaseScore(_gap_extend, _gap_open)
    {}

    Score(Score const &) = default;
    Score(Score &&) = default;

    Score & operator=(Score const &) = default;
    Score & operator=(Score &&) = default;

    ~Score() = default;
};

using SelectableAminoAcidMatrix =  Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> >;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

namespace impl
{
namespace score
{

template <typename TCurTag,
          typename TRunnable>
inline void
matrixTagDispatch(TagList<TCurTag, void> const &,
                  AminoAcidScoreMatrixID const m,
                  TRunnable const & runnable)
{
    using TUndType = std::underlying_type_t<AminoAcidScoreMatrixID>;
    if (Find<impl::score::MatrixTags, TCurTag>::VALUE == static_cast<TUndType>(m))
        runnable(TCurTag());
    else
        SEQAN_FAIL("ERROR: Recursing the ScoreMatrixTags failed, please report this as a BUG!");
}

template <typename TCurTag,
          typename TRestList,
          typename TRunnable>
inline void
matrixTagDispatch(TagList<TCurTag, TRestList> const &,
                  AminoAcidScoreMatrixID const m,
                  TRunnable const & runnable)
{
    using TUndType = std::underlying_type_t<AminoAcidScoreMatrixID>;
    if (Find<impl::score::MatrixTags, TCurTag>::VALUE == static_cast<TUndType>(m))
        runnable(TCurTag());
    else
        matrixTagDispatch(TRestList(), m, runnable);
}

} // namespace impl::score
} // namespace impl

// ----------------------------------------------------------------------------
// Function setScoreMatrixById()
// ----------------------------------------------------------------------------

/*!
 * @fn SelectableAminoAcidMatrix#setScoreMatrixById
 * @headerfile <seqan/score.h>
 * @brief Set the substitution score matrix
 *
 * @signature void setScoreMatrixById(SelectableAminoAcidMatrix & sc, AminoAcidScoreMatrixID const id)
 *
 * @param[in,out] sc    The score object to be modified.
 * @param[in]     id    The ID of the matrix (@link AminoAcidScoreMatrixID @endlink).
 */

template <typename TValue>
inline void
setScoreMatrixById(Score<TValue, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > & sc,
                   AminoAcidScoreMatrixID const m)
{
    sc._ident = m;
    impl::score::matrixTagDispatch(impl::score::MatrixTags(), m, [&sc] (auto const & tag)
    {
        setDefaultScoreMatrix(sc, tag);
    });
}

// ----------------------------------------------------------------------------
// Function getScoreMatrixId()
// ----------------------------------------------------------------------------

/*!
 * @fn SelectableAminoAcidMatrix#getScoreMatrixId
 * @headerfile <seqan/score.h>
 * @brief Get the matrix ID from the dynamic matrix object.
 *
 * @signature AminoAcidScoreMatrixID getScoreMatrixId(SelectableAminoAcidMatrix const & sc)
 *
 * @param[in] sc    The score object to be modified.
 */

template <typename TValue>
inline AminoAcidScoreMatrixID
getScoreMatrixId(Score<TValue, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > const & sc)
{
    return sc._ident;
}

}

#endif // SEQAN_SCORE_MATRIX_DYN_H_
