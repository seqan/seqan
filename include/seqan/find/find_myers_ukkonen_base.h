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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_BASE_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_BASE_H

#include <seqan/sequence.h>

#include <seqan/score/score_edit.h>
#include <seqan/find/find_base.h>
#include <seqan/find/find_pattern_base.h>
#include <seqan/find/find_begin.h>

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// MyersUkkonen
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class MyersPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Provides fast approximate searching of one string in another using Myer's fast bit-parallel algorithm with
 *        application of the Ukkonen- trick.
 *
 * @signature template <typename TNeedle[, typename TSpec[, typename TFindBeginPatternSpec]]>
 *            class Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >;
 *
 * @tparam TSpec   Specialization tag.  This is @link ApproximateFinderSearchTypeTags#FindInfix @endlink for
 *                 infix search or @link ApproximateFinderSearchTypeTags#FindPrefix @endlink for prefix search.
 *                 Defaults to @linkApproximateFinderSearchTypeTags#FindInfix @endlink.
 * @tparam TFindBeginPatternSpec
 *               Specialization of @link Pattern @endlink used to find the begin of matches.This must be a finder for
 *               prefix search, e.g. @link DPSearchPattern <tt>DPSearch&lt;TScore, FindPrefix&gt;</tt> @endlink or @link
 *               MyersPattern <tt>Myers&lt;FindPrefix&gt;</tt> @endlink. Specify <tt>void</tt> to suppress prefix
 *               searching. Default: @link DefaultFindBeginPatternSpec @endlink
 * @tparam TNeedle The needle type. Types: String
 *
 * The needle-length must be smaller than the highest number that can be stored in an unsigned int.
 */

template <typename TSpec = FindInfix,
          typename THasState = True,
          typename TFindBeginPatternSpec = typename DefaultFindBeginPatternSpec<EditDistanceScore, THasState>::Type>
struct Myers {};

struct NMatchesNone_;
struct NMatchesN_;
struct NMatchesAll_;

//FindInfix and FindPrefix are defined int find_base.h
template <typename TSpec, typename TFinderCharSetPolicy = NMatchesN_, typename TPatternCharSetPolicy = NMatchesN_>
struct AlignTextBanded; // search query in a parallelogram

// TODO(holtgrew): Really deprecated?
//deprecated shortcuts:

/*!
 * @typedef MyersUkkonen
 * @headerfile <seqan/find.h>
 * @brief Semi-global (query-global, text-local) pattern matching without
 *        findBegin() support.
 *
 * @signature typedef Myers<FindInfix, True, void> MyersUkkonen;
 *
 * @deprecated Use <tt>Myers&lt;FindInfix&gt;</tt> instead.
 */

typedef Myers<FindInfix, True, void> MyersUkkonen;

/*!
 * @typedef MyersUkkonenGlobal
 * @headerfile <seqan/find.h>
 * @brief Global (query-global, text-global) pattern matching without findBegin() support.
 *
 * @signature typedef Myers<FindInfix, True, void> MyersUkkonenGlobal;
 */

typedef Myers<FindPrefix, True, void> MyersUkkonenGlobal;

/*!
 * @typedef MyersUkkonenBanded
 * @headerfile <seqan/find.h>
 * @brief Global (query-global, text-local) pattern matching without findBegin() support.
 *
 * @signature Myers<AlignTextBanded<FindInfix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenBanded;
 */

/*!
 * @typedef MyersUkkonenGlobalBanded
 * @headerfile <seqan/find.h>
 * @brief Global (query-global, text-global) pattern matching without findBegin() support.
 *
 * @signature Myers<AlignTextBanded<FindPrefix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenGlobalBanded;
 */

typedef Myers<AlignTextBanded<FindInfix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenBanded;
typedef Myers<AlignTextBanded<FindPrefix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenGlobalBanded;

template <typename TValue>
struct MyersSmallAlphabet_:
    public Eval<ValueSize<TValue>::VALUE <= 8> {};

//////////////////////////////////////////////////////////////////////////////
// State Data
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
struct MyersSmallState_;

template <typename TNeedle, typename TSpec>
struct MyersLargeState_;

//////////////////////////////////////////////////////////////////////////////
// Pattern Data
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
struct MyersSmallPattern_;

template <typename TNeedle, typename TSpec>
struct MyersLargePattern_;


//____________________________________________________________________________
// bit 0 of the HP bit-vector
// 0 for begin-gap-free haystack
// 1 for global alignments of haystack

template <typename T>
struct MyersUkkonenHP0_ {
    enum { VALUE = 0 };
};

template <>
struct MyersUkkonenHP0_<FindPrefix> {
    enum { VALUE = 1 };
};

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_BASE_H
