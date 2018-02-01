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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the context that can be passed to the dp functions in order
// to reuse memory blocks in mutliple calls of the same function.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_DP_CONTEXT_H_
#define INCLUDE_SEQAN_ALIGN_DP_CONTEXT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename T>
struct GetDPScoreMatrix
{};

template <typename T>
struct GetDPTraceMatrix
{};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TScoreValue, typename TTraceValue,
          typename TScoreMatrixHost = String<TScoreValue>,
          typename TTraceMatrixHost = String<TTraceValue> >
struct DPContext
{
    TScoreMatrixHost _scoreMatrix;
    TTraceMatrixHost _traceMatrix;

    DPContext() : _scoreMatrix(), _traceMatrix()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function dpScoreMatrix()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts, typename TScoreMatHost, typename TTraceMatHost>
inline TScoreMatHost &
getDpScoreMatrix(DPContext<TScoreValue, TGapCosts, TScoreMatHost, TTraceMatHost> & dpContext)
{
    return dpContext._scoreMatrix;
}

template <typename TScoreValue, typename TGapCosts, typename TScoreMatHost, typename TTraceMatHost>
inline TScoreMatHost const &
getDpScoreMatrix(DPContext<TScoreValue, TGapCosts, TScoreMatHost, TTraceMatHost> const & dpContext)
{
    return dpContext._scoreMatrix;
}

// ----------------------------------------------------------------------------
// Function dpTraceMatrix()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts, typename TScoreMatHost, typename TTraceMatHost>
inline TTraceMatHost &
getDpTraceMatrix(DPContext<TScoreValue, TGapCosts, TScoreMatHost, TTraceMatHost> & dpContext)
{
    return dpContext._traceMatrix;
}

template <typename TScoreValue, typename TGapCosts, typename TScoreMatHost, typename TTraceMatHost>
inline TTraceMatHost const &
getDpTraceMatrix(DPContext<TScoreValue, TGapCosts, TScoreMatHost, TTraceMatHost> const & dpContext)
{
    return dpContext._traceMatrix;
}

// ----------------------------------------------------------------------------
// Function setDpScoreMatrix()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts, typename TScoreMatHost, typename TTraceMatHost>
inline void
setDpTraceMatrix(DPContext<TScoreValue, TGapCosts, TScoreMatHost, TTraceMatHost> & dpContext,
                 TScoreMatHost const & scoreMatrix)
{
    dpContext._scoreMatrix = scoreMatrix;
}

// ----------------------------------------------------------------------------
// Function setDpTraceMatrix()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapCosts, typename TScoreMatHost, typename TTraceMatHost>
inline void
setDpTraceMatrix(DPContext<TScoreValue, TGapCosts, TScoreMatHost, TTraceMatHost> & dpContext,
                 TTraceMatHost const & traceMatrix)
{
    dpContext._tarceMatrix = traceMatrix;
}

}

#endif // INCLUDE_SEQAN_ALIGN_DP_CONTEXT_H_
