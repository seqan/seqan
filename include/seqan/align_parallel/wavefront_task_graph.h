// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_GRAPH_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_GRAPH_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TTask, template <typename> typename TAllocator = std::allocator>
class WavefrontTaskGraph
{
    //-------------------------------------------------------------------------
    // Helper

    template <typename TDimType>
    struct Dimension
    {
        TDimType mDimH;
        TDimType mDimV;

        void setDimensionH(TDimType const dimH)
        {
            mDimH = dimH;
        }

        void setDimensionV(TDimType const dimV)
        {
            mDimV = dimV;
        }
    };

    //-------------------------------------------------------------------------
    // Typedefs

    using TTaskAllocator = TAllocator<TTask>;
    using TPointer       = typename TTaskAllocator<>
    using TDimension     = Dimension<TDimType>;

    //-------------------------------------------------------------------------
    // Private Member Variables

    // we could use a bare array, but this does not make sense.
    // Unless, we cannot use the 
    std::vector<std::vector<TTask>>     _mTaskDag;


    //-------------------------------------------------------------------------
    // Member Functions

    inline auto&
    operator[](TDimension const & dim)
    {
        return
    }
};

template <typename TGraph>
struct TaskGraphTraits;

template <typename ...TArgs>
struct TaskGraphTraits<WavefrontTaskGraph<TArgs...>>
{
    using TDimension = typename WavefrontTaskGraph<TArgs...>::TDimension;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename ...TArgs>
inline auto&
root(WavefrontTaskGraph<TArgs...> & me)
{
    return typename WavefrontTaskGraph<TArgs...>::TAllocator::reference(me._mTaskDag[0][0]);
}

template <typename ...TArgs>
inline auto&
sink(WavefrontTaskGraph<TArgs...> & me)
{
    using WavefrontTaskGraph<TArgs...>::Dimensions;
    return typename WavefrontTaskGraph<TArgs...>::TAllocator::reference(me._mTaskDag[length(me, HORIZONTAL) - 1][length(me, VERTICAL) - 1]);
}


}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_TASK_GRAPH_H_