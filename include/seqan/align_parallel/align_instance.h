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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INSTANCE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INSTANCE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
    
template <typename TScore, typename TDPTraits, typename TThreadContext>
struct DPTaskTraits
{
    // ----------------------------------------------------------------------------
    // Typedefs.
    // ----------------------------------------------------------------------------
    
    // DPTrait type forwarding.
    using TScoreValue       = typename Value<TScore>::Type;
    using TAlgorithmType    = typename TDPTraits::TAlgorithmType;
    using TTracebackType    = typename TDPTraits::TTracebackType;
    using TGapsType         = typename TDPTraits::TGapsType;
    
    // Types needed for the buffer.
    using TDPCell           = DPCell_<TScoreValue, TGapsType>;
    using TBufferValue      = Pair<TDPCell, typename TraceBitMap_<>::Type>;
    using TBuffer           = String<TBufferValue>;
    using TBlockBuffer      = DPTileBuffer<TBuffer>;
    
    // DPProfile type.
    using TDPProfile        = DPProfile_<TAlgorithmType, TGapsType, TTracebackType, Parallel>;
    
    // Parallel type.
    //    using TIsVectorized     = typename IsVectorized<typename TThreadContext::TExecutionPolicy>::Type;
    
    // ----------------------------------------------------------------------------
    // Static member functions.
    // ----------------------------------------------------------------------------
    
    // TODO(rrahn) Implement me!
    template <...>
    static auto createTaskGraph()
    {
        return createGraph(taskContext, tls, execPolicy);
    }
    
    // ----------------------------------------------------------------------------
    // Function createBlockBuffer()
    // ----------------------------------------------------------------------------
    
    template <typename TSeqH, typename TSeqV>
    static auto createBlockBuffer(TSeqH const & seqH, TSeqV const & seqV, TScore const & score, size_t const blockSize)
    {
        TBlockBuffer buffer;
        resize(buffer.horizontalBuffer, length(seqH), Exact());
        resize(buffer.verticalBuffer, length(seqV), Exact());
        
        TBufferValue tmp;
        tmp.i2 = _doComputeScore(tmp.i1, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), score, RecursionDirectionZero(), TDPProfile());
        for (auto itH = begin(buffer.horizontalBuffer, Standard()); itH != end(buffer.horizontalBuffer, Standard()); ++itH)
        {
            resize(*itH, length(front(seqH)), Exact());
            for (auto it = begin(*itH, Standard()); it != end(*itH, Standard()); ++it)
            {
                it->i2 = _doComputeScore(it->i1, TDPCell(), tmp.i1, TDPCell(), Nothing(), Nothing(), score, RecursionDirectionHorizontal(), TDPProfile());
                tmp.i1 = it->i1;
            }
        }
        tmp.i1 = decltype(tmp.i1){};
        tmp.i2 = _doComputeScore(tmp.i1, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), score, RecursionDirectionZero(), TDPProfile());
        
        for (auto itV = begin(buffer.verticalBuffer, Standard()); itV != end(buffer.verticalBuffer, Standard()); ++itV)
        {
            resize(*itV, length(front(seqV)) + 1, Exact());
            auto it = begin(*itV, Standard());
            it->i2 = tmp.i2;
            it->i1 = tmp.i1;
            ++it;
            for (; it != end(*itV, Standard()); ++it)
            {
                it->i2 = _doComputeScore(it->i1, TDPCell(), TDPCell(), tmp.i1, Nothing(), Nothing(), score, RecursionDirectionVertical(), TDPProfile());
                tmp.i1 = it->i1;
                tmp.i2 = it->i2;  // TODO(rrahn): Move out of loop.
            }
        }
        return buffer;
    }
    
    // ----------------------------------------------------------------------------
    // Function createBlocks()
    // ----------------------------------------------------------------------------
    
    template <typename TSeq>
    static auto createBlocks(TSeq const & seq, size_t const blockSize)
    {
        using TIter = typename Iterator<typename Infix<TSeq const>::Type, Standard>::Type;
        String<Range<TIter>> blocks;
        resize(blocks, (length(seq) + blockSize - 1) / blockSize, Exact());
        
        for (unsigned id = 0; id < length(blocks); ++id)
            blocks[id] = toRange(infix(seq, id * blockSize, _min(length(seq),(id + 1) * blockSize)));
        return blocks;
    }
    
};

template <typename TSeqH,
          typename TSeqV,
          typename TScore,
          typename TDPConfig,
          typename TDelegate,
          typename TThreadContext,
          typename TDPTaskTraits = DPTaskTraits<TScore, typename Traits<TDPConfig>::Type, TThreadContext>>
class AlignmentInstance
{
public:

    // ----------------------------------------------------------------------------
    // Member Variables.
    // ----------------------------------------------------------------------------

    TSeqH     const* mSeqHPtr     = nullptr;
    TSeqV     const* mSeqVPtr     = nullptr;
    TScore    const* mScorePtr    = nullptr;
    TDPConfig const* mDPConfigPtr = nullptr;

    TDelegate &     mDelegate;
    TThreadContext& mThreadContext;
    
    // ----------------------------------------------------------------------------
    // Constructors.
    // ----------------------------------------------------------------------------
    
    AlignmentInstance(TSeqH const & seqH,
                      TSeqV const & seqV,
                      TScore const & score,
                      TDPConfig const & config,
                      TDelegate & delegate,
                      TThreadContext & threadContext) :
        mSeqHPtr(&seqH),
        mSeqVPtr(&seqV),
        mScorePtr(&score),
        mDPConfigPtr(&config),
        mDelegate(delegate),
        mThreadContext(threadContext)
    {}

    // ----------------------------------------------------------------------------
    // Member Functions.
    // ----------------------------------------------------------------------------

    inline void
    operator()(uint16_t const mInstanceId)
    {
        // Initialize the strings.
        auto seqHBlocks = TDPTaskTraits::createBlocks(*mSeqHPtr, mDPConfigPtr->blockSize);
        auto seqVBlocks = TDPTaskTraits::createBlocks(*mSeqVPtr, mDPConfigPtr->blockSize);
        
        auto buffer = TDPTaskTraits::createBlockBuffer(seqHBlocks, seqVBlocks, *mScorePtr, mDPConfigPtr->blockSize);
        
        

    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
    
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INSTANCE_H_
