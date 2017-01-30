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
//  DP base task.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_BASE_2_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_BASE_2_H_

namespace seqan
{
namespace impl
{
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TTask>
class DPTaskBase
{
public:
    std::array<TTask*, 2>  mSuccessor;
    size_t                 mCol{0};
    size_t                 mRow{0};
    std::atomic<unsigned>  mRefCount{0};
    bool                   mLastTileH{false};
    bool                   mLastTileV{false};

    // ============================================================================
    // Member functions.

    // TODO(rrahn): Make global functions.
    inline void
    setRefCount(unsigned const n)
    {
        refCount.store(n, std::memory_order_relaxed);
    }

    inline unsigned
    decrementRefCount()
    {
        return --refCount;
    }

    inline unsigned
    incrementRefCount()
    {
        return ++refCount;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TTask>
auto column(TTask const & task)
{
    return task.mCol;
}

template <typename TTask>
auto row(TTask const & task)
{
    return task.mRow;
}

template <typename TTask>
bool inLastColumn(TTask const & task)
{
    return task.mLastTileH;
}

template <typename TTask>
auto inLastRow(TTask const & task)
{
    return task.mLastTileV;
}

template <typename TTask>
inline void
executeScalar(TTask & task)
{
    using TTaskTraits = typename TTask::TTraits;
    using TBuffer = typename TTaskTraits::TTileBuffer;
    using TDPScoutState = DPScoutState_<DPTiled<TBuffer>>;

    using TResultState = typename 
    using TDPTraits = typename Traits<typename TTaskTraits::TDPConfig>::Type;

    // Indirection:
    auto & tls = local(threadStorage(context(task)));
    auto & dpCache = cache(tls);

    auto & buffer = tileBuffer(context(task));

    // Something I really want to get refined.
    TDPScoutState scoutState(buffer.horizontalBuffer[column(task)],
                             buffer.verticalBuffer[row(task)]);  // Task local
    DPScout_<TDPCell, TDPScoutState> scout(scouteState);
    // We rather get the DPScout in here.

    // DEBUG: Remove!
//        auto bufHBegin = _taskContext.getTileBuffer().horizontalBuffer[_col];
//        auto bufVBegin = _taskContext.getTileBuffer().verticalBuffer[_row];

    computeTile(dpCache, scout,
                seqBlocksH(context(task))[column(task)],
                seqBlocksV(context(task))[row(task)],
                config(context(task)));

    if (isTrackingEnabled<typename TDPTraits::TAlgorithm>(inLastColumn(task), inLastRow(task)))
    {
        update(intermediate(tls, instance(context(task))),
               TResultState{maxScore(scout), maxHostPosition(scout)},
               column(task),
               row(task));
    }

    
    // DEBUG: Remove!
//        _taskContext.getDebugBuffer().matrix[_col][_row] = {bufHBegin, bufVBegin,
//                                                            _taskContext.getTileBuffer().horizontalBuffer[_col],
//                                                            _taskContext.getTileBuffer().verticalBuffer[_row],
//                                                            _col, _row, false};
    // TODO(rrahn): Add traceback later.
//    if (IsTracebackEnabled_<typename TTaskConfig::TDPConfig>::VALUE)
//    {
//        swap(getDpTraceMatrix(dpContext), pTls.mLocalTraceStore.localScalarTraceMatrix());
//        _taskContext.getTraceProxy().insert(std::make_pair(_col, _row),
//                                            TTraceProxyValue{&pTls.mLocalTraceStore,
//                                                             {0, static_cast<uint16_t>(length(pTls.mLocalTraceStore.mScalarTraceVec) - 1)},
//                                                              0});
//    }
    #ifdef DP_ALIGN_STATS
        ++serialCounter;
    #endif
}

}  // namespace impl
}  // namespace seqan
#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TASK_BASE_2_H_
