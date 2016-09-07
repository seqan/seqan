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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DAG_TASK_BASE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DAG_TASK_BASE_H_

namespace seqan
{

namespace debug
{

std::mutex m;

template <typename TStream, typename TString>
void printBuffer(TStream & stream, TString const & str)
{
    for (auto& val : str)
    {
        stream << val.i1 << "\t";
    }
    stream << "\n";
}

template <typename TTasks,
          typename TScoreValue, typename TSimdVec,
          typename TBufferH,
          typename TBufferV,
          typename TSimdBufferH,
          typename TSimdBufferV,
          typename TSimdDPContext,
          typename TDPConfig>
inline void
compareToScalar(TTasks const & pTasks,
                impl::dp::parallel::DPLocalStorage<TScoreValue, TSimdVec> & pTls,
                TBufferH const & pBufferH,
                TBufferV const & pBufferV,
                TSimdBufferH const & pSimdBufferH,
                TSimdBufferV const & pSimdBufferV,
                TSimdDPContext & pSimdDPContext,
                TDPConfig const & /*unused*/)
{
    using TTask = typename std::remove_pointer<typename Value<TTasks>::Type>::type;
    using TDPLocalStorage = impl::dp::parallel::DPLocalStorage<TScoreValue, TSimdVec>;
    using TStateThreadContext = impl::dp::parallel::StateThreadContext<TTask const, TDPLocalStorage>;
    using TBuffer = typename Value<TBufferH>::Type;
    using TDPCell = typename Value<typename Value<TBuffer>::Type, 1>::Type;
    using TDPScoutState = DPScoutState_<DPTiled<TBuffer, TStateThreadContext> >;


    String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.

    auto& taskContext = pTasks[0]->_taskContext;
    unsigned simdLane = 0;
    for (auto& task : pTasks)
    {
        TBuffer bufH = pBufferH[task->_col];
        TBuffer bufV = pBufferV[task->_row];
        TDPLocalStorage localStore;

        DPContext<TDPCell, typename TraceBitMap_<>::Type> dpContext;
        TDPScoutState scoutState(bufH, bufV, TStateThreadContext(*task, localStore));
        _computeAlignment(dpContext, traceSegments, scoutState,
                          taskContext.getSeqH()[task->_col], taskContext.getSeqV()[task->_row],
                          taskContext.getScore(), taskContext.getBand(), TDPConfig(),
                          false, false);

        // Check horizontal buffer!
        if (length(taskContext.getSeqV()[task->_row]) == length(taskContext.getSeqV()[0]))
        {
            for (unsigned i = 0; i < length(taskContext.getSeqH()[task->_col]); ++i)
            {
                if (bufH[i].i1._score != pSimdBufferH[i].i1._score[simdLane])
                {
                    for (auto& t : pTasks)
                    {
                        std::cout << "(" << t->_col << ", " << t->_row << ") ";
                    }
                    std::cout << "\n";
                }
                SEQAN_ASSERT_EQ_MSG(bufH[i].i1._score, pSimdBufferH[i].i1._score[simdLane], "i = %d, col = %d, row = %d, lc = %d, l0 = %d", i, task->_col, task->_row, length(taskContext.getSeqH()[task->_col]), length(taskContext.getSeqH()[0]));
                SEQAN_ASSERT_EQ(bufH[i].i1._horizontalScore, pSimdBufferH[i].i1._horizontalScore[simdLane]);
                SEQAN_ASSERT_EQ(bufH[i].i1._verticalScore, pSimdBufferH[i].i1._verticalScore[simdLane]);
                SEQAN_ASSERT_EQ(bufH[i].i2, pSimdBufferH[i].i2[simdLane]);
            }
        }

        // Check vertical buffer!
        if (length(taskContext.getSeqH()[task->_col]) == length(taskContext.getSeqH()[0]))
        {
            for (unsigned i = 0; i < length(taskContext.getSeqV()[task->_row]) + 1; ++i)
            {
                if (bufV[i].i1._score != pSimdBufferV[i].i1._score[simdLane])
                {
                    for (auto& t : pTasks)
                    {
                        std::cout << "(" << t->_col << ", " << t->_row << ") ";
                    }
                    std::cout << "\n";
                }
                SEQAN_ASSERT_EQ_MSG(bufV[i].i1._score, pSimdBufferV[i].i1._score[simdLane], "i = %d, col = %d, row = %d, lc = %d, l0 = %d", i, task->_col, task->_row, length(taskContext.getSeqV()[task->_row]), length(taskContext.getSeqV()[0]));
                SEQAN_ASSERT_EQ(bufV[i].i1._horizontalScore, pSimdBufferV[i].i1._horizontalScore[simdLane]);
                SEQAN_ASSERT_EQ(bufV[i].i1._verticalScore, pSimdBufferV[i].i1._verticalScore[simdLane]);
                SEQAN_ASSERT_EQ(bufV[i].i2, pSimdBufferV[i].i2[simdLane]);
            }
        }
        // Now we want to compare the trace matrix.
//        auto& traceMatScalar = getDpTraceMatrix(dpContext);
//        auto& traceMatSimd = getDpTraceMatrix(pSimdDPContext);
//        auto matSize = (length(taskContext.getSeqH()[task->_col]) + 1) * (length(taskContext.getSeqV()[task->_row]) + 1);
//        for (unsigned i = 0; i < matSize; ++i)
//        {
//            SEQAN_ASSERT_EQ(static_cast<uint8_t>(traceMatScalar[i]), static_cast<uint8_t>(traceMatSimd[i][simdLane]));
//        }

        ++simdLane;
    }
}

    // Call in simd block!
    //        debug::compareToScalar(tasks, pTls,
    //                               _taskContext.getTileBuffer().horizontalBuffer, _taskContext.getTileBuffer().verticalBuffer,
    //                               simdBufferH, simdBufferV, dpContext,
    //                               typename TTaskConfig::TDPConfig());

}

// ============================================================================
// Forwards
// ============================================================================

template <typename TTask>
struct IsDPTask;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// We require a fixed order of passed arguments.
// 1. parameter: seqH
// 2. parameter: seqV
// 3. parameter: score
// 4. parameter: band
// 5. parameter: TileBuffer
// 6. parameter: TraceProxy
// 7. parameter: TDPConfig
// 8. parameter: TDPScoutState
template <typename... TArgs>
struct DPTaskContext
{
    using TData         = std::tuple<TArgs...>;
    using TSize         = typename Size<typename std::decay<std::tuple_element_t<0, TData> >::type>::Type;

    using TTraceProxy_  = typename std::remove_pointer<std::tuple_element_t<5, TData> >::type;
    using TLTraceStore_ = typename TTraceProxy_::TLocalTraceStore;
    using TSimdVec      = typename TLTraceStore_::TSimdTraceValue;

    using TDPConfig     = typename std::decay<std::tuple_element_t<6, TData> >::type;
    using TBlockBuffer  = typename std::decay<std::tuple_element_t<7, TData> >::type;

    constexpr static const size_t VECTOR_SIZE = LENGTH<TSimdVec>::VALUE;

    TData                   _data;

    std::mutex              mLock;
    std::mutex              mLockEvent;
    std::condition_variable mReadyEvent;
    bool       mReady = false;

    DPTaskContext() = default;

    template <typename... TSubArgs>
    DPTaskContext(TSubArgs&& ...args)
    {
        static_assert(sizeof...(TArgs) == 8, "Requires 8 arguments.");
        static_assert(sizeof...(TSubArgs) <= sizeof...(TArgs), "Invalid number of parameters.");
        _fill(std::forward_as_tuple(args...), std::make_index_sequence<sizeof...(TSubArgs)>());
    }

    template <typename... TSubArgs, size_t... I>
    void _fill(std::tuple<TSubArgs...> src, std::index_sequence<I...>)
    {
        SEQAN_UNPACK_FUNC(std::get<I>(_data) = std::get<I>(src));
    }

    // Define TileBuffer
    // Define TraceBlock

    inline auto& getSeqH()
    {
        return *std::get<0>(_data);
    }

    inline auto& getSeqV()
    {
        return *std::get<1>(_data);
    }

    inline auto& getScore()
    {
        return *std::get<2>(_data);
    }

    inline auto& getBand()
    {
        return *std::get<3>(_data);
    }

    inline auto& getTileBuffer()
    {
        return *std::get<4>(_data);
    }

    inline auto& getTraceProxy()
    {
        return *std::get<5>(_data);
    }

    inline auto& getDebugBuffer()
    {
        return *std::get<6>(_data);
    }
};

template <typename TTaskConfig, typename TThreadLocalStorage, typename TVecExecPolicy, typename TParExecPolicy>
class DPTaskImpl;

template <typename TTask>
class DPTaskBase;

// Base class.
template <typename TTaskConfig, typename TThreadLocalStorage, typename TVecExecPolicy, typename TParExecPolicy>
class DPTaskBase<DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, TParExecPolicy> >
{
    using TDerivedTask = DPTaskImpl<TTaskConfig, TThreadLocalStorage, TVecExecPolicy, TParExecPolicy>;
    using TGapMode = typename GapTraits<typename TTaskConfig::TDPConfig>::Type;
    using TBuffer = typename TTaskConfig::TBlockBuffer;
    using TSize = typename TTaskConfig::TSize;
public:

    std::array<TDerivedTask*, 2>    successor;
    TTaskConfig&                    _taskContext;
    TDerivedTask&                   _derivedTask;
    typename TTaskConfig::TSize     _col = 0;
    typename TTaskConfig::TSize     _row = 0;
    bool                            _lastHBlock = false;
    bool                            _lastVBlock = false;

    DPTaskBase(TSize pCol, TSize pRow, TTaskConfig& pContext, TDerivedTask& pTask) :
        _taskContext(pContext),
        _derivedTask(pTask),
        _col(pCol),
        _row(pRow)
    {}

    template <typename TScoreValue, typename TSimdVec>
    inline void runScalar(impl::dp::parallel::DPLocalStorage<TScoreValue, TSimdVec> & pTls)
    {
        using TDPLocalStorage = impl::dp::parallel::DPLocalStorage<TScoreValue, TSimdVec>;
        using TTraceProxy = typename std::decay<decltype(_taskContext.getTraceProxy())>::type;
        using TTraceProxyValue = typename TTraceProxy::TTraceMatrixIdentifier;
        using TDPCell = DPCell_<TScoreValue, TGapMode>;
        using TStateThreadContext = impl::dp::parallel::StateThreadContext<DPTaskBase, TDPLocalStorage>;
        using TDPScoutState = DPScoutState_<DPTiled<TBuffer, TStateThreadContext> >;

        DPContext<TDPCell, typename TraceBitMap_<>::Type> dpContext;

        TDPScoutState scoutState(_taskContext.getTileBuffer().horizontalBuffer[_col],
                                 _taskContext.getTileBuffer().verticalBuffer[_row],
                                 TStateThreadContext(*this, pTls));  // Task local

        // DEBUG: Remove!
//        auto bufHBegin = _taskContext.getTileBuffer().horizontalBuffer[_col];
//        auto bufVBegin = _taskContext.getTileBuffer().verticalBuffer[_row];

        String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.
        _computeAlignment(dpContext, traceSegments, scoutState, _taskContext.getSeqH()[_col], _taskContext.getSeqV()[_row],
                          _taskContext.getScore(), _taskContext.getBand(), typename TTaskConfig::TDPConfig(),
                          _col == length(_taskContext.getSeqH()) - 1,
                          _row == length(_taskContext.getSeqV()) - 1);

        
        // DEBUG: Remove!
//        _taskContext.getDebugBuffer().matrix[_col][_row] = {bufHBegin, bufVBegin,
//                                                            _taskContext.getTileBuffer().horizontalBuffer[_col],
//                                                            _taskContext.getTileBuffer().verticalBuffer[_row],
//                                                            _col, _row, false};

        swap(getDpTraceMatrix(dpContext), pTls.mLocalTraceStore.localScalarTraceMatrix());
        _taskContext.getTraceProxy().insert(std::make_pair(_col, _row),
                                            TTraceProxyValue{&pTls.mLocalTraceStore,
                                                             {0, static_cast<uint16_t>(length(pTls.mLocalTraceStore.mScalarTraceVec) - 1)},
                                                             0});
        #ifdef DP_ALIGN_STATS
            ++serialCounter;
        #endif
    }

    template <typename TTasks, typename TScoreValue, typename TSimdVec>
    inline void runSimd(TTasks const & tasks, impl::dp::parallel::DPLocalStorage<TScoreValue, TSimdVec> & pTls)
    {
        using TDPLocalStorage = impl::dp::parallel::DPLocalStorage<TScoreValue, TSimdVec>;
        using TLocalTraceStore = typename TDPLocalStorage::TLocalTraceStore;
        using TSimdTraceMatrix = typename TLocalTraceStore::TSimdTraceMatrix;
        using TTraceProxy = typename std::decay<decltype(_taskContext.getTraceProxy())>::type;
        using TTraceProxyValue = typename TTraceProxy::TTraceMatrixIdentifier;

        using TStateThreadContext = impl::dp::parallel::StateThreadContext<TTasks const, TDPLocalStorage>;

        // Prepare scout state.

        // DEBUG: Remove!
//        for (auto& task : tasks)
//        {
//            auto& val = _taskContext.getDebugBuffer().matrix[task->_col][task->_row];
//            val.hBegin = _taskContext.getTileBuffer().horizontalBuffer[task->_col];
//            val.vBegin = _taskContext.getTileBuffer().verticalBuffer[task->_row];
//            val.col = task->_col;
//            val.row = task->_row;
//            val.isSimd = true;
//        }

        auto simdBufferH = impl::gatherSimdBuffer<TSimdVec>(tasks, _taskContext.getTileBuffer().horizontalBuffer, [](auto& task){ return task->_col; });
        auto simdBufferV = impl::gatherSimdBuffer<TSimdVec>(tasks, _taskContext.getTileBuffer().verticalBuffer, [](auto& task){ return task->_row; });

        // Prepare dpContext.
        DPContext<DPCell_<TSimdVec, TGapMode>,
                  typename Value<TSimdTraceMatrix>::Type,
                  String<DPCell_<TSimdVec, TGapMode>, Alloc<OverAligned> >,
                  TSimdTraceMatrix> dpContext;

        TStateThreadContext stateLocalContext(tasks, pTls);

        // Run alignment.
        impl::computeSimdBatch(dpContext, stateLocalContext, simdBufferH, simdBufferV,
                               _taskContext.getSeqH(), _taskContext.getSeqV(),
                               _taskContext.getScore(), _taskContext.getBand(),
                               typename TTaskConfig::TDPConfig());

        // Call in simd block!
//        debug::compareToScalar(tasks, pTls,
//                               _taskContext.getTileBuffer().horizontalBuffer, _taskContext.getTileBuffer().verticalBuffer,
//                               simdBufferH, simdBufferV, dpContext,
//                               typename TTaskConfig::TDPConfig());

        // Swap trace matrix into local thread store.
        swap(getDpTraceMatrix(dpContext), pTls.mLocalTraceStore.localSimdTraceMatrix());

        uint8_t simdLane = 0;
        for (const auto& task : tasks)
        {
            _taskContext.getTraceProxy().insert(std::make_pair(task->_col, task->_row),
                                                TTraceProxyValue{&pTls.mLocalTraceStore,
                                                                {1, static_cast<uint16_t>(length(pTls.mLocalTraceStore.mSimdTraceVec) - 1)},
                                                                simdLane++});
        }

        // Write back into buffer.
        impl::scatterSimdBuffer(_taskContext.getTileBuffer().horizontalBuffer, tasks, simdBufferH, [](auto& task){ return task->_col; });
        impl::scatterSimdBuffer(_taskContext.getTileBuffer().verticalBuffer, tasks, simdBufferV, [](auto& task){ return task->_row; });

        #ifdef DP_ALIGN_STATS
                ++simdCounter;
        #endif
        // DEBUG: Remove!
//        for (auto& task : tasks)
//        {
//            auto& val = _taskContext.getDebugBuffer().matrix[task->_col][task->_row];
//            val.hEnd = _taskContext.getTileBuffer().horizontalBuffer[task->_col];
//            val.vEnd = _taskContext.getTileBuffer().verticalBuffer[task->_row];
//        }
    }
};

// ----------------------------------------------------------------------------
// Class DPTaskGraph
// ----------------------------------------------------------------------------

template <typename TTask, typename TSpec = void>
class DPTaskGraph
{
public:
    using TTaskPtr_ = typename Value<DPTaskGraph>::Type;

    std::vector<std::vector<TTaskPtr_> > data;

    inline auto&
    get()
    {
        return data;
    }

    template <typename TPos>
    inline auto&
    operator[](TPos const pos)
    {
        return data[pos];
    }

    template <typename TPos>
    inline auto&
    operator[](TPos const pos) const
    {
        return data[pos];
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TTask>
struct IsDPTask : False
{};

template <typename TTask>
struct IsDPTask<DPTaskBase<TTask>> : True
{};

template <typename TTask, typename TSpec>
struct Value<DPTaskGraph<TTask, TSpec> >
{
    using Type = typename Pointer_<TTask>::Type;
};

template <typename TTask, typename TSpec>
struct Value<DPTaskGraph<TTask, TSpec> const>
{
    using Type = typename Value<DPTaskGraph<TTask, TSpec> >::Type const;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TTask, typename TSpec>
inline typename Reference<DPTaskGraph<TTask, TSpec> >::Type
firstTask(DPTaskGraph<TTask, TSpec> & graph)
{
    return front(front(graph.get()));
}

template <typename TTask, typename TSpec>
inline typename Reference<DPTaskGraph<TTask, TSpec> >::Type
lastTask(DPTaskGraph<TTask, TSpec> & graph)
{
    return back(back(graph.get()));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DAG_TASK_BASE_H_
