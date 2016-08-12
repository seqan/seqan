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
// 6. parameter: TraceBlock
// 7. parameter: TDpContext
// 8. parameter: TDPConfig
// 9. parameter: TDPScoutState
template <typename... TArgs>
struct DPTaskContext
{
    using TData         = std::tuple<TArgs...>;
    using TSize         = typename Size<typename std::decay<std::tuple_element_t<0, TData> >::type>::Type;
    using TDPContext    = typename std::decay<std::tuple_element_t<6, TData> >::type;
    using TDPConfig     = typename std::decay<std::tuple_element_t<7, TData> >::type;
    using TDPScoutState = typename std::decay<std::tuple_element_t<8, TData> >::type;

    TData _data;
    std::mutex mLock;
    std::mutex mLockEvent;
    std::condition_variable mReadyEvent;
    bool       mReady = false;
    unsigned   mSimdLength = 8;

    DPTaskContext() = default;

    template <typename... TSubArgs>
    DPTaskContext(TSubArgs&& ...args)
    {
        static_assert(sizeof...(TArgs) == 9, "Requires 9 arguments.");
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

    inline auto& getTraceBlock()
    {
        return *std::get<5>(_data);
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
    using TDPScoutState = typename TTaskConfig::TDPScoutState;
    using TSize = typename TTaskConfig::TSize;
public:

    typename TTaskConfig::TSize     _col = 0;
    typename TTaskConfig::TSize     _row = 0;
    TTaskConfig&                    _taskContext;
    TDerivedTask&                   _derivedTask;
    std::array<TDerivedTask*, 2>    successor;

    DPTaskBase(TSize pCol, TSize pRow, TTaskConfig& pContext, TDerivedTask& pTask) :
        _col(pCol),
        _row(pRow),
        _taskContext(pContext),
        _derivedTask(pTask)
    {}

    template <typename TDPContext>
    inline void runScalar(TDPContext & dpContext)
    {
        getDpTraceMatrix(dpContext) = _taskContext.getTraceBlock()[(_col * length(_taskContext.getSeqV())) + _row]; // Thread local.
        TDPScoutState scoutState(_taskContext.getTileBuffer().horizontalBuffer[_col], _taskContext.getTileBuffer().verticalBuffer[_row]);  // Task local
        String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.
        _computeAlignment(dpContext, traceSegments, scoutState, _taskContext.getSeqH()[_col], _taskContext.getSeqV()[_row],
                          _taskContext.getScore(), _taskContext.getBand(), typename TTaskConfig::TDPConfig(),
                          _col == length(_taskContext.getSeqH()) - 1,
                          _row == length(_taskContext.getSeqV()) - 1);
    }

    template <typename TTasks, typename TDPContext>
    inline void runSimd(TTasks const & tasks, TDPContext & /*dpContext*/)
    {
//        getDpTraceMatrix(dpContext) = _taskContext.getTraceBlock()[(_col * length(_taskContext.getSeqV())) + _row]; // Thread local.
//        TDPScoutState scoutState(_taskContext.getTileBuffer().horizontalBuffer[_col], _taskContext.getTileBuffer().verticalBuffer[_row]);  // Task local
//        String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.
//        _computeAlignment(dpContext, traceSegments, scoutState, _taskContext.getSeqH()[_col], _taskContext.getSeqV()[_row],
//                          _taskContext.getScore(), _taskContext.getBand(), typename TTaskConfig::TDPConfig(),
//                          _col == length(_taskContext.getSeqH()) - 1,
//                          _row == length(_taskContext.getSeqV()) - 1);
        using TSeqH = typename std::decay<decltype(_taskContext.getSeqH()[_col])>::type;
        using TSeqV = typename std::decay<decltype(_taskContext.getSeqV()[_row])>::type;
        StringSet<Gaps<TSeqH> > setH;
        StringSet<Gaps<TSeqV> > setV;
        for (auto& task : tasks)
        {
            Gaps<TSeqH> gapsH;
            assignSource(gapsH, _taskContext.getSeqH()[task->_col]);
            Gaps<TSeqV> gapsV;
            assignSource(gapsV, _taskContext.getSeqV()[task->_row]);
            appendValue(setH, gapsH);
            appendValue(setV, gapsV);
        }
        globalAlignment(setH, setV, _taskContext.getScore());
    }

    // TODO(rrahn): Remove after refactoring omp dag.
    template <typename TDPContext>
    inline void execute(TDPContext & dpContext)
    {
        runScalar(dpContext);
        _derivedTask.updateAndSpawnScalar();
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
