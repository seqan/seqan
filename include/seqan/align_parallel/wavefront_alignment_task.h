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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_TASK_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_TASK_H_

namespace seqan
{
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Config structre for the execution of one alignment using the wave-front model.
template <typename TDPSettings>
struct WavefrontAlignmentTaskConfig
{
    // ----------------------------------------------------------------------------
    // Member Typedefs.

    // DPTrait type forwarding.
    using TDPTraits       = typename TDPSettings::TTraits;
    using TScoreValue     = typename Value<typename TDPSettings::TScoringScheme>::Type;
    using TAlgorithmType  = typename TDPTraits::TAlgorithmType;
    using TTracebackType  = typename TDPTraits::TTracebackType;
    using TGapType        = typename TDPTraits::TGapType;

    // Wavefront Alignment Context.
    using TDPCell         = DPCell_<TScoreValue, TGapType>;
    using TBufferValue    = Pair<TDPCell, typename TraceBitMap_<>::Type>;
    using TBuffer         = String<TBufferValue>;
    using TBlockBuffer    = DPTileBuffer<TBuffer>;

    // DP Execution Context.
    using TDPProfile      = DPProfile_<TAlgorithmType, TGapType, TTracebackType, Parallel>;
    using TDPCache        = DPContext<TDPCell, typename TraceBitMap_<>::Type>;
    using TDPScout        = DPScout_<TDPCell, Default>;

    // Parallel Context.
    struct IntermediateTraits_
    {
        using TScoreValue   = decltype(maxScore(std::declval<TDPScout>()));
        using THostPosition = decltype(maxHostPosition(std::declval<TDPScout>()));
    };

    using TDPIntermediate = WavefrontAlignmentResult<IntermediateTraits_>;

    struct AlignThreadLocalConfig_
    {
        using TIntermediate = TDPIntermediate;
        using TCache        = TDPCache;

        using TLocalHost    = std::tuple<TIntermediate, TCache>;
    };

    using TThreadLocal = WavefrontAlignmentThreadLocalStorage<AlignThreadLocalConfig_>;
    using TAlignEvent  = WavefrontTaskEvent;
};

#ifdef SEQAN_SIMD_ENABLED
template <typename TDPSettings>
struct WavefrontAlignmentSimdTaskConfig : public WavefrontAlignmentTaskConfig<TDPSettings>
{
    // ----------------------------------------------------------------------------
    // Member Typedefs.

    using TBase_ = WavefrontAlignmentTaskConfig<TDPSettings>;

    using TDPSimdCell        = DPCell_<typename TDPSettings::TScoreValueSimd, typename TBase_::TGapType>;
    using TDPSimdTraceValue  = typename TraceBitMap_<typename TDPSettings::TScoreValueSimd>::Type;

    using TDPSimdScoreMatrix = String<TDPSimdCell, Alloc<OverAligned>>;
    using TDPSimdTraceMatrix = String<TDPSimdTraceValue, Alloc<OverAligned>>;
    using TDPSimdCache       = DPContext<TDPSimdCell, TDPSimdTraceValue, TDPSimdScoreMatrix, TDPSimdTraceMatrix>;

    using TDPScout_          = DPScout_<TDPSimdCell, SimdAlignmentScout<> >;
    using TDPIntermediate    = WavefrontAlignmentResult<typename TBase_::IntermediateTraits_>;

    // Parallel Context.
    struct SimdAlignThreadLocalConfig_
    {

        using TIntermediate = TDPIntermediate;
        using TCache        = typename TBase_::TDPCache;
        using TSimdCache    = TDPSimdCache;

        using TLocalHost    = std::tuple<TIntermediate, TCache, TSimdCache>;
    };

    using TThreadLocal = WavefrontAlignmentThreadLocalStorage<SimdAlignThreadLocalConfig_>;
    using TAlignEvent  = WavefrontTaskEvent;
};
#endif

// Incubator to setup the alignment job.
template <typename WavefrontAlignmentTaskConfigConcept>
struct WavefrontAlignmentTaskIncubator
{
    using TWatc = WavefrontAlignmentTaskConfigConcept;

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

    // ----------------------------------------------------------------------------
    // Function createBlockBuffer()
    // ----------------------------------------------------------------------------

    template <typename TSeqHBlocks, typename TSeqVBlovcks, typename TScore>
    static auto createBlockBuffer(TSeqHBlocks const & seqHBlocks, TSeqVBlovcks const & seqVBlocks, TScore const & score)
    {
        using TDPCell = typename TWatc::TDPCell;
        typename TWatc::TBlockBuffer buffer;
        resize(buffer.horizontalBuffer, length(seqHBlocks), Exact());
        resize(buffer.verticalBuffer, length(seqVBlocks), Exact());

        typename TWatc::TBufferValue tmp;

        using TDPMetaColH = DPMetaColumn_<typename TWatc::TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn>>;
        using TDPMetaColV = DPMetaColumn_<typename TWatc::TDPProfile, MetaColumnDescriptor<DPInitialColumn, FullColumn>>;

        TDPCell dummyCellD;
        TDPCell dummyCellH;
        TDPCell dummyCellV;
        tmp.i2 = _computeScore(tmp.i1, dummyCellD, dummyCellH, dummyCellV,  Nothing(), Nothing(), score,
                               RecursionDirectionZero(), typename TWatc::TDPProfile());
        for (auto itH = begin(buffer.horizontalBuffer, Standard());
             itH != end(buffer.horizontalBuffer, Standard());
             ++itH)
        {
            resize(*itH, length(front(seqHBlocks)), Exact());
            for (auto it = begin(*itH, Standard()); it != end(*itH, Standard()); ++it)
            {
                it->i2 = _computeScore(it->i1, dummyCellD, tmp.i1, dummyCellV, Nothing(), Nothing(), score,
                                       typename RecursionDirection_<TDPMetaColH, FirstCell>::Type(),
                                       typename TWatc::TDPProfile());
                tmp.i1 = it->i1;
            }
        }
        tmp.i1 = decltype(tmp.i1){};
        tmp.i2 = _computeScore(tmp.i1, dummyCellD, dummyCellH, dummyCellV, Nothing(), Nothing(), score,
                               RecursionDirectionZero(), typename TWatc::TDPProfile());

        for (auto itV = begin(buffer.verticalBuffer, Standard()); itV != end(buffer.verticalBuffer, Standard()); ++itV)
        {
            resize(*itV, length(front(seqVBlocks)) + 1, Exact());
            auto it = begin(*itV, Standard());
            it->i2 = tmp.i2;
            it->i1 = tmp.i1;
            ++it;
            for (; it != end(*itV, Standard()); ++it)
            {
                it->i2 = _computeScore(it->i1, dummyCellD, dummyCellH, dummyCellV, Nothing(), Nothing(), score,
                                       typename RecursionDirection_<TDPMetaColV, InnerCell>::Type(),
                                       typename TWatc::TDPProfile());
                _setVerticalScoreOfCell(it->i1, _verticalScoreOfCell(dummyCellV));
                tmp.i1 = it->i1;
                tmp.i2 = it->i2;  // TODO(rrahn): Move out of loop.
            }
        }
        return buffer;
    }

    // ----------------------------------------------------------------------------
    // Function createTaskGraph()
    // ----------------------------------------------------------------------------

    template <typename TWavefrontTaskContext>
    static auto createTaskGraph(TWavefrontTaskContext & taskContext)
    {
        using TDagTask = WavefrontTask<TWavefrontTaskContext>;

        std::vector<std::vector<std::shared_ptr<TDagTask>>> graph;

        resize(graph, length(taskContext.seqHBlocks));
        for (int i = length(taskContext.seqHBlocks); --i >= 0;)
        {
            resize(graph[i], length(taskContext.seqVBlocks));
            for (int j = length(taskContext.seqVBlocks); --j >= 0;)
            {
                using TSize = decltype(length(taskContext.seqHBlocks));
                TDagTask * successorRight = (static_cast<TSize>(i + 1) < length(taskContext.seqHBlocks))
                                                ?  graph[i+1][j].get()
                                                : nullptr;
                TDagTask * successorDown  = (static_cast<TSize>(j + 1) < length(taskContext.seqVBlocks))
                                                ? graph[i][j+1].get()
                                                : nullptr;
                graph[i][j] = std::make_shared<TDagTask>(taskContext,
                                                         std::array<TDagTask*, 2>{{successorRight, successorDown}},
                                                         static_cast<size_t>(i), static_cast<size_t>(j),
                                                         static_cast<size_t>(((i > 0) ? 1 : 0) + ((j > 0) ? 1 : 0)),
                                                         (static_cast<TSize>(i + 1) == length(taskContext.seqHBlocks)),
                                                         (static_cast<TSize>(j + 1) == length(taskContext.seqVBlocks)));
            }
        }
        return graph;
    }
};

// The actual alignment task that is executed by the wave-front model.
template <typename TSeqH,
          typename TSeqV,
          typename TDPSettings,
          typename TConfig = WavefrontAlignmentTaskConfig<TDPSettings>>
class WavefrontAlignmentTask
{
public:

    using TIncubator    = WavefrontAlignmentTaskIncubator<TConfig>;

    using TSeqHBlocks   = decltype(TIncubator::createBlocks(std::declval<TSeqH>(), std::declval<size_t>()));
    using TSeqVBlocks   = decltype(TIncubator::createBlocks(std::declval<TSeqV>(), std::declval<size_t>()));
    using TTileBuffer   = decltype(TIncubator::createBlockBuffer(std::declval<TSeqHBlocks>(),
                                                                 std::declval<TSeqVBlocks>(),
                                                                 std::declval<typename TDPSettings::TScoringScheme>()));

    using TTaskContext  = WavefrontAlignmentContext<TSeqHBlocks, TSeqVBlocks, TTileBuffer, TDPSettings>;

    // ----------------------------------------------------------------------------
    // Member Variables.
    // ----------------------------------------------------------------------------

    size_t              alignmentId{0};
    TSeqH const &       seqH;
    TSeqV const &       seqV;
    TDPSettings const & dpSettings;
    size_t              blockSize;

    // ----------------------------------------------------------------------------
    // Constructors.
    // ----------------------------------------------------------------------------

    WavefrontAlignmentTask() = delete;

    WavefrontAlignmentTask(TSeqH const & seqH,
                           TSeqV const & seqV,
                           TDPSettings const & dpSetting,
                           size_t const & blockSize) :
        seqH(seqH),
        seqV(seqV),
        dpSettings(dpSetting),
        blockSize(blockSize)
    {}


    WavefrontAlignmentTask(size_t const id,
                           TSeqH const & seqH,
                           TSeqV const & seqV,
                           TDPSettings const & dpSetting,
                           size_t const & blockSize) :
        alignmentId(id),
        seqH(seqH),
        seqV(seqV),
        dpSettings(dpSetting),
        blockSize(blockSize)
    {}

    // ----------------------------------------------------------------------------
    // Member Functions.
    // ----------------------------------------------------------------------------

    // This function now run's in a separate thread.
    template <typename TWavefrontExecutor,
              typename TCallback>
    inline void
    operator()(uint16_t const instanceId,
               TWavefrontExecutor & executor,
               TCallback && callback)
    {
        // Initialize the strings.
        auto seqHBlocks = TIncubator::createBlocks(seqH, blockSize);
        auto seqVBlocks = TIncubator::createBlocks(seqV, blockSize);

        // Create the buffer for the matrix.
        auto buffer = TIncubator::createBlockBuffer(seqHBlocks, seqVBlocks, dpSettings.scoringScheme);

        // Setup the task context and create task graph.
        TTaskContext taskContext{instanceId, seqHBlocks, seqVBlocks, buffer, dpSettings};
        auto taskGraph = TIncubator::createTaskGraph(taskContext);

        // Prepare event.
        WavefrontTaskEvent event;
        context(*taskGraph.back().back()).ptrEvent = &event;

        // Kick off the execution.
        using TWavefrontTaskExec = WavefrontTaskExecutor<std::decay_t<decltype(*taskGraph[0][0])>, TWavefrontExecutor>;
        spawn(executor, TWavefrontTaskExec{taskGraph[0][0].get(), &executor});

        // Wait for alignment to finish.
        wait(event);

        // Reduce.
        typename TConfig::TDPIntermediate interMax{};
        auto collectAndReset = [&](auto & threadLocalStorage)
        {
            updateMax(interMax, intermediate(threadLocalStorage, instanceId));
            clear(intermediate(threadLocalStorage, instanceId));
        };
        combineEach(*executor.ptrThreadLocal, collectAndReset);
        // Continue execution.
        callback(alignmentId, interMax._maxState.first);
    }

    template <typename TWavefrontExecutor,
              typename TSimdTaskQueue,
              typename TCallback>
    inline void
    operator()(uint16_t const instanceId,
               TWavefrontExecutor & executor,
               TSimdTaskQueue & taskQueue,
               TCallback && callback)
    {
        // Initialize the strings.
        auto seqHBlocks = TIncubator::createBlocks(seqH, blockSize);
        auto seqVBlocks = TIncubator::createBlocks(seqV, blockSize);

        // Create the buffer for the matrix.
        auto buffer = TIncubator::createBlockBuffer(seqHBlocks, seqVBlocks, dpSettings.scoringScheme);

        // Setup the task context and create task graph.
        TTaskContext taskContext{instanceId, seqHBlocks, seqVBlocks, buffer, dpSettings};
        auto taskGraph = TIncubator::createTaskGraph(taskContext);

        // Prepare event.
        WavefrontTaskEvent event;
        context(*taskGraph.back().back()).ptrEvent = &event;

        // Kick off the execution.
        using TWavefrontTaskExec = WavefrontTaskExecutor<TSimdTaskQueue, TWavefrontExecutor>;
        appendValue(taskQueue, *taskGraph[0][0]);
        spawn(executor, TWavefrontTaskExec{&taskQueue, &executor});

        // Wait for alignment to finish.
        wait(event);

        // Reduce.
        typename TConfig::TDPIntermediate interMax{};
        auto collectAndReset = [&](auto & threadLocalStorage)
        {
            updateMax(interMax, intermediate(threadLocalStorage, instanceId));
            clear(intermediate(threadLocalStorage, instanceId));
        };
        combineEach(*executor.ptrThreadLocal, collectAndReset);
        callback(alignmentId, interMax._maxState.first);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_TASK_H_
