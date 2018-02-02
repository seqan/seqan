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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

namespace impl
{

/*
 * Executor class that implements the correct execution mode.
 */
struct ParallelAlignmentExecutor
{
    template <typename TKernel,
              typename TSetH,
              typename TSetV,
              typename ...TArgs>
    auto operator()(Sequential const & /*execPolicy*/,
                    TKernel && kernel,
                    TSetH const & setH,
                    TSetV const & setV,
                    TArgs && ...args)
    {
        SEQAN_ASSERT_EQ(length(setH), length(setV));

        using TResult = decltype(kernel(setH, setV, std::forward<TArgs>(args)...));

        TResult superSet;
        resize(superSet, length(setH));

        auto zipCont = makeZipView(setH, setV, superSet);
#ifdef DP_PARALLEL_SHOW_PROGRESS
        ::impl::dp_parallel_progress::show_progress(length(setH));
#endif  // DP_PARALLEL_SHOW_PROGRESS
        for (auto && pwInst : zipCont)
        {
            std::get<2>(pwInst) = kernel(std::get<0>(pwInst), std::get<1>(pwInst), std::forward<TArgs>(args)...);
        }
        return superSet;
    }

    template <typename TKernel,
              typename TSetH,
              typename ...TArgs>
    auto operator()(ExecutionPolicy<Serial, Vectorial> const & /*execPolicy*/,
                    TKernel && kernel,
                    TSetH const & setH,
                    TArgs && ...args)
    {
#ifdef DP_PARALLEL_SHOW_PROGRESS
        ::impl::dp_parallel_progress::show_progress(length(setH));
#endif  // DP_PARALLEL_SHOW_PROGRESS
        // Automaically chooses vectorized code, or falls back to sequential code.
        return kernel(setH, std::forward<TArgs>(args)...);
    }

    template <typename TKernel,
              typename TSetH,
              typename TSetV,
              typename ...TArgs>
    auto operator()(SEQAN_UNUSED ExecutionPolicy<Parallel, Vectorial> const & execPolicy,  // maybe unused due to missing OMP support in clang.
                    TKernel && kernel,
                    TSetH const & setH,
                    TSetV const & setV,
                    TArgs && ...args)

    {
        SEQAN_ASSERT_EQ(length(setH), length(setV));

        using TPos = std::make_signed_t<decltype(length(setH))>;
        using TResult = decltype(kernel(setH, setV, std::forward<TArgs>(args)...));

        TPos chunkSize = _min(static_cast<TPos>(length(setH)), static_cast<TPos>(256));
        String<TPos> splitter;
        computeSplitters(splitter, length(setH), static_cast<TPos>(length(setH)/chunkSize));

        std::vector<TResult> superSet;
        superSet.resize(length(splitter));

#ifdef DP_PARALLEL_SHOW_PROGRESS
        ::impl::dp_parallel_progress::show_progress(length(setH));
#endif  // DP_PARALLEL_SHOW_PROGRESS

        SEQAN_OMP_PRAGMA(parallel for num_threads(numThreads(execPolicy)) schedule(guided))
        for (TPos job = 0; job < static_cast<TPos>(length(splitter)) - 1; ++job)  // TODO(rrahn): Why -1; Is there a bug in computeSplitters?
        {
            auto infSetH = infix(setH, splitter[job], splitter[job + 1]);
            auto infSetV = infix(setV, splitter[job], splitter[job + 1]);

            superSet[job] = kernel(infSetH, infSetV, std::forward<TArgs>(args)...);
        }
        // Reduce the result.
        TResult res;
        resize(res, length(setH));
        auto it = begin(res, Standard());
        for (auto && set : superSet)
        {
            arrayMoveForward(begin(set, Standard()), end(set, Standard()), it);
            it += length(set);
        }
        return res;
    }

    template <typename TKernel,
              typename TSetH,
              typename TSetV,
              typename ...TArgs>
    auto operator()(ExecutionPolicy<Parallel, Serial> const & execPolicy,
                    TKernel && kernel,
                    TSetH const & setH,
                    TSetV const & setV,
                    TArgs && ...args)

    {
        SEQAN_ASSERT_EQ(length(setH), length(setV));

        using TPos = std::make_signed_t<decltype(length(setH))>;
        using TResult = decltype(kernel(setH, setV, std::forward<TArgs>(args)...));

        Splitter<TPos> splitter(0, length(setH), numThreads(execPolicy));

        TResult superSet;
        resize(superSet, length(setH));

        auto zipCont = makeZipView(setH, setV, superSet);

#ifdef DP_PARALLEL_SHOW_PROGRESS
        ::impl::dp_parallel_progress::show_progress(length(setH));
#endif  // DP_PARALLEL_SHOW_PROGRESS

        SEQAN_OMP_PRAGMA(parallel for num_threads(length(splitter)))
        for (TPos job = 0; job < static_cast<TPos>(length(splitter)); ++job)
        {
            auto it = begin(zipCont, Standard()) + splitter[job];
            auto itEnd = begin(zipCont, Standard()) + splitter[job + 1];

            // NOTE(marehr): auto && seqPair does not work, thus declaring the
            // type explicitly, s.t. <=icpc 18.0.1 can compile the code (ticket
            // #03204483)
            using TSeqPair = decltype(*it);
            std::for_each(it, itEnd, [&](TSeqPair && seqPair)
            {
                std::get<2>(seqPair) = kernel(std::get<0>(seqPair), std::get<1>(seqPair), std::forward<TArgs>(args)...);
            });
        }
        return superSet;
    }
};

template <typename TWaveSpec, typename TVectorizationPolicy,
          typename TAlgorithmSpec,
          typename TSetH,
          typename TSetV,
          typename TScore,
          typename ...TArgs,
          std::enable_if_t<!std::is_same<WavefrontAlignment<TWaveSpec>, Serial>::value &&
                           !std::is_same<WavefrontAlignment<TWaveSpec>, Parallel>::value,
                           int> = 0>
inline auto
doWaveAlignment(ExecutionPolicy<WavefrontAlignment<TWaveSpec>, TVectorizationPolicy> const & execPolicy,
                TAlgorithmSpec const & /*tag*/,
                TSetH const & setH,
                TSetV const & setV,
                TScore const & scoringScheme,
                TArgs && .../*args*/)
{
    using TScoreValue = typename Value<TScore>::Type;

    // The vector containing the scores.
    std::vector<TScoreValue> res;
    res.resize(length(setH));

    auto dispatcher = [&res](auto && ...args)
    {
        alignExecBatch(std::forward<decltype(args)>(args)...,
                             [&res](auto const id, auto const score)
                             {
                                 res[id] = score;
                             });
    };

    // Differentiate between affine and linear gap costs.
    // TODO(rrahn): Setup configuration cascade.
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
    {
        struct DPConfigTraits
        {
            using TAlgorithmType SEQAN_UNUSED = TAlgorithmSpec;
            using TGapType       SEQAN_UNUSED = LinearGaps;
            using TBandType      SEQAN_UNUSED = BandOff;
            using TTracebackType SEQAN_UNUSED = TracebackOff;
            using TFormat        SEQAN_UNUSED = ArrayGaps;
        };

        using TDPSettings = seqan::DPSettings<TScore, DPConfigTraits>;

        TDPSettings settings;
        settings.scoringScheme = scoringScheme;
        dispatcher(execPolicy, setH, setV, settings);
    }
    else
    {
        struct DPConfigTraits
        {
            using TAlgorithmType SEQAN_UNUSED = TAlgorithmSpec;
            using TGapType       SEQAN_UNUSED = AffineGaps;
            using TBandType      SEQAN_UNUSED = BandOff;
            using TTracebackType SEQAN_UNUSED = TracebackOff;
            using TFormat        SEQAN_UNUSED = ArrayGaps;
        };

        using TDPSettings = seqan::DPSettings<TScore, DPConfigTraits>;

        TDPSettings settings;
        settings.scoringScheme = scoringScheme;
        dispatcher(execPolicy, setH, setV, settings);
    }
    return res;
}

} // namespace impl

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*
 * Wrapper functions for calling globalAlignmentScore and localAlignmentScore with an ExecutionPolicy.
 * Note the parallel interfaces are documented as part of the standard documentation in seqan/align module.
 */
template <typename TParallelPolicy, typename TVectorizationPolicy,
          typename ...TArgs,
          std::enable_if_t<std::is_same<TParallelPolicy, Serial>::value ||
                           std::is_same<TParallelPolicy, Parallel>::value,
                           int> = 0>
inline auto
globalAlignmentScore(ExecutionPolicy<TParallelPolicy, TVectorizationPolicy> const & execPolicy,
                     TArgs && ...args)
{
    auto kernel = [](auto && ...args)
    {
        return globalAlignmentScore(std::forward<decltype(args)>(args)...);
    };
    return impl::ParallelAlignmentExecutor{}(execPolicy, kernel, std::forward<TArgs>(args)...);
}

template <typename TParallelPolicy, typename TVectorizationPolicy,
          typename ...TArgs,
          std::enable_if_t<std::is_same<TParallelPolicy, Serial>::value ||
                           std::is_same<TParallelPolicy, Parallel>::value,
                           int> = 0>
inline auto
localAlignmentScore(ExecutionPolicy<TParallelPolicy, TVectorizationPolicy> const & execPolicy,
                    TArgs && ...args)
{
    auto kernel = [](auto && ...args)
    {
        return localAlignmentScore(std::forward<decltype(args)>(args)...);
    };
    return impl::ParallelAlignmentExecutor{}(execPolicy, kernel, std::forward<TArgs>(args)...);
}

// Wavefront execution of globalAlignmentScore w/ config.
template <typename TWaveSpec, typename TVectorizationPolicy,
          typename TSetH,
          typename TSetV,
          typename TScore,
          typename TConfig,
          std::enable_if_t<!std::is_same<WavefrontAlignment<TWaveSpec>, Serial>::value &&
                           !std::is_same<WavefrontAlignment<TWaveSpec>, Parallel>::value,
                           int> = 0>

inline auto
globalAlignmentScore(ExecutionPolicy<WavefrontAlignment<TWaveSpec>, TVectorizationPolicy> const & execPolicy,
                     TSetH const & setH,
                     TSetV const & setV,
                     TScore const & scoringScheme,
                     TConfig const & /*config*/)
{
    return impl::doWaveAlignment(execPolicy,
                                 GlobalAlignment_<typename SubstituteAlignConfig_<TConfig>::Type>{},
                                 setH,
                                 setV,
                                 scoringScheme);
}

// Wavefront execution of globalAlignmentScore w/o config.
template <typename TWaveSpec, typename TVectorizationPolicy,
          typename TSetH,
          typename TSetV,
          typename TScore,
          std::enable_if_t<!std::is_same<WavefrontAlignment<TWaveSpec>, Serial>::value &&
                           !std::is_same<WavefrontAlignment<TWaveSpec>, Parallel>::value,
                           int> = 0>

inline auto
globalAlignmentScore(ExecutionPolicy<WavefrontAlignment<TWaveSpec>, TVectorizationPolicy> const & execPolicy,
                     TSetH const & setH,
                     TSetV const & setV,
                     TScore const & scoringScheme)
{
    return globalAlignmentScore(execPolicy, setH, setV, scoringScheme, AlignConfig<>{});
}

template <typename TWaveSpec, typename TVectorizationPolicy,
          typename ...TArgs,
          std::enable_if_t<!std::is_same<WavefrontAlignment<TWaveSpec>, Serial>::value &&
                           !std::is_same<WavefrontAlignment<TWaveSpec>, Parallel>::value,
                           int> = 0>
inline auto
localAlignmentScore(ExecutionPolicy<WavefrontAlignment<TWaveSpec>, TVectorizationPolicy> const & execPolicy,
                    TArgs && ...args)
{
    return impl::doWaveAlignment(execPolicy, LocalAlignment_<>{}, std::forward<TArgs>(args)...);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_INTERFACE_H_
