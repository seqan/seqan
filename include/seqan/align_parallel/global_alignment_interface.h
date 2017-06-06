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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_GLOBAL_ALIGNMENT_INTERFACE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_GLOBAL_ALIGNMENT_INTERFACE_H_

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

struct ParallelAlignmentExecutor
{
    template <typename TSetH,
              typename TSetV,
              typename ...TArgs>
    auto operator()(Sequential const & /*execPolicy*/,
                    TSetH const & setH,
                    TSetV const & setV,
                    TArgs && ...args)
    {
        SEQAN_ASSERT_EQ(length(setH), length(setV));

        using TResult = decltype(globalAlignmentScore(setH, setV, std::forward<TArgs>(args)...));

        TResult superSet;
        resize(superSet, length(setH));

        auto zipCont = makeZipView(setH, setV, superSet);

        for (auto && pwInst : zipCont)
        {
            std::get<2>(pwInst) = globalAlignmentScore(std::get<0>(pwInst), std::get<1>(pwInst),
                                                       std::forward<TArgs>(args)...);
        }
        return superSet;
    }

    template <typename ...TArgs>
    auto operator()(ExecutionPolicy<Serial, Vectorial> const & /*execPolicy*/,
                    TArgs && ...args)
    {
        // Automaically chooses vectorized code, or falls back to sequential code.
        return globalAlignmentScore(std::forward<TArgs>(args)...);
    }

    template <typename TSetH,
              typename TSetV,
              typename ...TArgs>
    auto operator()(ExecutionPolicy<Parallel, Vectorial> const & execPolicy,
                    TSetH const & setH,
                    TSetV const & setV,
                    TArgs && ...args)

    {
        SEQAN_ASSERT_EQ(length(setH), length(setV));

        using TPos = decltype(length(setH));
        using TResult = decltype(globalAlignmentScore(setH, setV, std::forward<TArgs>(args)...));

        Splitter<TPos> splitter(0, length(setH), numThreads(execPolicy));

        std::vector<TResult> superSet;
        superSet.resize(length(splitter));

        SEQAN_OMP_PRAGMA(parallel for num_threads(length(splitter)))
        for (TPos job = 0; job < length(splitter); ++job)
        {
            auto infSetH = infix(setH, splitter[job], splitter[job + 1]);
            auto infSetV = infix(setV, splitter[job], splitter[job + 1]);

            superSet[job] = globalAlignmentScore(infSetH, infSetV, std::forward<TArgs>(args)...);
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

    template <typename TSetH,
              typename TSetV,
              typename ...TArgs>
    auto operator()(ExecutionPolicy<Parallel, Serial> const & execPolicy,
                    TSetH const & setH,
                    TSetV const & setV,
                    TArgs && ...args)

    {
        SEQAN_ASSERT_EQ(length(setH), length(setV));

        using TPos = decltype(length(setH));
        using TResult = decltype(globalAlignmentScore(setH, setV, std::forward<TArgs>(args)...));

        Splitter<TPos> splitter(0, length(setH), numThreads(execPolicy));

        TResult superSet;
        resize(superSet, length(setH));

        auto zipCont = makeZipView(setH, setV, superSet);

        SEQAN_OMP_PRAGMA(parallel for num_threads(length(splitter)))
        for (TPos job = 0; job < length(splitter); ++job)
        {
            auto it = begin(zipCont, Standard()) + splitter[job];
            auto itEnd = begin(zipCont, Standard()) + splitter[job + 1];

            std::for_each(it, itEnd, [&](auto && seqPair)
            {
                std::get<2>(seqPair) = globalAlignmentScore(std::get<0>(seqPair), std::get<1>(seqPair),
                                                            std::forward<TArgs>(args)...);
            });
        }
        return superSet;
    }
};

} // namespace impl

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TParallelPolicy, typename TVectorizationPolicy,
          typename ...TArgs,
          std::enable_if_t<std::is_same<TParallelPolicy, Serial>::value ||
                           std::is_same<TParallelPolicy, Parallel>::value,
                           int> = 0>
auto globalAlignmentScore(ExecutionPolicy<TParallelPolicy, TVectorizationPolicy> const & execPolicy,
                          TArgs && ...args)
{
    return impl::ParallelAlignmentExecutor{}(execPolicy, std::forward<TArgs>(args)...);
}

template <typename TWaveSpec, typename TVectorizationPolicy,
          typename TSetH,
          typename TSetV,
          typename TScore,
          typename ...TArgs,
          std::enable_if_t<!std::is_same<WavefrontAlignment<TWaveSpec>, Serial>::value &&
                           !std::is_same<WavefrontAlignment<TWaveSpec>, Parallel>::value,
                           int> = 0>
auto globalAlignmentScore(ExecutionPolicy<WavefrontAlignment<TWaveSpec>, TVectorizationPolicy> const & execPolicy,
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
        impl::alignExecBatch(std::forward<decltype(args)>(args)...,
                             [&res](auto const id, auto const score)
                             {
                                 res[id] = score;
                             });
    };

    // Differentiate between affine and linear gap costs.
    // TODO(rrahn): Setup configuration cascade.
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
    {
        struct DPConfigTraits : public DPTraits::GlobalLinear
        {
            using TTracebackType [[gnu::unused]] = seqan::TracebackOff;
        };

        using TDPSettings = seqan::DPSettings<TScore, DPConfigTraits>;

        TDPSettings settings;
        settings.mScoringScheme = scoringScheme;
        dispatcher(execPolicy, setH, setV, settings);
    }
    else
    {
        struct DPConfigTraits : public DPTraits::GlobalAffine
        {
            using TTracebackType [[gnu::unused]] = seqan::TracebackOff;
        };

        using TDPSettings = seqan::DPSettings<TScore, DPConfigTraits>;

        TDPSettings settings;
        settings.mScoringScheme = scoringScheme;
        dispatcher(execPolicy, setH, setV, settings);
    }
    return res;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_GLOBAL_ALIGNMENT_INTERFACE_H_
