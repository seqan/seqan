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

// #define DP_PARALLEL_SHOW_PROGRESS // Enable when debugging.

#include <seqan/align_parallel.h>

#include "../align/test_mock.h"

namespace test_align_parallel
{

template <typename TSets, typename TResults, typename ...TParams>
inline void
validateGlobal(TSets const & sets,
               TResults const & res,
               TParams && ...params)
{
    auto z = makeZipView(std::get<0>(sets), std::get<1>(sets), res);

    for (auto && inst : z)
    {
        auto tmp = globalAlignmentScore(std::get<0>(inst), std::get<1>(inst), std::forward<TParams>(params)...);
        SEQAN_ASSERT_EQ(tmp, std::get<2>(inst));
    }
}

template <typename TSets, typename TResults, typename ...TParams>
inline void
validateLocal(TSets const & sets,
              TResults const & res,
              TParams && ...params)
{
    auto z = makeZipView(std::get<0>(sets), std::get<1>(sets), res);

    for (auto && inst : z)
    {
        auto tmp = localAlignmentScore(std::get<0>(inst), std::get<1>(inst), std::forward<TParams>(params)...);
        SEQAN_ASSERT_EQ(tmp, std::get<2>(inst));
    }
}

}  // namespace test_align_parallel

// ----------------------------------------------------------------------------
// Class SimdAlignTest
// ----------------------------------------------------------------------------

// Common test class instance, which stores the types to be accessed.
template <typename TTuple>
class ParallelAlignInterfaceTest : public seqan::Test
{
public:
    using TExecPolicy = std::tuple_element_t<0, TTuple>;
};

// ----------------------------------------------------------------------------
// Configuration of typed tests for global alignment.
// ----------------------------------------------------------------------------

template <typename T>
class ParallelAlignInterfaceTestCommon : public ParallelAlignInterfaceTest<T>
{};

typedef
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Serial,                                             seqan::Serial>>,
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Parallel,                                           seqan::Serial>>,
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<>,                               seqan::Serial>>,
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>, seqan::Serial>>
#ifdef SEQAN_SIMD_ENABLED
        ,
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Serial,                                             seqan::Vectorial>>,
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Parallel,                                           seqan::Vectorial>>,
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<>,                               seqan::Vectorial>>,
        seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>, seqan::Vectorial>>
        > > > >
#endif // SEQAN_SIMD_ENABLED
        > > > > ParallelAlignInterfaceTestCommonTypes;

SEQAN_TYPED_TEST_CASE(ParallelAlignInterfaceTestCommon, ParallelAlignInterfaceTestCommonTypes);

SEQAN_TYPED_TEST(ParallelAlignInterfaceTestCommon, Global_Score)
{
    using namespace seqan;
    using TExecPolicy = typename TestFixture::TExecPolicy;

    auto sets = ::impl::test_align_mock::TestSequences_<Dna, ::impl::test_align_mock::EqualLengthSimd>::getSequences();

    Score<int, Simple> scoreLinear(4, -2, -4);
    TExecPolicy execPolicy;
    setNumThreads(execPolicy, 4);
    auto score = globalAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreLinear);

    test_align_parallel::validateGlobal(sets, score, scoreLinear);

    Score<int, Simple> scoreAffine(4, -2, -4, -10);
    score = globalAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreAffine);
    test_align_parallel::validateGlobal(sets, score, scoreAffine);
}

SEQAN_TYPED_TEST(ParallelAlignInterfaceTestCommon, Semi_Global_Score)
{
    using namespace seqan;
    using TExecPolicy = typename TestFixture::TExecPolicy;

    auto sets = ::impl::test_align_mock::TestSequences_<Dna, ::impl::test_align_mock::EqualLengthSimd>::getSequences();

    Score<int, Simple> scoreLinear(4, -2, -4);
    TExecPolicy execPolicy;
    setNumThreads(execPolicy, 4);
    auto score = globalAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreLinear, AlignConfig<true, false, false, true>());

    test_align_parallel::validateGlobal(sets, score, scoreLinear, AlignConfig<true, false, false, true>());

    Score<int, Simple> scoreAffine(4, -2, -4, -10);
    score = globalAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreAffine, AlignConfig<true, false, false, true>());
    test_align_parallel::validateGlobal(sets, score, scoreAffine, AlignConfig<true, false, false, true>());
}

SEQAN_TYPED_TEST(ParallelAlignInterfaceTestCommon, Local_Score)
{
    using namespace seqan;
    using TExecPolicy = typename TestFixture::TExecPolicy;

    auto sets = ::impl::test_align_mock::TestSequences_<Dna, ::impl::test_align_mock::EqualLengthSimd>::getSequences();

    Score<int, Simple> scoreLinear(4, -2, -4);
    TExecPolicy execPolicy;
    setNumThreads(execPolicy, 4);
    auto score = localAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreLinear);

    test_align_parallel::validateLocal(sets, score, scoreLinear);

    Score<int, Simple> scoreAffine(4, -2, -4, -10);
    score = localAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreAffine);
    test_align_parallel::validateLocal(sets, score, scoreAffine);
}
