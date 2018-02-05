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

#include <type_traits>

#include <seqan/basic.h>
#include <seqan/align_parallel.h>

// namespace test_align_parallel
// {
//
// template <typename TSets, typename TResults, typename ...TParams>
// inline void
// validateGlobal(TSets const & sets,
//                TResults const & res,
//                TParams && ...params)
// {
//     auto z = makeZipView(std::get<0>(sets), std::get<1>(sets), res);
//
//     for (auto && inst : z)
//     {
//         auto tmp = globalAlignmentScore(std::get<0>(inst), std::get<1>(inst), std::forward<TParams>(params)...);
//         SEQAN_ASSERT_EQ(tmp, std::get<2>(inst));
//     }
// }

// template <typename TSets, typename TResults, typename ...TParams>
// inline void
// validateLocal(TSets const & sets,
//               TResults const & res,
//               TParams && ...params)
// {
//     auto z = makeZipView(std::get<0>(sets), std::get<1>(sets), res);
//
//     for (auto && inst : z)
//     {
//         auto tmp = localAlignmentScore(std::get<0>(inst), std::get<1>(inst), std::forward<TParams>(params)...);
//         SEQAN_ASSERT_EQ(tmp, std::get<2>(inst));
//     }
// }

// }  // namespace test_align_parallel

// ----------------------------------------------------------------------------
// Class SimdAlignTest
// ----------------------------------------------------------------------------

// Common test class instance, which stores the types to be accessed.
// template <typename TTuple>
// class ParallelAlignBandedTest : public seqan::Test
// {
// public:
//     using TExecPolicy = std::tuple_element_t<0, TTuple>;
// };

// ----------------------------------------------------------------------------
// Configuration of typed tests for global alignment.
// ----------------------------------------------------------------------------

// template <typename T>
// class ParallelAlignInterfaceTestCommon : public ParallelAlignInterfaceTest<T>
// {};
//
// typedef
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Serial,                                             seqan::Serial>>,
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Parallel,                                           seqan::Serial>>,
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<>,                               seqan::Serial>>,
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>, seqan::Serial>>
// #ifdef SEQAN_SIMD_ENABLED
//         ,
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Serial,                                             seqan::Vectorial>>,
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::Parallel,                                           seqan::Vectorial>>,
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<>,                               seqan::Vectorial>>,
//         seqan::TagList<std::tuple<seqan::ExecutionPolicy<seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>, seqan::Vectorial>>
//         > > > >
// #endif // SEQAN_SIMD_ENABLED
//         > > > > ParallelAlignBandedTestCommonTypes;
//
// SEQAN_TYPED_TEST_CASE(ParallelAlignBandedTestCommon, ParallelAlignBandedTestCommonTypes);

// SEQAN_TYPED_TEST(ParallelAlignBandedTestCommon, Global_Score)
// {
//     using namespace seqan;
//     using TExecPolicy = typename TestFixture::TExecPolicy;
//
//     auto sets = ::impl::test_align_mock::TestSequences_<Dna, ::impl::test_align_mock::EqualLengthSimd>::getSequences();
//
//     Score<int, Simple> scoreLinear(4, -2, -4);
//     TExecPolicy execPolicy;
//     setNumThreads(execPolicy, 4);
//     setNumAl
//     auto score = globalAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreLinear, -100, 100);
//
//     test_align_parallel::validateGlobal(sets, score, scoreLinear);
//     //
//     // Score<int, Simple> scoreAffine(4, -2, -4, -10);
//     // score = globalAlignmentScore(execPolicy, std::get<0>(sets), std::get<1>(sets), scoreAffine);
//     // test_align_parallel::validateGlobal(sets, score, scoreAffine);
// }

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_construction)
{
    using namespace seqan;

    { // Custom constructor;
        impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                       impl::GridColumnSize{4}, impl::GridRowSize{5}};

        SEQAN_ASSERT_EQ(length(waveBand._grid), 20u);
        SEQAN_ASSERT_EQ(waveBand._hBandPos, 7u);
        SEQAN_ASSERT_EQ(waveBand._vBandPos, 4u);
        SEQAN_ASSERT_EQ(waveBand._colSize, 4u);
        SEQAN_ASSERT_EQ(waveBand._rowSize, 5u);
    }

    bool res = std::is_default_constructible<impl::DPWavefrontBand>::value;
    SEQAN_ASSERT_NOT(std::is_default_constructible<impl::DPWavefrontBand>::value);
    res = std::is_copy_constructible<impl::DPWavefrontBand>::value;
    SEQAN_ASSERT(res);
    res = std::is_move_constructible<impl::DPWavefrontBand>::value;
    SEQAN_ASSERT(res);
    res = std::is_copy_assignable<impl::DPWavefrontBand>::value;
    SEQAN_ASSERT(res);
    res = std::is_move_assignable<impl::DPWavefrontBand>::value;
    SEQAN_ASSERT(res);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_begin)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[0] = true;

    { // non-const container
        auto it = begin(waveBand, Standard());
        SEQAN_ASSERT(*it == waveBand._grid[0]);
    }

    { // const container
        impl::DPWavefrontBand const waveBandConst{waveBand};
        auto it = begin(waveBandConst, Standard());
        SEQAN_ASSERT(*it == waveBandConst._grid[0]);
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_end)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    back(waveBand._grid) = true;

    { // non-const container
        auto it = end(waveBand, Standard());
        SEQAN_ASSERT(*(it - 1) == back(waveBand._grid));
    }

    { // const container
        impl::DPWavefrontBand const waveBandConst{waveBand};
        auto it = end(waveBandConst, Standard());
        SEQAN_ASSERT(*(it - 1) == back(waveBandConst._grid));
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator)
{
    using namespace seqan;

    using TIterator      = Iter<impl::DPWavefrontBand, impl::DPWavefrontBandIterSpec>;
    using TIteratorConst = Iter<impl::DPWavefrontBand const, impl::DPWavefrontBandIterSpec>;

    SEQAN_ASSERT(std::is_default_constructible<TIterator>::value);

    {  // copy-constructor
        TIterator itA{};
        TIterator it{itA};

        SEQAN_ASSERT(it == itA);
    }

    {  // copy-constructor with const_iterator from iterator
        // TODO(rrahn): Fix adaptor iterator to make iterator comparable with different const qualified containers.
        // TIterator itA{};
        // TIteratorConst it{itA};
        // SEQAN_ASSERT(it == itA);
    }

    {  // copy-constructor with const_iterator from const_iterator
        TIteratorConst itA{};
        TIteratorConst it{itA};
        SEQAN_ASSERT(it == itA);
    }

    {  // move-constructor
        TIterator it{TIterator{}};
        SEQAN_ASSERT(it == TIterator{});
    }

    {  // move-constructor with const_iterator from iterator
        // TODO(rrahn): Fix adaptor iterator to make iterator comparable with different const qualified containers.
        // TIteratorConst it{TIterator{}};
        //SEQAN_ASSERT(it, TIterator{});
    }

    {  // move-constructor with const_iterator from const_iterator
        TIteratorConst it{TIteratorConst{}};
        SEQAN_ASSERT(it == TIteratorConst{});
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_increment_pre)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;

    auto it = begin(waveBand, Standard());

    SEQAN_ASSERT(*it == false);
    SEQAN_ASSERT(*++it == true);
    SEQAN_ASSERT(*++it == false);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_increment_post)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;

    auto it = begin(waveBand, Standard());

    SEQAN_ASSERT(*it++ == false);
    SEQAN_ASSERT(*it++ == true);
    SEQAN_ASSERT(*it == false);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_increment_offset)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;
    waveBand._grid[10] = true;

    auto it = begin(waveBand, Standard());

    SEQAN_ASSERT(*it == false);
    SEQAN_ASSERT(*(it + 1) == true);
    SEQAN_ASSERT(*(it + 10) == true);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_decrement_pre)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;

    auto it = begin(waveBand, Standard()) + 2;

    SEQAN_ASSERT(*it == false);
    SEQAN_ASSERT(*--it == true);
    SEQAN_ASSERT(*--it == false);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_decrement_post)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;

    auto it = begin(waveBand, Standard()) + 2;

    SEQAN_ASSERT(*it-- == false);
    SEQAN_ASSERT(*it-- == true);
    SEQAN_ASSERT(*it == false);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_decrement_offset)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;
    waveBand._grid[10] = true;

    auto it = begin(waveBand, Standard()) + 10;

    SEQAN_ASSERT(*it == true);
    SEQAN_ASSERT(*(it - 9) == true);
    SEQAN_ASSERT(*(it - 10) == false);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_difference)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;
    waveBand._grid[10] = true;

    auto itL = begin(waveBand, Standard());
    auto itR = begin(waveBand, Standard()) + 10;

    SEQAN_ASSERT(itR - itL == 10u);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_comparison)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};
    waveBand._grid[1] = true;
    waveBand._grid[10] = true;

    {
        auto itL = begin(waveBand, Standard());
        auto itR = begin(waveBand, Standard());

        SEQAN_ASSERT(itR == itL);
    }

    {
        auto itL = begin(waveBand, Standard());
        auto itR = begin(waveBand, Standard()) + 1;

        SEQAN_ASSERT(itR != itL);
    }
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_dereferencing)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());
    SEQAN_ASSERT(*it == false);

    waveBand._grid[0] = true;
    SEQAN_ASSERT(*it == true);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_container)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());
    SEQAN_ASSERT(std::addressof(container(it)) == std::addressof(waveBand));
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_position)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());
    SEQAN_ASSERT(position(it) == 0u);
    it += 10;
    SEQAN_ASSERT(position(it) == 10u);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_coordinate)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());

    impl::GridColumnIndex c;
    impl::GridRowIndex r;

    std::tie(c, r) = coordinates(it);
    SEQAN_ASSERT(c.get() == 0u);
    SEQAN_ASSERT(r.get() == 0u);

    it += 3;
    std::tie(c, r) = coordinates(it);
    SEQAN_ASSERT(c.get() == 0u);
    SEQAN_ASSERT(r.get() == 3u);

    it += 2;
    std::tie(c, r) = coordinates(it);
    SEQAN_ASSERT(c.get() == 1u);
    SEQAN_ASSERT(r.get() == 1u);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_column_index)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());

    impl::GridColumnIndex c = columnIndex(it);
    SEQAN_ASSERT(c.get() == 0u);

    it += 3;
    c = columnIndex(it);
    SEQAN_ASSERT(c.get() == 0u);

    it += 2;
    c = columnIndex(it);
    SEQAN_ASSERT(c.get() == 1u);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_row_index)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());

    impl::GridRowIndex r = rowIndex(it);
    SEQAN_ASSERT(r.get() == 0u);

    it += 3;
    r = rowIndex(it);
    SEQAN_ASSERT(r.get() == 3u);

    it += 2;
    r = rowIndex(it);
    SEQAN_ASSERT(r.get() == 1u);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_predecessor_left)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());
    *it = true;

    SEQAN_ASSERT(!hasPredecessorLeft(it));

    it += 3;
    SEQAN_ASSERT(!hasPredecessorLeft(it));

    it += 1;
    SEQAN_ASSERT(hasPredecessorLeft(it));

    it += 1;
    SEQAN_ASSERT(!hasPredecessorLeft(it));
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_predecessor_above)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = begin(waveBand, Standard());
    *it = true;

    SEQAN_ASSERT(!hasPredecessorAbove(it));

    it += 1;
    SEQAN_ASSERT(hasPredecessorAbove(it));

    it += 2;
    SEQAN_ASSERT(!hasPredecessorAbove(it));

    it += 1;
    SEQAN_ASSERT(!hasPredecessorAbove(it));
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_successor_right)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = end(waveBand, Standard()) - 1;
    *it = true;

    SEQAN_ASSERT(!hasSuccessorRight(it));
    it -= 3;
    *it = true;

    SEQAN_ASSERT(!hasSuccessorRight(it));
    --it;
    SEQAN_ASSERT(hasSuccessorRight(it));

    it -= 3;
    SEQAN_ASSERT(hasSuccessorRight(it));
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_has_successor_below)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = end(waveBand, Standard()) - 1;
    *it = true;

    SEQAN_ASSERT(!hasSuccessorBelow(it));
    it -= 1;

    SEQAN_ASSERT(hasSuccessorBelow(it));
    --it;
    SEQAN_ASSERT(!hasSuccessorBelow(it));
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_successor_right)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = end(waveBand, Standard()) - 1;
    *it = true;

    it -= 4;
    SEQAN_ASSERT(*successorRight(it) == true);
    --it;
    SEQAN_ASSERT(*successorRight(it) == false);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_iterator_successor_below)
{
    using namespace seqan;

    impl::DPWavefrontBand waveBand{impl::HorizontalBandPos{7}, impl::VerticalBandPos{4},
                                   impl::GridColumnSize{4}, impl::GridRowSize{5}};

    auto it = end(waveBand, Standard()) - 1;
    *it = true;
    --it;
    SEQAN_ASSERT(*successorBelow(it) == true);
    --it;
    SEQAN_ASSERT(*successorBelow(it) == false);
}

SEQAN_DEFINE_TEST(test_align_parallel_algorithm_banded_dp_wavefront_band_transform_to_grid)
{
    using namespace seqan;
                     //0  1  2  3  4  5  6
    CharString seqH = "Garfield the cat";
    CharString seqV = "Garfield the fat cat";

    { // w/o Band
        auto band = transformToGrid(seqH, seqV, 3, DPBandConfig<BandOff>{});
        for (auto it = begin(band, Standard{}); it != end(band, Standard{}); ++it)
            SEQAN_ASSERT(*it);
    }

    {  // w/ Band
        auto band = transformToGrid(seqH, seqV, 3, DPBandConfig<BandOn>{-4, 5});
        // Column 0
        SEQAN_ASSERT(band._grid[0]);
        SEQAN_ASSERT(band._grid[1]);
        SEQAN_ASSERT(band._grid[2]);
        SEQAN_ASSERT(!band._grid[3]);
        SEQAN_ASSERT(!band._grid[4]);
        SEQAN_ASSERT(!band._grid[5]);
        SEQAN_ASSERT(!band._grid[6]);

        // Column 1
        SEQAN_ASSERT(band._grid[7]);
        SEQAN_ASSERT(band._grid[8]);
        SEQAN_ASSERT(band._grid[9]);
        SEQAN_ASSERT(band._grid[10]);
        SEQAN_ASSERT(!band._grid[11]);
        SEQAN_ASSERT(!band._grid[12]);
        SEQAN_ASSERT(!band._grid[13]);

        //Column 2
        SEQAN_ASSERT(band._grid[14]);
        SEQAN_ASSERT(band._grid[15]);
        SEQAN_ASSERT(band._grid[16]);
        SEQAN_ASSERT(band._grid[17]);
        SEQAN_ASSERT(band._grid[18]);
        SEQAN_ASSERT(!band._grid[19]);
        SEQAN_ASSERT(!band._grid[20]);

        // Column 3
        SEQAN_ASSERT(!band._grid[21]);
        SEQAN_ASSERT(band._grid[22]);
        SEQAN_ASSERT(band._grid[23]);
        SEQAN_ASSERT(band._grid[24]);
        SEQAN_ASSERT(band._grid[25]);
        SEQAN_ASSERT(band._grid[26]);
        SEQAN_ASSERT(!band._grid[27]);

        // Column 4
        SEQAN_ASSERT(!band._grid[28]);
        SEQAN_ASSERT(!band._grid[29]);
        SEQAN_ASSERT(band._grid[30]);
        SEQAN_ASSERT(band._grid[31]);
        SEQAN_ASSERT(band._grid[32]);
        SEQAN_ASSERT(band._grid[33]);
        SEQAN_ASSERT(band._grid[34]);

        // Column 5
        SEQAN_ASSERT(!band._grid[35]);
        SEQAN_ASSERT(!band._grid[36]);
        SEQAN_ASSERT(!band._grid[37]);
        SEQAN_ASSERT(band._grid[38]);
        SEQAN_ASSERT(band._grid[39]);
        SEQAN_ASSERT(band._grid[40]);
        SEQAN_ASSERT(band._grid[41]);
    }
}
