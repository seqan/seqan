// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/index/index_fm_device.h>

#include "test_cuda_common.h"

using namespace seqan;

// TODO(esiragusa): move this into metaprogramming algebra
namespace seqan {
template <typename T1, typename T2>
struct Pair<T1, T2, Tag<void> > {};
}

// ============================================================================
// Types
// ============================================================================

typedef TagList<FibreRawText,
        TagList<FibreLF
        > >
    FMIndexFibres;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class CudaIndexTest
// ----------------------------------------------------------------------------

template <typename TType>
class CudaIndexTest : public Test
{
public:
    typedef TType                                       TIndex;
    typedef typename Host<TIndex>::Type                 TText;
    typedef typename Device<TIndex>::Type               TCudaIndex;

    TText    text;
    TIndex   index;

    CudaIndexTest() :
        text(),
        index(text)
    {
        // TODO(esiragusa): init generic text.
        appendValue(text, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

        // TODO(esiragusa): reverse text on FM-index only.
        reverse(text);
        indexCreate(index);
        reverse(text);
    }
};

typedef TagList<DnaStringSetFMIndex> CudaIndexTestTypes;

SEQAN_TYPED_TEST_CASE(CudaIndexTest, CudaIndexTestTypes);

// ----------------------------------------------------------------------------
// Class CudaIndexFibreTest
// ----------------------------------------------------------------------------

template <typename TTypes>
class CudaIndexFibreTest : public CudaIndexTest<typename Value<TTypes, 1>::Type> {};

// TODO(esiragusa): use metaprogramming algebra.
//typedef Product<DnaStringSetFMIndex, FMIndexFibres>::Type CudaIndexFibreTestTypes;

typedef TagList<Pair<DnaStringSetFMIndex, FibreRawText, Tag<void> >,
        TagList<Pair<DnaStringSetFMIndex, FibreLF, Tag<void> >
        > >
    CudaIndexFibreTestTypes;

SEQAN_TYPED_TEST_CASE(CudaIndexFibreTest, CudaIndexFibreTestTypes);

// ----------------------------------------------------------------------------
// Class CudaIndexCountTest
// ----------------------------------------------------------------------------

template <typename TTypes>
class CudaIndexCountTest : public CudaIndexTest<typename Value<TTypes, 1>::Type>
{
public:
    typedef typename Value<TTypes, 1>::Type             TIndex;
    typedef typename Value<TTypes, 2>::Type             TNeedles;
    typedef CudaIndexTest<TIndex>                       TBase;
    typedef typename Size<TIndex>::Type                 TSize;
    typedef typename Device<TNeedles>::Type             TCudaNeedles;

    TNeedles needles;
    TSize    occurrences;

    CudaIndexCountTest() :
        TBase()
    {
        // TODO(esiragusa): append generic needles.
        appendValue(needles, "ACGT");
        appendValue(needles, "CGT");
        appendValue(needles, "GTA");

        occurrences = countOccurrences(this->index, needles);
    }
};

// TODO(esiragusa): use metaprogramming algebra.
//typedef Product<DnaStringSetFMIndex, DnaStringSet>::Type CudaIndexCountTestTypes;

typedef TagList<Pair<DnaStringSetFMIndex, DnaStringSet, Tag<void> > > CudaIndexCountTestTypes;

SEQAN_TYPED_TEST_CASE(CudaIndexCountTest, CudaIndexCountTestTypes);

// ============================================================================
// Tests
// ============================================================================

// ----------------------------------------------------------------------------
// Test assign()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(CudaIndexTest, Assign)
{
    typedef typename TestFixture::TIndex        TIndex;
    typedef typename TestFixture::TCudaIndex    TCudaIndex;

    cudaDeviceReset();

    TCudaIndex cudaIndex;
    assign(cudaIndex, this->index);
    SEQAN_ASSERT_EQ(length(cudaIndex), length(this->index));

//    TIndex        index;
//    assign(index, cudaIndex);
//    SEQAN_ASSERT(index == this->index);
}

// ----------------------------------------------------------------------------
// Test value() on Index Fibres
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(CudaIndexFibreTest, Values)
{
    typedef FibreLF                                 TTag;
    typedef typename TestFixture::TIndex            TIndex;
    typedef typename TestFixture::TCudaIndex        TCudaIndex;
    typedef typename Fibre<TIndex, TTag>::Type      TFibre;
    typedef typename Fibre<TCudaIndex, TTag>::Type  TCudaFibre;
    typedef typename View<TCudaFibre>::Type         TCudaFibreView;
    typedef typename Size<TFibre>::Type             TSize;

    cudaDeviceReset();

    TCudaIndex cudaIndex;
    assign(cudaIndex, this->index);

    TFibre & fibre = getFibre(this->index, TTag());
    TCudaFibre & cudaFibre = getFibre(cudaIndex, TTag());
    SEQAN_ASSERT_EQ(length(fibre), length(cudaFibre));

    TCudaFibreView cudaFibreView = view(cudaFibre);
    for (TSize pos = 0; pos < length(fibre); pos++)
    {
        testGetValue<<<1,1>>>(cudaFibreView, pos, fibre[pos]);
        cudaDeviceSynchronize();
        SEQAN_ASSERT_EQ(cudaGetLastError(), cudaSuccess);
    }
}

// ----------------------------------------------------------------------------
// Test countOccurrences()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(CudaIndexCountTest, Count)
{
    typedef typename TestFixture::TCudaIndex    TCudaIndex;
    typedef typename TestFixture::TCudaNeedles  TCudaNeedles;

    cudaDeviceReset();

    TCudaIndex   cudaIndex;
    TCudaNeedles cudaNeedles;

    assign(cudaIndex, this->index);
    assign(cudaNeedles, this->needles);

    SEQAN_ASSERT_EQ(countOccurrences(cudaIndex, cudaNeedles), this->occurrences);
}

// ============================================================================
// Register Tests
// ============================================================================

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
