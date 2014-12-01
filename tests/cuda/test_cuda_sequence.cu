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
#include <thrust/count.h>

#include "test_cuda_common.h"

using namespace seqan;

// ============================================================================
// Types
// ============================================================================

typedef TagList<String<char, Alloc<> >,
        TagList<String<Dna, Alloc<> >,
        TagList<String<Dna5, Alloc<> >
//        TagList<String<Dna, Packed<> >
        > > > //>
        StringTypes;

// TODO(esiragusa): test StringSets.
//typedef TagList<StringSet<CharString, Owner<ConcatDirect<> > >,
//        TagList<StringSet<DnaString, Owner<ConcatDirect<> > >
//        > >
//    TStringSetTypes;

// TODO(esiragusa): use metaprogramming algebra.
//typedef Product<StringTypes, Owner<ConcatDirect<> > >::Type TStringSetTypes;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class CudaSequenceTest
// ----------------------------------------------------------------------------

template <typename TType>
class CudaSequenceTest : public Test
{
public:
    typedef TType                             TString;
    typedef typename Device<TString>::Type    TCudaString;
    typedef typename Value<TString>::Type     TAlphabet;

    TString str;

    CudaSequenceTest() :
        str("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
    {}
};

SEQAN_TYPED_TEST_CASE(CudaSequenceTest, StringTypes);


// ============================================================================
// Tests
// ============================================================================

// ----------------------------------------------------------------------------
// Test assign()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(CudaSequenceTest, Assign)
{
    typedef typename TestFixture::TString       TString;
    typedef typename TestFixture::TCudaString   TCudaString;
    typedef typename TestFixture::TAlphabet     TAlphabet;

    cudaDeviceReset();

    TCudaString cudaStr;
    assign(cudaStr, this->str);
    SEQAN_ASSERT_EQ(length(cudaStr), length(this->str));
    SEQAN_ASSERT_EQ(thrust::count(begin(cudaStr, Standard()), end(cudaStr, Standard()), TAlphabet('A')), 10u);
    SEQAN_ASSERT_EQ(thrust::count(begin(cudaStr, Standard()), end(cudaStr, Standard()), TAlphabet('C')), 10u);
    SEQAN_ASSERT_EQ(thrust::count(begin(cudaStr, Standard()), end(cudaStr, Standard()), TAlphabet('G')), 10u);
    SEQAN_ASSERT_EQ(thrust::count(begin(cudaStr, Standard()), end(cudaStr, Standard()), TAlphabet('T')), 10u);

//    TString str;
//    assign(cudaStr, str);
//    SEQAN_ASSERT_EQ(str, this->str);
}

// ----------------------------------------------------------------------------
// Test getValue()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(CudaSequenceTest, GetValue)
{
    typedef typename TestFixture::TString       TString;
    typedef typename TestFixture::TCudaString   TCudaString;
    typedef typename View<TCudaString>::Type    TCudaStringView;
    typedef typename Size<TString>::Type        TSize;

    cudaDeviceReset();

    TCudaString cudaStr;
    assign(cudaStr, this->str);
    TCudaStringView cudaStrView = view(cudaStr);

    for (TSize pos = 0; pos < length(this->str); pos++)
    {
        testGetValue<<<1,1>>>(cudaStrView, pos, getValue(this->str, pos));
        cudaDeviceSynchronize();
        SEQAN_ASSERT_EQ(cudaGetLastError(), cudaSuccess);
    }
}

// ============================================================================
// Register Tests
// ============================================================================

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
