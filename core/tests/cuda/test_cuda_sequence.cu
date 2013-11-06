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

using namespace seqan;

// ============================================================================
// Types
// ============================================================================

typedef TagList<String<char, Alloc<> >,
        TagList<String<Dna, Alloc<> >,
        TagList<String<Dna5, Alloc<> >,
        TagList<String<Dna, Packed<> >
        > > > >
        StringTypes;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class CudaSequenceTest
// ----------------------------------------------------------------------------

template <typename TString>
class CudaSequenceTest : public Test
{
public:
    typedef TString                           Type;
    typedef typename Device<TString>::Type    TDeviceType;

    TString str;

    CudaSequenceTest() :
        str("ACGT")
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
    typedef typename TestFixture::Type          TString;
    typedef typename TestFixture::TDeviceType   TDeviceString;

    TString       str;
    TDeviceString deviceStr;

    assign(deviceStr, this->str);

    SEQAN_ASSERT_EQ(length(deviceStr), length(this->str));

//    assign(str, deviceStr);
//    SEQAN_ASSERT(str == this->str);
}

// ============================================================================
// Register Tests
// ============================================================================

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
