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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements test for the find methods.
// ==========================================================================
#ifndef EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_TRAVERSE_H_
#define EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_TRAVERSE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/journaled_string_tree.h>

#include "test_journaled_string_tree_mock_generator.h"

// Global variable to read the config file.
seqan::DataParallelTestConfig<char> testConfig;

namespace seqan
{

template <typename TContainer>
struct DummyCaller_
{
    TContainer * _containerPtr;

    DummyCaller_(TContainer & other) : _containerPtr(&other)
    {}
};

//template <typename T> SEQAN_CONCEPT_IMPL(DummyCaller_<T>,       (JstTraversalConcept));
//template <typename T> SEQAN_CONCEPT_IMPL(DummyCaller_<T> const, (JstTraversalConcept));

template <typename TValue>
struct DummyDelegator_
{

    StringSet<String<TValue> > _processedSeq;

    DummyDelegator_(unsigned size)
    {
        resize(_processedSeq, size, Exact());
    }

    template <typename TTraverser>
    void operator()(TTraverser & traverser)
    {
        typedef typename Positions<TTraverser>::Type TPosVec;
#ifdef TEST_DEBUG_OUTPUT
        std::cerr << "position(traverser): ";
#endif
        TPosVec posVec = positions(traverser);
        for (unsigned i = 0; i < length(posVec); ++i)
        {
#ifdef TEST_DEBUG_OUTPUT
            std::cerr << "("<< posVec[i].i1 << ", " <<  posVec[i].i2 << ")" << "; ";
#endif
            appendValue(_processedSeq[posVec[i].i1],
                        value(stringSet(container(traverser)), posVec[i].i1)[posVec[i].i2]);
        }
#ifdef TEST_DEBUG_OUTPUT
        std::cerr << std::endl;
#endif
    }
};

template <typename TContainer>
struct GetState<DummyCaller_<TContainer> >
{
    typedef Nothing Type;
};

template <typename TContainer>
struct GetState<DummyCaller_<TContainer> const>
{
    typedef Nothing const Type;
};

template <typename TContainer>
struct GetJstTraverser<DummyCaller_<TContainer> >
{
    typedef DummyCaller_<TContainer> TDummy_;
    typedef typename GetState<TDummy_>::Type TState;
    typedef JstTraverser<TContainer, TState, JstTraverserSpec<> > Type;
};

template <typename TDumyContainer, typename TValue, typename TContainer, typename TState, typename TSpec, typename TTag>
inline typename Size<JstTraverser<TContainer, TState, TSpec> >::Type
deliverContext(DummyCaller_<TDumyContainer> & /*dummy*/,
               DummyDelegator_<TValue> & delegator,
               JstTraverser<TContainer, TState, TSpec> & traverser,
               TTag const & /*tag*/)
{
    delegator(traverser);
    return 1;
}

template <typename TContainer>
inline Nothing
getState(DummyCaller_<TContainer> & /*dummy*/)
{
    return Nothing();
}

template <typename TContainer>
inline Nothing const
getState(DummyCaller_<TContainer const> & /*dummy*/)
{
    return Nothing();
}

template <typename TContainer, typename TState>
inline void
setState(DummyCaller_<TContainer> & /*dummy*/,
         TState const & /*state*/)
{
    // no-op.
}

}

template <typename TTestSeq, typename TCompareSeq, typename TSize>
bool compareResults(TTestSeq const & testSeq, TCompareSeq const & compSeq, TSize const & windowSize)
{
    SEQAN_ASSERT_EQ(length(testSeq), length(compSeq));
    for (unsigned i = 0; i < length(testSeq); ++i)
        if (isNotEqual(testSeq[i], prefix(compSeq[i], length(compSeq[i]) - (windowSize -1))))
            return false;
    return true;
}


template <typename TMock, typename TTester, typename TSize>
void _printDebugInfo(TMock const & mockGen, TTester const & dpTester, TSize const & windowSize)
{
    std::cerr << "Host: " << host(mockGen._seqData) << std::endl;
    for (unsigned i = 0; i < length(dpTester._processedSeq); ++i)
        std::cerr << "Compare 1: " << dpTester._processedSeq[i] << "\nCompare 2: "
        << prefix(mockGen._seqData[i], length(mockGen._seqData[i]) - (windowSize - 1)) << "\n" << std::endl;
}

bool _runTestForConfiguration(unsigned posConf,
                              unsigned varConf,
                              unsigned covConf,
                              unsigned refLength,
                              unsigned windowSize,
                              seqan::StringTreeDefault const & /*stringTreeTag*/)
{
    using namespace seqan;

    typedef DummyDelegator_<char> TSequenceAppender;

    typedef String<MockVariantData<char> > TVarData;
    typedef String<String<bool, Packed<> > > TCovData;
    typedef MockGenerator_<unsigned, char> TMockGenerator;

    typedef TMockGenerator::TStringTree TStringTree;
    typedef JstTraverser<TStringTree, Nothing, JstTraverserSpec<> > TTraverser;
    typedef DummyCaller_<TStringTree> TDummyCaller;

    TVarData varData;
    TCovData covData;
    testConfig.getTestConfiguration(varData, covData, posConf, varConf, covConf);

    // Initialize the mock generator.
    TMockGenerator mockGen;
    // Generate the mock for the current configuration.
    mockGen.generate(varData, covData, refLength);

    TStringTree jst(host(mockGen._seqData), mockGen._varStore);
    journalNextBlock(jst, windowSize);

    TSequenceAppender seqAppender(length(mockGen._seqData));
    TTraverser traverser(jst, windowSize);
    TDummyCaller dummyCaller(jst);

    traverse(dummyCaller, seqAppender, traverser);

#ifdef TEST_DEBUG_OUTPUT
    _printDebugInfo(mockGen, seqAppender, windowSize);
#endif

    return compareResults(seqAppender._processedSeq, mockGen._seqData, windowSize);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_jst_traversal_concept)
{
    using namespace seqan;

    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    typedef DummyCaller_<TJst> TDummy;

    bool res =  Is<JstTraversalConcept<TDummy> >::VALUE;
    SEQAN_ASSERT_EQ(res, true);
}

// ----------------------------------------------------------------------------
// Test all SNPs.
// ----------------------------------------------------------------------------

// Test all at position 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 0, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 0, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 0, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 0, 5, 101, 30, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_0_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 0, 5, 101, 30, seqan::StringTreeDefault()));
}

// ----------------------------------------------------------------------------
// Test all deletions.
// ----------------------------------------------------------------------------

// Test all beginning at position 0, all dels, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_1_5_journaled_string_tree)
{
//    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_1_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 1, 5, 101, 50, seqan::StringTreeDefault()));
}

// ----------------------------------------------------------------------------
// Test all Insertions.
// ----------------------------------------------------------------------------

// Test all at position 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_0_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(0, 2, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_1_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(1, 2, 5, 101, 10, seqan::StringTreeDefault()));
}

// Test all at position 30, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 0, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 1, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 2, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 3, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 4, 101, 50, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_2_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(2, 2, 5, 101, 50, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_3_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(3, 2, 5, 101, 30, seqan::StringTreeDefault()));
}

// Test different positions including 0, all snps, different coverages.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_2_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 0, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 0, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_2_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 1, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 1, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_2_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 2, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 2, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_2_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 3, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 3, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_2_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 4, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 4, 101, 30, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_4_2_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 5, 101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(4, 2, 5, 101, 30, seqan::StringTreeDefault()));
}

// ----------------------------------------------------------------------------
// Test special variant comobinations.
// ----------------------------------------------------------------------------

// Test all deletions with different size.

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_0_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 0, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 0, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_2_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 2, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 2, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_5_3_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(5, 3, 5, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_6_4_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(6, 4, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(6, 4, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_7_4_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(7, 4, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(7, 4, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_1_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 1, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 1, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_3_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 3, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 3, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_4_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 4, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 4, 101, 10, seqan::StringTreeDefault()));
}

SEQAN_DEFINE_TEST(test_journaled_journaled_string_tree_traverse_config_9_5_5_journaled_string_tree)
{
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 5, 101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(9, 5, 5, 101, 10, seqan::StringTreeDefault()));
}

#endif  // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_TRAVERSE_H_
