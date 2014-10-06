// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Test cases for the JST traversal.
// ==========================================================================

#define DISABLE_COMMON_2
//#define TEST_DEBUG_OUTPUT

#ifndef EXTRAS_TESTS_TEST_JOURNALED_STRING_TREE_TRAVERSAL_2_H_
#define EXTRAS_TESTS_TEST_JOURNALED_STRING_TREE_TRAVERSAL_2_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/journaled_string_tree.h>

#include "test_journaled_string_tree_mock_generator.h"

// ============================================================================
// Some test specific classes and helpers.
// ============================================================================

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

template <typename T> SEQAN_CONCEPT_IMPL((JstTraversalConcept), DummyCaller_<T>);
template <typename T> SEQAN_CONCEPT_IMPL((JstTraversalConcept), DummyCaller_<T> const);

template <typename TPair>
struct CompareLessFunctor_
{
    CompareLessFunctor_()
    {}

    inline bool
    operator()(TPair const & lhs, TPair const & rhs)
    {
        return lhs.i1 < rhs.i1;
    }
};

template <typename TValue, typename TTestSeq>
struct DummyDelegator_
{
    typedef Pair<unsigned, TValue> TPair;

    StringSet<String<TPair> > _processedSeq;
    TTestSeq                  _testSeq;

    DummyDelegator_(unsigned size, TTestSeq const & testSeq)
    {
        resize(_processedSeq, size, Exact());
        _testSeq = testSeq;
    }

    template <typename TTraverser>
    void operator()(TTraverser & traverser)
    {
        typedef typename Position<TTraverser>::Type TPosVec;
        SEQAN_OMP_PRAGMA(critical(insert))
        {
#ifdef TEST_DEBUG_OUTPUT
            printf("Thread: %i position(traverser): ", omp_get_thread_num());
#endif
            TPosVec posVec = position(traverser);
            for (unsigned i = 0; i < length(posVec); ++i)
            {
#ifdef TEST_DEBUG_OUTPUT
                printf("(%lu, %lu) ",posVec[i].i1, posVec[i].i2);
#endif
                appendValue(_processedSeq[posVec[i].i1],
                            TPair(posVec[i].i2, value(_testSeq, posVec[i].i1)[posVec[i].i2]));
            }

#ifdef TEST_DEBUG_OUTPUT
            printf("\n");
#endif
        }
    }

    void postProcess()
    {
        for (unsigned i = 0; i < length(_processedSeq); ++i)
            sort(_processedSeq[i], CompareLessFunctor_<TPair>());
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
    typedef JstTraverser<TContainer, TState, JstTraverserConfig<> > Type;
};

template <typename TDumyContainer, typename TDelegator, typename TContainer, typename TState, typename TSpec, typename TTag>
inline typename Size<JstTraverser<TContainer, TState, TSpec> >::Type
notify(DummyCaller_<TDumyContainer> & /*dummy*/,
       TDelegator & delegator,
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

template <typename TContainer>
inline void
initState(DummyCaller_<TContainer> & /*dummy*/)
{
    // no-op.
}

}

// ============================================================================
// Defines the typed test classes.
// ============================================================================

using namespace seqan;

template <typename TJstTraversalSpecConfig_, typename TConfig_>
struct GlobalJstConfig
{
    typedef typename TJstTraversalSpecConfig_::TContextPosition TContextPosition;
    typedef typename TJstTraversalSpecConfig_::TFullContext TFullContext;

    static const unsigned DELTA_CONFIG = TConfig_::DELTA_CONFIG;
    static const unsigned POSITION_CONFIG = TConfig_::POSITION_CONFIG;
    static const unsigned COVERAGE_CONFIG = TConfig_::COVERAGE_CONFIG;
};


template <typename TContextPosition_, typename TFullContext_>
struct JstTraversalSpecConfig
{
    typedef TContextPosition_ TContextPosition;
    typedef TFullContext_ TFullContext;
};

// General configuration struct.
template <unsigned _DELTA_CONFIG, unsigned _POSITION_CONFIG, unsigned _COVERAGE_CONFIG>
struct JstTraversalConfig
{
    static const unsigned DELTA_CONFIG = _DELTA_CONFIG;
    static const unsigned POSITION_CONFIG = _POSITION_CONFIG;
    static const unsigned COVERAGE_CONFIG = _COVERAGE_CONFIG;
};

// Some nice helper to save a lot of boiler plate code.
template <unsigned I1, unsigned I2, unsigned I3, unsigned V1, unsigned V2, unsigned V3>
struct JstConfigGenRecursive_
{
    typedef TagList<JstTraversalConfig<V1, V2, V3>,
            typename JstConfigGenRecursive_<I1, I2, I3, V1, V2, V3 - 1>::Type > Type;
};

template <unsigned I1, unsigned I2, unsigned I3>
struct JstConfigGenRecursive_<I1, I2, I3, 0, 0, 0>
{
    typedef TagList<JstTraversalConfig<0, 0, 0> > Type;
};

template <unsigned I1, unsigned I2, unsigned I3, unsigned V1>
struct JstConfigGenRecursive_<I1, I2, I3, V1, 0, 0>
{
    typedef TagList<JstTraversalConfig<V1, 0, 0>,
            typename JstConfigGenRecursive_<I1, I2, I3, V1 - 1, I2, I3>::Type > Type;
};

template <unsigned I1, unsigned I2, unsigned I3, unsigned V1, unsigned V2>
struct JstConfigGenRecursive_<I1, I2, I3, V1, V2, 0>
{
    typedef TagList<JstTraversalConfig<V1, V2, 0>,
            typename JstConfigGenRecursive_<I1, I2, I3, V1, V2 - 1, I3>::Type > Type;
};

template <unsigned I1, unsigned I2, unsigned I3>
struct JstConfigGenerator_
{
    typedef typename JstConfigGenRecursive_<I1, I2, I3, I1, I2, I3>::Type Type;
};

// Recursively build test types.
template <typename TTagList1, typename TTagList2, typename TOriginalList2>
struct TestTypeBuilderHelperRecursive_;

// Recursion anchor.
template <typename TTag1, typename TTag2, typename TOriginal>
struct TestTypeBuilderHelperRecursive_<TagList<TTag1, void>, TagList<TTag2, void>, TOriginal>
{
    typedef TagList<GlobalJstConfig<TTag1, TTag2> > Type;
};

// Recursion anchor.
template <typename TTag1, typename TSubList1, typename TTag2, typename TOriginal>
struct TestTypeBuilderHelperRecursive_<TagList<TTag1, TSubList1>, TagList<TTag2, void>, TOriginal>
{
    typedef TagList<GlobalJstConfig<TTag1, TTag2>,
                    typename TestTypeBuilderHelperRecursive_<TSubList1, TOriginal, TOriginal>::Type
                    > Type;
};

// Recursion.
template <typename TTag1, typename TSubList1, typename TTag2, typename TSubList2, typename TOriginal>
struct TestTypeBuilderHelperRecursive_<TagList<TTag1, TSubList1>, TagList<TTag2, TSubList2>, TOriginal>
{
    typedef TagList<GlobalJstConfig<TTag1, TTag2>,
                    typename TestTypeBuilderHelperRecursive_<TagList<TTag1, TSubList1>, TSubList2, TOriginal>::Type
                    > Type;
};

// Test Type Builder.
template <typename TTagList1, typename TTagList2>
struct TestTypeBuilderHelper_
{
    typedef typename TestTypeBuilderHelperRecursive_<TTagList1, TTagList2, TTagList2>::Type Type;
};


// Define test class fixture substitution.
template <typename TConfig_>
class JstTraversalTest : public Test
{
public:
    typedef TConfig_ TConfig;
};

// Define common test class.
template <typename T>
class JstTraversalTestCommon : public JstTraversalTest<T>
{};

template <typename T>
class JstTraversalTestCommon2 : public JstTraversalTest<T>
{};

typedef
    TagList<JstTraversalConfig<5, 4, 0>,
    TagList<JstTraversalConfig<5, 4, 1>,
    TagList<JstTraversalConfig<5, 4, 2>,
    TagList<JstTraversalConfig<5, 4, 3>,
    TagList<JstTraversalConfig<5, 4, 4>,
    TagList<JstTraversalConfig<5, 4, 5>,
    TagList<JstTraversalConfig<6, 5, 3> ,
    TagList<JstTraversalConfig<7, 5, 3>,
    TagList<JstTraversalConfig<9, 6, 1>,
    TagList<JstTraversalConfig<9, 6, 3>,
    TagList<JstTraversalConfig<9, 6, 4>,
    TagList<JstTraversalConfig<9, 6, 5>,
    JstConfigGenerator_<4, 3, 5>::Type
    > > > > > > > > > > > >
    JstTestCases;

typedef TagList<JstTraversalSpecConfig<ContextPositionLeft, True> > JstTraversalSpec1;
typedef TestTypeBuilderHelper_<JstTraversalSpec1, JstTestCases>::Type JstConfigTypes1;

typedef TagList<JstTraversalSpecConfig<ContextPositionRight, False> > JstTraversalSpec2;
typedef TestTypeBuilderHelper_<JstTraversalSpec2, JstTestCases>::Type JstConfigTypes2;

// Define typed test specific configuration types.
#ifndef DISABLE_COMMON_1
SEQAN_TYPED_TEST_CASE(JstTraversalTestCommon, JstConfigTypes1);
#endif

#ifndef DISABLE_COMMON_2
SEQAN_TYPED_TEST_CASE(JstTraversalTestCommon2, JstConfigTypes2);
#endif

// ============================================================================
// Functions to run the actual tests.
// ============================================================================

template <typename TTestSeq, typename TCompareSeq, typename TSize>
bool compareResults(TTestSeq const & testSeq, TCompareSeq const & compSeq, TSize const & /*windowSize*/)
{
    SEQAN_ASSERT_EQ(length(testSeq), length(compSeq));
    for (unsigned i = 0; i < length(testSeq); ++i)
        for (unsigned j = 0; j < length(testSeq[i]); ++j)
        {
            if (testSeq[i][j].i2 != compSeq[i][j])
                return false;
        }
//        if (isNotEqual(testSeq[i].i2, prefix(compSeq[i], length(compSeq[i]) - (windowSize -1))))
//            return false;
    return true;
}

template <typename TMock, typename TTester, typename TSize>
void _printDebugInfo(TMock const & mockGen, TTester const & dpTester, TSize const & windowSize)
{
    std::cerr << "Host: " << host(mockGen._seqData) << std::endl;
    for (unsigned i = 0; i < length(dpTester._processedSeq); ++i)
    {
        std::cerr << "Traversed: ";
        for (unsigned j = 0; j < length(dpTester._processedSeq[i]); ++j)
           std::cerr << dpTester._processedSeq[i][j].i2;
        std::cerr << "\nGenerated: " << prefix(mockGen._seqData[i], length(mockGen._seqData[i]) - (windowSize - 1)) <<
                     "\n" << std::endl;
    }
}

// With block size and parallel tag.
template <typename TContextPos, typename TFullContext, typename TParallel>
bool _runTestForConfiguration(unsigned posConf,
                              unsigned varConf,
                              unsigned covConf,
                              TContextPos,
                              TFullContext,
                              unsigned refLength,
                              unsigned windowSize,
                              unsigned blockSize,
                              seqan::StringTreeDefault const & /*stringTreeTag*/,
                              seqan::Tag<TParallel>  const & parallelTag)
{
    using namespace seqan;

    typedef String<MockVariantData<char> > TVarData;
    typedef String<String<bool, Packed<> > > TCovData;
    typedef MockGenerator_<unsigned, char> TMockGenerator;

    typedef typename TMockGenerator::TStringTree TStringTree;
    typedef JstTraverser<TStringTree, Nothing, JstTraverserConfig<TContextPos, TFullContext> > TTraverser;
    typedef DummyCaller_<TStringTree> TDummyCaller;
    typedef GetStringSet<TStringTree>::Type TJournalSet;
    typedef DummyDelegator_<char, TJournalSet> TSequenceAppender;

    TVarData varData;
    TCovData covData;
    testConfig.getTestConfiguration(varData, covData, posConf, varConf, covConf);

    // Initialize the mock generator.
    TMockGenerator mockGen;
    // Generate the mock for the current configuration.
    mockGen.generate(varData, covData, refLength);

    TStringTree jst(host(mockGen._seqData), mockGen._varStore);
    if (blockSize > 0)
        setBlockSize(jst, blockSize);

//    TSequenceAppender seqAppender(length(mockGen._seqData));
    TSequenceAppender seqAppender(length(mockGen._seqData), mockGen._seqData);
    TTraverser traverser(jst, windowSize);
    TDummyCaller dummyCaller(jst);

    traverse(dummyCaller, seqAppender, traverser, parallelTag);

    seqAppender.postProcess();
    bool res = compareResults(seqAppender._processedSeq, mockGen._seqData, windowSize);

#ifdef TEST_DEBUG_OUTPUT_RES
    if (!res)
        _printDebugInfo(mockGen, seqAppender, windowSize);
#endif
    return res;
}

// with block size + no parallel tag.
template <typename TContextPos, typename TFullContext>
inline bool _runTestForConfiguration(unsigned posConf,
                                     unsigned varConf,
                                     unsigned covConf,
                                     TContextPos contextPos,
                                     TFullContext fullContext,
                                     unsigned refLength,
                                     unsigned windowSize,
                                     unsigned blockSize,
                                     seqan::StringTreeDefault const & stringTreeTag)
{
    return _runTestForConfiguration(posConf, varConf, covConf, contextPos, fullContext, refLength, windowSize,
                                    blockSize, stringTreeTag, seqan::Serial());
}

// no block size + parallel tag.
template <typename TContextPos, typename TFullContext, typename TParallel>
inline bool _runTestForConfiguration(unsigned posConf,
                                     unsigned varConf,
                                     unsigned covConf,
                                     TContextPos contextPos,
                                     TFullContext fullContext,
                                     unsigned refLength,
                                     unsigned windowSize,
                                     seqan::StringTreeDefault const & stringTreeTag,
                                     seqan::Tag<TParallel> const & parallelTag)
{
    return _runTestForConfiguration(posConf, varConf, covConf, contextPos, fullContext, refLength, windowSize, 0,
                                    stringTreeTag, parallelTag);
}

// no block size + no parallel tag.
template <typename TContextPos, typename TFullContext>
inline bool _runTestForConfiguration(unsigned posConf,
                                     unsigned varConf,
                                     unsigned covConf,
                                     TContextPos contextPos,
                                     TFullContext fullContext,
                                     unsigned refLength,
                                     unsigned windowSize,
                                     seqan::StringTreeDefault const & stringTreeTag)
{
    return _runTestForConfiguration(posConf, varConf, covConf, contextPos, fullContext, refLength, windowSize, 0,
                                    stringTreeTag);
}

// ============================================================================
// Test full JST generation in serial mode.
// ============================================================================

template <typename TConfig>
inline void _testJstTraversalCommonFullSerial(TConfig /*config*/)
{
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 3, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 20, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 30, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 50, seqan::StringTreeDefault()));
}

#ifndef DISABLE_COMMON_1
SEQAN_TYPED_TEST(JstTraversalTestCommon, FullSerial)
{
    _testJstTraversalCommonFullSerial(typename TestFixture::TConfig());
}
#endif

#ifndef DISABLE_COMMON_2
SEQAN_TYPED_TEST(JstTraversalTestCommon2, FullSerial)
{
    _testJstTraversalCommonFullSerial(typename TestFixture::TConfig());
}
#endif

// ============================================================================
// Test block-wise JST generation in serial mode.
// ============================================================================

template <typename TConfig>
inline void _testJstTraversalCommonBlockSerial(TConfig /*config*/)
{
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 3, 2, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 20, 2, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 30, 2, seqan::StringTreeDefault()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 50, 2, seqan::StringTreeDefault()));
}

#ifndef DISABLE_COMMON_1
SEQAN_TYPED_TEST(JstTraversalTestCommon, BlockSerial)
{
    _testJstTraversalCommonBlockSerial(typename TestFixture::TConfig());
}
#endif

#ifndef DISABLE_COMMON_2
SEQAN_TYPED_TEST(JstTraversalTestCommon2, BlockSerial)
{
    _testJstTraversalCommonBlockSerial(typename TestFixture::TConfig());
}
#endif

// ============================================================================
// Test full JST generation in parallel mode.
// ============================================================================

template <typename TConfig>
inline void _testJstTraversalCommonFullParallel(TConfig /*config*/)
{
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 3, seqan::StringTreeDefault(), seqan::Parallel()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 20, seqan::StringTreeDefault(), seqan::Parallel()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 30, seqan::StringTreeDefault(), seqan::Parallel()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 50, seqan::StringTreeDefault(), seqan::Parallel()));
}

#ifndef DISABLE_COMMON_1
SEQAN_TYPED_TEST(JstTraversalTestCommon, FullParallel)
{
    _testJstTraversalCommonFullParallel(typename TestFixture::TConfig());
}
#endif

#ifndef DISABLE_COMMON_2
SEQAN_TYPED_TEST(JstTraversalTestCommon2, FullParallel)
{
    _testJstTraversalCommonFullParallel(typename TestFixture::TConfig());
}
#endif

// ============================================================================
// Test block-wise JST generation in parallel mode.
// ============================================================================

template <typename TConfig>
inline void _testJstTraversalCommonBlockParallel(TConfig /*config*/)
{
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 3, 2, seqan::StringTreeDefault(), seqan::Parallel()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 20, 2, seqan::StringTreeDefault(), seqan::Parallel()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 30, 2, seqan::StringTreeDefault(), seqan::Parallel()));
    SEQAN_ASSERT(_runTestForConfiguration(TConfig::DELTA_CONFIG, TConfig::POSITION_CONFIG, TConfig::COVERAGE_CONFIG,
                                          typename TConfig::TContextPosition(), typename TConfig::TFullContext(),
                                          101, 50, 2, seqan::StringTreeDefault(), seqan::Parallel()));
}

#ifndef DISABLE_COMMON_1
SEQAN_TYPED_TEST(JstTraversalTestCommon, BlockParallel)
{
    _testJstTraversalCommonBlockParallel(typename TestFixture::TConfig());
}
#endif

#ifndef DISABLE_COMMON_2
SEQAN_TYPED_TEST(JstTraversalTestCommon2, BlockParallel)
{
    _testJstTraversalCommonBlockParallel(typename TestFixture::TConfig());
}
#endif

#endif  // EXTRAS_TESTS_TEST_JOURNALED_STRING_TREE_TRAVERSAL_2_H_
