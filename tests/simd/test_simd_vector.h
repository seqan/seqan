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
// Author: David Weese <david.weese@fu-berlin.de>
//         Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// ==========================================================================
// Tests for SIMD vectors.
// ==========================================================================

#ifndef SEQAN_CORE_TESTS_SIMD_TEST_SIMD_VECTOR_H_
#define SEQAN_CORE_TESTS_SIMD_TEST_SIMD_VECTOR_H_

#include <seqan/simd.h>
#include <random>
#include <seqan/sequence.h>
#include <seqan/misc/bit_twiddling.h>

#if defined(SEQAN_SIMD_ENABLED)
namespace seqan {

template <int ROWS, typename TVector>
inline void test_matrix_transpose()
{
    typedef typename Value<TVector>::Type TValue;
    typedef TVector TMatrix[LENGTH<TVector>::VALUE];
    const int COLS = LENGTH<TVector>::VALUE;

    String<TValue> random;
    resize(random, ROWS * COLS);

    std::mt19937 rng;
    // http://stackoverflow.com/questions/31460733/why-arent-stduniform-int-distributionuint8-t-and-stduniform-int-distri
    std::uniform_int_distribution<uint64_t> pdf(0, MaxValue<TValue>::VALUE);
    for (unsigned i = 0; i < length(random); ++i)
        random[i] = static_cast<TValue>(pdf(rng));

    TMatrix tmp;
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            tmp[i][j] = random[i * COLS + j];

//    for(int i=0;i<ROWS;++i)
//        print(std::cout, tmp[i]) << std::endl;

    transpose<ROWS>(tmp);

//    std::cout << std::endl;
//    std::cout << std::endl;
//    for(int i=0;i<DIM;++i)
//        print(std::cout, tmp[i]) << std::endl;
#if defined(__x86_64__) || defined(__amd64__)
    _mm_empty();  // Fixes icpc warning #13203: No EMMS instruction before call to function
#endif // defined(__x86_64__) || defined(__amd64__)
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            SEQAN_ASSERT_EQ(tmp[i][j], random[j * ROWS + i]);
}

template <typename TSimdVector>
void fillVectors(TSimdVector & a, TSimdVector & b)
{
    using namespace seqan;
    constexpr auto length = LENGTH<TSimdVector>::VALUE;

    for (auto i = 0; i < length; ++i)
    {
        a[i] = (i - 1) * 3;
        b[i] = length - i;
    }
}

template <typename TSimdVector, typename TSize>
void reverseIndexSequence(TSimdVector & idx, TSize length)
{
    for (auto i = 0; i < length; ++i)
    {
        // note: umesimd swizzle interface has no a[i] = i; support.
        assignValue(idx, i, length - i - 1);
    }
}

template <typename TSimdVector>
constexpr auto trueValue()
{
    using TSimdMaskVector = typename SimdMaskVector<TSimdVector>::Type;
    using TValue = typename Value<TSimdMaskVector>::Type;
    return static_cast<TValue>(-1);
}

} // namespace seqan

// ----------------------------------------------------------------------------
// Configuration of typed tests for simd vectors.
// ----------------------------------------------------------------------------

template <typename TSimdVector_>
class SimdVectorTestCommon : public seqan::Test
{
public:
    using TValue = typename seqan::Value<TSimdVector_>::Type;
    constexpr static auto const LENGTH = seqan::LENGTH<TSimdVector_>::VALUE;
    using TSimdVector = TSimdVector_;
};

template <typename TSimdVector_>
class SimdVectorTestGather : public seqan::Test
{
public:
    using TValue = typename seqan::Value<TSimdVector_>::Type;
    constexpr static auto const LENGTH = seqan::LENGTH<TSimdVector_>::VALUE;
    using TSimdVector = TSimdVector_;
};

typedef
        seqan::TagList<seqan::SimdVector<int8_t, 16>::Type,
        seqan::TagList<seqan::SimdVector<int16_t, 8>::Type,
        seqan::TagList<seqan::SimdVector<int32_t, 4>::Type,
        seqan::TagList<seqan::SimdVector<int64_t, 2>::Type,
        seqan::TagList<seqan::SimdVector<uint8_t, 16>::Type,
        seqan::TagList<seqan::SimdVector<uint16_t, 8>::Type,
        seqan::TagList<seqan::SimdVector<uint32_t, 4>::Type,
        seqan::TagList<seqan::SimdVector<uint64_t, 2>::Type
        #if SEQAN_SIZEOF_MAX_VECTOR >= 32
        , // Extension of the list above
        seqan::TagList<seqan::SimdVector<int8_t,  32>::Type,
        seqan::TagList<seqan::SimdVector<int16_t, 16>::Type,
        seqan::TagList<seqan::SimdVector<int32_t,  8>::Type,
        seqan::TagList<seqan::SimdVector<int64_t,  4>::Type,
        seqan::TagList<seqan::SimdVector<uint8_t,  32>::Type,
        seqan::TagList<seqan::SimdVector<uint16_t, 16>::Type,
        seqan::TagList<seqan::SimdVector<uint32_t,  8>::Type,
        seqan::TagList<seqan::SimdVector<uint64_t,  4>::Type
        #if SEQAN_SIZEOF_MAX_VECTOR >= 64
        , // Extension of the list above
        seqan::TagList<seqan::SimdVector<int8_t,  64>::Type,
        seqan::TagList<seqan::SimdVector<int16_t, 32>::Type,
        seqan::TagList<seqan::SimdVector<int32_t, 16>::Type,
        seqan::TagList<seqan::SimdVector<int64_t,  8>::Type,
        seqan::TagList<seqan::SimdVector<uint8_t,  64>::Type,
        seqan::TagList<seqan::SimdVector<uint16_t, 32>::Type,
        seqan::TagList<seqan::SimdVector<uint32_t, 16>::Type,
        seqan::TagList<seqan::SimdVector<uint64_t,  8>::Type
        > > > > > > > >
        #endif
        > > > > > > > >
        #endif
        > > > > > > > >
        SimdVectorCommonCommonTypes;

SEQAN_TYPED_TEST_CASE(SimdVectorTestCommon, SimdVectorCommonCommonTypes);
SEQAN_TYPED_TEST_CASE(SimdVectorTestGather, SimdVectorCommonCommonTypes);

SEQAN_DEFINE_TEST(test_simd_types)
{
    using namespace seqan;

    // SimdVector16Char
    static_assert(std::is_same<SimdVector<int8_t, 16>::Type,  SimdVector16SChar>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int16_t, 8>::Type,  SimdVector8Short>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int32_t, 4>::Type,  SimdVector4Int>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int64_t, 2>::Type,  SimdVector2Int64>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint8_t, 16>::Type, SimdVector16UChar>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint16_t, 8>::Type, SimdVector8UShort>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint32_t, 4>::Type, SimdVector4UInt>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint64_t, 2>::Type, SimdVector2UInt64>::value, "should be the same type");

    static_assert(LENGTH<SimdVector4UInt>::VALUE == 4, "128bit register fits 4 int's");
    static_assert(LENGTH<SimdVector<uint32_t, 4>::Type>::VALUE == 4, "128bit register fits 4 int's");
    SimdVector<uint32_t, 4>::Type a128 = {0, 1, 2, 3};
    for (uint32_t i = 0; i < 4; ++i) {
        // std::cout << "DEBUG: " << i << ": " << a128[i] << " = " << i << std::endl;
        SEQAN_ASSERT_EQ(a128[i], i);
    }

    // SimdVector32Char
#if SEQAN_SIZEOF_MAX_VECTOR >= 32
    static_assert(std::is_same<SimdVector<int8_t,  32>::Type,  SimdVector32SChar>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int16_t, 16>::Type,  SimdVector16Short>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int32_t,  8>::Type,  SimdVector8Int>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int64_t,  4>::Type,  SimdVector4Int64>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint8_t,  32>::Type, SimdVector32UChar>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint16_t, 16>::Type, SimdVector16UShort>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint32_t,  8>::Type, SimdVector8UInt>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint64_t,  4>::Type, SimdVector4UInt64>::value, "should be the same type");

    static_assert(LENGTH<SimdVector8UInt>::VALUE == 8, "256bit register fits 8 int's");
    static_assert(LENGTH<SimdVector<uint32_t, 8>::Type>::VALUE == 8, "256bit register fits 8 int's");
    SimdVector<uint32_t, 8>::Type a256 = {0, 1, 2, 3, 4, 5, 6, 7};
    for (uint32_t i = 0; i < 8; ++i) {
        // std::cout << "DEBUG: " << i << ": " << a256[i] << " = " << i << std::endl;
        SEQAN_ASSERT_EQ(a256[i], i);
    }
#endif

#if SEQAN_SIZEOF_MAX_VECTOR >= 64
    static_assert(std::is_same<SimdVector<int8_t,  64>::Type,  SimdVector64SChar>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int16_t, 32>::Type,  SimdVector32Short>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int32_t, 16>::Type,  SimdVector16Int>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<int64_t,  8>::Type,  SimdVector8Int64>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint8_t,  64>::Type, SimdVector64UChar>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint16_t, 32>::Type, SimdVector32UShort>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint32_t, 16>::Type, SimdVector16UInt>::value, "should be the same type");
    static_assert(std::is_same<SimdVector<uint64_t,  8>::Type, SimdVector8UInt64>::value, "should be the same type");

    static_assert(LENGTH<SimdVector16UInt>::VALUE == 16, "512bit register fits 16 int's");
    static_assert(LENGTH<SimdVector<uint32_t, 16>::Type>::VALUE == 16, "512bit register fits 16 int's");
    SimdVector<uint32_t, 16>::Type a512 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    for (uint32_t i = 0; i < 16; ++i) {
        // std::cout << "DEBUG: " << i << ": " << a512[i] << " = " << i << std::endl;
        SEQAN_ASSERT_EQ(a512[i], i);
    }
#endif
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, MetaFunctions)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;

    // NOTE(marehr): explicit namespace is necessary for msvc 2015:
    // error C2039: 'VALUE': is not a member of '`global namespace''
    constexpr auto length = seqan::LENGTH<TSimdVector>::VALUE;
    using TValue = typename Value<TSimdVector>::Type;
    typedef typename SimdVector<TValue, length>::Type TSimdVectorNew;

    bool sameType = IsSameType<TSimdVector, TSimdVectorNew>::VALUE;
    SEQAN_ASSERT(sameType);
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, SizeOf)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u};

    // only on windows are the values unequal
    SEQAN_ASSERT_GEQ(sizeof(a), sizeof(TValue) * length);

    // on linux we assume that the sizes are equal
#ifndef STDLIB_VS
    SEQAN_ASSERT_EQ(sizeof(a), sizeof(TValue) * length);
#endif
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, SubscriptType)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);
    TValue c = a[0];

    bool sameType = IsSameType<TValue, decltype(c)>::VALUE;
    SEQAN_ASSERT(sameType);
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, ClearVector)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    auto zero = static_cast<TValue>(0);
    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    clearVector(a);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)a[i] << " = " << 0 << std::endl;
        SEQAN_ASSERT_EQ(a[i], zero);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, CreateVector)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    auto scalar = static_cast<TValue>(23);
    auto a = createVector<TSimdVector>(scalar);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)a[i] << " = " << 23 << std::endl;
        SEQAN_ASSERT_EQ(a[i], scalar);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, FillVectorConstant)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u};

    fillVector(a, 5);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)a[i] << " = " << i << std::endl;
        SEQAN_ASSERT_EQ(a[i], static_cast<TValue>(5));
    }
}

template <typename TSimdVector, std::size_t... index >
inline void
call_fill_vector(TSimdVector & a, std::index_sequence<index...>)
{
    seqan::fillVector(a, index...);
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, FillVector)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u};

    // calls seqan::fillVector(a, 0, 1, 2, 3, ..., length-1);
    call_fill_vector(a, std::make_index_sequence<length>{});

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)a[i] << " = " << i << std::endl;
        SEQAN_ASSERT_EQ(a[i], static_cast<TValue>(i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, CmpEqual)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TSimdMaskVector = typename SimdMaskVector<TSimdVector>::Type;
    using TValue = typename TestFixture::TValue;
    using TBoolValue = decltype(trueValue<TSimdVector>());
    constexpr auto length = TestFixture::LENGTH;

    TBoolValue false_ = 0,
               true_ = trueValue<TSimdVector>();

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    // There is never a match for the most instantiations of this test.
    a[1] = 23;
    b[1] = 23;

    TSimdMaskVector c = cmpEq(a, b);

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = (i == 1) ? 23 : (-3 + i * 3);
        TValue b_i = (i == 1) ? 23 : (length - i);
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " == " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] == b[i] ? true_ : false_);
        SEQAN_ASSERT_EQ(c[i], a_i == b_i ? true_ : false_);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, CmpGt)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TSimdMaskVector = typename SimdMaskVector<TSimdVector>::Type;
    using TValue = typename TestFixture::TValue;
    using TBoolValue = decltype(trueValue<TSimdVector>());
    constexpr auto length = TestFixture::LENGTH;

    TBoolValue false_ = 0,
               true_ = trueValue<TSimdVector>();

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    TSimdMaskVector c = cmpGt(a, b);

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " > " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] > b[i] ? true_ : false_);
        SEQAN_ASSERT_EQ(c[i], a_i > b_i ? true_ : false_);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Max)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = max(a, b);

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = max (" << (int)a[i] << ", " << (int)b[i] << ")" << std::endl;
        SEQAN_ASSERT_EQ(c[i], std::max(a[i], b[i]));
        SEQAN_ASSERT_EQ(c[i], std::max(a_i, b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Min)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = min(a, b);

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = min (" << (int)a[i] << ", " << (int)b[i] << ")" << std::endl;
        SEQAN_ASSERT_EQ(c[i], std::min(a[i], b[i]));
        SEQAN_ASSERT_EQ(c[i], std::min(a_i, b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, BitwiseOr)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = a | b;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " | " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] | b[i]);
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i | b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, BitwiseOrAssign)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u}, c{0u};
    fillVectors(a, b);

    c = a;
    c |= b;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " | " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] | b[i]);
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i | b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, BitwiseAnd)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = a & b;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " & " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] & b[i]);
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i & b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, BitwiseAndAssign)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u}, c{0u};
    fillVectors(a, b);

    c = a;
    c &= b;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " & " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] & b[i]);
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i & b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, BitwiseNot)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = ~a;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = ~" << (int)a[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(~a[i]));
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(~a_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Addition)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = a + b;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " + " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a[i] + b[i]));
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i + b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Subtraction)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = a - b;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " - " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a[i] - b[i]));
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i - b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Multiplication)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = a * b;

    for (size_t i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " * " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a[i] * b[i]));
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i * b_i));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Division)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = a / b;

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " / " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] / b[i]);
        SEQAN_ASSERT_EQ(c[i], a_i / b_i);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, BitwiseAndNot)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = andNot(a, b);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = (~" << (int)a[i] << ") & " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(~a[i] & b[i]));
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(~(-3 + i * 3) & (length - i)));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, ShiftRightLogical)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    // ensure that a >= 0, because (-3) >> 2 has undefined behavior according to
    // C++ 11 standard.
    a = a + createVector<TSimdVector>(3);

    auto c = shiftRightLogical(a, 2);

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = i * 3;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " >> " << (int)2 << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] >> 2);
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a_i >> 2));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Blend)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    auto c = blend(b, a, cmpGt(a, b));

    for (auto i = 0; i < length; ++i)
    {
        TValue a_i = -3 + i * 3, b_i = length - i;
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << " > " << (int)b[i] << " ? " << (int)a[i] << " : " << (int)b[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i] > b[i] ? (TValue)a[i] : (TValue)b[i]);
        SEQAN_ASSERT_EQ(c[i], a_i > b_i ? a_i : b_i);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Storeu)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, b{0u};
    fillVectors(a, b);

    TValue c[length];
    storeu(c, a);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i]);
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(-3 + i * 3));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Load)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, c{0u};
    fillVectors(a, c);

    alignas(TSimdVector) TValue b[length];
    storeu(b, a);
    c = load<TSimdVector>(b);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[i]);
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(-3 + i * 3));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Gather)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    using TValue = typename TestFixture::TValue;
    constexpr auto length = TestFixture::LENGTH;

    TSimdVector a{0u}, idx{0u};
    fillVectors(a, idx);
    reverseIndexSequence(idx, length);

    TValue b[length];
    storeu(b, a);
    auto c = gather(b, idx);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[idx[i]] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[idx[i]]);
        SEQAN_ASSERT_EQ(c[i], a[length - i - 1]);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, ShuffleConstant1)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    constexpr auto length = TestFixture::LENGTH;
    typedef typename SimdSwizzleVector<TSimdVector>::Type TSimdSwizzleVector;

    TSimdVector a{0u}, b{0u};
    auto idx = createVector<TSimdSwizzleVector>(1);
    fillVectors(a, b);

    auto c = shuffleVector(a, idx);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[idx[i]] << ", idx: " << (int)idx[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[idx[i]]);
        SEQAN_ASSERT_EQ(c[i], a[1]);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, ShuffleConstant2)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    constexpr auto length = TestFixture::LENGTH;
    typedef typename SimdSwizzleVector<TSimdVector>::Type TSimdSwizzleVector;

    TSimdVector a{0u}, b{0u};
    auto idx = createVector<TSimdSwizzleVector>(length - 2);
    fillVectors(a, b);

    auto c = shuffleVector(a, idx);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[idx[i]] << ", idx: " << (int)idx[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[idx[i]]);
        SEQAN_ASSERT_EQ(c[i], a[length-2]);
    }
}

SEQAN_TYPED_TEST(SimdVectorTestCommon, Shuffle)
{
    using namespace seqan;
    using TSimdVector = typename TestFixture::TSimdVector;
    constexpr auto length = TestFixture::LENGTH;
    typedef typename SimdSwizzleVector<TSimdVector>::Type TSimdSwizzleVector;

    TSimdVector a{0u}, b{0u};
    TSimdSwizzleVector idx{0u};
    fillVectors(a, b);
    reverseIndexSequence(idx, length);

    auto c = shuffleVector(a, idx);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (int)c[i] << " = " << (int)a[idx[i]] << ", idx: " << (int)idx[i] << std::endl;
        SEQAN_ASSERT_EQ(c[i], a[idx[i]]);
        SEQAN_ASSERT_EQ(c[i], a[length - i - 1]);
    }
}

template <typename TSimdVector, typename TValue, typename TArrayValue>
inline void test_gather_array()
{
    using namespace seqan;
    constexpr auto length = LENGTH<TSimdVector>::VALUE;

    TSimdVector idx{0u};
    reverseIndexSequence(idx, length);

    TArrayValue a[2*length];

    // fill gather array
    for (auto i = 0; i < 2*length; ++i)
    {
        a[i] = (i-1)*3;
    }

    auto c = gather(a, idx);

    for (auto i = 0; i < length; ++i)
    {
        // std::cout << "DEBUG: " << i << " / " << length << ": " << (TValue)c[i] << " = " << (TValue)a[idx[i]] << std::endl;
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a[idx[i]]));
        SEQAN_ASSERT_EQ(c[i], static_cast<TValue>(a[length - i - 1]));
    }
}

SEQAN_TYPED_TEST(SimdVectorTestGather, CharArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        int8_t
    >();
}

SEQAN_TYPED_TEST(SimdVectorTestGather, ShortArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        int16_t
    >();
}

SEQAN_TYPED_TEST(SimdVectorTestGather, IntArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        int32_t
    >();
}

SEQAN_TYPED_TEST(SimdVectorTestGather, LongArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        int64_t
    >();
}

SEQAN_TYPED_TEST(SimdVectorTestGather, UCharArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        uint8_t
    >();
}

SEQAN_TYPED_TEST(SimdVectorTestGather, UShortArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        uint16_t
    >();
}

SEQAN_TYPED_TEST(SimdVectorTestGather, UIntArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        uint32_t
    >();
}

SEQAN_TYPED_TEST(SimdVectorTestGather, ULongArray)
{
    test_gather_array<
        typename TestFixture::TSimdVector,
        typename TestFixture::TValue,
        uint64_t
    >();
}

#ifdef __SSE4_1__

SEQAN_DEFINE_TEST(test_simd_transpose_8x8)
{
    seqan::test_matrix_transpose<8, seqan::SimdVector<unsigned char, 8>::Type>();
}

SEQAN_DEFINE_TEST(test_simd_transpose_16x16)
{
    seqan::test_matrix_transpose<16, seqan::SimdVector<unsigned char, 16>::Type>();
}

#endif  // #ifdef __SSE4_1__
#ifdef __AVX2__

SEQAN_DEFINE_TEST(test_simd_transpose_32x32)
{
    seqan::test_matrix_transpose<32, seqan::SimdVector<unsigned char, 32>::Type >();
}

#endif  // #ifdef __AVX2__
#endif  // SEQAN_SIMD_ENABLED

#endif  // #ifndef SEQAN_CORE_TESTS_SIMD_TEST_SIMD_VECTOR_H_
