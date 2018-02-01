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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_ALPHABET_H_
#define TESTS_BASIC_TEST_BASIC_ALPHABET_H_

using namespace seqan;

//Helper class that counts ctors, dtors and copys
struct Test1
{
    static int m_ctor_count;
    static int m_dtor_count;
    static int m_copy_count;
    static int m_move_count;

    Test1(): x(0xfade)
    {
        ++m_ctor_count;
    }
    Test1(Test1 const & obj): x(obj.x)
    {
        ++m_ctor_count;
        ++m_copy_count;
    }
    Test1(Test1 const & obj, Move): x(obj.x)
    {
        ++m_move_count;
        obj.x = 0x105e;
    }
    ~Test1()
    {
        x = 0xdead;
        ++m_dtor_count;
    }
    Test1 & operator = (Test1 const & obj)
    {
        x = obj.x;
        ++m_copy_count;
        return *this;
    }
    mutable int x;
};

int Test1::m_ctor_count = 0;
int Test1::m_dtor_count = 0;
int Test1::m_copy_count = 0;
int Test1::m_move_count = 0;


inline void
move(Test1 & target, Test1 & source)
{
    ++Test1::m_move_count;
    target.x = source.x;
    source.x = 0x105e;
}
inline void
move(Test1 const & target, Test1 & source)
{
    ++Test1::m_move_count;
    target.x = source.x;
    source.x = 0x105e;
}
inline void
move(Test1 & target, Test1 const & source)
{
    ++Test1::m_move_count;
    target.x = source.x;
    source.x = 0x105e;
}
inline void
move(Test1 const & target, Test1 const & source)
{
    ++Test1::m_move_count;
    target.x = source.x;
    source.x = 0x105e;
}

//////////////////////////////////////////////////////////////////////////////
//Test value array function for a class type

SEQAN_DEFINE_TEST(test_basic_alphabet_interface)
{
    {
        Test1 a; //1 ctor
        a.x = 0xbeef;

        char c_buf1[200 * sizeof(Test1)];
        Test1 * a_buf1 = (Test1 *) c_buf1;

        char c_buf2[200 * sizeof(Test1)];
        Test1 * a_buf2 = (Test1 *) c_buf2;

        arrayConstruct(a_buf1, a_buf1 + 100); //100 ctor
        for (int i=0; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i].x, 0xfade);
        }

        arrayConstruct(a_buf2, a_buf2 + 100, a); //100 ctor, 100 copy
        for (int i=0; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf2[i].x, 0xbeef);
        }

        arrayConstruct(a_buf1, a_buf1+ 100); //100 ctor
        arrayDestruct(a_buf1, a_buf1 + 100); //100 dtor
        for (int i=0; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i].x, 0xdead);
        }

        arrayConstructCopy(a_buf1, a_buf1 + 100, a_buf2); //100 ctor, 100 copy
        arrayDestruct(a_buf1, a_buf1 + 100); //100 dtor
        for (int i=0; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf2[i].x, 0xdead);
        }

        arrayFill(a_buf1, a_buf1 + 100, a); //100 copy
        for (int i=0; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i].x, 0xbeef);
        }

        arrayCopyForward(a_buf1, a_buf1 + 100, a_buf2); //100 copy

        for (int i=0; i < 100; ++i) a_buf1[i].x = i;

        arrayCopy(a_buf1, a_buf1 + 50, a_buf1 + 20); //50 copy
        for (int i=0; i < 50; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i+20].x, i);
        }

        arrayCopy(a_buf1 + 80, a_buf1 + 100, a_buf1 + 75); //20 copy
        for (int i=80; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i-5].x, i);
        }

        for (int i=0; i < 100; ++i) a_buf1[i].x = i;

        arrayClearSpace(a_buf1, 100, 50, 70); //20 ctor, 70 dtor, 50 copy
        for (int i=50; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i+20].x, i);
        }
        arrayConstruct(a_buf1, a_buf1 + 70, a); //70 ctor, 70 copy

        arrayClearSpace(a_buf1, 120, 70, 50); //70 dtor, 50 copy
        for (int i=50; i < 100; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i].x, i);
        }

        arrayConstruct(a_buf1, a_buf1 + 50, a); //70 ctor, 70 copy

        for (int i=0; i < 100; ++i) a_buf1[i].x = i;

        arrayClearSpace(a_buf1, 100, 90, 110); //10 ctor, 90 dtor, 10 copy
        for (int i=110; i < 120; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i].x, i-20);
        }

        arrayDestruct(a_buf1 + 110, a_buf1 + 120);    //10 dtor
        arrayDestruct(a_buf2, a_buf2 + 100); //100 dtor


        /* Commented out until move construction is fixed.
           See http://trac.mi.fu-berlin.de/seqan/ticket/380 for more information.
        // TODO(holtgrew): Fix move construction
        arrayConstruct(a_buf2, a_buf2 + 23); //23 ctor
        arrayConstructMove(a_buf2, a_buf2 + 23, a_buf1); // 23 move
        for (int i = 0; i < 23; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i].x, 0xfade);
            SEQAN_ASSERT_EQ(a_buf2[i].x, 0x105e);
        }

        // TODO(holtgrew): Fix moving of values in arrays
        arrayMove(a_buf1, a_buf1 + 23, a_buf1 + 5); // 23 move
        for (int i = 0; i < 23; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i + 5].x, 0xfade);
        }
        for (int i = 0; i < 5; ++i)
        {
            SEQAN_ASSERT_EQ(a_buf1[i].x, 0x105e);
        }

        arrayMove(a_buf1 + 5, a_buf1 + 28, a_buf1); // 23 move

        arrayDestruct(a_buf1, a_buf1 + 23); //23 dtor
        */

        //1 dtor for a
    }

    /* Commented out until http://trac.mi.fu-berlin.de/seqan/ticket/380 is fixed.
    // TODO(holtgrew): Fix moving of values in arrays
    SEQAN_ASSERT_EQ(Test1::m_ctor_count, 574);
    SEQAN_ASSERT_EQ(Test1::m_dtor_count, 574);
    SEQAN_ASSERT_EQ(Test1::m_copy_count, 700);
    SEQAN_ASSERT_EQ(Test1::m_move_count, 69);


    SEQAN_ASSERT_EQ(gapValue<char>(), '-');
    SEQAN_ASSERT_EQ(gapValue<int>(), int());
    */
}

//////////////////////////////////////////////////////////////////////////////
//Test value array functions for some types

template <typename T, typename T_>
void TestArrayFunctions(T_ const _val1, T_ const _val2)
{
    T val1 = (T)_val1;
    T val2 = (T)_val2;

    T a = val1;

    T a_buf1[200];

    T a_buf2[200];

    arrayConstruct(a_buf1, a_buf1 + 100); //nothing happens

    arrayConstruct(a_buf2, a_buf2 + 100, a);
    for (int i=0; i < 100; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf2[i], val1);
    }

    arrayDestruct(a_buf1, a_buf1 + 100); //nothing happens

    arrayConstructCopy(a_buf2, a_buf2 + 100, a_buf1);
    for (int i=0; i < 100; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf1[i], val1);
    }

    a = val2;

    arrayFill(a_buf1, a_buf1 + 100, a);
    for (int i=0; i < 100; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf1[i], val2);
    }

    arrayCopyForward(a_buf1, a_buf1 + 100, a_buf2);
    for (int i=0; i < 100; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf2[i], val2);
    }

    for (int i=0; i < 100; ++i) a_buf1[i] = (T) i;

    arrayCopy(a_buf1, a_buf1 + 50, a_buf1 + 20);
    for (int i=0; i < 50; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf1[i+20], (T)i);
    }

    arrayCopy(a_buf1 + 80, a_buf1 + 100, a_buf1 + 75);
    for (int i=80; i < 100; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf1[i-5], (T)i);
    }

    for (int i=0; i < 100; ++i) a_buf1[i] = (T) i;

    arrayClearSpace(a_buf1, 100, 50, 70);
    for (int i=50; i < 100; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf1[i+20], (T)i);
    }

    arrayClearSpace(a_buf1, 120, 70, 50);
    for (int i=50; i < 100; ++i)
    {
        SEQAN_ASSERT_EQ(a_buf1[i], (T)i);
    }

}

//////////////////////////////////////////////////////////////////////////////
//Test SimpleType instances

template <typename T>
void TestSimpleType()
{
    T a = T();
    T b(a);
    b = a;

    typename Value<T>::Type val;
    val = a;
    a = val;

    T a_buf1[200];
    T a_buf2[200];
    T a_buf3[200];

    arrayConstruct(a_buf1, a_buf1 + 100);
    arrayConstruct(a_buf1, a_buf1 + 100, a);
    arrayDestruct(a_buf1, a_buf1 + 100);
    arrayConstructCopy(a_buf2, a_buf2 + 100, a_buf1);
    arrayFill(a_buf1, a_buf1 + 100, a);
    arrayCopyForward(a_buf1, a_buf1 + 100, a_buf2);
    arrayCopy(a_buf1, a_buf1 + 50, a_buf1 + 20);
    arrayCopy(a_buf1 + 80, a_buf1 + 100, a_buf1 + 75);
    arrayMoveForward(a_buf1, a_buf1 + 10, a_buf3);
    arrayMoveBackward(a_buf3, a_buf3 + 10, a_buf1);
    arrayClearSpace(a_buf1, 100, 50, 70);
    arrayClearSpace(a_buf1, 120, 70, 50);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource>
void TestConversion()
{
    static TSource source;
    TTarget target = source;

    SEQAN_ASSERT_EQ(target, source);
    SEQAN_ASSERT_EQ(source, target);
    SEQAN_ASSERT(!(target != source));
    SEQAN_ASSERT(!(source != target));
    SEQAN_ASSERT(!(target < source));
    SEQAN_ASSERT(!(source < target));
    SEQAN_ASSERT(target <= source);
    SEQAN_ASSERT(source <= target);
    SEQAN_ASSERT(!(target > source));
    SEQAN_ASSERT(!(source > target));
    SEQAN_ASSERT(target >= source);
    SEQAN_ASSERT(source >= target);

    TSource const c_source = TSource();
    target = c_source;

    target = TSource();
    assign(target, source);

    TSource a_source_1[200];
    TTarget a_target[200];

    arrayCopyForward(a_source_1, a_source_1 + 100, a_target);
    arrayCopy(a_source_1, a_source_1 + 50, a_target + 20);
    arrayCopyBackward(a_source_1, a_source_1 + 100, a_target);
    arrayMoveForward(a_source_1, a_source_1 + 100, a_target);
    arrayMove(a_source_1, a_source_1 + 50, a_target + 20);
    arrayMoveBackward(a_source_1, a_source_1 + 100, a_target);
}

SEQAN_DEFINE_TEST(test_basic_conversions)
{
    TestConversion<char, Dna>();
    TestConversion<char, Dna5>();
    TestConversion<char, Rna>();
    TestConversion<char, Rna5>();
    TestConversion<char, Iupac>();
    TestConversion<char, AminoAcid>();

    TestConversion<Dna, char>();
    TestConversion<Dna, uint8_t>();
    TestConversion<Dna, Unicode>();
    TestConversion<Dna, int>();
    TestConversion<Dna, Dna5>();
    TestConversion<Dna, Iupac>();

    TestConversion<Dna5, char>();
    TestConversion<Dna5, uint8_t>();
    TestConversion<Dna5, Unicode>();
    TestConversion<Dna5, Dna>();
    TestConversion<Dna5, Iupac>();

    TestConversion<Rna, char>();
    TestConversion<Rna, uint8_t>();
    TestConversion<Rna, Unicode>();
    TestConversion<Rna, int>();
    TestConversion<Rna, Rna5>();

    TestConversion<Rna5, char>();
    TestConversion<Rna5, uint8_t>();
    TestConversion<Rna5, Unicode>();
    TestConversion<Rna5, Rna>();

    TestConversion<Iupac, char>();
    TestConversion<Iupac, uint8_t>();
    TestConversion<Iupac, Unicode>();
    TestConversion<Iupac, Dna>();
    TestConversion<Iupac, Dna5>();

    TestConversion<AminoAcid, char>();
    TestConversion<AminoAcid, uint8_t>();
    TestConversion<AminoAcid, Unicode>();

    typedef SimpleType<int, void> ST;

    TestConversion<int, ST>();
    TestConversion<unsigned int, ST>();
    TestConversion<short, ST>();
    TestConversion<unsigned short, ST>();
    TestConversion<char, ST>();
    TestConversion<signed char, ST>();
    TestConversion<unsigned char, ST>();
    TestConversion<ST, ST>();

}

//////////////////////////////////////////////////////////////////////////////
//Test infimum / supremum values

template <typename T>
void TestExtremeValuesSigned()
{
    long double minVal = -1;
    for(unsigned e = 1; e < BitsPerValue<T>::VALUE; ++e)
        minVal = 2*minVal;

    long double maxVal = -minVal - 1;

    bool isSigned = IsSameType< typename MakeSigned_<T>::Type, T >::VALUE;
    SEQAN_ASSERT(isSigned);

    long double maxDelta = maxVal - std::numeric_limits<T>::max();
    long double minDelta = minVal - (long double)std::numeric_limits<T>::min();
    SEQAN_ASSERT(maxDelta <= maxVal/1000);
    SEQAN_ASSERT(-maxVal/1000 <= maxDelta);
    SEQAN_ASSERT(minDelta <= maxVal/1000);
    SEQAN_ASSERT(-maxVal/1000 <= minDelta);
}

template <typename T>
void TestExtremeValuesUnsigned()
{
    long double maxVal = 1;
    for(unsigned e = 0; e < BitsPerValue<T>::VALUE; ++e)
        maxVal = 2*maxVal;
    maxVal = maxVal - 1;

    bool isUnsigned = IsSameType< typename MakeUnsigned_<T>::Type, T >::VALUE;
    SEQAN_ASSERT(isUnsigned);

    long double maxDelta = maxVal - std::numeric_limits<T>::max();
    SEQAN_ASSERT_LEQ(maxDelta, maxVal/1000);
    SEQAN_ASSERT_LEQ(-maxVal/1000, maxDelta);
    SEQAN_ASSERT_EQ((T)0, std::numeric_limits<T>::min());
}

SEQAN_DEFINE_TEST(test_basic_alphabet_extreme_values)
{
    TestExtremeValuesSigned<signed char>();
    TestExtremeValuesSigned<signed short>();
    TestExtremeValuesSigned<signed int>();
    TestExtremeValuesSigned<signed long>();
    TestExtremeValuesUnsigned<unsigned char>();
    TestExtremeValuesUnsigned<unsigned short>();
    TestExtremeValuesUnsigned<unsigned int>();
    TestExtremeValuesUnsigned<unsigned long>();
    TestExtremeValuesSigned<int64_t>();
/*    TestExtremeValues<float>();
    TestExtremeValues<double>();
    TestExtremeValues<long double>();*/
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_basic_simple_types) {
    TestSimpleType<Dna>();
    TestSimpleType<Dna5>();
    TestSimpleType<Rna>();
    TestSimpleType<Rna5>();
    TestSimpleType<Iupac>();
    TestSimpleType<AminoAcid>();
    TestSimpleType<bool>();
}

SEQAN_DEFINE_TEST(test_basic_array_functions)
{
    TestArrayFunctions<char>(0xde, 0xad);
    TestArrayFunctions<signed char>(0xde, 0xad);
    TestArrayFunctions<unsigned char>(0xde, 0xad);
    TestArrayFunctions<short>(0xdead, 0xbeef);
    TestArrayFunctions<unsigned short>(0xdead, 0xbeef);
    TestArrayFunctions<int>(0xdead, 0xbeef);
    TestArrayFunctions<unsigned int>(0xdead, 0xbeef);
    TestArrayFunctions<float>(3.1, 1.2);
    TestArrayFunctions<double>(3.1, 1.2);
    TestArrayFunctions<long double>(3.1, 1.2);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_value_size)
{
    SEQAN_ASSERT_EQ(+ValueSize<bool>::VALUE, 2u);
    SEQAN_ASSERT_EQ(+ValueSize<int8_t>::VALUE, 256u);
    SEQAN_ASSERT_EQ(+ValueSize<uint8_t>::VALUE, 256u);
    SEQAN_ASSERT_EQ(+ValueSize<int16_t>::VALUE, 65536u);
    SEQAN_ASSERT_EQ(+ValueSize<uint16_t>::VALUE, 65536u);
    SEQAN_ASSERT_EQ(+ValueSize<int32_t>::VALUE, (uint64_t)4294967296ll);
    SEQAN_ASSERT_EQ(+ValueSize<uint32_t>::VALUE, (uint64_t)4294967296ll);
    SEQAN_ASSERT_EQ(+ValueSize<int64_t>::VALUE, 0u);
    SEQAN_ASSERT_EQ(+ValueSize<uint64_t>::VALUE, 0u);

    SEQAN_ASSERT_EQ(valueSize<bool>(), 2u);
    SEQAN_ASSERT_EQ(valueSize<int8_t>(), 256u);
    SEQAN_ASSERT_EQ(valueSize<uint8_t>(), 256u);
    SEQAN_ASSERT_EQ(valueSize<int16_t>(), 65536u);
    SEQAN_ASSERT_EQ(valueSize<uint16_t>(), 65536u);
    SEQAN_ASSERT_EQ(valueSize<int32_t>(), (uint64_t)4294967296ll);
    SEQAN_ASSERT_EQ(valueSize<uint32_t>(), (uint64_t)4294967296ll);
    SEQAN_ASSERT_EQ(valueSize<int64_t>(), 0u);
    SEQAN_ASSERT_EQ(valueSize<uint64_t>(), 0u);
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_ALPHABET_H_
