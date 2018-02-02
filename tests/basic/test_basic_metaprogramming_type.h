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
// Tests for the type query / modification part of the metaprogramming
// library.
// ==========================================================================

// Note that we use +Or<T>::VALUE and +false below.  The plus before the metafunction is such that no const-ref is
// generated for the static const bool member which causes linker errors.  The plus before the constant is there
// to suppress bool/integer comparison warnings with MSVC++.

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_TYPE_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_TYPE_H_

SEQAN_DEFINE_TEST(test_basic_metaprogramming_type_same_type)
{
    using namespace seqan;

    // Test for the values of the VALUE members.
    SEQAN_ASSERT_EQ((+IsSameType<bool, int>::VALUE),  +false);
    SEQAN_ASSERT_EQ((+IsSameType<bool, bool>::VALUE), +true);

    // Test for type of the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<IsSameType<bool, int>::Type,  False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<IsSameType<bool, bool>::Type, True>::VALUE),  true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_type_make_signed)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<unsigned char>::Type,   signed char>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<signed   char>::Type,   signed char>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<         char>::Type,   signed char>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<unsigned short>::Type,  signed short>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<signed   short>::Type,  signed short>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<         short>::Type,  signed short>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<unsigned int>::Type,    signed int>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<signed   int>::Type,    signed int>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<         int>::Type,    signed int>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<unsigned long>::Type,   signed long>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<signed   long>::Type,   signed long>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<         long>::Type,   signed long>::VALUE), true);

    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<int8_t  >::Type, int8_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<int16_t >::Type, int16_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<uint16_t>::Type, int16_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<int32_t >::Type, int32_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<uint32_t>::Type, int32_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<int64_t >::Type, int64_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeSigned<uint64_t>::Type, int64_t>::VALUE), true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_type_make_unsigned)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<unsigned char>::Type,   unsigned char>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<signed   char>::Type,   unsigned char>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<         char>::Type,   unsigned char>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<unsigned short>::Type,  unsigned short>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<signed   short>::Type,  unsigned short>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<         short>::Type,  unsigned short>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<unsigned int>::Type,    unsigned int>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<signed   int>::Type,    unsigned int>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<         int>::Type,    unsigned int>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<unsigned long>::Type,   unsigned long>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<signed   long>::Type,   unsigned long>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<         long>::Type,   unsigned long>::VALUE), true);

    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<int8_t  >::Type, uint8_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<uint8_t >::Type, uint8_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<int16_t >::Type, uint16_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<uint16_t>::Type, uint16_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<int32_t >::Type, uint32_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<uint32_t>::Type, uint32_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<int64_t >::Type, uint64_t>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<MakeUnsigned<uint64_t>::Type, uint64_t>::VALUE), true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_type_remove_reference)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ((TestTypeEq<typename RemoveReference<unsigned &>::Type, unsigned>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<typename RemoveReference<unsigned  >::Type, unsigned>::VALUE),      true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_type_remove_const)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ((TestTypeEq<typename RemoveConst<unsigned const>::Type, unsigned>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<typename RemoveConst<unsigned      >::Type, unsigned>::VALUE), true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_type_is_const)
{
    using namespace seqan;

    // Test for the values of the VALUE members.
    SEQAN_ASSERT_EQ((+IsConst_<bool>::VALUE),       +false);
    SEQAN_ASSERT_EQ((+IsConst_<bool const>::VALUE), +true);

    // Test for type of the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<IsConst_<bool>::Type,       False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<IsConst_<bool const>::Type, True>::VALUE),  true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_type_class_identifier)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(ClassIdentifier_<int>::getID(), ClassIdentifier_<int>::getID());
    SEQAN_ASSERT_NEQ(ClassIdentifier_<int>::getID(), ClassIdentifier_<bool>::getID());
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_TYPE_H_
