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
// Tests for the SeqAn type header.
// ==========================================================================

#ifndef TEST_BASIC_TEST_BASIC_TYPE_H_
#define TEST_BASIC_TEST_BASIC_TYPE_H_

#include <seqan/basic.h>

// Test structs used in this header.  We will use TestStruct2_ for
// testing default metafunction implementations and TestStruct1_ for
// our own.
struct TestStruct1_ {};
struct TestStruct2_ {};

namespace seqan {

template <>
struct Value<TestStruct1_>
{
    typedef int Type;
};

template <>
struct Value<TestStruct1_ const>
{
    typedef int const Type;
};

}  // namespace seqan

// Test our Value metafunction.  The tested point here is that the
// mechanisms works the way it should.
SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_value)
{
    bool b = IsSameType<int, Value<TestStruct1_>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const, Value<TestStruct1_ const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);

    // TODO(holtgrew): Return-self should probably be removed in the future, this is part of the elements-are-containers issue.
    b = IsSameType<TestStruct2_, Value<TestStruct2_>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<TestStruct2_ const, Value<TestStruct2_ const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct GetValue<TestStruct1_>
{
    typedef int Type;
};

template <>
struct GetValue<TestStruct1_ const>
{
    typedef int Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_get_value)
{
    bool b = IsSameType<int, GetValue<TestStruct1_>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int, GetValue<TestStruct1_ const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);

    // TODO(holtgrew): Return-self should probably be removed in the future, this is part of the elements-are-containers issue.
    b = IsSameType<TestStruct2_ const &, GetValue<TestStruct2_>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<TestStruct2_ const &, GetValue<TestStruct2_ const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Reference<TestStruct1_>
{
    typedef int & Type;
};

template <>
struct Reference<TestStruct1_ const>
{
    typedef int const & Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_reference)
{
    bool b = IsSameType<int &, Reference<TestStruct1_>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const &, Reference<TestStruct1_ const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);

    // TODO(holtgrew): Return-self should probably be removed in the future, this is part of the elements-are-containers issue.
    b = IsSameType<TestStruct2_ &, Reference<TestStruct2_>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<TestStruct2_ const &, Reference<TestStruct2_ const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Size<TestStruct1_>
{
    typedef unsigned short Type;
};

template <>
struct Size<TestStruct1_ const>
{
    typedef unsigned short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_size)
{
    bool b = IsSameType<Size<TestStruct1_>::Type, unsigned short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Size<TestStruct1_ const>::Type, unsigned short>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<Size<TestStruct2_>::Type, size_t>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Size<TestStruct2_ const>::Type, size_t>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Difference<TestStruct1_>
{
    typedef short Type;
};

template <>
struct Difference<TestStruct1_ const>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_difference)
{
    bool b = IsSameType<Difference<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Difference<TestStruct1_ const>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<Difference<TestStruct2_>::Type, ptrdiff_t>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Difference<TestStruct2_ const>::Type, ptrdiff_t>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Position<TestStruct1_>
{
    typedef short Type;
};

template <>
struct Position<TestStruct1_ const>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_position)
{
    bool b = IsSameType<Position<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Position<TestStruct1_ const>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<Position<TestStruct2_>::Type, size_t>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Position<TestStruct2_ const>::Type, size_t>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Host<TestStruct1_>
{
    typedef int Type;
};

template <>
struct Host<TestStruct1_ const>
{
    typedef int const Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_host)
{
}

template <typename T>
struct MyClass_ {};

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_spec)
{
    bool b = IsSameType<Spec<MyClass_<int> >::Type, int>::Type::VALUE;
    SEQAN_ASSERT(b);
}

template <typename T>
struct MyClass2_ {};

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_deepest_spec)
{
    bool b = IsSameType<DeepestSpec<MyClass2_<int> >::Type, int>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<DeepestSpec<MyClass2_<MyClass_<int> > >::Type, int>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Cargo<TestStruct1_>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_cargo)
{
    bool b = IsSameType<Cargo<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Cargo<TestStruct1_ const>::Type, short const>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<Cargo<TestStruct2_>::Type, Nothing>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Cargo<TestStruct2_ const>::Type, Nothing const>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct VertexDescriptor<TestStruct1_>
{
    typedef short Type;
};

template <>
struct VertexDescriptor<TestStruct1_ const>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_vertex_descriptor)
{
    bool b = IsSameType<VertexDescriptor<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<VertexDescriptor<TestStruct1_ const>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<VertexDescriptor<TestStruct2_>::Type, void *>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<VertexDescriptor<TestStruct2_ const>::Type, void *>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Id<TestStruct1_>
{
    typedef short Type;
};

template <>
struct Id<TestStruct1_ const>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_id)
{
    bool b = IsSameType<Id<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Id<TestStruct1_ const>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<Id<TestStruct2_>::Type, unsigned int>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Id<TestStruct2_ const>::Type, unsigned int>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Key<TestStruct1_>
{
    typedef short Type;
};

template <>
struct Key<TestStruct1_ const>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_key)
{
    bool b = IsSameType<Key<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Key<TestStruct1_ const>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<Key<TestStruct2_>::Type, TestStruct2_>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Key<TestStruct2_ const>::Type, TestStruct2_>::Type::VALUE;
    SEQAN_ASSERT(b);
}

namespace seqan {

template <>
struct Object<TestStruct1_>
{
    typedef short Type;
};

template <>
struct Object<TestStruct1_ const>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_object)
{
    bool b = IsSameType<Object<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Object<TestStruct1_ const>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);

    // Has no default implementation.
}

namespace seqan {

template <>
struct Source<TestStruct1_>
{
    typedef short Type;
};

template <>
struct Source<TestStruct1_ const>
{
    typedef short Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_source)
{
    bool b = IsSameType<Source<TestStruct1_>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Source<TestStruct1_ const>::Type, short>::Type::VALUE;
    SEQAN_ASSERT(b);

    b = IsSameType<Source<TestStruct2_>::Type, TestStruct2_>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Source<TestStruct2_ const>::Type, TestStruct2_>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_parameter)
{
    bool b = IsSameType<Parameter_<int>::Type, int &>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Parameter_<int const>::Type, int const &>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Parameter_<int *>::Type, int *>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Parameter_<int const *>::Type, int const *>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Parameter_<int [5]>::Type, int *>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<Parameter_<int const [5]>::Type, int const *>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(seqan_basic_type_to_parameter)
{
    // TODO(holtgrew): _toParameter() necessary/good?

    int i;
    int const ci = 99;
    int a[10];

    _toParameter<int>(& i) = 10;
    SEQAN_ASSERT_EQ(i, 10);

    *_toParameter<int *>(& i) = 20;
    SEQAN_ASSERT_EQ(i, 20);

    Pointer_<int>::Type p1 = _toPointer(i);
    *p1 = 30;
    SEQAN_ASSERT_EQ(i, 30);

    Pointer_<int *>::Type p2 = _toPointer(p1);
    *p2 = 40;
    SEQAN_ASSERT_EQ(i, 40);

    Pointer_<int[10]>::Type p3 = _toPointer(a);
    p3[1] = 50;
    SEQAN_ASSERT_EQ(a[1], 50);

    Pointer_<int const *>::Type p4 = _toPointer(ci);
    SEQAN_ASSERT_EQ(*p4, 99);
}

namespace seqan {

template <>
struct LENGTH<TestStruct1_>
{
    enum { VALUE = 10 };
};

template <>
struct LENGTH<TestStruct1_ const>
{
    enum { VALUE = 10 };
};

}  // namespace seqan

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_length)
{
    unsigned len = LENGTH<TestStruct1_>::VALUE;
    SEQAN_ASSERT_EQ(len, 10u);
    len = LENGTH<TestStruct1_>::VALUE;
    SEQAN_ASSERT_EQ(len, 10u);

    // TODO(holtgrew): elements-are-containers should go away!
    len = LENGTH<TestStruct2_>::VALUE;
    SEQAN_ASSERT_EQ(len, 1u);
    len = LENGTH<TestStruct2_>::VALUE;
    SEQAN_ASSERT_EQ(len, 1u);
}

SEQAN_DEFINE_TEST(seqan_basic_type_metafunction_is_integral)
{
    bool b = IsInteger<char>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsInteger<short int>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsInteger<int>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsInteger<long int>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsInteger<int64_t>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsInteger<TestStruct1_>::Type::VALUE;
    SEQAN_ASSERT_NOT(b);
}

#endif  // TEST_BASIC_TEST_BASIC_TYPE_H_
