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
// Tests for the Pair, Triple and Tuple hierarchy in module basic.
// ==========================================================================

// TODO(holtgrew): Test/Implement Metafunction.GetValue.

#include <sstream>

#ifndef TESTS_BASIC_TEST_BASIC_AGGREGATES_H_
#define TESTS_BASIC_TEST_BASIC_AGGREGATES_H_

using namespace seqan;

// ==========================================================================
// Helper Code
// ==========================================================================

// Helper type for testing transports.  Each constructed objects gets an id
// and increases the static nextId member.  Assignment only updates the id,
// setting sets the id as well and movement copies id and value and sets the
// id of the other to -1 to signal emptiness.

struct Transportable_
{
    static int nextId;
    int id;
    int value;

    Transportable_()
    {
        id = nextId++;
    }

    explicit
    Transportable_(int i)
    {
        id = nextId++;
        value = i;
    }

    Transportable_(const Transportable_ & other)
    {
        id = nextId++;
        value = other.value;
    }

    Transportable_(Transportable_ & other, Move const & /*tag*/)
    {
        id = other.id;
        value = other.value;
        other.id = -1;
    }

    Transportable_ & operator=(Transportable_ const & other)
    {
        value = other.value;
        return *this;
    }
};

int Transportable_::nextId = 0;

void assign(Transportable_ & target, Transportable_ & source)
{
    target.value = source.value;
}

void assign(Transportable_ & target, Transportable_ const & source)
{
    target.value = source.value;
}

void move(Transportable_ & target, Transportable_ & source)
{
    target.value = source.value;
    target.id = source.id;
    source.id = -1;
}

void set(Transportable_ & target, Transportable_ & source)
{
    target.value = source.value;
    target.id = source.id;
}

void set(Transportable_ & target, Transportable_ const & source)
{
    target.value = source.value;
    target.id = source.id;
}

// ==========================================================================
// Generic Test Helpers.
// ==========================================================================

// ==========================================================================
// Actual Tests Definitions.
// ==========================================================================

// --------------------------------------------------------------------------
// Tests for Base Pair.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_metafunctions)
{
    typedef Pair<int, unsigned> TPair;
    // Metafunction LENGTH
    unsigned l = LENGTH<TPair>::VALUE;
    SEQAN_ASSERT_EQ(l, 2u);
    // Metafunction Value
    bool b = IsSameType<typename Value<TPair, 1>::Type, int>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = IsSameType<typename Value<TPair, 2>::Type, unsigned>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // Metafunction Spec
    b = IsSameType<typename Spec<TPair>::Type, void>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // // Metafunction Key
    // b = IsSameType<typename Key<TPair>::Type, int>::VALUE;
    // SEQAN_ASSERT_EQ(b, true);
    // // Metafunction Cargo
    // b = IsSameType<typename Cargo<TPair>::Type, unsigned>::VALUE;
    // SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_constructors)
{
    typedef Pair<int, unsigned> TPair;

    // Default constructor.
    {
        TPair p;
    }

    // Normal Constructor Pair(i1, i2).
    {
        TPair p(1, 2);

        SEQAN_ASSERT_EQ(p.i1, 1);
        SEQAN_ASSERT_EQ(p.i2, 2u);
    }

    // Copy Constructor Pair(p).
    {
        TPair p(1, 2);
        TPair p2(p);

        SEQAN_ASSERT_EQ(p2.i1, 1);
        SEQAN_ASSERT_EQ(p2.i2, 2u);
    }

    // Conversion constructor from other pair.
    {
        Pair<int, int, Pack> p2(1, 2);
        TPair p(p2);

        SEQAN_ASSERT_EQ(p.i1, 1);
        SEQAN_ASSERT_EQ(p.i2, 2u);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_assign)
{
    // Test with ints.
    {
        Pair<int, int> p1;
        Pair<int, int> p2(1, 2);
        assign(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }

    // Test with Transportable_ type.
    {
        Pair<Transportable_, Transportable_> p1;
        Pair<Transportable_, Transportable_> p2(Transportable_(1), Transportable_(2));
        assign(p1, p2);
        SEQAN_ASSERT_NEQ(p1.i1.id, p2.i1.id);
        SEQAN_ASSERT_EQ(p1.i1.value, p2.i1.value);
        SEQAN_ASSERT_NEQ(p1.i2.id, p2.i2.id);
        SEQAN_ASSERT_EQ(p1.i2.value, p2.i2.value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_set)
{
    // Test with ints.
    {
        Pair<int, int> p1;
        Pair<int, int> p2(1, 2);
        set(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }

    // Test with Transportable_ type.
    {
        Pair<Transportable_, Transportable_> p1;
        Pair<Transportable_, Transportable_> p2(Transportable_(1), Transportable_(2));
        set(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1.id, p2.i1.id);
        SEQAN_ASSERT_EQ(p1.i1.value, p2.i1.value);
        SEQAN_ASSERT_EQ(p1.i2.id, p2.i2.id);
        SEQAN_ASSERT_EQ(p1.i2.value, p2.i2.value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_move)
{
    // Test with ints.
    {
        Pair<int, int> p1;
        Pair<int, int> p2(1, 2);
        move(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }

    // Test with Transportable_ type.
    {
        Pair<Transportable_, Transportable_> p1;
        Pair<Transportable_, Transportable_> p2(Transportable_(1), Transportable_(2));
        int p2I1Id = p2.i1.id;
        int p2I1Value = p2.i1.value;
        int p2I2Id = p2.i2.id;
        int p2I2Value = p2.i2.value;
        move(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1.id, p2I1Id);
        SEQAN_ASSERT_EQ(p1.i1.value, p2I1Value);
        SEQAN_ASSERT_EQ(p1.i2.id, p2I2Id);
        SEQAN_ASSERT_EQ(p1.i2.value, p2I2Value);
        SEQAN_ASSERT_EQ(p2.i1.id, -1);
        SEQAN_ASSERT_EQ(p2.i2.id, -1);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_value)
{
    SEQAN_ASSERT_FAIL("No reference in Pair<>");
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_get_value)
{
    Pair<int, int> p(-1, 1);
    SEQAN_ASSERT_EQ(getValueI1(p), -1);
    SEQAN_ASSERT_EQ(getValueI2(p), 1);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_assign_value)
{
    // Test with ints.
    {
        Pair<int, int> p1;
        Pair<int, int> p2(1, 2);

        assignValueI1(p1, p2.i1);
        assignValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }

    // Test with Transportable_ type.
    {
        Pair<Transportable_, Transportable_> p1;
        Pair<Transportable_, Transportable_> p2(Transportable_(1), Transportable_(2));
        int p1I1Id = p1.i1.id;
        int p1I2Id = p1.i2.id;
        int p2I1Id = p2.i1.id;
        int p2I1Value = p2.i1.value;
        int p2I2Id = p2.i2.id;
        int p2I2Value = p2.i2.value;

        assignValueI1(p1, p2.i1);
        assignValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1.id, p1I1Id);
        SEQAN_ASSERT_EQ(p1.i1.value, p2I1Value);
        SEQAN_ASSERT_EQ(p1.i2.id, p1I2Id);
        SEQAN_ASSERT_EQ(p1.i2.value, p2I2Value);
        SEQAN_ASSERT_EQ(p2.i1.id, p2I1Id);
        SEQAN_ASSERT_EQ(p2.i1.value, p2I1Value);
        SEQAN_ASSERT_EQ(p2.i2.id, p2I2Id);
        SEQAN_ASSERT_EQ(p2.i2.value, p2I2Value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_set_value)
{
    // Test with ints.
    {
        Pair<int, int> p1;
        Pair<int, int> p2(1, 2);

        setValueI1(p1, p2.i1);
        setValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }

    // Test with Transportable_ type.
    {
        Pair<Transportable_, Transportable_> p1;
        Pair<Transportable_, Transportable_> p2(Transportable_(1), Transportable_(2));
        int p2I1Id = p2.i1.id;
        int p2I1Value = p2.i1.value;
        int p2I2Id = p2.i2.id;
        int p2I2Value = p2.i2.value;

        setValueI1(p1, p2.i1);
        setValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1.id, p2I1Id);
        SEQAN_ASSERT_EQ(p1.i1.value, p2I1Value);
        SEQAN_ASSERT_EQ(p1.i2.id, p2I2Id);
        SEQAN_ASSERT_EQ(p1.i2.value, p2I2Value);
        SEQAN_ASSERT_EQ(p2.i1.id, p2I1Id);
        SEQAN_ASSERT_EQ(p2.i1.value, p2I1Value);
        SEQAN_ASSERT_EQ(p2.i2.id, p2I2Id);
        SEQAN_ASSERT_EQ(p2.i2.value, p2I2Value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_move_value)
{
    // Test with ints.
    {
        Pair<int, int> p1;
        Pair<int, int> p2(1, 2);

        moveValueI1(p1, p2.i1);
        moveValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }

    // Test with Transportable_ type.
    {
        Pair<Transportable_, Transportable_> p1;
        Pair<Transportable_, Transportable_> p2(Transportable_(1), Transportable_(2));
        int p2I1Id = p2.i1.id;
        int p2I1Value = p2.i1.value;
        int p2I2Id = p2.i2.id;
        int p2I2Value = p2.i2.value;

        moveValueI1(p1, p2.i1);
        moveValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1.id, p2I1Id);
        SEQAN_ASSERT_EQ(p1.i1.value, p2I1Value);
        SEQAN_ASSERT_EQ(p1.i2.id, p2I2Id);
        SEQAN_ASSERT_EQ(p1.i2.value, p2I2Value);
        SEQAN_ASSERT_EQ(p2.i1.id, -1);
        SEQAN_ASSERT_EQ(p2.i1.value, p2I1Value);
        SEQAN_ASSERT_EQ(p2.i2.id, -1);
        SEQAN_ASSERT_EQ(p2.i2.value, p2I2Value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_comparison_same_spec)
{
    Pair<int, int> p00(0, 0);
    Pair<int, int> p01(0, 1);
    Pair<int, int> p10(1, 0);

    SEQAN_ASSERT(p00 == p00);
    SEQAN_ASSERT_NOT(p01 == p10);

    SEQAN_ASSERT(p00 != p10);
    SEQAN_ASSERT_NOT(p01 != p01);

    SEQAN_ASSERT(p01 > p00);
    SEQAN_ASSERT_NOT(p01 > p10);

    SEQAN_ASSERT(p01 >= p01);
    SEQAN_ASSERT(p01 >= p00);
    SEQAN_ASSERT_NOT(p01 >= p10);

    SEQAN_ASSERT(p00 <= p01);
    SEQAN_ASSERT(p00 <= p00);
    SEQAN_ASSERT_NOT(p10 <= p01);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_comparison_different_spec)
{
    Pair<int, int64_t> p00(0, 0);
    Pair<int, int64_t, Pack> p01(0, 1);
    Pair<int, short int, BitPacked<20, 12> > p10(1, 0);

    SEQAN_ASSERT(p00 == p00);
    SEQAN_ASSERT_NOT(p01 == p10);

    SEQAN_ASSERT(p00 != p10);
    SEQAN_ASSERT_NOT(p01 != p01);

    SEQAN_ASSERT(p01 > p00);
    SEQAN_ASSERT_NOT(p01 > p10);

    SEQAN_ASSERT(p01 >= p01);
    SEQAN_ASSERT(p01 >= p00);
    SEQAN_ASSERT_NOT(p01 >= p10);

    SEQAN_ASSERT(p00 <= p01);
    SEQAN_ASSERT(p00 <= p00);
    SEQAN_ASSERT_NOT(p10 <= p01);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_base_stream_output)
{
    std::stringstream s;
    Pair<int, int> p(-1, 1);
    s << p;
    SEQAN_ASSERT_EQ(s.str(), "< -1 , 1 >");
}

// --------------------------------------------------------------------------
// Tests for Pack Pair.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_metafunctions)
{
    typedef Pair<int, unsigned, Pack> TPair;
    // Metafunction LENGTH
    unsigned l = LENGTH<TPair>::VALUE;
    SEQAN_ASSERT_EQ(l, 2u);
    // Metafunction Value
    bool b = IsSameType<typename Value<TPair, 1>::Type, int>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = IsSameType<typename Value<TPair, 2>::Type, unsigned>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // Metafunction Spec
    b = IsSameType<typename Spec<TPair>::Type, Pack>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // // Metafunction Key
    // b = IsSameType<typename Key<TPair>::Type, int>::VALUE;
    // SEQAN_ASSERT_EQ(b, true);
    // // Metafunction Cargo
    // b = IsSameType<typename Cargo<TPair>::Type, unsigned>::VALUE;
    // SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_constructors)
{
    typedef Pair<int, unsigned, Pack> TPair;

    // Default constructor.
    {
        TPair p{};

        SEQAN_ASSERT_EQ(p.i1, 0);
        SEQAN_ASSERT_EQ(p.i2, 0u);
    }

    // Normal Constructor Pair(i1, i2).
    {
        TPair p(1, 2);

        SEQAN_ASSERT_EQ(p.i1, 1);
        SEQAN_ASSERT_EQ(p.i2, 2u);
    }

    // Copy Constructor Pair(p).
    {
        TPair p(1, 2);
        TPair p2(p);

        SEQAN_ASSERT_EQ(p2.i1, 1);
        SEQAN_ASSERT_EQ(p2.i2, 2u);
    }

    // Conversion constructor from other pair.
    {
        Pair<int, int> p2(1, 2);
        TPair p(p2);

        SEQAN_ASSERT_EQ(p.i1, 1);
        SEQAN_ASSERT_EQ(p.i2, 2u);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_assign)
{
    // Test with ints.
    {
        Pair<int, int, Pack> p1;
        Pair<int, int, Pack> p2(1, 2);
        assign(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_set)
{
    // Test with ints.
    {
        Pair<int, int, Pack> p1;
        Pair<int, int, Pack> p2(1, 2);
        set(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_move)
{
    // Test with ints.
    {
        Pair<int, int, Pack> p1;
        Pair<int, int, Pack> p2(1, 2);
        move(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_stream_output)
{
    std::stringstream s;
    Pair<int, int, Pack> p(-1, 1);
    s << p;
    SEQAN_ASSERT_EQ(s.str(), "< -1 , 1 >");

}

// SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_value)

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_get_value)
{
    Pair<int, int, Pack> p(-1, 1);
    SEQAN_ASSERT_EQ(getValueI1(p), -1);
    SEQAN_ASSERT_EQ(getValueI2(p), 1);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_assign_value)

{
    // Test with ints.
    {
        Pair<int, int, Pack> p1;
        Pair<int, int, Pack> p2(1, 2);

        assignValueI1(p1, p2.i1);
        assignValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_set_value)
{
    // Test with ints.
    {
        Pair<int, int, Pack> p1;
        Pair<int, int, Pack> p2(1, 2);

        setValueI1(p1, p2.i1);
        setValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_move_value)
{
    // Test with ints.
    {
        Pair<int, int, Pack> p1;
        Pair<int, int, Pack> p2(1, 2);

        int i1 = p2.i1;
        int i2 = p2.i2;

        moveValueI1(p1, i1);
        moveValueI2(p1, i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_comparison_same_spec)
{
    Pair<int, int, Pack> p00(0, 0);
    Pair<int, int, Pack> p01(0, 1);
    Pair<int, int, Pack> p10(1, 0);

    SEQAN_ASSERT(p00 == p00);
    SEQAN_ASSERT_NOT(p01 == p10);

    SEQAN_ASSERT(p00 != p10);
    SEQAN_ASSERT_NOT(p01 != p01);

    SEQAN_ASSERT(p01 > p00);
    SEQAN_ASSERT_NOT(p01 > p10);

    SEQAN_ASSERT(p01 >= p01);
    SEQAN_ASSERT(p01 >= p00);
    SEQAN_ASSERT_NOT(p01 >= p10);

    SEQAN_ASSERT(p00 <= p01);
    SEQAN_ASSERT(p00 <= p00);
    SEQAN_ASSERT_NOT(p10 <= p01);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_packed_comparison_different_spec)
{
    Pair<int, int64_t> p00(0, 0);
    Pair<int, int64_t, Pack> p01(0, 1);
    Pair<int, short int, BitPacked<20, 12> > p10(1, 0);

    SEQAN_ASSERT(p00 == p00);
    SEQAN_ASSERT_NOT(p01 == p10);

    SEQAN_ASSERT(p00 != p10);
    SEQAN_ASSERT_NOT(p01 != p01);

    SEQAN_ASSERT(p01 > p00);
    SEQAN_ASSERT_NOT(p01 > p10);

    SEQAN_ASSERT(p01 >= p01);
    SEQAN_ASSERT(p01 >= p00);
    SEQAN_ASSERT_NOT(p01 >= p10);

    SEQAN_ASSERT(p00 <= p01);
    SEQAN_ASSERT(p00 <= p00);
    SEQAN_ASSERT_NOT(p10 <= p01);
}

// --------------------------------------------------------------------------
// Tests for Bit Pack Pair.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_metafunctions)
{
    typedef BitPacked<16, 16> TSpec;
    typedef Pair<int, unsigned, TSpec> TPair;
    // Metafunction LENGTH
    unsigned l = LENGTH<TPair>::VALUE;
    SEQAN_ASSERT_EQ(l, 2u);
    // Metafunction Value
    bool b = IsSameType<typename Value<TPair, 1>::Type, int>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = IsSameType<typename Value<TPair, 2>::Type, unsigned>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // Metafunction Spec
    b = IsSameType<typename Spec<TPair>::Type, TSpec>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // // Metafunction Key
    // b = IsSameType<typename Key<TPair>::Type, int>::VALUE;
    // SEQAN_ASSERT_EQ(b, true);
    // // Metafunction Cargo
    // b = IsSameType<typename Cargo<TPair>::Type, unsigned>::VALUE;
    // SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_constructors)
{
    typedef Pair<int, unsigned, BitPacked<16, 16> > TPair;

    // Default constructor.
    {
        TPair p;
    }

    // Normal Constructor Pair(i1, i2).
    {
        TPair p(1, 2);

        SEQAN_ASSERT_EQ(p.i1, 1);
        SEQAN_ASSERT_EQ(p.i2, 2u);
    }

    // Copy Constructor Pair(p).
    {
        TPair p(1, 2);
        TPair p2(p);

        SEQAN_ASSERT_EQ(p2.i1, 1);
        SEQAN_ASSERT_EQ(p2.i2, 2u);
    }

    // Conversion constructor from other pair.
    {
        Pair<int, int, Pack> p2(1, 2);
        TPair p(p2);

        SEQAN_ASSERT_EQ(p.i1, 1);
        SEQAN_ASSERT_EQ(p.i2, 2u);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_assign)
{
    // Test with ints.
    {
        Pair<int, int, BitPacked<16, 16> > p1;
        Pair<int, int, BitPacked<16, 16> > p2(1, 2);
        assign(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_set)
{
    // Test with ints.
    {
        Pair<int, int, BitPacked<16, 16> > p1;
        Pair<int, int, BitPacked<16, 16> > p2(1, 2);
        set(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }

}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_move)
{
    // Test with ints.
    {
        Pair<int, int, BitPacked<16, 16> > p1;
        Pair<int, int, BitPacked<16, 16> > p2(1, 2);
        move(p1, p2);
        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

// SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_value)

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_get_value)
{
    Pair<int, int, BitPacked<16, 16> > p(-1, 1);
    SEQAN_ASSERT_EQ(getValueI1(p), -1);
    SEQAN_ASSERT_EQ(getValueI2(p), 1);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_assign_value)
{
    // Test with ints.
    {
        Pair<int, int, BitPacked<16, 16> > p1;
        Pair<int, int, BitPacked<16, 16> > p2(1, 2);

        assignValueI1(p1, p2.i1);
        assignValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_set_value)
{
    // Test with ints.
    {
        Pair<int, int, BitPacked<16, 16> > p1;
        Pair<int, int, BitPacked<16, 16> > p2(1, 2);

        setValueI1(p1, p2.i1);
        setValueI2(p1, p2.i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_move_value)
{
    // Test with ints.
    {
        Pair<int, int, BitPacked<16, 16> > p1;
        Pair<int, int, BitPacked<16, 16> > p2(1, 2);

        int i1 = p2.i1;
        int i2 = p2.i2;

        moveValueI1(p1, i1);
        moveValueI2(p1, i2);

        SEQAN_ASSERT_EQ(p1.i1, p2.i1);
        SEQAN_ASSERT_EQ(p1.i2, p2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_comparison_same_spec)
{
    Pair<int, int, BitPacked<16, 16> > p00(0, 0);
    Pair<int, int, BitPacked<16, 16> > p01(0, 1);
    Pair<int, int, BitPacked<16, 16> > p10(1, 0);

    SEQAN_ASSERT(p00 == p00);
    SEQAN_ASSERT_NOT(p01 == p10);

    SEQAN_ASSERT(p00 != p10);
    SEQAN_ASSERT_NOT(p01 != p01);

    SEQAN_ASSERT(p01 > p00);
    SEQAN_ASSERT_NOT(p01 > p10);

    SEQAN_ASSERT(p01 >= p01);
    SEQAN_ASSERT(p01 >= p00);
    SEQAN_ASSERT_NOT(p01 >= p10);

    SEQAN_ASSERT(p00 <= p01);
    SEQAN_ASSERT(p00 <= p00);
    SEQAN_ASSERT_NOT(p10 <= p01);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_comparison_different_spec)
{
    Pair<int, int64_t> p00(0, 0);
    Pair<int, int64_t, Pack> p01(0, 1);
    Pair<int, short int, BitPacked<20, 12> > p10(1, 0);

    SEQAN_ASSERT(p00 == p00);
    SEQAN_ASSERT_NOT(p01 == p10);

    SEQAN_ASSERT(p00 != p10);
    SEQAN_ASSERT_NOT(p01 != p01);

    SEQAN_ASSERT(p01 > p00);
    SEQAN_ASSERT_NOT(p01 > p10);

    SEQAN_ASSERT(p01 >= p01);
    SEQAN_ASSERT(p01 >= p00);
    SEQAN_ASSERT_NOT(p01 >= p10);

    SEQAN_ASSERT(p00 <= p01);
    SEQAN_ASSERT(p00 <= p00);
    SEQAN_ASSERT_NOT(p10 <= p01);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_pair_bit_packed_stream_output)
{
    std::stringstream s;
    Pair<int, int> p(-1, 1);
    s << p;
    SEQAN_ASSERT_EQ(s.str(), "< -1 , 1 >");
}

// --------------------------------------------------------------------------
// Tests for Base Triple.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_metafunctions)
{
    typedef Triple<int, unsigned, double> TTriple;
    // Metafunction LENGTH
    unsigned l = LENGTH<TTriple>::VALUE;
    SEQAN_ASSERT_EQ(l, 3u);
    // Metafunction Value
    bool b = IsSameType<typename Value<TTriple, 1>::Type, int>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = IsSameType<typename Value<TTriple, 2>::Type, unsigned>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = IsSameType<typename Value<TTriple, 3>::Type, double>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // Metafunction Spec
    b = IsSameType<typename Spec<TTriple>::Type, void>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_constructors)
{
    typedef Triple<int, unsigned, double, Pack> TTriple;

    // Default constructor.
    {
        TTriple t;
    }

    // Normal Constructor Triple(i1, i2).
    {
        TTriple t(1, 2, 3);

        SEQAN_ASSERT_EQ(t.i1, 1);
        SEQAN_ASSERT_EQ(t.i2, 2u);
        SEQAN_ASSERT_EQ(t.i3, 3.);
    }

    // Copy Constructor Triple(p).
    {
        TTriple t(1, 2, 3);
        TTriple t2(t);

        SEQAN_ASSERT_EQ(t2.i1, 1);
        SEQAN_ASSERT_EQ(t2.i2, 2u);
        SEQAN_ASSERT_EQ(t2.i3, 3.);
    }

    // Conversion constructor from other triple.
    {
        Triple<int, int, int> t2(1, 2, 3);
        TTriple t(t2);

        SEQAN_ASSERT_EQ(t.i1, 1);
        SEQAN_ASSERT_EQ(t.i2, 2u);
        SEQAN_ASSERT_EQ(t.i3, 3.);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_assign)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);
        assign(t1, t2);
        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_set)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);
        set(t1, t2);
        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_move)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);
        move(t1, t2);
        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_get_value)
{
    Triple<int, int, int> t(-1, 1, 0);
    SEQAN_ASSERT_EQ(getValueI1(t), -1);
    SEQAN_ASSERT_EQ(getValueI2(t), 1);
    SEQAN_ASSERT_EQ(getValueI3(t), 0);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_assign_value)
{
    // Test with ints.
    {
        Triple<int, int, int> t1;
        Triple<int, int, int> t2(1, 2, 3);

        assignValueI1(t1, t2.i1);
        assignValueI2(t1, t2.i2);
        assignValueI3(t1, t2.i3);

        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_set_value)
{
    // Test with ints.
    {
        Triple<int, int, int> t1;
        Triple<int, int, int> t2(1, 2, 3);

        setValueI1(t1, t2.i1);
        setValueI2(t1, t2.i2);
        setValueI3(t1, t2.i3);

        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_move_value)
{
    // Test with ints.
    {
        Triple<int, int, int> t1;
        Triple<int, int, int> t2(1, 2, 3);

        moveValueI1(t1, t2.i1);
        moveValueI2(t1, t2.i2);
        moveValueI3(t1, t2.i3);

        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_comparison_same_spec)
{
    Triple<int, int, int> t001(0, 0, 1);
    Triple<int, int, int> t010(0, 1, 0);
    Triple<int, int, int> t100(1, 0, 0);

    SEQAN_ASSERT(t001 == t001);
    SEQAN_ASSERT_NOT(t010 == t100);

    SEQAN_ASSERT(t001 != t100);
    SEQAN_ASSERT_NOT(t010 != t010);

    SEQAN_ASSERT(t010 > t001);
    SEQAN_ASSERT_NOT(t010 > t100);

    SEQAN_ASSERT(t010 >= t010);
    SEQAN_ASSERT(t010 >= t001);
    SEQAN_ASSERT_NOT(t010 >= t100);

    SEQAN_ASSERT(t001 <= t010);
    SEQAN_ASSERT(t001 <= t001);
    SEQAN_ASSERT_NOT(t100 <= t010);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_comparison_different_spec)
{
    Triple<int, int64_t, char> t001(0, 0, 1);
    Triple<int, int64_t, char, Pack> t010(0, 1, 0);
    Triple<int, short int, char > t100(1, 0, 0);

    SEQAN_ASSERT(t001 == t001);
    SEQAN_ASSERT_NOT(t010 == t100);

    SEQAN_ASSERT(t001 != t100);
    SEQAN_ASSERT_NOT(t010 != t010);

    SEQAN_ASSERT(t010 > t001);
    SEQAN_ASSERT_NOT(t010 > t100);

    SEQAN_ASSERT(t010 >= t010);
    SEQAN_ASSERT(t010 >= t001);
    SEQAN_ASSERT_NOT(t010 >= t100);

    SEQAN_ASSERT(t001 <= t010);
    SEQAN_ASSERT(t001 <= t001);
    SEQAN_ASSERT_NOT(t100 <= t010);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_base_stream_output)
{
    std::stringstream s;
    Triple<int, int, int> t(-1, 1, 0);
    s << t;
    SEQAN_ASSERT_EQ(s.str(), "< -1 , 1 , 0 >");
}

// --------------------------------------------------------------------------
// Tests for Pack Triple.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_metafunctions)
{
    typedef Triple<int, unsigned, double, Pack> TTriple;
    // Metafunction LENGTH
    unsigned l = LENGTH<TTriple>::VALUE;
    SEQAN_ASSERT_EQ(l, 3u);
    // Metafunction Value
    bool b = IsSameType<typename Value<TTriple, 1>::Type, int>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = IsSameType<typename Value<TTriple, 2>::Type, unsigned>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = IsSameType<typename Value<TTriple, 3>::Type, double>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // Metafunction Spec
    b = IsSameType<typename Spec<TTriple>::Type, Pack>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_constructors)
{
    typedef Triple<int, unsigned, double, Pack> TTriple;

    // Default constructor.
    {
        TTriple p;
    }

    // Normal Constructor Triple(i1, i2).
    {
        TTriple t(1, 2, 3);

        SEQAN_ASSERT_EQ(t.i1, 1);
        SEQAN_ASSERT_EQ(t.i2, 2u);
        SEQAN_ASSERT_EQ(t.i3, 3.);
    }

    // Copy Constructor Triple(p).
    {
        TTriple t(1, 2, 3);
        TTriple t2(t);

        SEQAN_ASSERT_EQ(t2.i1, 1);
        SEQAN_ASSERT_EQ(t2.i2, 2u);
        SEQAN_ASSERT_EQ(t2.i3, 3.);
    }

    // Conversion constructor from other triple.
    {
        Triple<int, int, int> t2(1, 2, 3);
        TTriple t(t2);

        SEQAN_ASSERT_EQ(t.i1, 1);
        SEQAN_ASSERT_EQ(t.i2, 2u);
        SEQAN_ASSERT_EQ(t.i3, 3.);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_assign)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);
        assign(t1, t2);
        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_set)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);
        set(t1, t2);
        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_move)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);
        move(t1, t2);
        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_stream_output)
{
    std::stringstream s;
    Triple<int, int, int, Pack> p(-1, 1, 0);
    s << p;
    SEQAN_ASSERT_EQ(s.str(), "< -1 , 1 , 0 >");

}

// SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_value)

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_get_value)
{
    Triple<int, int, int, Pack> t(-1, 1, 0);
    SEQAN_ASSERT_EQ(getValueI1(t), -1);
    SEQAN_ASSERT_EQ(getValueI2(t), 1);
    SEQAN_ASSERT_EQ(getValueI3(t), 0);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_assign_value)

{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);

        assignValueI1(t1, t2.i1);
        assignValueI2(t1, t2.i2);
        assignValueI3(t1, t2.i3);

        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_set_value)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);

        setValueI1(t1, t2.i1);
        setValueI2(t1, t2.i2);
        setValueI3(t1, t2.i3);

        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
        SEQAN_ASSERT_EQ(t1.i3, t2.i3);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_move_value)
{
    // Test with ints.
    {
        Triple<int, int, int, Pack> t1;
        Triple<int, int, int, Pack> t2(1, 2, 3);

        int i1 = t2.i1;
        int i2 = t2.i2;
        int i3 = t2.i3;

        moveValueI1(t1, i1);
        moveValueI2(t1, i2);
        moveValueI3(t1, i3);

        SEQAN_ASSERT_EQ(t1.i1, t2.i1);
        SEQAN_ASSERT_EQ(t1.i2, t2.i2);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_comparison_same_spec)
{
    Triple<int, int, int, Pack> t001(0, 0, 1);
    Triple<int, int, int, Pack> t010(0, 1, 0);
    Triple<int, int, int, Pack> t100(1, 0, 0);

    SEQAN_ASSERT(t001 == t001);
    SEQAN_ASSERT_NOT(t010 == t100);

    SEQAN_ASSERT(t001 != t100);
    SEQAN_ASSERT_NOT(t010 != t010);

    SEQAN_ASSERT(t010 > t001);
    SEQAN_ASSERT_NOT(t010 > t100);

    SEQAN_ASSERT(t010 >= t010);
    SEQAN_ASSERT(t010 >= t001);
    SEQAN_ASSERT_NOT(t010 >= t100);

    SEQAN_ASSERT(t001 <= t010);
    SEQAN_ASSERT(t001 <= t001);
    SEQAN_ASSERT_NOT(t100 <= t010);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_triple_packed_comparison_different_spec)
{
    Triple<int, int64_t, char> t001(0, 0, 1);
    Triple<int, int64_t, char, Pack> t010(0, 1, 0);
    Triple<int, short int, char, BitPacked<20, 12> > t100(1, 0, 0);

    SEQAN_ASSERT(t001 == t001);
    SEQAN_ASSERT_NOT(t010 == t100);

    SEQAN_ASSERT(t001 != t100);
    SEQAN_ASSERT_NOT(t010 != t010);

    SEQAN_ASSERT(t010 > t001);
    SEQAN_ASSERT_NOT(t010 > t100);

    SEQAN_ASSERT(t010 >= t010);
    SEQAN_ASSERT(t010 >= t001);
    SEQAN_ASSERT_NOT(t010 >= t100);

    SEQAN_ASSERT(t001 <= t010);
    SEQAN_ASSERT(t001 <= t001);
    SEQAN_ASSERT_NOT(t100 <= t010);
}

// --------------------------------------------------------------------------
// Tests for Base Tuple.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_metafunctions)
{
    typedef Tuple<int, 2> TTuple;
    // Metafunction LENGTH
    unsigned l = LENGTH<TTuple>::VALUE;
    SEQAN_ASSERT_EQ(l, 2u);
    // Metafunction Value
    bool b = IsSameType<typename Value<TTuple>::Type, int>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // Metafunction Spec
    b = IsSameType<typename Spec<TTuple>::Type, void>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_constructors)
{
    // Default constructor.
    Tuple<int, 2> t;
    (void)t;
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_assign)
{
    typedef Tuple<int, 2> TTuple;
    TTuple t1;
    TTuple t2;

    t1.i[0] = 1;
    t1.i[1] = -1;

    assign(t2, t1);

    SEQAN_ASSERT_EQ(t1.i[0], 1);
    SEQAN_ASSERT_EQ(t1.i[1], -1);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_set)
{
    // Test with ints.
    {
        Tuple<int, 2> t1;
        Tuple<int, 2> t2;
        t2.i[0] = 1;
        t2.i[1] = 2;

        set(t1, t2);

        SEQAN_ASSERT_EQ(t1.i[0], t2.i[0]);
        SEQAN_ASSERT_EQ(t1.i[1], t2.i[1]);
    }

    // Test with Transportable_ type.
    {
        Tuple<Transportable_, 2> t1;
        Tuple<Transportable_, 2> t2;
        t2.i[0] = Transportable_(1);
        t2.i[1] = Transportable_(2);

        set(t1, t2);

        SEQAN_ASSERT_EQ(t1.i[0].id, t2.i[0].id);
        SEQAN_ASSERT_EQ(t1.i[0].value, t2.i[0].value);
        SEQAN_ASSERT_EQ(t1.i[1].id, t2.i[1].id);
        SEQAN_ASSERT_EQ(t1.i[1].value, t2.i[1].value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_move)
{
    // Test with ints.
    {
        Tuple<int, 2> t1;
        Tuple<int, 2> t2;
        t2.i[0] = 1;
        t2.i[1] = 2;

        move(t1, t2);

        SEQAN_ASSERT_EQ(t1.i[0], t2.i[0]);
        SEQAN_ASSERT_EQ(t1.i[1], t2.i[1]);
    }

    // Test with Transportable_ type.
    {
        Tuple<Transportable_, 2> t1;
        Tuple<Transportable_, 2> t2;
        t2.i[0] = Transportable_(1);
        t2.i[1] = Transportable_(2);
        int t2I1Id = t2.i[0].id;
        int t2I1Value = t2.i[0].value;
        int t2I2Id = t2.i[1].id;
        int t2I2Value = t2.i[1].value;

        move(t1, t2);

        SEQAN_ASSERT_EQ(t1.i[0].id, t2I1Id);
        SEQAN_ASSERT_EQ(t1.i[0].value, t2I1Value);
        SEQAN_ASSERT_EQ(t1.i[1].id, t2I2Id);
        SEQAN_ASSERT_EQ(t1.i[1].value, t2I2Value);
        SEQAN_ASSERT_EQ(t2.i[0].id, -1);
        SEQAN_ASSERT_EQ(t2.i[1].id, -1);
    }
}

// SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_value);  // TODO(holtgrew): Need proxy for this.

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_get_value)
{
    Tuple<int, 2> t;
    t.i[0] = 1;
    t.i[1] = -1;
    SEQAN_ASSERT_EQ(getValue(t, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t, 1), -1);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_assign_value)
{
    // Test with ints.
    {
        Tuple<int, 2> t1;
        Tuple<int, 2> t2;
        t2.i[0] = 1;
        t2.i[1] = -1;

        assignValue(t1, 0, t2.i[0]);
        assignValue(t1, 1, t2.i[1]);

        SEQAN_ASSERT_EQ(t1.i[0], t2.i[0]);
        SEQAN_ASSERT_EQ(t1.i[1], t2.i[1]);
    }

    // Test with Transportable_ type.
    {
        Tuple<Transportable_, 2> t1;
        Tuple<Transportable_, 2> t2;
        t2.i[0] = Transportable_(1);
        t2.i[1] = Transportable_(-1);
        int t1I1Id = t1.i[0].id;
        int t1I2Id = t1.i[1].id;
        int t2I1Id = t2.i[0].id;
        int t2I1Value = t2.i[0].value;
        int t2I2Id = t2.i[1].id;
        int t2I2Value = t2.i[1].value;

        assignValue(t1, 0, t2.i[0]);
        assignValue(t1, 1, t2.i[1]);

        SEQAN_ASSERT_EQ(t1.i[0].id, t1I1Id);
        SEQAN_ASSERT_EQ(t1.i[0].value, t2I1Value);
        SEQAN_ASSERT_EQ(t1.i[1].id, t1I2Id);
        SEQAN_ASSERT_EQ(t1.i[1].value, t2I2Value);
        SEQAN_ASSERT_EQ(t2.i[0].id, t2I1Id);
        SEQAN_ASSERT_EQ(t2.i[0].value, t2I1Value);
        SEQAN_ASSERT_EQ(t2.i[1].id, t2I2Id);
        SEQAN_ASSERT_EQ(t2.i[1].value, t2I2Value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_set_value)
{
    // Test with ints.
    {
        Tuple<int, 2> t1;
        Tuple<int, 2> t2;
        t2.i[0] = 1;
        t2.i[1] = 2;

        setValue(t1, 0, t2.i[0]);
        setValue(t1, 1, t2.i[1]);

        SEQAN_ASSERT_EQ(t1.i[0], t2.i[0]);
        SEQAN_ASSERT_EQ(t1.i[1], t2.i[1]);
    }

    // Test with Transportable_ type.
    {
        Tuple<Transportable_, 2> t1;
        Tuple<Transportable_, 2> t2;
        t2.i[0] = Transportable_(1);
        t2.i[1] = Transportable_(2);
        int t2I1Id = t2.i[0].id;
        int t2I1Value = t2.i[0].value;
        int t2I2Id = t2.i[1].id;
        int t2I2Value = t2.i[1].value;

        setValue(t1, 0, t2.i[0]);
        setValue(t1, 1, t2.i[1]);

        SEQAN_ASSERT_EQ(t1.i[0].id, t2I1Id);
        SEQAN_ASSERT_EQ(t1.i[0].value, t2I1Value);
        SEQAN_ASSERT_EQ(t1.i[1].id, t2I2Id);
        SEQAN_ASSERT_EQ(t1.i[1].value, t2I2Value);
        SEQAN_ASSERT_EQ(t2.i[0].id, t2I1Id);
        SEQAN_ASSERT_EQ(t2.i[0].value, t2I1Value);
        SEQAN_ASSERT_EQ(t2.i[1].id, t2I2Id);
        SEQAN_ASSERT_EQ(t2.i[1].value, t2I2Value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_move_value)
{
    // Test with ints.
    {
        Tuple<int, 2> t1;
        Tuple<int, 2> t2;
        t2.i[0] = 1;
        t2.i[1] = 2;

        setValue(t1, 0, t2.i[0]);
        setValue(t1, 1, t2.i[1]);

        SEQAN_ASSERT_EQ(t1.i[0], t2.i[0]);
        SEQAN_ASSERT_EQ(t1.i[1], t2.i[1]);
    }

    // Test with Transportable_ type.
    {
        Tuple<Transportable_, 2> t1;
        Tuple<Transportable_, 2> t2;
        t2.i[0] = Transportable_(1);
        t2.i[1] = Transportable_(2);
        int t2I1Id = t2.i[0].id;
        int t2I1Value = t2.i[0].value;
        int t2I2Id = t2.i[1].id;
        int t2I2Value = t2.i[1].value;

        moveValue(t1, 0, t2.i[0]);
        moveValue(t1, 1, t2.i[1]);

        SEQAN_ASSERT_EQ(t1.i[0].id, t2I1Id);
        SEQAN_ASSERT_EQ(t1.i[0].value, t2I1Value);
        SEQAN_ASSERT_EQ(t1.i[1].id, t2I2Id);
        SEQAN_ASSERT_EQ(t1.i[1].value, t2I2Value);
        SEQAN_ASSERT_EQ(t2.i[0].id, -1);
        SEQAN_ASSERT_EQ(t2.i[0].value, t2I1Value);
        SEQAN_ASSERT_EQ(t2.i[1].id, -1);
        SEQAN_ASSERT_EQ(t2.i[1].value, t2I2Value);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_shift_left)
{
    // Test with ints.
    {
        Tuple<int, 3> t1;
        t1.i[0] = 1;
        t1.i[1] = 2;
        t1.i[2] = 3;

        shiftLeft(t1);

        SEQAN_ASSERT_EQ(t1.i[0], 2);
        SEQAN_ASSERT_EQ(t1.i[1], 3);
        SEQAN_ASSERT_EQ(t1.i[1], 3);
    }

    // TODO(holtgrew): Test with transportable after fixing set/move/assign here.
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_shift_right)
{
    // Test with ints.
    {
        Tuple<int, 3> t1;
        t1.i[0] = 1;
        t1.i[1] = 2;
        t1.i[2] = 3;

        shiftRight(t1);

        SEQAN_ASSERT_EQ(t1.i[0], 1);
        SEQAN_ASSERT_EQ(t1.i[1], 1);
        SEQAN_ASSERT_EQ(t1.i[2], 2);
    }

    // TODO(holtgrew): Test with transportable after fixing set/move/assign here.
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_clear)
{
    // Test with ints.
    {
        Tuple<int, 3> t1;
        t1.i[0] = 1;
        t1.i[1] = 2;
        t1.i[2] = 3;

        clear(t1);

        SEQAN_ASSERT_EQ(t1.i[0], 0);
        SEQAN_ASSERT_EQ(t1.i[1], 0);
        SEQAN_ASSERT_EQ(t1.i[2], 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_length)
{
    // Test with ints.
    {
        Tuple<int, 3> t1;
        SEQAN_ASSERT_EQ(length(t1), 3u);
    }
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_comparison_same_spec)
{
    Tuple<int, 2> p00; p00.i[0] = 0; p00.i[1] = 0;
    Tuple<int, 2> p01; p01.i[0] = 0; p01.i[1] = 1;
    Tuple<int, 2> p10; p10.i[0] = 1; p10.i[1] = 0;

    SEQAN_ASSERT(p00 == p00);
    SEQAN_ASSERT_NOT(p01 == p10);

    SEQAN_ASSERT(p00 != p10);
    SEQAN_ASSERT_NOT(p01 != p01);

    SEQAN_ASSERT(p01 > p00);
    SEQAN_ASSERT_NOT(p01 > p10);

    SEQAN_ASSERT(p01 >= p01);
    SEQAN_ASSERT(p01 >= p00);
    SEQAN_ASSERT_NOT(p01 >= p10);

    SEQAN_ASSERT(p00 <= p01);
    SEQAN_ASSERT(p00 <= p00);
    SEQAN_ASSERT_NOT(p10 <= p01);
}

// SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_comparison_different_spec)  // Not worth the bother, probably.

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_base_stream_output)
{
    Tuple<int, 3> t1;
    t1.i[0] = 1;
    t1.i[1] = 2;
    t1.i[2] = 3;

    std::stringstream s;
    s << t1;
    SEQAN_ASSERT_EQ(s.str(), "[1 2 3]");  // TODO(holtgrew): Why not "< 1 , 2 , 3 >", as for triples?
}

// --------------------------------------------------------------------------
// Tests for Bit-Pack Tuple.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_metafunctions)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    // Metafunction LENGTH
    unsigned l = LENGTH<TTuple>::VALUE;
    SEQAN_ASSERT_EQ(l, 2u);
    // Metafunction Value
    bool b = IsSameType<typename Value<TTuple>::Type, char>::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    // Metafunction Spec
    b = IsSameType<typename Spec<TTuple>::Type, BitPacked<> >::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_constructors)
{
    // Default constructor.
    Tuple<char, 2, BitPacked<> > t;
    (void)t;
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_assign)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    TTuple t1;
    TTuple t2;

    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    assign(t2, t1);

    SEQAN_ASSERT_EQ(getValue(t2, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t2, 1), 2);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_set)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    TTuple t1;
    TTuple t2;

    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    set(t2, t1);

    SEQAN_ASSERT_EQ(getValue(t2, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t2, 1), 2);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_move)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    TTuple t1;
    TTuple t2;

    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    move(t2, t1);

    SEQAN_ASSERT_EQ(getValue(t2, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t2, 1), 2);
}

// SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_value);  // TODO(holtgrew): Need proxy for this.

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_get_value)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    TTuple t1;

    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    SEQAN_ASSERT_EQ(getValue(t1, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t1, 1), 2);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_assign_value)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    TTuple t1;

    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    SEQAN_ASSERT_EQ(getValue(t1, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t1, 1), 2);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_set_value)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    TTuple t1;

    clear(t1);
    setValue(t1, 0, 1);
    setValue(t1, 1, 2);

    SEQAN_ASSERT_EQ(getValue(t1, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t1, 1), 2);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_move_value)
{
    typedef Tuple<char, 2, BitPacked<> > TTuple;
    TTuple t1;

    clear(t1);
    moveValue(t1, 0, 1);
    moveValue(t1, 1, 2);

    SEQAN_ASSERT_EQ(getValue(t1, 0), 1);
    SEQAN_ASSERT_EQ(getValue(t1, 1), 2);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_shift_left)
{
    Tuple<char, 2, BitPacked<> > t1;
    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    shiftLeft(t1);

    SEQAN_ASSERT_EQ(getValue(t1, 0), 2);
    SEQAN_ASSERT_EQ(getValue(t1, 1), 0);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_shift_right)
{
    Tuple<char, 2, BitPacked<> > t1;
    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    shiftRight(t1);

    SEQAN_ASSERT_EQ(getValue(t1, 0), 0);
    SEQAN_ASSERT_EQ(getValue(t1, 1), 1);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_clear)
{
    Tuple<char, 2, BitPacked<> > t1;
    clear(t1);
    assignValue(t1, 0, 1);
    assignValue(t1, 1, 2);

    clear(t1);

    SEQAN_ASSERT_EQ(getValue(t1, 0), 0);
    SEQAN_ASSERT_EQ(getValue(t1, 1), 0);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_length)
{
    Tuple<char, 2, BitPacked<> > t1;
    SEQAN_ASSERT_EQ(length(t1), 2u);
}

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_comparison_same_spec)
{
    Tuple<char, 2, BitPacked<> > t00;
    clear(t00);
       assignValue(t00, 0, 0);
    assignValue(t00, 1, 0);

    Tuple<char, 2, BitPacked<> > t01;
    clear(t01);
    assignValue(t01, 0, 0);
    assignValue(t01, 1, 1);

    Tuple<char, 2, BitPacked<> > t10;
    clear(t10);
    assignValue(t10, 0, 1);
    assignValue(t10, 1, 0);

    SEQAN_ASSERT(t00 == t00);
    SEQAN_ASSERT_NOT(t01 == t10);

    SEQAN_ASSERT(t00 != t10);
    SEQAN_ASSERT_NOT(t01 != t01);

    SEQAN_ASSERT(t01 > t00);
    SEQAN_ASSERT_NOT(t01 > t10);

    SEQAN_ASSERT(t01 >= t01);
    SEQAN_ASSERT(t01 >= t00);
    SEQAN_ASSERT_NOT(t01 >= t10);

    SEQAN_ASSERT(t00 <= t01);
    SEQAN_ASSERT(t00 <= t00);
    SEQAN_ASSERT_NOT(t10 <= t01);
}

// SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_comparison_different_spec)  // Not worth the bother, probably.

SEQAN_DEFINE_TEST(test_basic_aggregates_tuple_bit_packed_stream_output)
{
    Tuple<char, 2, BitPacked<> > t1;
    clear(t1);
    assignValue(t1, 0, 'a');
    assignValue(t1, 1, 'b');

    std::stringstream s;
    s << t1;
    SEQAN_ASSERT_EQ(s.str(), "[a b]");

}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_AGGREGATES_H_
