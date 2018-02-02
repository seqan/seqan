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

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_SMART_POINTER_HOLDER_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_SMART_POINTER_HOLDER_H_

// ==========================================================================
// Helper Code
// ==========================================================================

struct CDStruct
{
    static int nextId;
    static int defaultConstructions;
    static int copyConstructions;
    static int moveConstructions;
    static int moves;
    static int destructions;
    static int assignments;
    static int sets;
    static CDStruct const * lastOther;

    // Note that we use construct this type in char arrays below.  This is no
    // problem since it only contains ints and there should be no alignment
    // problems.

    int id;
    int copiedFrom;
    int movedFrom;
    int assignedFrom;
    int setFrom;

    CDStruct() : copiedFrom(-1), movedFrom(-1), assignedFrom(-1), setFrom(-1)
    {
        id = nextId++;
        defaultConstructions += 1;
    }

    CDStruct(CDStruct const & other)
            : copiedFrom(other.id), movedFrom(-1), assignedFrom(-1), setFrom(-1)
    {
        id = nextId++;
        lastOther = &other;
        copyConstructions += 1;
    }

    CDStruct(CDStruct & other, seqan::Move const & /*tag*/)
            : copiedFrom(-1), movedFrom(other.id), assignedFrom(-1), setFrom(-1)
    {
        lastOther = &other;
        moveConstructions += 1;
    }

    CDStruct & operator=(CDStruct const & other)
    {
        lastOther = &other;
        assignments += 1;
        copiedFrom = -1;
        movedFrom = -1;
        assignedFrom = other.id;
        setFrom = -1;
        return *this;
    }

    ~CDStruct()
    {
        destructions += 1;
    }
};

void move(CDStruct & target, CDStruct const & source)
{
    CDStruct::lastOther = &source;
    CDStruct::moves += 1;
    target.copiedFrom = -1;
    target.assignedFrom = -1;
    target.movedFrom = source.id;
    target.setFrom = -1;
}


void move(CDStruct & target, CDStruct & source)
{
    move(target, const_cast<CDStruct const &>(source));
}

void set(CDStruct & target, CDStruct const & source)
{
    CDStruct::lastOther = &source;
    CDStruct::sets += 1;
    target.copiedFrom = -1;
    target.assignedFrom = -1;
    target.movedFrom = -1;
    target.setFrom = source.id;
}


void set(CDStruct & target, CDStruct & source)
{
    set(target, const_cast<CDStruct const &>(source));
}

int CDStruct::nextId = 0;
int CDStruct::defaultConstructions = 0;
int CDStruct::copyConstructions = 0;
int CDStruct::moveConstructions = 0;
int CDStruct::moves = 0;
int CDStruct::destructions = 0;
int CDStruct::assignments = 0;
int CDStruct::sets = 0;
CDStruct const * CDStruct::lastOther = 0;

void resetCDStructStatics()
{
    CDStruct::lastOther = 0;
    CDStruct::defaultConstructions = 0;
    CDStruct::copyConstructions = 0;
    CDStruct::moveConstructions = 0;
    CDStruct::moves = 0;
    CDStruct::destructions = 0;
    CDStruct::assignments = 0;
    CDStruct::sets = 0;
}

// ==========================================================================
// Tests
// ==========================================================================

// --------------------------------------------------------------------------
// Test for Simple Holder
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_simple_metafunctions)
{
    using namespace seqan;

    // Tests for non-const holder.
    {
        typedef Holder<int, Simple> THolder;

        typedef typename Value<THolder>::Type TValue;
        bool b = IsSameType<int, TValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename GetValue<THolder>::Type TGetValue;
        b = IsSameType<int const &, TGetValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Reference<THolder>::Type TReference;
        b = IsSameType<int &, TReference>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Spec<THolder>::Type TSpec;
        b = IsSameType<Simple, TSpec>::Type::VALUE;
        SEQAN_ASSERT(b);
    }

    // Tests for const holder.
    {
        typedef Holder<int, Simple> const THolder;

        typedef typename Value<THolder>::Type TValue;
        bool b = IsSameType<int, TValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename GetValue<THolder>::Type TGetValue;
        b = IsSameType<int const &, TGetValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Reference<THolder>::Type TReference;
        b = IsSameType<int &, TReference>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Spec<THolder>::Type TSpec;
        b = IsSameType<Simple, TSpec>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_simple_constructors)
{
    using namespace seqan;

    // Simple Holder always copies, wo we do not really have to test much
    // here.

    // Construct from value.
    {
        typedef Holder<int, Simple> THolder;

        int x = 10;
        THolder holder(x);

        SEQAN_ASSERT_EQ(value(holder), 10);
        SEQAN_ASSERT_NEQ(&value(holder), &x);
        SEQAN_ASSERT_NOT(empty(holder));
        SEQAN_ASSERT_NOT(dependent(holder));
    }

    // Construct from other holder.
    {
        typedef Holder<int, Simple> THolderSource;
        typedef Holder<int, Simple> THolderTarget;

        int x = 10;
        THolderSource holderSource(x);
        THolderTarget holderTarget(holderSource);

        SEQAN_ASSERT_EQ(value(holderTarget), 10);
        SEQAN_ASSERT_NEQ(&value(holderTarget), &value(holderSource));
        SEQAN_ASSERT_NOT(empty(holderTarget));
        SEQAN_ASSERT_NOT(dependent(holderTarget));
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_simple_transport)
{
    using namespace seqan;

    // assign()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder s(s1);
        THolder t(s2);
        assign(t, s);

        SEQAN_ASSERT_EQ(value(t).assignedFrom, value(s).id);
    }

    // move()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder s(s1);
        THolder t(s2);
        move(t, s);

        SEQAN_ASSERT_EQ(value(t).assignedFrom, value(s).id);
    }

    // set()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder s(s1);
        THolder t(s2);
        set(t, s);

        SEQAN_ASSERT_EQ(value(t).assignedFrom, value(s).id);
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_simple_transport_value)
{
    using namespace seqan;

    // assignValue()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder h(s1);
        assignValue(h, s2);

        SEQAN_ASSERT_EQ(value(h).assignedFrom, s2.id);
    }

    // moveValue()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder h(s1);
        moveValue(h, s2);

        SEQAN_ASSERT_EQ(value(h).movedFrom, s2.id);
    }

    // setValue()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder h(s1);
        setValue(h, s2);

        SEQAN_ASSERT_EQ(value(h).setFrom, s2.id);
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_simple_dependencies)
{
    using namespace seqan;

    // Since Simple Holder always copies, some simple tests are enough.

    // Construct from value;  Clearing should have no effect.
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct x;
        THolder holder(x);

        SEQAN_ASSERT_NOT(empty(holder));
        SEQAN_ASSERT_NOT(dependent(holder));
        clear(holder);
        SEQAN_ASSERT_NOT(empty(holder));

        SEQAN_ASSERT_EQ(value(holder).copiedFrom, x.id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &x);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // Test detach - empty implementation.
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct x;
        THolder holder(x);

        SEQAN_ASSERT_NOT(empty(holder));
        SEQAN_ASSERT_NOT(dependent(holder));
        clear(holder);
        SEQAN_ASSERT_NOT(empty(holder));

        detach(holder);

        SEQAN_ASSERT_EQ(value(holder).copiedFrom, x.id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &x);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // Test create - trivial implementation.
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct x;
        THolder holder(x);
        SEQAN_ASSERT_EQ(value(holder).copiedFrom, x.id);
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &x);

        create(holder);  // Does nothing.

        CDStruct y;
        create(holder, y);  // Assigns.
        SEQAN_ASSERT_EQ(value(holder).assignedFrom, y.id);
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &y);

        CDStruct z;
        create(holder, z, Move());  // Moves
        SEQAN_ASSERT_EQ(value(holder).assignedFrom, z.id);
        // SEQAN_ASSERT_EQ(value(holder).movedFrom, z.id);
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &z);

        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 3);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        // SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        // SEQAN_ASSERT_EQ(CDStruct::assignments, 1);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 2);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // Copy construct from same type.  Should copy all the same and not be
    // dependent.
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct x;
        THolder holderSource(x);
        THolder holderTarget(holderSource);

        SEQAN_ASSERT_EQ(value(holderTarget).copiedFrom, value(holderSource).id);
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &value(holderSource));

        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 2);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_simple_accessor)
{
    using namespace seqan;

    typedef Holder<int, Simple> THolder;

    THolder holder(10);

    SEQAN_ASSERT_EQ(value(holder), 10);
    SEQAN_ASSERT_EQ(getValue(holder), 10);

    THolder const constHolder(10);

    SEQAN_ASSERT_EQ(value(constHolder), 10);
    SEQAN_ASSERT_EQ(getValue(constHolder), 10);
}

// --------------------------------------------------------------------------
// Test for Tristate Holder
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_tristate_metafunctions)
{
    using namespace seqan;

    // Tests for non-const holder.
    {
        typedef Holder<int, Tristate> THolder;

        typedef typename Value<THolder>::Type TValue;
        bool b = IsSameType<int, TValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename GetValue<THolder>::Type TGetValue;
        b = IsSameType<int const &, TGetValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Reference<THolder>::Type TReference;
        b = IsSameType<int &, TReference>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Spec<THolder>::Type TSpec;
        b = IsSameType<Tristate, TSpec>::Type::VALUE;
        SEQAN_ASSERT(b);
    }

    // Tests for const holder.
    {
        typedef Holder<int, Tristate> const THolder;

        typedef typename Value<THolder>::Type TValue;
        bool b = IsSameType<int, TValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename GetValue<THolder>::Type TGetValue;
        b = IsSameType<int const &, TGetValue>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Reference<THolder>::Type TReference;
        b = IsSameType<int &, TReference>::Type::VALUE;
        SEQAN_ASSERT(b);

        typedef typename Spec<THolder>::Type TSpec;
        b = IsSameType<Tristate, TSpec>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_tristate_constructors)
{
    using namespace seqan;

    typedef Holder<CDStruct, Tristate> THolder;
    typedef Holder<CDStruct, Tristate> TConstHolder;

    // Default constructor
    {
        THolder holder;
        TConstHolder constHolder;
    }

    // Constructor with reference
    {
        resetCDStructStatics();
        CDStruct s;

        THolder holder(s);
        TConstHolder constHolder(s);

        SEQAN_ASSERT_NOT(empty(holder));
        SEQAN_ASSERT(dependent(holder));
        SEQAN_ASSERT_EQ(&value(holder), &s);

        SEQAN_ASSERT_NOT(empty(constHolder));
        SEQAN_ASSERT(dependent(constHolder));
        SEQAN_ASSERT_EQ(&value(constHolder), &s);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    {
        resetCDStructStatics();
        CDStruct const s;

        THolder holder(s);
        TConstHolder constHolder(s);

        SEQAN_ASSERT_NOT(empty(holder));
        SEQAN_ASSERT_NOT(dependent(holder));
        SEQAN_ASSERT_NEQ(&value(holder), &s);

        // TODO(holtgrew): Change tests once Holder is fixed for const values.

        SEQAN_ASSERT_NOT(empty(holder));
        // SEQAN_ASSERT(dependent(constHolder));
        SEQAN_ASSERT_NOT(dependent(constHolder));
        // SEQAN_ASSERT_EQ(&value(constHolder), &s);
        SEQAN_ASSERT_NEQ(&value(constHolder), &s);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &s);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        // SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 2);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // TODO(holtgrew): Complete the following once holders of const values work!

    // Copy constructor const to const
    // Copy constructor const to non-const
    // Copy constructor non-const to const
    // Copy constructor non-const to non-const
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_tristate_transport)
{
    using namespace seqan;

    // TODO(holtgrew): Complete the following once holders of const values work!

    // Also check dependencies on the following.
    // FAST
    //   const <-- set/move <-- non-const
    //   const <-- set/move <-- const
    //   non-const <-- set/move <-- non-const
    // COPY
    //   non-const <-- set/move <-- const
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_tristate_transport_value)
{
    using namespace seqan;

    // assignValue()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder h(s1);
        assignValue(h, s2);

        SEQAN_ASSERT_EQ(value(h).assignedFrom, s2.id);
    }

    // moveValue()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder h(s1);
        moveValue(h, s2);

        SEQAN_ASSERT_EQ(value(h).movedFrom, s2.id);
    }

    // setValue()
    {
        typedef Holder<CDStruct, Simple> THolder;

        resetCDStructStatics();

        CDStruct s1;
        CDStruct s2;
        THolder h(s1);
        setValue(h, s2);

        SEQAN_ASSERT_EQ(value(h).setFrom, s2.id);
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_tristate_dependencies)
{
    using namespace seqan;

    typedef Holder<CDStruct, Tristate> THolder;

    // Default construction -> assignment.
    {
        resetCDStructStatics();

        THolder holder;
        CDStruct s;

        SEQAN_ASSERT(empty(holder));
        assignValue(holder, s);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &s);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // Default construction -> create, default construction.
    {
        resetCDStructStatics();

        THolder holder;

        SEQAN_ASSERT(empty(holder));
        create(holder);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    // Default construction -> create, copy construction.
    {
        resetCDStructStatics();

        THolder holder;
        CDStruct s;

        SEQAN_ASSERT(empty(holder));
        create(holder, s);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &s);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // Default construction -> create -> create.
    {
        resetCDStructStatics();

        THolder holder;

        SEQAN_ASSERT(empty(holder));
        create(holder);
        create(holder);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // Construction with value -> detach.
    {
        resetCDStructStatics();

        CDStruct s;
        THolder holder(s);

        SEQAN_ASSERT(dependent(holder));

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);

        detach(holder);

        SEQAN_ASSERT_NEQ(&value(holder), &s);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &s);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }

    // Copy construction.
    {
        resetCDStructStatics();

        CDStruct s;
        THolder holder(s);
        THolder holder2(s);

        SEQAN_ASSERT(dependent(holder));
        SEQAN_ASSERT(dependent(holder2));
        SEQAN_ASSERT_EQ(&value(holder), &s);
        SEQAN_ASSERT_EQ(&value(holder2), &s);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_smart_pointer_holder_tristate_accessor)
{
    using namespace seqan;

    typedef Holder<int, Tristate> THolder;

    THolder holder(10);

    SEQAN_ASSERT_EQ(value(holder), 10);
    SEQAN_ASSERT_EQ(getValue(holder), 10);

    THolder const constHolder(10);

    SEQAN_ASSERT_EQ(value(constHolder), 10);
    SEQAN_ASSERT_EQ(getValue(constHolder), 10);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_SMART_POINTER_HOLDER_H_
