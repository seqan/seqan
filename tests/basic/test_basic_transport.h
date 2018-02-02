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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_TRANSPORT_H_
#define TESTS_BASIC_TEST_BASIC_TRANSPORT_H_

// ==========================================================================
// Helper Code
// ==========================================================================

struct MoveYes {};
struct MoveNo {};

template <typename TSpec>
struct TransportObj_
{
    static int assignments;
    static int sets;
    static int moves;
    static int nextId;

    int id;
    int assignedFrom;
    int setFrom;
    int movedFrom;

    TransportObj_()
            : id(nextId++), assignedFrom(-1), setFrom(-1), movedFrom(-1)
    {}

    TransportObj_(TransportObj_ & other)
            : id(nextId++), assignedFrom(other.id), setFrom(-1), movedFrom(-1)
    {
    }

    TransportObj_(TransportObj_ & other, seqan::Move)
            : id(nextId++), assignedFrom(-1), setFrom(-1), movedFrom(other.id)
    {
        bool b = IsSameType<TSpec, MoveYes>::Type::VALUE;
        SEQAN_ASSERT(b);
    }

    TransportObj_(TransportObj_ && other)
            : id(nextId++), assignedFrom(-1), setFrom(-1), movedFrom(other.id)
    {
        bool b = IsSameType<TSpec, MoveYes>::Type::VALUE;
        SEQAN_ASSERT(b);
    }

    TransportObj_ &
    operator=(TransportObj_ const & other)
    {
        assignments += 1;
        assignedFrom = other.id;
        return *this;
    }
};

void set(TransportObj_<MoveYes> & target, TransportObj_<MoveYes> const & other)
{
    target.setFrom = other.id;
    TransportObj_<MoveYes>::sets += 1;
}

void set(TransportObj_<MoveYes> & target, TransportObj_<MoveYes> & other)
{
    set(target, static_cast<TransportObj_<MoveYes> const &>(other));
}

void move(TransportObj_<MoveYes> & target, TransportObj_<MoveYes> const & other)
{
    target.movedFrom = other.id;
    TransportObj_<MoveYes>::moves += 1;
}

void move(TransportObj_<MoveYes> & target, TransportObj_<MoveYes> & other)
{
    move(target, static_cast<TransportObj_<MoveYes> const &>(other));
}

template <typename TSpec>
int TransportObj_<TSpec>::assignments = 0;
template <typename TSpec>
int TransportObj_<TSpec>::sets = 0;
template <typename TSpec>
int TransportObj_<TSpec>::moves = 0;
template <typename TSpec>
int TransportObj_<TSpec>::nextId = 0;

template <typename TSpec>
void resetTransportObjStatics(TSpec const &)
{
    TransportObj_<TSpec>::assignments = 0;
    TransportObj_<TSpec>::sets = 0;
    TransportObj_<TSpec>::moves = 0;
}

namespace seqan {

template <>
struct HasMoveConstructor<TransportObj_<MoveYes> >
{
    typedef True Type;
};

template <>
struct HasMoveConstructor<TransportObj_<MoveYes> const>
{
    typedef True Type;
};

}  // namespace seqan

// ==========================================================================
// Actual Tests
// ==========================================================================

SEQAN_DEFINE_TEST(test_basic_transport_has_move_constructor)
{
    // Test HasMoveConstructor metafunction, default and specialization from above.
    bool b = HasMoveConstructor<int>::Type::VALUE;
    SEQAN_ASSERT_NOT(b);
    b = HasMoveConstructor<TransportObj_<MoveNo> >::Type::VALUE;
    SEQAN_ASSERT_NOT(b);

    b = HasMoveConstructor<TransportObj_<MoveYes> >::Type::VALUE;
    SEQAN_ASSERT(b);
    b = HasMoveConstructor<TransportObj_<MoveYes> const>::Type::VALUE;
    SEQAN_ASSERT(b);
}

// Test default overloads of set/assign/set
SEQAN_DEFINE_TEST(test_basic_transport_default_overloads)
{
    // assign()
    {
        resetTransportObjStatics(MoveNo());
        TransportObj_<MoveNo> obj1, obj2;
        assign(obj1, obj2);
        SEQAN_ASSERT_EQ(obj1.assignedFrom, obj2.id);
        SEQAN_ASSERT_EQ(obj1.setFrom, -1);
        SEQAN_ASSERT_EQ(obj1.movedFrom, -1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::assignments, 1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::sets, 0);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::moves, 0);
    }
    // set()
    {
        resetTransportObjStatics(MoveNo());
        TransportObj_<MoveNo> obj1, obj2;
        set(obj1, obj2);
        SEQAN_ASSERT_EQ(obj1.assignedFrom, obj2.id);
        SEQAN_ASSERT_EQ(obj1.setFrom, -1);
        SEQAN_ASSERT_EQ(obj1.movedFrom, -1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::assignments, 1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::sets, 0);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::moves, 0);
    }
    // move()
    {
        resetTransportObjStatics(MoveNo());
        TransportObj_<MoveNo> obj1, obj2;
        move(obj1, obj2);
        SEQAN_ASSERT_EQ(obj1.assignedFrom, obj2.id);
        SEQAN_ASSERT_EQ(obj1.setFrom, -1);
        SEQAN_ASSERT_EQ(obj1.movedFrom, -1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::assignments, 1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::sets, 0);
        SEQAN_ASSERT_EQ(TransportObj_<MoveNo>::moves, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_transport_assign_move_set)
{
    // assign()
    {
        resetTransportObjStatics(MoveYes());
        TransportObj_<MoveYes> obj1, obj2;
        assign(obj1, obj2);
        SEQAN_ASSERT_EQ(obj1.assignedFrom, obj2.id);
        SEQAN_ASSERT_EQ(obj1.setFrom, -1);
        SEQAN_ASSERT_EQ(obj1.movedFrom, -1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::assignments, 1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::sets, 0);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::moves, 0);
    }
    // set()
    {
        resetTransportObjStatics(MoveYes());
        TransportObj_<MoveYes> obj1, obj2;
        set(obj1, obj2);
        SEQAN_ASSERT_EQ(obj1.assignedFrom, -1);
        SEQAN_ASSERT_EQ(obj1.setFrom, obj2.id);
        SEQAN_ASSERT_EQ(obj1.movedFrom, -1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::assignments, 0);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::sets, 1);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::moves, 0);
    }
    // move()
    {
        resetTransportObjStatics(MoveYes());
        TransportObj_<MoveYes> obj1, obj2;
        move(obj1, obj2);
        SEQAN_ASSERT_EQ(obj1.assignedFrom, -1);
        SEQAN_ASSERT_EQ(obj1.setFrom, -1);
        SEQAN_ASSERT_EQ(obj1.movedFrom, obj2.id);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::assignments, 0);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::sets, 0);
        SEQAN_ASSERT_EQ(TransportObj_<MoveYes>::moves, 1);
    }
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_TRANSPORT_H_
