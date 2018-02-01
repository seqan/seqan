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
// Tests for the SeqAn value construction / destruction.
// ==========================================================================

#ifndef TEST_BASIC_TEST_BASIC_CONSTRUCT_DESTRUCT_H_
#define TEST_BASIC_TEST_BASIC_CONSTRUCT_DESTRUCT_H_

#include <seqan/basic.h>

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

SEQAN_DEFINE_TEST(test_basic_construct_destruct_metafunction_is_simple)
{
    bool b = IsSimple<bool>::VALUE;
    SEQAN_ASSERT(b);
    b = IsSimple<int>::VALUE;
    SEQAN_ASSERT(b);
    b = IsSimple<CDStruct>::VALUE;
    SEQAN_ASSERT_NOT(b);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_construct_value_pointer)
{
    // Default construction.
    {
        char space[sizeof(CDStruct)];
        resetCDStructStatics();

        CDStruct * ptr = reinterpret_cast<CDStruct * >(&space[0]);
        valueConstruct(ptr);

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Copy construction.
    {
        char space[sizeof(CDStruct)];
        CDStruct other;
        resetCDStructStatics();

        CDStruct * ptr = reinterpret_cast<CDStruct * >(&space[0]);
        valueConstruct(ptr, other);

        SEQAN_ASSERT_EQ(ptr->copiedFrom, other.id);
        SEQAN_ASSERT_EQ(ptr->movedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &other);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Move construction.
    {
        char space[sizeof(CDStruct)];
        CDStruct other;
        resetCDStructStatics();

        CDStruct * ptr = reinterpret_cast<CDStruct * >(&space[0]);
        valueConstruct(ptr, other, seqan::Move());

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, other.id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &other);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_destruct_value_pointer)
{
    // Destruct default-constructed object.
    {
        resetCDStructStatics();
        CDStruct obj;

        valueDestruct(&obj);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Destruct object created by valueConstruct()
    {
        char space[sizeof(CDStruct)];
        resetCDStructStatics();

        valueConstruct(reinterpret_cast<CDStruct * >(&space[0]));
        valueDestruct(reinterpret_cast<CDStruct * >(&space[0]));

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_construct_pointer)
{
    // arrayConstruct() calling default constructors
    {
        char space[2 * sizeof(CDStruct)];
        resetCDStructStatics();

        CDStruct * itBeg = reinterpret_cast<CDStruct *>(&space[0]);
        CDStruct * itEnd = itBeg + 2;
        arrayConstruct(itBeg, itEnd);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 2);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // arrayConstruct() calling copy constructors
    {
        CDStruct prototype;
        char space[2 * sizeof(CDStruct)];
        resetCDStructStatics();

        CDStruct * itBeg = reinterpret_cast<CDStruct *>(&space[0]);
        CDStruct * itEnd = itBeg + 2;
        arrayConstruct(itBeg, itEnd, prototype);

        SEQAN_ASSERT_EQ(itBeg[0].copiedFrom, prototype.id);
        SEQAN_ASSERT_EQ(itBeg[0].movedFrom, -1);
        SEQAN_ASSERT_EQ(itBeg[0].assignedFrom, -1);
        SEQAN_ASSERT_EQ(itBeg[1].copiedFrom, prototype.id);
        SEQAN_ASSERT_EQ(itBeg[1].movedFrom, -1);
        SEQAN_ASSERT_EQ(itBeg[1].assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &prototype);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 2);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_construct_copy_pointer)
{
    // Construct two objects.
    char space[2 * sizeof(CDStruct)];
    CDStruct * itBeg = reinterpret_cast<CDStruct *>(&space[0]);
    CDStruct * itEnd = itBeg + 2;
    arrayConstruct(itBeg, itEnd);
    resetCDStructStatics();

    // Use arrayConstructCopy() to copy them.
    char space2[2 * sizeof(CDStruct)];
    CDStruct * targetBeg = reinterpret_cast<CDStruct *>(&space2[0]);
    arrayConstructCopy(itBeg, itEnd, targetBeg);

    SEQAN_ASSERT_EQ(targetBeg[0].copiedFrom, itBeg[0].id);
    SEQAN_ASSERT_EQ(targetBeg[0].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetBeg[0].assignedFrom, -1);
    SEQAN_ASSERT_EQ(targetBeg[1].copiedFrom, itBeg[1].id);
    SEQAN_ASSERT_EQ(targetBeg[1].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetBeg[1].assignedFrom, -1);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, &(itBeg[1]));
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 2);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 0);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_construct_move_pointer)
{
    // Construct two objects.
    char space[2 * sizeof(CDStruct)];
    CDStruct * itBeg = reinterpret_cast<CDStruct *>(&space[0]);
    CDStruct * itEnd = itBeg + 2;
    arrayConstruct(itBeg, itEnd);
    resetCDStructStatics();

    // Use arrayConstructCopy() to move them.
    char space2[2 * sizeof(CDStruct)];
    CDStruct * targetBeg = reinterpret_cast<CDStruct *>(&space2[0]);
    arrayConstructMove(itBeg, itEnd, targetBeg);

    // TODO(holtgrew): Move-construction disabled, copying now. Should be changed back later.
    // SEQAN_ASSERT_EQ(targetBeg[0].copiedFrom, -1);
    // SEQAN_ASSERT_EQ(targetBeg[0].movedFrom, itBeg[0].id);
    SEQAN_ASSERT_EQ(targetBeg[0].copiedFrom, itBeg[0].id);
    SEQAN_ASSERT_EQ(targetBeg[0].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetBeg[0].assignedFrom, -1);
    // SEQAN_ASSERT_EQ(targetBeg[1].copiedFrom, -1);
    // SEQAN_ASSERT_EQ(targetBeg[1].movedFrom, itBeg[1].id);
    SEQAN_ASSERT_EQ(targetBeg[1].copiedFrom, itBeg[1].id);
    SEQAN_ASSERT_EQ(targetBeg[1].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetBeg[1].assignedFrom, -1);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, &(itBeg[1]));
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    // SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    // SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 2);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 2);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 0);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_destruct_pointer)
{
    // arrayDestruct() calling default constructors
    char space[2 * sizeof(CDStruct)];
    resetCDStructStatics();

    CDStruct * itBeg = reinterpret_cast<CDStruct *>(&space[0]);
    CDStruct * itEnd = itBeg + 2;
    arrayConstruct(itBeg, itEnd);
    arrayDestruct(itBeg, itEnd);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 2);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 0);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 2);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_fill_pointer)
{
    CDStruct prototype;
    char space[2 * sizeof(CDStruct)];
    resetCDStructStatics();
    // TODO(holtgrew): Need to initialize arrays before assigning with fill_n() through arrayFill()?  Also see http://stackoverflow.com/questions/5490756.

    CDStruct * itBeg = reinterpret_cast<CDStruct *>(&space[0]);
    CDStruct * itEnd = itBeg + 2;
    arrayFill(itBeg, itEnd, prototype);

    SEQAN_ASSERT_EQ(itBeg[0].copiedFrom, -1);
    SEQAN_ASSERT_EQ(itBeg[0].movedFrom, -1);
    SEQAN_ASSERT_EQ(itBeg[0].assignedFrom, prototype.id);
    SEQAN_ASSERT_EQ(itBeg[1].copiedFrom, -1);
    SEQAN_ASSERT_EQ(itBeg[1].movedFrom, -1);
    SEQAN_ASSERT_EQ(itBeg[1].assignedFrom, prototype.id);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, &prototype);
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 0);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 2);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_copy_forward_pointer)
{
    // TODO(holtgrew): A more refined version of this test would copy forward inside the same array.
    char source[2 * sizeof(CDStruct)];
    CDStruct * sourcePtr = reinterpret_cast<CDStruct *>(&source[0]);
    arrayConstruct(sourcePtr, sourcePtr + 2);

    char target[2 * sizeof(CDStruct)];
    CDStruct * targetPtr = reinterpret_cast<CDStruct *>(&target[0]);
    arrayConstruct(targetPtr, targetPtr + 2);
    resetCDStructStatics();

    arrayCopyForward(sourcePtr, sourcePtr + 2, targetPtr);

    SEQAN_ASSERT_EQ(targetPtr[0].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].assignedFrom, sourcePtr[0].id);
    SEQAN_ASSERT_EQ(targetPtr[1].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].assignedFrom, sourcePtr[1].id);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[1]);
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 0);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 2);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_copy_backward_pointer)
{
    // TODO(holtgrew): A more refined version of this test would copy backward inside the same array.
    char source[2 * sizeof(CDStruct)];
    CDStruct * sourcePtr = reinterpret_cast<CDStruct *>(&source[0]);
    arrayConstruct(sourcePtr, sourcePtr + 2);

    char target[2 * sizeof(CDStruct)];
    CDStruct * targetPtr = reinterpret_cast<CDStruct *>(&target[0]);
    arrayConstruct(targetPtr, targetPtr + 2);
    resetCDStructStatics();

    arrayCopyBackward(sourcePtr, sourcePtr + 2, targetPtr);

    SEQAN_ASSERT_EQ(targetPtr[0].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].assignedFrom, sourcePtr[0].id);
    SEQAN_ASSERT_EQ(targetPtr[1].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].assignedFrom, sourcePtr[1].id);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[0]);
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 0);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 2);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_copy_pointer)
{
    // TODO(holtgrew): A more refined version of this test would copy (1) forward and (2) backward inside the same array.
    char source[2 * sizeof(CDStruct)];
    CDStruct * sourcePtr = reinterpret_cast<CDStruct *>(&source[0]);
    arrayConstruct(sourcePtr, sourcePtr + 2);

    char target[2 * sizeof(CDStruct)];
    CDStruct * targetPtr = reinterpret_cast<CDStruct *>(&target[0]);
    arrayConstruct(targetPtr, targetPtr + 2);
    resetCDStructStatics();

    arrayCopy(sourcePtr, sourcePtr + 2, targetPtr);

    SEQAN_ASSERT_EQ(targetPtr[0].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].assignedFrom, sourcePtr[0].id);
    SEQAN_ASSERT_EQ(targetPtr[1].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].movedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].assignedFrom, sourcePtr[1].id);

    // Stack order is different between GCC versions.
    if (sourcePtr >= targetPtr)
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[1]);
    else
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[0]);
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 0);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 2);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_move_forward_pointer)
{
    // TODO(holtgrew): A more refined version of this test would move forward inside the same array.
    char source[2 * sizeof(CDStruct)];
    CDStruct * sourcePtr = reinterpret_cast<CDStruct *>(&source[0]);
    arrayConstruct(sourcePtr, sourcePtr + 2);

    char target[2 * sizeof(CDStruct)];
    CDStruct * targetPtr = reinterpret_cast<CDStruct *>(&target[0]);
    arrayConstruct(targetPtr, targetPtr + 2);
    resetCDStructStatics();

    arrayMoveForward(sourcePtr, sourcePtr + 2, targetPtr);

    SEQAN_ASSERT_EQ(targetPtr[0].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].movedFrom, sourcePtr[0].id);
    SEQAN_ASSERT_EQ(targetPtr[0].assignedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].movedFrom, sourcePtr[1].id);
    SEQAN_ASSERT_EQ(targetPtr[1].assignedFrom, -1);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[1]);
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 2);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_move_backward_pointer)
{
    // TODO(holtgrew): A more refined version of this test would move backward inside the same array.
    char source[2 * sizeof(CDStruct)];
    CDStruct * sourcePtr = reinterpret_cast<CDStruct *>(&source[0]);
    arrayConstruct(sourcePtr, sourcePtr + 2);

    char target[2 * sizeof(CDStruct)];
    CDStruct * targetPtr = reinterpret_cast<CDStruct *>(&target[0]);
    arrayConstruct(targetPtr, targetPtr + 2);
    resetCDStructStatics();

    arrayMoveBackward(sourcePtr, sourcePtr + 2, targetPtr);

    SEQAN_ASSERT_EQ(targetPtr[0].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].movedFrom, sourcePtr[0].id);
    SEQAN_ASSERT_EQ(targetPtr[0].assignedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].movedFrom, sourcePtr[1].id);
    SEQAN_ASSERT_EQ(targetPtr[1].assignedFrom, -1);

    SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[0]);
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 2);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_move_pointer)
{
    // TODO(holtgrew): A more refined version of this test would move (1) forward and (2) backward inside the same array.
    char source[2 * sizeof(CDStruct)];
    CDStruct * sourcePtr = reinterpret_cast<CDStruct *>(&source[0]);
    arrayConstruct(sourcePtr, sourcePtr + 2);

    char target[2 * sizeof(CDStruct)];
    CDStruct * targetPtr = reinterpret_cast<CDStruct *>(&target[0]);
    arrayConstruct(targetPtr, targetPtr + 2);
    resetCDStructStatics();

    arrayMove(sourcePtr, sourcePtr + 2, targetPtr);

    SEQAN_ASSERT_EQ(targetPtr[0].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[0].movedFrom, sourcePtr[0].id);
    SEQAN_ASSERT_EQ(targetPtr[0].assignedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].copiedFrom, -1);
    SEQAN_ASSERT_EQ(targetPtr[1].movedFrom, sourcePtr[1].id);
    SEQAN_ASSERT_EQ(targetPtr[1].assignedFrom, -1);

    // Stack order is different between GCC versions.
    if (sourcePtr >= targetPtr)
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[1]);
    else
        SEQAN_ASSERT_EQ(CDStruct::lastOther, &sourcePtr[0]);
    SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::moves, 2);
    SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
    SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
}

SEQAN_DEFINE_TEST(test_basic_construct_destruct_array_clear_space_pointer)
{
    // Clear only.
    {
        char arr[6 * sizeof(CDStruct)];
        CDStruct * arrPtr = reinterpret_cast<CDStruct *>(&arr[0]);
        arrayConstruct(arrPtr, arrPtr + 4);
        resetCDStructStatics();

        arrayClearSpace(arrPtr, 4, 4, 0);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 4);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Do not move, clear first entry.
    {
        char arr[6 * sizeof(CDStruct)];
        CDStruct * arrPtr = reinterpret_cast<CDStruct *>(&arr[0]);
        arrayConstruct(arrPtr, arrPtr + 4);
        resetCDStructStatics();

        arrayClearSpace(arrPtr, 4, 1, 1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Move left, overlapping.
    {
        char arr[6 * sizeof(CDStruct)];
        CDStruct * arrPtr = reinterpret_cast<CDStruct *>(&arr[0]);
        arrayConstruct(arrPtr, arrPtr + 4);
        resetCDStructStatics();

        int oldId1 = arrPtr[1].id;
        int oldId2 = arrPtr[2].id;
        int oldId3 = arrPtr[3].id;
        arrayClearSpace(arrPtr, 4, 1, 0);
        SEQAN_ASSERT_EQ(arrPtr[0].movedFrom, oldId1);
        SEQAN_ASSERT_EQ(arrPtr[1].movedFrom, oldId2);
        SEQAN_ASSERT_EQ(arrPtr[2].movedFrom, oldId3);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(&arrPtr[3]));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 3);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Move left, non-overlapping.
    {
        char arr[6 * sizeof(CDStruct)];
        CDStruct * arrPtr = reinterpret_cast<CDStruct *>(&arr[0]);
        arrayConstruct(arrPtr, arrPtr + 4);
        resetCDStructStatics();

        int oldId2 = arrPtr[2].id;
        int oldId3 = arrPtr[3].id;
        arrayClearSpace(arrPtr, 4, 2, 0);
        SEQAN_ASSERT_EQ(arrPtr[0].movedFrom, oldId2);
        SEQAN_ASSERT_EQ(arrPtr[1].movedFrom, oldId3);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(&arrPtr[3]));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 2);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 2);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Move right, overlapping.
    {
        char arr[6 * sizeof(CDStruct)];
        CDStruct * arrPtr = reinterpret_cast<CDStruct *>(&arr[0]);
        arrayConstruct(arrPtr, arrPtr + 4);
        resetCDStructStatics();

        int oldId2 = arrPtr[2].id;
        int oldId3 = arrPtr[3].id;
        arrayClearSpace(arrPtr, 4, 2, 3);
        SEQAN_ASSERT_EQ(arrPtr[3].movedFrom, oldId2);
        // SEQAN_ASSERT_EQ(arrPtr[4].movedFrom, oldId3);
        SEQAN_ASSERT_EQ(arrPtr[4].movedFrom, -1);
        SEQAN_ASSERT_EQ(arrPtr[4].copiedFrom, oldId3);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(&arrPtr[2]));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        // SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        // SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 1);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 3);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }

    // Move right, non-overlapping.
    {
        char arr[6 * sizeof(CDStruct)];
        CDStruct * arrPtr = reinterpret_cast<CDStruct *>(&arr[0]);
        arrayConstruct(arrPtr, arrPtr + 4);
        resetCDStructStatics();

        int oldId2 = arrPtr[2].id;
        int oldId3 = arrPtr[3].id;
        arrayClearSpace(arrPtr, 4, 2, 4);
        // SEQAN_ASSERT_EQ(arrPtr[4].movedFrom, oldId2);
        // SEQAN_ASSERT_EQ(arrPtr[5].movedFrom, oldId3);
        SEQAN_ASSERT_EQ(arrPtr[4].copiedFrom, oldId2);
        SEQAN_ASSERT_EQ(arrPtr[5].copiedFrom, oldId3);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(&arrPtr[3]));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        // SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        // SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 2);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 2);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 4);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
    }
}

#endif  // TEST_BASIC_TEST_BASIC_CONSTRUCT_DESTRUCT_H_
