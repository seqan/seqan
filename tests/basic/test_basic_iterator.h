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

#ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_H_
#define TESTS_BASIC_TEST_BASIC_ITERATOR_H_

#include <vector>

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

// Disable for windows vs as the two-phase template lookup is broken. See the
// following link http://stackoverflow.com/questions/6273176/what-exactly-is-broken-with-microsoft-visual-cs-two-phase-template-instanti
#if !defined(STDLIB_VS)
// --------------------------------------------------------------------------
// Pointer adaptions to test the positional iterator.
// --------------------------------------------------------------------------

template <typename TValue, typename TPos>
inline TValue &
value(TValue * me,
      TPos pos)
{
    return me[pos];
}

template <typename TValue, typename TPos, typename TValue2>
inline void
assignValue(TValue * me,
            TPos pos,
            TValue2 const & _value)
{
    seqan::assign(value(me, pos), _value);
}

template <typename TValue, typename TValue2, typename TPos>
inline void
moveValue(TValue * me,
          TPos pos,
          TValue2 const & _value)
{
    move(value(me, pos), _value);
}
#endif  // !defined(STDLIB_VS) 
// --------------------------------------------------------------------------
// Tests for Pointer Adaption to Iterator Concept
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_metafunctions)
{
    using namespace seqan;

    // Pointers.
    {
        typedef int * TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Value<TIterator>::Type, int>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Const-Pointers.
    {
        typedef int const * TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        // TODO(holtgrew): This is inconsistent
        b = IsSameType<typename Value<TIterator>::Type, int const>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_transport)
{
    using namespace seqan;

    // assign()
    {
        int x = 1, y = 2;
        int * ptr = &x;
        assign(ptr, &y);
        SEQAN_ASSERT_EQ(ptr, &y);
    }
    // move()
    {
        int x = 1, y = 2;
        int * ptr = &x;
        move(ptr, &y);
        SEQAN_ASSERT_EQ(ptr, &y);
    }
    // set()
    {
        int x = 1, y = 2;
        int * ptr = &x;
        set(ptr, &y);
        SEQAN_ASSERT_EQ(ptr, &y);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_transport_value)
{
    using namespace seqan;

    // assignValue()
    {
        CDStruct cs1, cs2;
        resetCDStructStatics();

        CDStruct * ptr = &cs1;
        *ptr = cs2;

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->setFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, cs2.id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &cs2);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 1);
    }
    // moveValue()
    {
        CDStruct cs1, cs2;
        resetCDStructStatics();

        CDStruct * ptr = &cs1;
        moveValue(ptr, cs2);

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, cs2.id);
        SEQAN_ASSERT_EQ(ptr->setFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &cs2);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    // setValue()
    {
        CDStruct cs1, cs2;
        resetCDStructStatics();

        CDStruct * ptr = &cs1;
        setValue(ptr, cs2);

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->setFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, -1);

        SEQAN_ASSERT_EQ(ptr, &cs2);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_movement)
{
    using namespace seqan;

    // goNext/operator++
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        int * ptr = &arr[4];
        goNext(ptr);
        SEQAN_ASSERT_EQ(ptr, &arr[5]);

        ptr = &arr[4];
        ptr++;
        SEQAN_ASSERT_EQ(ptr, &arr[5]);

        ptr = &arr[4];
        ++ptr;
        SEQAN_ASSERT_EQ(ptr, &arr[5]);
    }
    // goPrevious/operator--
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        int * ptr = &arr[4];
        goPrevious(ptr);
        SEQAN_ASSERT_EQ(ptr, &arr[3]);

        ptr = &arr[4];
        ptr--;
        SEQAN_ASSERT_EQ(ptr, &arr[3]);

        ptr = &arr[4];
        --ptr;
        SEQAN_ASSERT_EQ(ptr, &arr[3]);
    }
    // goFurther/operator+=/operator-=
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        int * ptr = &arr[4];
        goFurther(ptr, 2);
        SEQAN_ASSERT_EQ(ptr, &arr[6]);

        ptr = &arr[4];
        goFurther(ptr, -2);
        SEQAN_ASSERT_EQ(ptr, &arr[2]);

        ptr = &arr[4];
        ptr += 2;
        SEQAN_ASSERT_EQ(ptr, &arr[6]);

        ptr = &arr[4];
        ptr -= 2;
        SEQAN_ASSERT_EQ(ptr, &arr[2]);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_arithmetics)
{
    using namespace seqan;

    int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    int * ptr = &arr[4];

    int * ptr2 = ptr + 2;
    SEQAN_ASSERT_EQ(ptr2, &arr[6]);

    ptr2 = ptr - 2;
    SEQAN_ASSERT_EQ(ptr2, &arr[2]);

    ptr2 = ptr + 2;
    SEQAN_ASSERT_EQ(ptr2 - ptr, 2);
}

// --------------------------------------------------------------------------
// Tests for STL Iterator Adaption to Iterator Concept
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_metafunctions)
{
    using namespace seqan;

    // Iterator.
    {
        typedef Iter<std::vector<int>, StdIteratorAdaptor> TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Value<TIterator>::Type, int>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int &>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, std::vector<int> >::VALUE;
        SEQAN_ASSERT(b);
    }
    // Const-Iterators.
    {
        typedef Iter<std::vector<int> const, StdIteratorAdaptor> TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        // TODO(holtgrew): This is inconsistent
        b = IsSameType<typename Value<TIterator>::Type, int>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, std::vector<int> const>::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_constructors)
{
    using namespace seqan;

    typedef Iter<std::vector<int>, StdIteratorAdaptor> TIterator;

    std::vector<int> vec;
    vec.push_back(0);
    vec.push_back(1);

    // Construct from STL iterator.
    TIterator it(vec.begin());

    // Copy constructor.
    TIterator it2(it);

    SEQAN_ASSERT_EQ(*it, *it2);
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_transport)
{
    using namespace seqan;

    typedef Iter<std::vector<int>, StdIteratorAdaptor> TIterator;

    std::vector<int> vec;
    vec.push_back(0);
    vec.push_back(1);

    TIterator it(vec.begin());

    // assign()
    {
        TIterator it2;
        assign(it2, it);
        SEQAN_ASSERT(it == it2);
    }
    // set()
    {
        TIterator it2;
        seqan::set(it2, it);
        SEQAN_ASSERT(it == it2);
    }
    // move()
    {
        TIterator it2;
        move(it2, it);
        SEQAN_ASSERT(it == it2);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_transport_value)
{
    using namespace seqan;

    typedef Iter<std::vector<CDStruct>, StdIteratorAdaptor> TIterator;

    // assignValue()
    {
        std::vector<CDStruct> vec;
        vec.resize(3);
        resetCDStructStatics();

        TIterator it = vec.begin();
        assignValue(it, vec[1]);

        SEQAN_ASSERT_EQ(it->copiedFrom, -1);
        SEQAN_ASSERT_EQ(it->movedFrom, -1);
        SEQAN_ASSERT_EQ(it->setFrom, -1);
        SEQAN_ASSERT_EQ(it->assignedFrom, vec[1].id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &vec[1]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 1);
    }
    // moveValue()
    {
        std::vector<CDStruct> vec;
        vec.resize(3);
        resetCDStructStatics();

        TIterator it = vec.begin();
        moveValue(it, vec[1]);

        SEQAN_ASSERT_EQ(it->copiedFrom, -1);
        SEQAN_ASSERT_EQ(it->movedFrom, vec[1].id);
        SEQAN_ASSERT_EQ(it->setFrom, -1);
        SEQAN_ASSERT_EQ(it->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &vec[1]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    // There is no setValue() with iterator semantics.
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_movement)
{
    using namespace seqan;

    typedef Iter<std::vector<int>, StdIteratorAdaptor> TIterator;

    // goNext/operator++
    {
        std::vector<int> vec;
        for (int i = 0; i < 10; ++i)
            vec.push_back(i);

        TIterator it(vec.begin() + 4);
        goNext(it);
        SEQAN_ASSERT_EQ(*it, vec[5]);

        it = vec.begin() + 4;
        it++;
        SEQAN_ASSERT_EQ(*it, vec[5]);

        it = vec.begin() + 4;
        ++it;
        SEQAN_ASSERT_EQ(*it, vec[5]);
    }
    // goPrevious/operator--
    {
        std::vector<int> vec;
        for (int i = 0; i < 10; ++i)
            vec.push_back(i);

        TIterator it(vec.begin() + 4);
        goPrevious(it);
        SEQAN_ASSERT_EQ(*it, vec[3]);

        it = vec.begin() + 4;
        it--;
        SEQAN_ASSERT_EQ(*it, vec[3]);

        it = vec.begin() + 4;
        --it;
        SEQAN_ASSERT_EQ(*it, vec[3]);
    }
    // goFurther/operator+=/operator-=
    {
        std::vector<int> vec;
        for (int i = 0; i < 10; ++i)
            vec.push_back(i);

        TIterator it(vec.begin() + 4);
        goFurther(it, 2);
        SEQAN_ASSERT_EQ(*it, vec[6]);

        it = vec.begin() + 4;
        goFurther(it, -2);
        SEQAN_ASSERT_EQ(*it, vec[2]);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_arithmetics)
{
    using namespace seqan;

    typedef Iter<std::vector<int>, StdIteratorAdaptor> TIterator;

    std::vector<int> vec;
    for (int i = 0; i < 10; ++i)
        vec.push_back(i);

    TIterator it = vec.begin() + 4;

    TIterator it2 = vec.begin() + 6;
    SEQAN_ASSERT_EQ(&*it2, &vec[6]);

    it2 = it - 2;
    SEQAN_ASSERT_EQ(&*it2, &vec[2]);

    it2 = it + 2;
    SEQAN_ASSERT_EQ(it2 - it, 2);
}

// --------------------------------------------------------------------------
// Tests for Adaptor Iterator
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_metafunctions)
{
    using namespace seqan;

    // Pointers.
    {
        typedef Iter<int *, AdaptorIterator<int *> > TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Value<TIterator>::Type, int>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int &>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, int *>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Const-Pointers.
    {
        typedef Iter<int const *, AdaptorIterator<int const *> > TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        // TODO(holtgrew): This is inconsistent
        b = IsSameType<typename Value<TIterator>::Type, int const>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, int const *>::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_constructors)
{
    using namespace seqan;

    typedef Iter<int *, AdaptorIterator<int *> > TIterator;

    int container[] = { 0, 1, 2, 3 };

    // Default constructor.
    TIterator it;
    // From container and iterator.
    TIterator it2(container, &container[2]);
    SEQAN_ASSERT_EQ(it2.data_container, &container[0]);
    SEQAN_ASSERT_EQ(it2.data_iterator, &container[2]);
    // Copy constructor.
    TIterator it3(it2);
    SEQAN_ASSERT_EQ(it3.data_container, it2.data_container);
    SEQAN_ASSERT_EQ(it3.data_iterator, it2.data_iterator);
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_transport)
{
    using namespace seqan;

    typedef Iter<int *, AdaptorIterator<int *> > TIterator;

    int container[] = { 0, 1, 2, 3 };

    TIterator it(&container[0], &container[0]);

    // assign()
    {
        TIterator it2;
        assign(it2, it);
        SEQAN_ASSERT_EQ(it.data_container, it2.data_container);
        SEQAN_ASSERT_EQ(it.data_iterator, it2.data_iterator);
    }
    // set()
    {
        TIterator it2;
        seqan::set(it2, it);
        SEQAN_ASSERT_EQ(it.data_container, it2.data_container);
        SEQAN_ASSERT_EQ(it.data_iterator, it2.data_iterator);
    }
    // move()
    {
        TIterator it2;
        seqan::move(it2, it);
        SEQAN_ASSERT_EQ(it.data_container, it2.data_container);
        SEQAN_ASSERT_EQ(it.data_iterator, it2.data_iterator);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_transport_value)
{
    using namespace seqan;

    typedef Iter<CDStruct *, AdaptorIterator<CDStruct *> > TIterator;

    // assignValue()
    {
        CDStruct values[5];
        resetCDStructStatics();

        TIterator it(&values[0], &values[0]);
        assignValue(it, values[2]);

        SEQAN_ASSERT_EQ(it->copiedFrom, -1);
        SEQAN_ASSERT_EQ(it->movedFrom, -1);
        SEQAN_ASSERT_EQ(it->setFrom, -1);
        SEQAN_ASSERT_EQ(it->assignedFrom, values[2].id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &values[2]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 1);
    }
    // moveValue()
    {
        CDStruct values[5];
        resetCDStructStatics();

        TIterator it(&values[0], &values[0]);
        moveValue(it, values[2]);

        SEQAN_ASSERT_EQ(it->copiedFrom, -1);
        SEQAN_ASSERT_EQ(it->movedFrom, values[2].id);
        SEQAN_ASSERT_EQ(it->setFrom, -1);
        SEQAN_ASSERT_EQ(it->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &values[2]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    // Implementation of setValue() cannot work for all cases.
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_movement)
{
    using namespace seqan;

    typedef Iter<int *, AdaptorIterator<int *> > TIterator;

    // goNext/operator++
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        {
            TIterator it(&arr[0], &arr[4]);
            goNext(it);
            SEQAN_ASSERT_EQ(it.data_iterator, &arr[5]);
        }
        {
            TIterator it(&arr[0], &arr[4]);
            it++;
            SEQAN_ASSERT_EQ(it.data_iterator, &arr[5]);
        }
        {
            TIterator it(&arr[0], &arr[4]);
            ++it;
            SEQAN_ASSERT_EQ(it.data_iterator, &arr[5]);
        }
    }
    // goPrevious/operator--
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        {
            TIterator it(&arr[0], &arr[4]);
            goPrevious(it);
            SEQAN_ASSERT_EQ(it.data_iterator, &arr[3]);
        }
        {
            TIterator it(&arr[0], &arr[4]);
            it--;
            SEQAN_ASSERT_EQ(it.data_iterator, &arr[3]);
        }
        {
            TIterator it(&arr[0], &arr[4]);
            --it;
            SEQAN_ASSERT_EQ(it.data_iterator, &arr[3]);
        }
    }
    // goFurther/operator+=/operator-=
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        {
            TIterator it(&arr[0], &arr[4]);
            goFurther(it, 2);
            SEQAN_ASSERT_EQ(it, &arr[6]);
        }
        {
            TIterator it(&arr[0], &arr[4]);
            goFurther(it, -2);
            SEQAN_ASSERT_EQ(it.data_iterator, &arr[2]);
        }
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_arithmetics)
{
    using namespace seqan;

    int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    typedef Iter<int *, AdaptorIterator<int *> > TIterator;

    TIterator it(&arr[0], &arr[4]);

    TIterator it2 = it + 2;
    SEQAN_ASSERT_EQ(it2, &arr[6]);

    it2 = it - 2;
    SEQAN_ASSERT_EQ(it2, &arr[2]);

    it2 = it + 2;
    SEQAN_ASSERT_EQ(it2 - it, 2);
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_rooted_functions)
{
    using namespace seqan;

    int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    typedef Iter<int *, AdaptorIterator<int *> > TIterator;

    TIterator it(&arr[0], &arr[2]);
    SEQAN_ASSERT_EQ(container(it), &arr[0]);
    // TODO(holtgrew): Tests for atEnd() and position() depend on an actual adaption which is in the sequence module and we do not want to have a dependency on this here in the test for basic.  Maybe add a test for this in the test for sequence module.
}

// --------------------------------------------------------------------------
// Tests for Positional Iterator
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_position_metafunctions)
{
    using namespace seqan;

    // Pointers.
    {
        typedef Iter<int *, PositionIterator> TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Value<TIterator>::Type, int>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int &>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, int *>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Const-Pointers.
    {
        typedef Iter<int const *, PositionIterator> TIterator;

        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);

        // TODO(holtgrew): This is inconsistent
        b = IsSameType<typename Value<TIterator>::Type, int const>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<typename Container<TIterator>::Type, int const *>::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_constructors)
{
    using namespace seqan;

    typedef Iter<int *, PositionIterator> TIterator;
    //typedef Iter<int *, AdaptorIterator<int *> > TIterator2;

    int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    // Default constructor.
    TIterator itD;
    // Construct from container and position.
    TIterator it(&arr[0], 2);
    SEQAN_ASSERT_EQ(it.data_container, &arr[0]);
    SEQAN_ASSERT_EQ(position(it), 2u);
    // TODO(holtgrew): Deactivated because of actually missing begin() function from adaption. Maybe test this in the sequence module?
    // Copy construct from iterator of different type to same container
    // that supports the position() function.
    // TIterator2 it2(&arr[0], &arr[2]);
    // TIterator it3(it2);
    // SEQAN_ASSERT_EQ(it3.data_container, &arr[0]);
    // SEQAN_ASSERT_EQ(position(it3), 2);
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_transport)
{
    using namespace seqan;

    typedef Iter<int *, PositionIterator> TIterator;

    int container[] = { 0, 1, 2, 3 };

    TIterator it(&container[0], 0);

    // assign()
    {
        TIterator it2;
        assign(it2, it);
        SEQAN_ASSERT_EQ(it.data_container, it2.data_container);
        SEQAN_ASSERT_EQ(it.data_position, it2.data_position);
    }
    // set()
    {
        TIterator it2;
        seqan::set(it2, it);
        SEQAN_ASSERT_EQ(it.data_container, it2.data_container);
        SEQAN_ASSERT_EQ(it.data_position, it2.data_position);
    }
    // move()
    {
        TIterator it2;
        seqan::move(it2, it);
        SEQAN_ASSERT_EQ(it.data_container, it2.data_container);
        SEQAN_ASSERT_EQ(it.data_position, it2.data_position);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_transport_value)
{
    using namespace seqan;

    typedef Iter<CDStruct *, PositionIterator> TIterator;

    // assignValue()
    {
        CDStruct values[5];
        resetCDStructStatics();

        TIterator it(&values[0], 0);
        
        assignValue(it, values[2]);

        SEQAN_ASSERT_EQ(it->copiedFrom, -1);
        SEQAN_ASSERT_EQ(it->movedFrom, -1);
        SEQAN_ASSERT_EQ(it->setFrom, -1);
        SEQAN_ASSERT_EQ(it->assignedFrom, values[2].id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &values[2]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 1);
    }
    // moveValue()
    {
        CDStruct values[5];
        resetCDStructStatics();

        TIterator it(&values[0], 0);
        moveValue(it, values[2]);

        SEQAN_ASSERT_EQ(it->copiedFrom, -1);
        SEQAN_ASSERT_EQ(it->movedFrom, values[2].id);
        SEQAN_ASSERT_EQ(it->setFrom, -1);
        SEQAN_ASSERT_EQ(it->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &values[2]);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    // Implementation of setValue() cannot work for all cases.
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_movement)
{
    using namespace seqan;

    typedef Iter<int *, PositionIterator> TIterator;

    // goNext/operator++
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        {
            TIterator it(&arr[0], 4);
            goNext(it);
            SEQAN_ASSERT_EQ(it.data_position, 5u);
        }
        {
            TIterator it(&arr[0], 4);
            it++;
            SEQAN_ASSERT_EQ(it.data_position, 5u);
        }
        {
            TIterator it(&arr[0], 4);
            ++it;
            SEQAN_ASSERT_EQ(it.data_position, 5u);
        }
    }
    // goPrevious/operator--
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        {
            TIterator it(&arr[0], 4);
            goPrevious(it);
            SEQAN_ASSERT_EQ(it.data_position, 3u);
        }
        {
            TIterator it(&arr[0], 4);
            it--;
            SEQAN_ASSERT_EQ(it.data_position, 3u);
        }
        {
            TIterator it(&arr[0], 4);
            --it;
            SEQAN_ASSERT_EQ(it.data_position, 3u);
        }
    }
    // goFurther/operator+=/operator-=
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        {
            TIterator it(&arr[0], 4);
            goFurther(it, 2);
            SEQAN_ASSERT_EQ(it.data_position, 6u);
        }
        {
            TIterator it(&arr[0], 4);
            goFurther(it, -2);
            SEQAN_ASSERT_EQ(it.data_position, 2u);
        }
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_arithmetics)
{
    using namespace seqan;

    int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    typedef Iter<int *, PositionIterator> TIterator;

    TIterator it(&arr[0], 4);

    TIterator it2 = it + 2;
    SEQAN_ASSERT_EQ(it2.data_position, 6u);

    it2 = it - 2;
    SEQAN_ASSERT_EQ(it2.data_position, 2u);

    it2 = it + 2;
    SEQAN_ASSERT_EQ(it2 - it, 2);
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_rooted_functions)
{
    using namespace seqan;

    int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    typedef Iter<int *, PositionIterator> TIterator;

    TIterator it(&arr[0], 2);
    SEQAN_ASSERT_EQ(container(it), &arr[0]);
    // TODO(holtgrew): Tests for atEnd() and position() depend on an actual adaption which is in the sequence module and we do not want to have a dependency on this here in the test for basic.  Maybe add a test for this in the test for sequence module.
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_H_
