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

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <ctime>

#include <map>
#include <set>

#include <seqan/map.h>
#include <seqan/misc/set.h>
#include <seqan/find.h>

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TMap>
void Test_STLSet()
{
    //typedef typename Value<TMap>::Type TValue;
    typedef typename Iterator<TMap>::Type TIterator;

    TMap map;

    SEQAN_ASSERT(length(map) == 0);

    insert(map, 2);
    SEQAN_ASSERT(hasKey(map,2));
    SEQAN_ASSERT(!hasKey(map,5));

    insert(map, 5);
    SEQAN_ASSERT(hasKey(map,5));

    insert(map, 1);
    SEQAN_ASSERT(hasKey(map, 1));

    int arr1[] = {1, 2, 5};
    TIterator it;
    it = begin(map);
    int i;
    for (i = 0; !atEnd(it); ++i)
    {
        SEQAN_ASSERT(i < 3);
        SEQAN_ASSERT(key(it) == arr1[i]);
        goNext(it);
    }
    SEQAN_ASSERT(i == 3);

    it = find(map, 2);

    TIterator it2(it);
    SEQAN_ASSERT(it == it2);
    SEQAN_ASSERT(key(it) == 2);
    SEQAN_ASSERT(key(it2) == 2);

    goNext(it2);
    SEQAN_ASSERT(it != it2);

    erase(map, it);
    SEQAN_ASSERT(!hasKey(map,2));

    erase(map, 5);
    SEQAN_ASSERT(!hasKey(map,5));

    clear(map);
    SEQAN_ASSERT(length(map) == 0);
}


template <typename TMap>
void Test_NoCargo_Single()
{
    //typedef typename Value<TMap>::Type TValue;
    typedef typename Iterator<TMap>::Type TIterator;

    TMap map;

    SEQAN_ASSERT(length(map) == 0);

    insert(map, 2);
    SEQAN_ASSERT(map[2]);
    SEQAN_ASSERT(!map[5]);

    TMap map2;
    map2 = map;
    SEQAN_ASSERT(map2[2]);
    SEQAN_ASSERT(!map2[5]);

    insert(map, 5);
    SEQAN_ASSERT(map[5]);

    insert(map, 1);
    SEQAN_ASSERT(hasKey(map, 1));

    int arr1[] = {1, 2, 5};
    TIterator it = begin(map);
    int i;
    for (i = 0; !atEnd(it); ++i)
    {
        SEQAN_ASSERT(i < 3);
        SEQAN_ASSERT(key(it) == arr1[i]);
        goNext(it);
    }
    SEQAN_ASSERT(i == 3);

    it = find(map, 2);
    erase(map, it);
    SEQAN_ASSERT(!map[2]);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TMap>
void Test_STLMap()
{
    typedef typename Value<TMap>::Type TValue;
    typedef typename Iterator<TMap>::Type TIterator;

    TMap map;

    SEQAN_ASSERT(length(map) == 0);

    insert(map, TValue(3, 1) );
    insert(map, TValue(1, 2) );
    insert(map, 10);
       map[10] = 3;

    insert(map, 2, 4);

    SEQAN_ASSERT(length(map) == 4);

    int arr1[] = {2, 4, 1, 3};
    int arr2[] = {1, 2, 3, 10};
    TIterator it;
    it = begin(map);
    int i;
    for (i = 0; !atEnd(it); ++i)
    {
        SEQAN_ASSERT(i < 4);
        SEQAN_ASSERT(cargo(it) == arr1[i]);
        SEQAN_ASSERT(key(it) == arr2[i]);
        //SEQAN_ASSERT(key(value(it)) == arr2[i]);
        goNext(it);
    }
    SEQAN_ASSERT(i == 4);

    map[8] = 5;
    map[2] = 6;

    SEQAN_ASSERT(mapValue(map, 8) == 5);

    int arr3[] = {2, 6, 1, 5, 3};
    int arr4[] = {1, 2, 3, 8, 10};
    i = 0;
    for (it = begin(map); it != end(map); ++it)
    {
        SEQAN_ASSERT(i < 5);
        SEQAN_ASSERT(cargo(it) == arr3[i]);
        SEQAN_ASSERT(key(it) == arr4[i]);
        ++i;
    }
    SEQAN_ASSERT(i == 5);

    it = find(map, 7);
    SEQAN_ASSERT(it);
    SEQAN_ASSERT(key(it) == 8);
    SEQAN_ASSERT(cargo(it) == 5);

    TIterator it2(it);
    SEQAN_ASSERT(it2);
    SEQAN_ASSERT(key(it2) == 8);
    SEQAN_ASSERT(cargo(it2) == 5);

    //SEQAN_ASSERT(key(value(it)) == 8);
    //SEQAN_ASSERT(cargo(*it) == 5);
    //cargo(value(it)) = 20;
    //SEQAN_ASSERT(cargo(it) == 20);

    SEQAN_ASSERT(it != begin(map));
    SEQAN_ASSERT(it == it);

    SEQAN_ASSERT(hasKey(map, 8));
    SEQAN_ASSERT(length(map) == 5);
    erase(map, it);
    SEQAN_ASSERT(!hasKey(map, 8));
    SEQAN_ASSERT(length(map) == 4);

    erase(map, 2);
    SEQAN_ASSERT(!hasKey(map, 2));
    SEQAN_ASSERT(length(map) == 3);

    erase(map, 8);
    SEQAN_ASSERT(length(map) == 3);

    clear(map);
    SEQAN_ASSERT(length(map) == 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TMap>
void Test_Cargo_Single()
{
    typedef typename Value<TMap>::Type TValue;
    typedef typename Iterator<TMap>::Type TIterator;

    TMap map;

    SEQAN_ASSERT(length(map) == 0);

    insert(map, TValue(3, 1) );
    insert(map, TValue(1, 2) );
    insert(map, 10, 3);
    insert(map, 2, 4);

    SEQAN_ASSERT(length(map) == 4);

    int arr1[] = {2, 4, 1, 3};
    int arr2[] = {1, 2, 3, 10};
    TIterator it = begin(map);
    int i;
    for (i = 0; !atEnd(it); ++i)
    {
        SEQAN_ASSERT(i < 4);
        SEQAN_ASSERT(cargo(it) == arr1[i]);
        SEQAN_ASSERT(key(it) == arr2[i]);
        SEQAN_ASSERT(key(value(it)) == arr2[i]);
        goNext(it);
    }
    SEQAN_ASSERT(i == 4);

    map[8] = 5;
    map[2] = 6;

    SEQAN_ASSERT(mapValue(map, 8) == 5);

    int arr3[] = {2, 6, 1, 5, 3};
    int arr4[] = {1, 2, 3, 8, 10};
    i = 0;
    for (it = begin(map); it != end(map); ++it)
    {
        SEQAN_ASSERT(i < 5);
        SEQAN_ASSERT(cargo(it) == arr3[i]);
        SEQAN_ASSERT(key(it) == arr4[i]);
        ++i;
    }
    SEQAN_ASSERT(i == 5);

    it = find(map, 7);
    SEQAN_ASSERT(it != 0);
    SEQAN_ASSERT(key(it) == 8);
    SEQAN_ASSERT(cargo(it) == 5);

    SEQAN_ASSERT(key(value(it)) == 8);
    SEQAN_ASSERT(cargo(*it) == 5);
    cargo(value(it)) = 20;
    SEQAN_ASSERT(cargo(it) == 20);

    SEQAN_ASSERT(it != begin(map));
    SEQAN_ASSERT(it == it);

    SEQAN_ASSERT(hasKey(map, 8));
    SEQAN_ASSERT(length(map) == 5);
    erase(map, it);
    SEQAN_ASSERT(!hasKey(map, 8));
    SEQAN_ASSERT(length(map) == 4);

    erase(map, 2);
    SEQAN_ASSERT(!hasKey(map, 2));
    SEQAN_ASSERT(length(map) == 3);

    erase(map, 8);
    SEQAN_ASSERT(length(map) == 3);

    clear(map);
    SEQAN_ASSERT(length(map) == 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TMap>
void Test_Cargo_Multiple()
{
    //typedef typename Value<TMap>::Type TValue;
    //typedef typename Iterator<TMap>::Type TIterator;

    TMap map;

    insert(map, 2, 4);
    SEQAN_ASSERT(length(map) == 1);

    insert(map, 2, 5);
    SEQAN_ASSERT(length(map) == 1);

    add(map, 2, 6);
    SEQAN_ASSERT(length(map) == 2);

    add(map, 3, 7);
    SEQAN_ASSERT(length(map) == 3);

    eraseAll(map, 2);
    SEQAN_ASSERT(length(map) == 1);
}

//////////////////////////////////////////////////////////////////////////////

void Test_Skiplist_Extra()
{
    typedef Pair<char, int> TValue;
    typedef Map<TValue, Skiplist< > > TMap;
    TMap map;

    map[8] = 5;
    SEQAN_ASSERT(value(map, 8) == TValue(8, 5) );

    map[2] = 3;
    TMap map2 = map;
    SEQAN_ASSERT(length(map2) == 2);
    SEQAN_ASSERT(map2[8] == 5);
    SEQAN_ASSERT(map2[2] == 3);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TMap1, typename TMap2>
void CompareMaps_(TMap1 & map1, TMap2 & map2)
{
    typedef typename Iterator<TMap1>::Type TIter1;
    typedef typename Iterator<TMap2>::Type TIter2;

    SEQAN_ASSERT(length(map1) == length(map2));

    TIter1 it1(map1);
    TIter2 it2(map2);

    while (!atEnd(it1))
    {
        SEQAN_ASSERT(!atEnd(it2));
        SEQAN_ASSERT(key(it1) == key(it2));
        SEQAN_ASSERT(cargo(it1) == cargo(it2));

        goNext(it1);
        goNext(it2);
    }

    SEQAN_ASSERT(atEnd(it2));
}

//____________________________________________________________________________

void Test_Skiplist_Stress()
{
    typedef Pair<int, int> TValue;
    typedef Map<TValue, Skiplist< > > TSkipList;
    typedef Iterator<TSkipList>::Type TSkipListIterator;

    typedef std::map<int, int> TSTDMap;
    typedef Iterator<TSTDMap>::Type TSTDMapIterator;

    TSkipList sl;
    TSkipList sl2;
    TSTDMap stdmap;

    TSkipListIterator sl_it(sl);
    TSTDMapIterator stdmap_it(stdmap);

    for (int j = 0; j < 100; ++j)
    {
        cout << ".";

        for (int i = 0; i < 500; ++i)
        {
            int _key = rand();
            int value = rand();

            insert(sl, _key, value);
            insert(stdmap, _key, value);
        }

        CompareMaps_(sl, stdmap);

        for (int i = 0; i < 100; ++i)
        {
            int _key = (rand() + (rand() << 16));

            sl_it = find(sl, _key);
            stdmap_it = find(stdmap, _key);

            SEQAN_ASSERT(atEnd(sl_it) == atEnd(stdmap_it));
            if (! atEnd(sl_it))
            {
                SEQAN_ASSERT(key(sl_it) == key(stdmap_it));
                SEQAN_ASSERT(cargo(sl_it) == cargo(stdmap_it));

                erase(sl, key(sl_it));
                erase(stdmap, key(stdmap_it));
            }
        }

        CompareMaps_(sl, stdmap);

        sl2 = sl;
        CompareMaps_(sl2, stdmap);
    }
    cout << "\n";

}

//////////////////////////////////////////////////////////////////////////////

void Main_Test_Map()
{
    Test_Cargo_Single< Map< Pair<char, int>, Skiplist< > > >();
    Test_Cargo_Single< Map< Pair<char, int>, VectorSet< > > >();

    Test_STLMap< std::map<char, int> >();

    Test_NoCargo_Single< Map< char, Skiplist< > > >();
    Test_NoCargo_Single< Map< char, VectorSet< > > >();

    Test_STLSet< std::set<char> >();

       Test_Cargo_Multiple< Map< Pair<char, int>, Skiplist< > > >();

    Test_Skiplist_Extra();

    Test_Skiplist_Stress();
}


SEQAN_DEFINE_TEST(test_map_map)
{
    Main_Test_Map();
}

