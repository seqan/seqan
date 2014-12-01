// ===========================================================================
//                 SGIP - Solution of Graph Isomorphism Problem
// ===========================================================================
// Copyright (C) 2012 by Jialu Hu
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ===========================================================================
// Author: Jialu Hu <Jialu.Hu@fu-berlin.de>
// ===========================================================================

#ifndef APPS_SGIP_SGIP_OUTPUT_H_
#define APPS_SGIP_SGIP_OUTPUT_H_

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>

#include <seqan/basic.h>
#include <seqan/sequence.h>


// ============================================================================
// Forward
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Matrix_;
typedef seqan::Tag<Matrix_> GMatrix;

struct General_;
typedef seqan::Tag<General_> General;

struct GeneralMap_;
typedef seqan::Tag<GeneralMap_>  GeneralMap;

struct DegreeMap_;
typedef seqan::Tag<DegreeMap_> DegreeMap;

struct ParityMap_;
typedef seqan::Tag<ParityMap_> ParityMap;

struct CheckMap_;
typedef seqan::Tag<CheckMap_> CheckMap;

struct FFFF_;
typedef seqan::Tag<FFFF_> FFFF;

// ============================================================================
// Function
// ============================================================================

// --------------------------------------------------------------------------
// Function outputLabel()
// --------------------------------------------------------------------------

template <typename TStr, typename TSpec>
void outputLabel(TStr &, seqan::Tag<TSpec>);

// Function outputLabel for String.
template <typename TStr>
void outputLabel(TStr & str, General const &)
{
    using namespace seqan;

    typedef typename Iterator<TStr>::Type TIterator;
    TIterator it = begin(str,Standard());
    TIterator itEnd = end(str,Standard());
    std::cout << "String: ";
    while (it != itEnd)
    {
        std::cout << getValue(it) << " ";
        goNext(it);
    }
    std::cout << std::endl;
    std::cout << "length: " << length(str) << std::endl;
}

template <typename TValue, typename TSpec, typename TTag>
void outputLabel(seqan::String<TValue, TSpec> &, TTag const &){}

template <typename TValue, typename TSpec>
void outputLabel(seqan::String<TValue, TSpec> & mat, General const &)
{
    using namespace seqan;

    typedef typename Iterator<seqan::String<TValue, TSpec> >::Type TIterator;
    TIterator it = begin(mat);
    std::cout << "String: ";
    while (!atEnd(it))
    {
        std::cout << getValue(it) << " ";
        goNext(it);
    }
    std::cout << std::endl;
    std::cout << "length: " << length(mat) << std::endl;
}

template <typename TValue, typename TSpec>
void outputLabel(seqan::String<TValue, TSpec> & mat, size_t n)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (mat[i + j * n])
                std::cout << j << " " << i << std::endl;
            if (mat[j + i * n])
                std::cout << i << " " << j << std::endl;
        }
    }
}

template <typename TValue, typename TSpec>
void outputLabel(seqan::String<TValue, TSpec> & mat, FFFF const &)
{
    size_t i = 0, j;
    size_t len = length(mat);
    while (i < len)
    {
        unsigned num = 0;
        for (j = 0; j < 16; j++, i++)
        {
            if (i < len)
            {
                if (mat[i])
                    num = num * 2 + 1;
                else
                    num = num * 2;
            }
            else
            {
                break;
            }
        }
        std::cout << std::setbase(16) << num << " ";
    }
    std::cout << std::setbase(10) << std::endl;
}

// Function outputLabel for StringSet.
template <typename TString>
void outputLabel(seqan::StringSet<TString> & mat, General const &)
{
    using namespace seqan;

    typedef typename Iterator<StringSet<TString>, Rooted>::Type TIterator;
    typedef typename Iterator<TString, Rooted>::Type TStrIterator;
    TIterator it = begin(mat);
    while (!atEnd(it))
    {
        TStrIterator ti = begin(value(it));
        while (!atEnd(ti))
        {
            std::cout << getValue(ti) << " ";
            goNext(ti);
        }
        std::cout << std::endl;
        goNext(it);
    }
}

template <typename TString>
void outputLabel(seqan::StringSet<TString> & mat, size_t n)
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            std::cout << mat[i][j] << " ";
        std::cout << std::endl;
    }
}

// --------------------------------------------------------------------------
// Function outputMap()
// --------------------------------------------------------------------------

// Output map structure type.
template <typename TMap>
void outputMap(TMap & map, GeneralMap const &)
{
    typedef typename TMap::iterator TIterator;
    TIterator it = map.begin();
    TIterator itEnd = map.end();
    for (; it < itEnd; it++)
        std::cout << it->first << "                   " << it->second << std::endl;
}

template <typename TKey, typename TValue>
void outputMap(std::unordered_map<TKey, TValue> & map, ParityMap const &)
{
    using namespace seqan;
    typedef typename std::unordered_map<TKey, TValue>::iterator TIterator;
    typedef typename Iterator<TValue, Rooted>::Type    TStrIterator;
    TIterator it = map.begin();
    TIterator itEnd = map.end();
    for (; it != itEnd; it++)
    {
        std::cout << it->first << " ";
        TValue str = it->second;
        TStrIterator ip = begin(str);
        TStrIterator ipEnd = end(str);
        while (ip != ipEnd)
        {
            std::cout << getValue(ip) << " ";
            goNext(ip);
        }
        std::cout << std::endl;
    }
}

template <typename TKey, typename TValue>
void outputMap(std::unordered_map<TKey, TValue> & map, CheckMap const &)
{
    using namespace seqan;
    typedef typename std::unordered_map<TKey, TValue>::iterator TIterator;
    //typedef typename Iterator<TValue, Rooted>::Type    TStrIterator;
    TIterator it = map.begin();
    TIterator itEnd = map.end();
    for (; it != itEnd; it++)
        std::cout << it->second << " ";
    std::cout << std::endl;
}

template <typename TKey, typename TValue>
void outputMap(std::unordered_map<TKey, TValue> & map, DegreeMap const &)
{
    using namespace seqan;
    typedef typename std::unordered_map<TKey, TValue>::iterator TIterator;
    //typedef typename Iterator<TValue, Rooted>::Type    TStrIterator;
    TIterator it = map.begin();
    TIterator itEnd = map.end();
    for (; it != itEnd; it++)
        std::cout << it->first << " " << it->second.inDegree << " " << it->second.outDegree << std::endl;
    std::cout << std::endl;
}

template <typename TKey, typename TValue, typename TComp>
void outputMap(std::multimap<TKey, TValue, TComp> & degreemap, DegreeMap const &)
{
    typedef typename std::multimap<TKey, TValue, TComp>::iterator TIterator;
    TIterator it = degreemap.begin();
    TIterator itEnd = degreemap.end();
    for (; it != itEnd; it++)
        std::cout << it->second << " " << it->first.inDegree << " " << it->first.outDegree << std::endl;
}

#endif  // #ifndef APPS_SGIP_SGIP_OUTPUT_H_
