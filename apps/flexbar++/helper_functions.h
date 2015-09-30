// ==========================================================================
//                             helper_functions.h
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================


#ifndef HELPERFUNCTIONS_H_
#define HELPERFUNCTIONS_H_

#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "read.h"


// seqan->std interface functions

std::string prefix(const std::string& str, unsigned int len) noexcept
{
    return str.substr(0, len);
}

inline std::string seqanToStd(const seqan::Dna5QString& rhs) noexcept
{
    std::string ret;
    ret.resize(length(rhs));
    char c;
    std::transform(begin(rhs), end(rhs), ret.begin(), [&c](const auto& element){
        seqan::assign(c, element);
        return c;
    });
    return ret;
}

void append(std::string& str1, const std::string& str2)
{
    str1 += str2;
}

unsigned int length(const std::string& str) noexcept
{
    return str.size();
}

template <typename T>
unsigned int length(const std::vector<T>& vec) noexcept
{
    return vec.size();
}

template <typename T>
void resize(std::vector<T>& vec, unsigned int len)
{
    vec.resize(len);
}


void insert(std::string& dest, unsigned int k, const std::string& token)
{
    dest.insert(k, token);
}

//

// always use the forward read for barcode detection
template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value>>
std::string getPrefix(const TRead<TSeq>& read, unsigned len) noexcept
{
    return static_cast<const std::string>(prefix(read.seq, len));
}

template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value>>
std::string getPrefix(const TRead<TSeq>& read, unsigned len, bool = false) noexcept
{
    (void)len;
    return seqanToStd(read.demultiplex);
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) noexcept
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
        std::back_inserter(result), std::plus<T>());
    return result;
}


template <typename TDest, typename TSource>
void insertAfterFirstToken(TDest& dest, TSource&& source)
{
    unsigned const len = length(dest);
    for (unsigned k = 0;k < len;++k)
        if (dest[k] == ' ' || k == (len - 1))
        {
            insert(dest, k, std::forward<TSource>(source));
            break;
        }
}

template <class F, class... Ts>
void for_each_argument(F f, Ts&&... a) {
    // destructor of temps blocks until all threads are finished
    std::make_tuple(f(std::forward<Ts>(a))...); 
}

template<typename Trem, typename TremVal, typename TRead>
auto _eraseSeqs(const Trem& rem, const TremVal remVal, std::vector<TRead>& reads) noexcept
{
    const auto oldSize = reads.size();
    auto it = rem.cbegin();
    reads.erase(std::remove_if(reads.begin(), reads.end(),
        [&rem, &it, remVal](const auto& read) {(void)read;return *(it++) == remVal;}), reads.end());
    return oldSize - reads.size();
}

template<typename Trem, typename TremVal, typename... TContainer>
auto _eraseSeqs(const Trem& rem, const TremVal remVal, TContainer&&... container) noexcept
{
    const auto numRemoveElements = std::count(begin(rem), end(rem), remVal);
    auto eraseElements = [&rem, numRemoveElements, remVal](auto& seq)  // erase Elements using the remove erase idiom
    {
        const auto beginAddr = &*seqan::begin(seq);
        std::remove_if(seqan::begin(seq), seqan::end(seq),
            [&rem, &beginAddr, remVal](const auto& element) {
            return rem[&element - beginAddr] == remVal;});
        resize(seq, length(seq) - numRemoveElements);
        return 0;
    };
    for_each_argument(eraseElements, std::forward<TContainer>(container)...);
    return numRemoveElements;
}

// parallel execution brings no speedup here, because memory io is the bottle neck
template<typename Trem, typename TremVal, typename... TContainer>
auto _eraseSeqsDisabled(const Trem& rem, const TremVal remVal, TContainer&&... container)
{
    const auto numRemoveElements = std::count(begin(rem), end(rem), remVal);
    // eraseElementsWrapper returns a future 
    auto eraseElementsWrapper = [&rem, numRemoveElements, remVal](auto& seq)
    {
        auto eraseElements = [&rem, numRemoveElements, remVal](auto* seq)  // erase Elements using the remove erase idiom
        {
            const auto beginAddr = &*begin(*seq);
            std::remove_if(begin(*seq), end(*seq),
                [&rem, &beginAddr, remVal](const auto& element) {
                return rem[&element - beginAddr] == remVal;});
            resize(*seq, length(*seq) - numRemoveElements);
        };
        return std::async(std::launch::async | std::launch::deferred, eraseElements, &seq);
    };
    // blocks until all futures are ready
    for_each_argument(eraseElementsWrapper, std::forward<TContainer>(container)...);
    return numRemoveElements;
}


#endif // HELPERFUNCTIONS_H_
