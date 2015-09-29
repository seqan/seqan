// ==========================================================================
//                             adapterTrimming.h
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
// This file provides the adapter trimming functionality of seqan-flexbar
// which is based in the implementation of the original flexbar program in 
// [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================

#ifndef HELPERFUNCTIONS_H_
#define HELPERFUNCTIONS_H_

#include <future>
#include <string>
#include <regex>

#include <seqan/basic.h>
#include <seqan/sequence.h>

struct MultiStringMatcher
{
    template <typename TContainer>
    MultiStringMatcher(const TContainer& patterns)
        : _patterns(patterns)
    {}
    template <typename TToken>
    int getMatchIndex(const TToken& token) const noexcept
    {
        unsigned int index = 0;
        for (const auto& pattern : _patterns)
        {
            if (pattern == token)   // likely
                return index;
            ++index;
        }
        return -1;
    }
private:
    const std::vector<std::string> _patterns;
};

// seqan->std interface functions

std::string prefix(const std::string& str, unsigned int len)
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

unsigned int length(const std::string& str)
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

template<typename R>
bool is_ready(std::future<R> const& f)
{
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
        std::back_inserter(result), std::plus<T>());
    return result;
}

template<typename TSeq>
struct ReadBase 
{
    using seqType = TSeq;

    TSeq seq;
    std::string id;
    int demuxResult;

    ReadBase() = default;
    ReadBase(const ReadBase& rhs) = default;
    ReadBase(ReadBase&& rhs) noexcept(std::is_nothrow_move_constructible<TSeq>::value)
    {
        seq = std::move(rhs.seq);
        id = std::move(rhs.id);
        demuxResult = rhs.demuxResult;
    }

    bool operator==(const ReadBase& rhs) const
    {
        return seq == rhs.seq && id == rhs.id && demuxResult == rhs.demuxResult;
    }
    ReadBase& operator=(const ReadBase& rhs) = default;
    ReadBase& operator=(const ReadBase&& rhs) noexcept(std::is_nothrow_move_assignable<TSeq>::value)
    {
        seq = std::move(rhs.seq);
        id = std::move(rhs.id);
        demuxResult = rhs.demuxResult;
        return *this;
    }
    inline unsigned int minSeqLen() const noexcept
    {
        return length(seq);
    }
};

template<typename TSeq>
struct Read : ReadBase<TSeq>
{
    //Read() = default;
    //Read(const Read& rhs) = default;
    //Read(Read&& rhs) noexcept
    //{
    //    ReadBase<TSeq>::ReadBase(std::move(rhs));
    //}
    //bool operator==(const Read& rhs) const
    //{
    //    return ReadBase<TSeq>::operator==(rhs);
    //}
    //Read& operator=(const Read& rhs) = default;
    //Read& operator=(const Read&& rhs) noexcept
    //{
    //    ReadBase<TSeq>::operator=(std::move(rhs));
    //    return *this;
    //}
};

template<typename TSeq>
struct ReadMultiplex : ReadBase<TSeq>
{
    TSeq demultiplex;

    ReadMultiplex() = default;
    ReadMultiplex(const ReadMultiplex& rhs) = default;
    ReadMultiplex(ReadMultiplex&& rhs) noexcept(std::is_nothrow_move_constructible<TSeq>::value)
    : ReadBase<TSeq>(std::move(rhs))
    {
        demultiplex = std::move(rhs.demultiplex);
    }

    bool operator==(const ReadMultiplex& rhs) const
    {
        return ReadBase<TSeq>::operator==(rhs) && demultiplex == rhs.demultiplex;
    }
    ReadMultiplex& operator=(const ReadMultiplex& rhs) = default;
    ReadMultiplex& operator=(const ReadMultiplex&& rhs)  noexcept(std::is_nothrow_move_assignable<TSeq>::value)
    {
        ReadBase<TSeq>::operator=(std::move(rhs));
        demultiplex = std::move(rhs.demultiplex);
        return *this;
    }
};

template<typename TSeq>
struct ReadPairedEnd : ReadBase<TSeq>
{
    TSeq seqRev;
    std::string idRev;

    ReadPairedEnd() = default;
    ReadPairedEnd(const ReadPairedEnd& rhs) = default;
    ReadPairedEnd(ReadPairedEnd&& rhs) noexcept(std::is_nothrow_move_constructible<TSeq>::value)
    : ReadBase<TSeq>(std::move(rhs))
    {
        seqRev = std::move(rhs.seqRev);
        idRev = std::move(rhs.idRev);
    }

    bool operator==(const ReadPairedEnd& rhs) const
    {
        return ReadBase<TSeq>::operator==(rhs) && seqRev == rhs.seqRev && idRev == rhs.idRev;
    }
    ReadPairedEnd& operator=(const ReadPairedEnd& rhs) = default;
    ReadPairedEnd& operator=(const ReadPairedEnd&& rhs)  noexcept(std::is_nothrow_move_assignable<TSeq>::value)
    {
        ReadBase<TSeq>::operator=(std::move(rhs));
        seqRev = std::move(rhs.seqRev);
        idRev = std::move(rhs.idRev);
        return *this;
    }
    inline unsigned int minSeqLen() const noexcept
    {
        return std::min(length(ReadBase<TSeq>::seq), length(seqRev));
    }
};

template<typename TSeq>
struct ReadMultiplexPairedEnd : ReadPairedEnd<TSeq>
{
    TSeq demultiplex;

    ReadMultiplexPairedEnd() = default;
    ReadMultiplexPairedEnd(const ReadMultiplexPairedEnd& rhs) = default;
    ReadMultiplexPairedEnd(ReadMultiplexPairedEnd&& rhs) noexcept(std::is_nothrow_move_constructible<TSeq>::value)
    : ReadPairedEnd<TSeq>(std::move(rhs))
    {  
        demultiplex = std::move(rhs.demultiplex);
    }

    bool operator==(const ReadMultiplexPairedEnd& rhs) const
    {
        return ReadPairedEnd<TSeq>::operator==(rhs) && demultiplex == rhs.demultiplex;
    }
    ReadMultiplexPairedEnd& operator=(const ReadMultiplexPairedEnd& rhs) = default;
    ReadMultiplexPairedEnd& operator=(const ReadMultiplexPairedEnd&& rhs)  noexcept(std::is_nothrow_move_assignable<TSeq>::value)
    {
        ReadPairedEnd<TSeq>::operator=(std::move(rhs));
        demultiplex = std::move(rhs.demultiplex);
        return *this;
    }
};

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
        [&rem, &it, remVal](const auto& element) {
        return *(it++) == remVal;}), reads.end());
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
