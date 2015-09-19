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

template<typename TSeq>
struct ReadBase 
{
    using seqType = TSeq;

    TSeq seq;
    std::string id;
    int demuxResult;

    ReadBase() = default;
    ReadBase(const ReadBase& rhs) = default;

    bool operator==(const ReadBase& rhs) const
    {
        return seq == rhs.seq && id == rhs.id && demuxResult == rhs.demuxResult;
    }

    ReadBase& operator=(const ReadBase&& rhs)
    {
        seq = std::move(rhs.seq);
        id = std::move(rhs.id);
        demuxResult = rhs.demuxResult;
        return *this;
    }
};

template<typename TSeq>
struct Read : ReadBase<TSeq>
{
};

template<typename TSeq>
struct ReadMultiplex : ReadBase<TSeq>
{
    std::string demultiplex;

    ReadMultiplex() = default;
    ReadMultiplex(const ReadMultiplex& rhs) = default;

    bool operator==(const ReadMultiplex& rhs) const
    {
        return ReadBase::operator==(rhs) && demultiplex == rhs.demultiplex;
    }
    ReadMultiplex& operator=(const ReadMultiplex&& rhs)
    {
        ReadBase::operator=(std::move(rhs));
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

    bool operator==(const ReadPairedEnd& rhs) const
    {
        return ReadBase::operator==(rhs) && seqRev == rhs.seqRev && idRev == rhs.idRev;
    }
    ReadPairedEnd& operator=(const ReadPairedEnd&& rhs)
    {
        ReadBase::operator=(std::move(rhs));
        seqRev = std::move(rhs.seqRev);
        idRev = std::move(rhs.idRev);
        return *this;
    }
};

template<typename TSeq>
struct ReadMultiplexPairedEnd : ReadPairedEnd<TSeq>
{
    std::string demultiplex;

    ReadMultiplexPairedEnd() = default;
    ReadMultiplexPairedEnd(const ReadMultiplexPairedEnd& rhs) = default;

    bool operator==(const ReadMultiplexPairedEnd& rhs) const
    {
        return ReadPairedEnd::operator==(rhs) && demultiplex == rhs.demultiplex;
    }
    ReadMultiplexPairedEnd& operator=(const ReadMultiplexPairedEnd&& rhs)
    {
        ReadPairedEnd::operator=(std::move(rhs));
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
    auto temp = std::make_tuple(f(std::forward<Ts>(a))...); 
}

template<typename Trem, typename TremVal, typename TRead>
auto _eraseSeqs(const Trem& rem, const TremVal remVal, std::vector<TRead>& reads)
{
    const auto numRemoveElements = std::count(rem.begin(), rem.end(), remVal);
    const auto beginAddr = (void*)&*(reads.begin());
    std::remove_if(reads.begin(), reads.end(),
        [&rem, &beginAddr, remVal](const auto& element) {
        return rem[&element - beginAddr] == remVal;});
    reads.resize(reads.size() - numRemoveElements);
    return numRemoveElements;
}

template<typename Trem, typename TremVal, typename... TContainer>
auto _eraseSeqs(const Trem& rem, const TremVal remVal, TContainer&&... container)
{
    const auto numRemoveElements = std::count(begin(rem), end(rem), remVal);
    auto eraseElements = [&rem, numRemoveElements, remVal](auto& seq)  // erase Elements using the remove erase idiom
    {
        const auto beginAddr = (void*)&*seqan::begin(seq);
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
            const auto beginAddr = (void*)&*begin(*seq);
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

std::string prefix2(const std::string& str, unsigned int len)
{
    return str.substr(0, len);
}

#endif // HELPERFUNCTIONS_H_