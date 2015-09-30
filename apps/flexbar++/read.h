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


#ifndef READ_H_
#define READ_H_

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

#endif