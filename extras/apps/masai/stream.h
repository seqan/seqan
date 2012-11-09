// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// This file contains the FileWriter Stream class specialization.
// ==========================================================================

#ifndef SANDBOX_ESIRAGUSA_APPS_MASAI_STREAM_H_
#define SANDBOX_ESIRAGUSA_APPS_MASAI_STREAM_H_

namespace seqan {

struct MMapWriter;
struct FileWriter;

template <typename TValue, typename TSpec>
struct FileStream;

template <typename TValue, typename TSpec>
struct Value<Stream<FileStream<TValue, TSpec> > >
{
    typedef TValue Type;
};

// FileStream with async. file I/O

template <typename TValue>
class Stream<FileStream<TValue, FileWriter> >
{
public:
    typedef File< Async<> >                             TFile;
    typedef SimpleBuffer<TValue>                        TBuffer;
    typedef PageFrame<TValue, TFile, Dynamic<> >        TPageFrame;
    typedef PageChain<TPageFrame>                       TPageChain;
    typedef	typename Iterator<TBuffer, Standard>::Type	TIterator;

    static const unsigned pageSize = 4*1024*1024;      // 1M entries per page
    static const unsigned numPages = 2;                // double-buffering

    TFile       file;
    TPageChain  chain;
    unsigned    nextPageNo;
    TIterator   it, itEnd;

    Stream():
        chain(numPages)
    {
        _initialize();
    }

    ~Stream()
    {
        close(*this);
    }

    inline void _initialize()
    {
        it = NULL;
        itEnd = NULL;
        nextPageNo = 0;
    }

    inline bool _advanceBuffer()
    {
        bool result = true;

        // write previously provided buffer to disk
        if (nextPageNo != 0)
            result &= _writeBuffer();

        // step one buffer ahead
        TPageFrame &pf = *chain.getReadyPage();

        // allocate page memory if not done already
        if (pf.begin == NULL)
            allocPage(pf, pageSize, file);

        // assign correct page number
        pf.pageNo = nextPageNo++;

        it = pf.begin;
        itEnd = pf.end;

        return result;
    }

    inline bool _stopWriting()
    {
        // has anything been written?
        if (nextPageNo == 0)
            return true;

        // write previously provided buffer to disk
        bool result = _writeBuffer();

        // wait for all outstanding operations to finish
        result &= _flush();

        _initialize();

        return result;
    }

private:

    inline bool _writeBuffer()
    {
        // shrink buffer size if it was not fully written (like the last buffer)
        if (it != itEnd)
            resize(*chain.last, it - chain.last->begin);
        return writePage(*chain.last, file);
    }

    inline bool _flush()
    {
        bool result = true;
        TPageFrame *p = chain.first;
        while (p != NULL)
        {
            result &= waitFor(*p);
            freePage(*p, file);
            p = p->next;
        }
        result &= flush(file);

        return result;
    }
};



// FileStream with memory mapping

template <typename TValue>
class Stream<FileStream<TValue, MMapWriter> >
{
public:
    typedef File< Async<> >                             TFile;
    typedef SimpleBuffer<TValue>                        TBuffer;
    typedef PageFrame<TValue, TFile, Dynamic<> >        TPageFrame;
    typedef PageChain<TPageFrame>                       TPageChain;
    //        typedef String<TPageFrame*>                         TPageChain;
    typedef	typename Iterator<TBuffer, Standard>::Type	TIterator;

    static const unsigned pageSize = 4*1024*1024;      // 1M entries per page
    static const unsigned numPages = 2;                // double-buffering

    TFile       file;
    TPageChain  chain;
    unsigned    nextPageNo;
    TIterator   it, itEnd;

    Stream():
    chain(numPages)
    {
        _initialize();
    }

    ~Stream()
    {
        close(*this);
    }

    inline void _initialize()
    {
        it = NULL;
        itEnd = NULL;
        nextPageNo = 0;
    }

    inline bool _advanceBuffer()
    {
        bool result = true;

        // write previously provided buffer to disk
        if (nextPageNo != 0)
            result &= _writeBuffer();

        // step one buffer ahead
        TPageFrame &pf = *chain.getReadyPage();

        // assign correct page number
        pf.pageSize = pageSize;
        pf.pageNo = nextPageNo++;

        // map next page into memory
        resize(file, nextPageNo * pageSize * sizeof(TValue));
        mapWritePage(pf, file, pageSize);

        it = pf.begin;
        itEnd = pf.end;

        return result;
    }

    inline bool _stopWriting()
    {
        // has anything been written?
        if (nextPageNo == 0)
            return true;

        // write previously provided buffer to disk
        bool result = _writeBuffer();

        // wait for all outstanding operations to finish
        result &= _flush();

        resize(file, (nextPageNo * pageSize - (itEnd - it)) * sizeof(TValue));

        _initialize();

        return result;
    }

private:

    inline bool _writeBuffer()
    {
        // shrink buffer size if it was not fully written (like the last buffer)
        if (it != itEnd)
            resize(*chain.last, it - chain.last->begin);
        //            std::cout << chain.last->pageNo << '\t';
        //            std::cout << *chain.last << std::endl;
        return unmapPage(*chain.last, file);
    }

    inline bool _flush()
    {
        bool result = true;
        TPageFrame *p = chain.first;
        while (p != NULL)
        {
            result &= waitFor(*p);
            freePage(*p, file);
            p = p->next;
        }
        result &= flush(file);

        return result;
    }
};




template <typename TValue, typename TSpec, typename TFilename, typename TFlags>
inline int
open(Stream<FileStream<TValue, TSpec> > & stream, TFilename const &filename, TFlags const &flags)
{
    stream._initialize();
    return open(stream.file, filename, flags);
}

template <typename TValue, typename TSpec>
inline int
streamWriteChar(Stream<FileStream<TValue, TSpec> > & stream, TValue const &value)
{
    if (stream.it == stream.itEnd)
    {
        if (!stream._advanceBuffer())
            return 1;
    }
    *(stream.it++) = value;
    return 0;
}

template <typename TValue, typename TSpec, typename TSourceIter, typename TCount>
inline size_t
streamWriteBlock(Stream<FileStream<TValue, TSpec> > & stream, TSourceIter srcIter, TCount count)
{
    typedef typename Size<Stream<FileStream<TValue, TSpec> > >::Type TSSize;

    TCount written = 0;

    while (written < count)
    {
        // do we need a new buffer?
        if (stream.it == stream.itEnd)
        {
            if (!stream._advanceBuffer())
                return written;
        }

        // how many values can we write into the buffer?
        TSSize cnt = _min((TSSize)(count - written), (TSSize)(stream.itEnd - stream.it));
        // write them
        arrayCopyForward(srcIter, srcIter + cnt, stream.it);

        stream.it += cnt;
        srcIter += cnt;
        written += cnt;
    }
    return written;
}

template <typename TValue, typename TSpec>
inline int
close(Stream<FileStream<TValue, TSpec> > & stream)
{
    stream._stopWriting();
    return close(stream.file);
}

}

#endif  // #ifndef SANDBOX_ESIRAGUSA_APPS_MASAI_STREAM_H_
