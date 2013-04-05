// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// fast block-based file reader/writer with stream interfacce
// supports async I/O / memory-mapping
// ==========================================================================

#ifndef SEQAN_FILE_FILE_STREAM_H_
#define SEQAN_FILE_FILE_STREAM_H_

//#define SEQAN_DEBUG_FILESTREAM

namespace seqan {

//////////////////////////////////////////////////////////////////////////////


// Stream Modes
struct WriteOnly
{
    typedef False   Readable;
    typedef True    Writable;
};

struct ReadOnly
{
    typedef True    Readable;
    typedef False   Writable;
};

struct ReadAndWrite
{
    typedef True    Readable;
    typedef True    Writable;
};

template <typename TMode, typename TSpec = MMap<>, typename TValue = char>
struct FileStream;

template <typename TMode, typename TSpec, typename TValue>
struct Value<Stream<FileStream<TMode, TSpec, TValue> > >
{
    typedef TValue Type;
};

//    template <typename TValue>
//    struct Size<Stream<FileStream<TSpec, TValue> > >:
//        public Size<typename Stream<FileStream<TSpec, TValue> >::TFile> {};
//
//    template <typename TValue>
//    struct Position<Stream<FileStream<TSpec, TValue> > >:
//        public Position<typename Stream<FileStream<TSpec, TValue> >::TFile> {};



    // FileStream with async. file I/O 
/*
template <typename TValue>
class Stream<FileStream<FileWriter, TValue> >
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
//            std::cout << chain.last->pageNo << '\t';
//            std::cout << *chain.last << std::endl;
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
*/

template <typename TMode, typename TConfig, typename TValue, typename TPageFrame>
inline bool
_readBuffer(Stream<FileStream<TMode, File<TConfig>, TValue> > &stream, TPageFrame &pf)
{
    // allocate page memory if not done already
    if (pf.begin == NULL)
        allocPage(pf, stream.pageSize, stream.file);

    // only read in read-mode
    if (!TMode::Readable::VALUE || (int)pf.pageNo > stream.lastPageNo)
    {
        resize(pf, 0);
        pf.status = READY;
        return true;
    }

    // shrink buffer size at the end of file
    if ((int)pf.pageNo < stream.lastPageNo)
        resize(pf, stream.pageSize);
    else
        resize(pf, stream.fileSize - stream.lastPageNo * stream.pageSize);

    // read next page into memory
    return readPage(pf, stream.file);
}

#ifndef PLATFORM_WINDOWS
template <typename TMode, typename TConfig, typename TValue, typename TPageFrame>
inline bool
_readBuffer(Stream<FileStream<TMode, MMap<TConfig>, TValue> > &stream, TPageFrame &pf)
{
    typedef Stream<FileStream<TMode, MMap<TConfig>, TValue> >   TFileStream;
    typedef typename TFileStream::TFile                         TFile;
    typedef typename Size<TFile>::Type                          TSize;

    _setCapacity(pf, stream.pageSize);
    TSize size = stream.pageSize;

    // how to handle pages crossing/beyond the end of file
    if ((int)pf.pageNo >= stream.lastPageNo)
    {
        if (TMode::Writable::VALUE)
        {
            // 1. stream is writable:
            //
            // increase file size to next page boundary
            // and map the whole page
            stream.lastPageNo = pf.pageNo;
            stream.fileSize = (TSize)(stream.lastPageNo + 1) * (TSize)stream.pageSize;
            resize(stream.file, stream.fileSize);
        }
        else
        {
            // 2. stream is read-only:
            //
            // adapt buffer to file size
            // map only the contents of the file
            if ((int)pf.pageNo == stream.lastPageNo)
            {
                size = stream.fileSize - (TSize)stream.lastPageNo * (TSize)stream.pageSize;
            }
            else
            {
                // don't try to read pages behind the end the file
                resize(pf, 0);
                pf.status = READY;
                return true;
            }
        }
    }

    if (TMode::Writable::VALUE)
        return mapWritePage(pf, stream.file, size);
    else
        return mapReadPage(pf, stream.file, size);
}
#endif


template <typename TMode, typename TValue, typename TSpec, typename TPageFrame>
inline bool
_preprocessBuffer(Stream<FileStream<TMode, TSpec, TValue> > &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = READY;
    return true;
}

template <typename TMode, typename TValue, typename TSpec, typename TPageFrame>
inline bool
_postprocessBuffer(Stream<FileStream<TMode, TSpec, TValue> > &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = POSTPROCESSED;
    return true;
}


template <typename TMode, typename TConfig, typename TValue, typename TPageFrame>
inline bool
_writeBuffer(Stream<FileStream<TMode, File<TConfig>, TValue> > &stream, TPageFrame &pf)
{
    // only write in write-mode
    if (!TMode::Writable::VALUE)
    {
        pf.status = UNUSED;
        return true;
    }

    // shrink buffer size if it was not fully written (like the last buffer)
    if (stream.it != stream.itEnd)
        resize(pf, stream.it - pf.begin);

    return writePage(pf, stream.file);
}

#ifndef PLATFORM_WINDOWS
template <typename TMode, typename TConfig, typename TValue, typename TPageFrame>
inline bool
_writeBuffer(Stream<FileStream<TMode, MMap<TConfig>, TValue> > &stream, TPageFrame &pf)
{
    typedef Stream<FileStream<TMode, MMap<TConfig>, TValue> >   TFileStream;
    typedef typename TFileStream::TFile                         TFile;
    typedef typename Size<TFile>::Type                          TSize;

    TSize pageLen = length(pf);
    pf.status = UNUSED;
    unmapPage(pf, stream.file);

    if (TMode::Writable::VALUE && (pageLen < stream.pageSize))
    {
        // shrink the file again, it was enlarged earlier (in _readBuffer)
        stream.lastPageNo = pf.pageNo;
        stream.fileSize = (TSize)pf.pageNo * (TSize)stream.pageSize + pageLen;
        resize(stream.file, stream.fileSize);
    }
    return true;
}
#endif

template <typename TValue, typename TMode, typename TSpec, typename TPageFrame>
inline bool
tryFetchBuffer(Stream<FileStream<TMode, TSpec, TValue> > &stream, TPageFrame &pf, bool &inProgress, bool doWait = false)
{
    inProgress = true;
    if (pf.status == UNUSED)
    {
        // buffer is idle, now read asynchronously
        pf.pageNo = stream._getNextFetchPageNo();
        if (!_readBuffer(stream, pf))
            return false;
    }

    if (pf.status == READING)
    {
        bool inProgress = false;
        if (doWait)
        {
            if (!waitFor(pf))
                return false;
        }
        else
        {
            if (!waitFor(pf, 0, inProgress))
                return false;
        }
        
        // still in operation?
        if (inProgress)
            return true;

        // was reading, now preprocessing
        pf.status = PREPROCESSING;
        if (!_preprocessBuffer(stream, pf, doWait))
            return false;
    }

    if (pf.status == PREPROCESSING)
    {
        // we get here only if doWait==false and the
        // asynchronous postprocessing is in progress
        SEQAN_ASSERT_NOT(doWait);
        return true;
    }

    if (pf.status == READY)
    {
        // if in write mode, extend half-read buffer to page size
        if (TMode::Writable::VALUE)
            resize(pf, capacity(pf));
        
        // all done
        inProgress = false;
        return true;
    }

    SEQAN_FAIL("PageFrame has inconsistent state.");
    return false;
}

template <typename TValue, typename TMode, typename TSpec, typename TPageFrame>
inline bool
tryReleaseBuffer(Stream<FileStream<TMode, TSpec, TValue> > &stream, TPageFrame &pf, bool &inProgress, bool doWait = false)
{
    inProgress = true;
    #ifdef SEQAN_DEBUG_FILESTREAM
    std::cout << "Release:" << std::endl;
    #endif
    if (pf.status == READY)
    {
        pf.status = POSTPROCESSING;
        _postprocessBuffer(stream, pf, doWait);
    }

    if (pf.status == POSTPROCESSING)
    {
        #ifdef SEQAN_DEBUG_FILESTREAM
        std::cout<<pf<<std::endl;
        #endif
        // we get here only if doWait==false and the
        // asynchronous postprocessing is in progress
        SEQAN_ASSERT_NOT(doWait);
        return true;
    }

    if (pf.status == POSTPROCESSED)
    {
        #ifdef SEQAN_DEBUG_FILESTREAM
        std::cout<<pf<<std::endl;
        #endif
        // buffer was postprocessed, now write asynchronously
        if (!_writeBuffer(stream, pf))
            return false;
    }

    if (pf.status == WRITING)
    {
        #ifdef SEQAN_DEBUG_FILESTREAM
        std::cout<<pf<<std::endl;
        #endif
        bool inProgress = false;
        if (doWait)
        {
            if (!waitFor(pf))
                return false;
        }
        else
        {
            if (!waitFor(pf, 0, inProgress))
                return false;
        }

        // still in operation?
        if (inProgress)
            return true;
    }

    // transform ready after write into unused
    if (pf.status == READY)
        pf.status = UNUSED;

    if (pf.status == UNUSED)
    {
        #ifdef SEQAN_DEBUG_FILESTREAM
        std::cout<<pf<<std::endl;
        #endif
        // all done
        inProgress = false;
        return true;
    }
    
    SEQAN_FAIL("PageFrame has inconsistent state.");
    return false;
}

// FileStream with memory mapping

template <typename TMode, typename TSpec, typename TValue>
class Stream<FileStream<TMode, TSpec, TValue> >
{
public:
    typedef File< Async<> >                             TFile;
    typedef typename Size<TFile>::Type                  TSize;
    typedef Buffer<TValue>                              TBuffer;
    typedef Buffer<TValue, PageFrame<TFile, Dynamic> >  TPageFrame;
    typedef PageChain<TPageFrame>                       TPageChain;
//        typedef String<TPageFrame*>                         TPageChain;
    typedef	typename Iterator<TBuffer, Standard>::Type	TIterator;

    static const unsigned pageSize = 4*1024*1024;      // 1M entries per page
//    static const unsigned pageSize = 4*1024;      // 1M entries per page
    static const unsigned numPages = 2;                // double-buffering

    TFile       file;
    TSize       fileSize;
    TPageFrame  *framePtr;
    TPageChain  inChain;
    TPageChain  outChain;

    unsigned    nextReadyPageNo;
    unsigned    nextFetchPageNo;
    int         lastPageNo;

    bool        inChainNeedsHousekeeping;
    bool        outChainNeedsHousekeeping;

    TIterator   it, itEnd;

    Stream()
    {
        _initialize();
    }

    ~Stream()
    {
        close(*this);
    }

    inline void _initialize()
    {
        fileSize = 0;
        framePtr = NULL;
        it = NULL;
        itEnd = NULL;
        nextReadyPageNo = 0;
        nextFetchPageNo = 0;
        lastPageNo = -1;
        inChainNeedsHousekeeping = false;
        outChainNeedsHousekeeping = false;
    }

    inline void _printStats()
    {
        std::cout << "in: " << inChain.size() << "\tout:" << outChain.size() << "\tframe:" << (__uint64)framePtr << std::endl;
    }

    inline int _getNextFetchPageNo()
    {
        return nextFetchPageNo++;
    }

    inline bool _injectUnusedBuffer(TPageFrame &frame)
    {
        bool inProgress;
        frame.status = UNUSED;
        if (!tryFetchBuffer(*this, frame, inProgress))
            return false;

        #ifdef SEQAN_DEBUG_FILESTREAM
        std::cout << "ADDED INTO IN-CHAIN (" << inChain.frames << "):" << std::endl;
        std::cout << frame << std::endl;
        #endif
        inChain.pushBack(frame);
        inChainNeedsHousekeeping |= inProgress;
        return true;
    }

    inline bool _houseKeeping()
    {
        bool inProgress;

        if (outChainNeedsHousekeeping)
        {
            bool someAreInProgress = false;
            for (TPageFrame *p = outChain.first; p != NULL;)
            {
                TPageFrame &page = *p;
                p = p->next;
                
                if (!tryReleaseBuffer(*this, page, inProgress))
                    return false;

                if (!inProgress)
                {
                    outChain.erase(page);
                    _injectUnusedBuffer(page);
                }
                else
                    someAreInProgress = true;
            }
            outChainNeedsHousekeeping = someAreInProgress;
        }

        if (inChainNeedsHousekeeping)
        {
            bool someAreInProgress = false;
            for (TPageFrame *p = inChain.first; p != NULL; p = p->next)
            {
                if (!tryFetchBuffer(*this, *p, inProgress))
                    return false;
                someAreInProgress |= inProgress;
            }
            inChainNeedsHousekeeping = someAreInProgress;
        }
        return true;
    }

    // ,---, ,---, ,---, ,---, ,---,
    // | 0 | | 1 | | 2 | | 3 | | 4 |
    // '---' '---' '---' '---' '---'
    //
    // 

    inline bool _advanceBuffer()
    {
        bool result = true;
        bool inProgress;

        // release old buffer
        if (framePtr != NULL)
        {
            // shrink current buffer size if it was not fully written (like the last buffer)
            TPageFrame &oldFrame = *framePtr;
            if (it != itEnd)
                resize(oldFrame, it - oldFrame.begin);

            // start to release the current buffer
            result &= tryReleaseBuffer(*this, oldFrame, inProgress);
            if (inProgress)
            {
                // writing is in progress
                // -> we put the frame into the out-chain
                #ifdef SEQAN_DEBUG_FILESTREAM
                std::cout << "ADDED INTO OUT-CHAIN (" << outChain.frames << "):" << std::endl;
                std::cout << oldFrame << std::endl;
                #endif
                outChain.pushBack(oldFrame);
                outChainNeedsHousekeeping |= inProgress;
            }
            else
                // writing finished immediately
                // -> we reinject the now unused buffer into the in-chain
                _injectUnusedBuffer(oldFrame);
        }

        // watch current progress
        // initiate next asynchronous task for each buffer
        // reinject written buffers for prefetching the next ones
        _houseKeeping();

        // fetch new buffer
        if (inChain.empty())
        {
            // we need more incoming buffers
            if (inChain.size() + outChain.size() < numPages)
            {
                // we can afford to create new buffers
                _injectUnusedBuffer(*(new TPageFrame()));
            }
            else
            {
                // all out-chain buffers are in progress (otherwise in-chain would not be empty)
                // so we take and wait for the first out-chain buffer
                TPageFrame &anyOutChainBuffer = outChain.popFront();
                tryReleaseBuffer(*this, anyOutChainBuffer, inProgress, true);
                _injectUnusedBuffer(anyOutChainBuffer);
            }
        }

        // we implicitly assume the front buffer
        // to be the next in sequence nextReadyPageNo, nextReadyPageNo+1, ...
        framePtr = &inChain.popFront();
        SEQAN_ASSERT_EQ(framePtr->pageNo, nextReadyPageNo);
        ++nextReadyPageNo;

        // now fetch-and-wait for new buffer
        result &= tryFetchBuffer(*this, *framePtr, inProgress, true);
        it = framePtr->begin;
        itEnd = framePtr->end;

        if (TMode::Writable::VALUE && (unsigned)(itEnd - it) < pageSize)
            itEnd = it + pageSize;

        return result;
    }

    inline bool _stop()
    {
        // has anything been written?
        if (nextReadyPageNo == 0)
            return true;

        // shrink current buffer size if it was not fully written (like the last buffer)
        TPageFrame &oldFrame = *framePtr;
        resize(oldFrame, it - oldFrame.begin);

        // start to release the current buffer
        bool dummy;
        bool result = tryReleaseBuffer(*this, oldFrame, dummy, true);
        freePage(oldFrame, file);
        delete &oldFrame;

        while (!outChain.empty())
        {
            result &= tryReleaseBuffer(*this, outChain.front(), dummy, true);
            freePage(outChain.front(), file);
            delete &outChain.popFront();
        }

        while (!inChain.empty())
        {
            // cancel(inChain.front(), file);
            freePage(inChain.front(), file);
            delete &inChain.popFront();
        }

        _initialize();

        return result;
    }
};


template <typename TMode, typename TSpec, typename TValue, typename TFilename, typename TFlags>
inline int
open(Stream<FileStream<TMode, TSpec, TValue> > & stream, TFilename const &filename, TFlags const &flags)
{
    stream._initialize();
    bool result = open(stream.file, filename, flags);
    stream.fileSize = (result)? size(stream.file) / sizeof(TValue) : 0ul;
    if (stream.fileSize != 0)
        stream.lastPageNo = (stream.fileSize - 1) / stream.pageSize;
    else
        stream.lastPageNo = -1;
    return result;
}

template <typename TMode, typename TSpec, typename TValue>
inline int
streamWriteChar(Stream<FileStream<TMode, TSpec, TValue> > & stream, TValue const &value)
{
    if (stream.it == stream.itEnd)
    {
        if (!stream._advanceBuffer())
            return 1;

        SEQAN_ASSERT_LT(stream.it, stream.itEnd);
    }
    *(stream.it++) = value;
    return 0;
}

template <typename TMode, typename TSpec, typename TValue, typename TSourceIter, typename TCount>
inline TCount
streamWriteBlock(Stream<FileStream<TMode, TSpec, TValue> > & stream, TSourceIter srcIter, TCount count)
{
    typedef typename Size<Stream<FileStream<TSpec, TValue> > >::Type TSSize;

    TSSize remaining = count;
    while (remaining != 0)
    {
        // do we need a new buffer?
        if (stream.it == stream.itEnd)
        {
            if (!stream._advanceBuffer())
                return count - remaining;
            
            SEQAN_ASSERT_LT(stream.it, stream.itEnd);
        }

        // how many values can we write into the buffer?
        TSSize cnt = _min(remaining, (TSSize)(stream.itEnd - stream.it));
        // write them
        arrayCopyForward(srcIter, srcIter + cnt, stream.it);

        stream.it += cnt;
        srcIter += cnt;
        remaining -= cnt;
    }
    return count;
}

template <typename TMode, typename TSpec, typename TValue>
inline int
streamReadChar(TValue &value, Stream<FileStream<TMode, TSpec, TValue> > & stream)
{
    if (stream.it == stream.itEnd)
    {
        // TODO(weese):
        // Stop this horrible return value inconsistency,
        // some stream functions return 1 on error and 0 on success
        // or the number of transferred values :-/
        //
        // until this is fixed we return EOF in case of an error
        
        if (!stream._advanceBuffer())
            return EOF;

        if (stream.it == stream.itEnd)
            return EOF;
    }
    value = *(stream.it++);
    return 0;
}

template <typename TDestIter, typename TMode, typename TSpec, typename TValue, typename TCount>
inline TCount
streamReadBlock(TDestIter dstIter, Stream<FileStream<TMode, TSpec, TValue> > & stream, TCount count)
{
    typedef typename Size<Stream<FileStream<TSpec, TValue> > >::Type TSSize;

    TSSize remaining = count;
    while (remaining != 0)
    {
        // do we need a new buffer?
        if (stream.it == stream.itEnd)
        {
            if (!stream._advanceBuffer())
                return count - remaining;

            if (stream.it == stream.itEnd)
                return EOF;
        }

        // how many values can we write into the buffer?
        TSSize cnt = _min(remaining, (TSSize)(stream.itEnd - stream.it));
        // write them
        arrayCopyForward(stream.it, stream.it + cnt, dstIter);

        stream.it += cnt;
        dstIter += cnt;
        remaining -= cnt;
    }
    return count;
}

template <typename TMode, typename TSpec, typename TValue>
inline int
close(Stream<FileStream<TMode, TSpec, TValue> > & stream)
{
    if (stream.file)
    {
        stream._stop();
        return close(stream.file);
    }
    return true;
}

} // namespace seqan

#endif  // #ifndef SEQAN_FILE_FILE_STREAM_H_
