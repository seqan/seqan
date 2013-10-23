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
// fast block-based file reader/writer with buffer interfacce
// supports async I/O / memory-mapping
// ==========================================================================

#ifndef SEQAN_FILE_FILE_STREAM_H_
#define SEQAN_FILE_FILE_STREAM_H_

//#define SEQAN_DEBUG_FILESTREAM

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FileStreamBuffer
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec = Async<> >
struct FileStreamBuffer : public std::basic_streambuf<TValue>
{
    typedef std::basic_streambuf<TValue>                TBase;
    typedef std::char_traits<TValue>                    TTraits;

    typedef File<Async<> >                              TFile;
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

    // Iterators to the current position in the page, the end of the page and the end of the page as it is on the disk.
    // itReadEnd is relevant if the file is opened in read+write mode on the last page.  itEnd can point behind the end
    // of the file on the disk whereas itReadEnd points to the end of the file on the disk.  itReadEnd will be
    // incremented if something is written on the last page.

    unsigned modeFlags;  // bitmask with open flags
    bool eofFlag;  // whether we are at the end of the file
    int errorFlag;  // error flag, 0 on no error

    FileStreamBuffer()
    {
        _initialize();
    }

    ~FileStreamBuffer()
    {
        close(*this);
    }
    
    inline void _initialize()
    {
        this->setp(NULL, NULL);
        this->setg(NULL, NULL, NULL);

        fileSize = 0;
        framePtr = NULL;
        nextReadyPageNo = 0;
        nextFetchPageNo = 0;
        lastPageNo = -1;
        inChainNeedsHousekeeping = false;
        outChainNeedsHousekeeping = false;
        modeFlags = 0;
        eofFlag = false;
        errorFlag = 0;
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
            if (this->pptr != this->epptr)
                resize(oldFrame, this->pptr - this->pbase);

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

        // framePtr->begin+pageSize always points to the end of the page, whereas framePtr->end can point in the middle.
        setg(framePtr->begin, framePtr->begin, framePtr->end);
        setp(framePtr->begin, framePtr->begin + pageSize);

        return result;
    }

    // Flush everything to disk.  If preserveFileSize then the file size is preserved.
    inline bool _stop(bool preserveFileSize = false)
    {
        // has anything been written?
        if (nextReadyPageNo == 0)
            return true;

        // shrink current buffer size if it was not fully written (like the last buffer)
        TPageFrame &oldFrame = *framePtr;
        resize(oldFrame, std::max(this->pptr, this->egptr()) - oldFrame.begin);

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

        TSize _fileSize = fileSize;
        _initialize();
        if (preserveFileSize)
            fileSize = _fileSize;

        return result;
    }


    typename TTraits::int_type
    virtual overflow(typename TTraits::int_type val)
    {
        if (val != TTraits::eof())
        {
            *this->pptr = TTraits::to_char_type(val);
            ++this->pptr;
        }
        if (this->pptr == this->epptr && !_advanceBuffer())
            return TTraits::eof();
        else
            return val;
    }

    typename TTraits::int_type
    virtual underflow()
    {
        if (_advanceBuffer())
            return TTraits::to_int_type(*this->gptr);
        else
            return TTraits::eof();
    }

    TValue* eback() const   { return TBase::eback(); }
    TValue* gptr()  const   { return TBase::gptr();  }
    TValue* egptr() const   { return TBase::egptr(); }

    TValue* pbase() const   { return TBase::pbase(); }
    TValue* pptr()  const   { return TBase::pptr();  }
    TValue* epptr() const   { return TBase::epptr(); }
};

template <typename TValue, typename TDirection, typename TSpec>
class FileStream : public BasicStream<TValue, TDirection>::Type
{
    typedef FileStreamBuffer<TValue, TDirection, TSpec> TStreamBuffer;

    TStreamBuffer buffer;

    FileStream(const char *fileName, int openMode)
    {
        open(*this, fileName, openMode);
        this->init(buffer);
    }

    ~FileStream()
    {
        close(*this);
    }
};

template <typename TValue, typename TDirection, typename TSpec>
inline bool
open(FileStream<TValue, TDirection, TSpec> &stream, const char *fileName, int openMode)
{
    return open(stream.buffer, fileName, openMode);
}

template <typename TValue, typename TDirection, typename TSpec>
inline bool
close(FileStream<TValue, TDirection, TSpec> &stream)
{
    return close(stream.buffer);
}

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
struct Value<FileStreamBuffer<TValue, TDirection, TSpec> >
{
    typedef TValue Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
struct Position<FileStreamBuffer<TValue, TDirection, TSpec> >
{
    typedef FileStreamBuffer<TValue, TDirection, TSpec> TStream_;
    typedef typename TStream_::TFile TFile_;

    typedef typename Position<TFile_>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
struct Size<FileStreamBuffer<TValue, TDirection, TSpec> >
{
    typedef FileStreamBuffer<TValue, TDirection, TSpec> TStream_;
    typedef typename TStream_::TFile TFile_;

    typedef typename Size<TFile_>::Type Type;
};


// ============================================================================
// Functions
// ============================================================================

// Variant for File<TFileSpec> (fallback if MMap<> not available as on Windows).
template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_readBuffer(FileStreamBuffer<TValue, TDirection, TSpec> &buffer, TPageFrame &pf)
{
    // allocate page memory if not done already
    if (pf.begin == NULL)
        if (!allocPage(pf, buffer.pageSize, buffer.file))
            SEQAN_FAIL("ERROR: Could not allocate page of size %u.", (unsigned)buffer.pageSize);

    // only read in read-mode
    if (IsSameType<TDirection, Output>::VALUE || (int)pf.pageNo > buffer.lastPageNo)
    {
        resize(pf, 0);
        pf.status = READY;
        return true;
    }

    // shrink buffer size at the end of file
    buffer.egptr() = pf.end;
    if ((int)pf.pageNo < buffer.lastPageNo)
        resize(pf, buffer.pageSize);
    else
        resize(pf, buffer.fileSize - buffer.lastPageNo * buffer.pageSize);

    // read next page into memory
    return readPage(pf, buffer.file);
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_readBuffer(FileStreamBuffer<TValue, TDirection, MMap<TConfig> > &buffer, TPageFrame &pf)
{
    typedef FileStreamBuffer<TValue, TDirection, MMap<TConfig> >   TFileStream;
    typedef typename TFileStream::TFile                         TFile;
    typedef typename Size<TFile>::Type                          TSize;

    _setCapacity(pf, buffer.pageSize);
    TSize size = buffer.pageSize;

    // Get the page size on the disk, adjusted for last page; required for setting itReadEnd.
    unsigned pageSizeOnDisk = size;

    // how to handle pages crossing/beyond the end of file
    if ((int)pf.pageNo >= buffer.lastPageNo)
    {
        if (!IsSameType<TDirection, Input>::VALUE)
        {
            pageSizeOnDisk = buffer.fileSize - pf.pageNo * buffer.pageSize;

            // 1. buffer is writable:
            //
            // increase file size to next page boundary
            // and map the whole page
            buffer.lastPageNo = pf.pageNo;
            buffer.fileSize = (TSize)(buffer.lastPageNo + 1) * (TSize)buffer.pageSize;
            resize(buffer.file, buffer.fileSize);
        }
        else
        {
            // 2. buffer is read-only:
            //
            // adapt buffer to file size
            // map only the contents of the file
            if ((int)pf.pageNo == buffer.lastPageNo)
            {
                pageSizeOnDisk = size;
                size = buffer.fileSize - (TSize)buffer.lastPageNo * (TSize)buffer.pageSize;
            }
            else
            {
                pageSizeOnDisk = 0;
                // don't try to read pages behind the end the file
                resize(pf, 0);
                pf.status = READY;
                return true;
            }
        }
    }

    bool res = false;
    if (!IsSameType<TDirection, Input>::VALUE)
        res = mapWritePage(pf, buffer.file, size);
    else
        res = mapReadPage(pf, buffer.file, size);
    buffer.egptr() = pf.begin + pageSizeOnDisk;
    return res;
}
#endif


template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_preprocessBuffer(FileStreamBuffer<TValue, TDirection, TSpec> &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = READY;
    return true;
}

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_postprocessBuffer(FileStreamBuffer<TValue, TDirection, TSpec> &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = POSTPROCESSED;
    return true;
}


template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_writeBuffer(FileStreamBuffer<TValue, TDirection, TSpec> &buffer, TPageFrame &pf)
{
    // only write in write-mode
    if (IsSameType<TDirection, Input>::VALUE)
    {
        pf.status = UNUSED;
        return true;
    }

    // shrink buffer size if it was not fully written (like the last buffer)
    if (buffer.pptr() != buffer.epptr())
        resize(pf, buffer.pptr() - pf.begin);

    return writePage(pf, buffer.file);
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_writeBuffer(FileStreamBuffer<TValue, TDirection, MMap<TConfig> > &buffer, TPageFrame &pf)
{
    typedef FileStreamBuffer<TValue, TDirection, MMap<TConfig> >   TFileStream;
    typedef typename TFileStream::TFile                         TFile;
    typedef typename Size<TFile>::Type                          TSize;

    TSize pageLen = length(pf);
    pf.status = UNUSED;
    unmapPage(pf, buffer.file);

    if (!IsSameType<TDirection, Input>::VALUE && (pageLen < buffer.pageSize))
    {
        // shrink the file again, it was enlarged earlier (in _readBuffer)
        buffer.lastPageNo = pf.pageNo;
        buffer.fileSize = (TSize)pf.pageNo * (TSize)buffer.pageSize + pageLen;
        resize(buffer.file, buffer.fileSize);
    }
    return true;
}
#endif

// ----------------------------------------------------------------------------
// Function tryFetchBuffer()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
tryFetchBuffer(FileStreamBuffer<TValue, TDirection, TSpec> &buffer, TPageFrame &pf, bool &inProgress, bool doWait = false)
{
    inProgress = true;
    if (pf.status == UNUSED)
    {
        // buffer is idle, now read asynchronously
        pf.pageNo = buffer._getNextFetchPageNo();
        if (!_readBuffer(buffer, pf))
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
        if (!_preprocessBuffer(buffer, pf, doWait))
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
        if (!IsSameType<TDirection, Input>::VALUE)
            resize(pf, capacity(pf));

        // all done
        inProgress = false;
        return true;
    }

    SEQAN_FAIL("PageFrame has inconsistent state.");
    return false;
}

// ----------------------------------------------------------------------------
// Function tryReleaseBuffer()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
tryReleaseBuffer(FileStreamBuffer<TValue, TDirection, TSpec> &buffer, TPageFrame &pf, bool &inProgress, bool doWait = false)
{
    inProgress = true;
    #ifdef SEQAN_DEBUG_FILESTREAM
    std::cout << "Release:" << std::endl;
    #endif
    if (pf.status == READY)
    {
        pf.status = POSTPROCESSING;
        _postprocessBuffer(buffer, pf, doWait);
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
        if (!_writeBuffer(buffer, pf))
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

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDirection, typename TValue, typename TFilename, typename TFlags>
inline bool
open(FileStreamBuffer<TValue, TDirection, TSpec> & buffer, TFilename const &filename, TFlags const &flags)
{
    buffer._initialize();
    buffer.modeFlags = flags;  // must be after _initialize()!
    bool result = open(buffer.file, filename, flags);
    buffer.fileSize = (result)? size(buffer.file) / sizeof(TValue) : 0ul;
    if (buffer.fileSize != 0)
        buffer.lastPageNo = (buffer.fileSize - 1) / buffer.pageSize;
    else
        buffer.lastPageNo = -1;
    return result;
}

} // namespace seqan

#endif  // #ifndef SEQAN_FILE_FILE_STREAM_H_
