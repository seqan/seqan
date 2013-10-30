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
// fast block-based file reader/writer with buffer interface
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
// Class FilePager
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec = Async<> >
struct FilePageTable
{
    typedef File<Async<> >                              TFile;
    typedef typename Size<TFile>::Type                  TSize;
    typedef Buffer<TValue>                              TBuffer;
    typedef Buffer<TValue, PageFrame<TFile, Dynamic> >  TPageFrame;
    typedef PageChain<TPageFrame>                       TPageChain;
//        typedef String<TPageFrame*>                         TPageChain;
    typedef	typename Iterator<TBuffer, Standard>::Type	TIterator;

    typedef short int                                   TPageId;
    typedef String<TPageId>                             TPageDir;

//    static const unsigned pageSize = 4*1024;      // 1M entries per page
    static const unsigned numPages = 2;                // double-buffering
    static const unsigned pageSize = 4*1024*1024;      // 1M entries per page

    TFile       file;
    TSize       fileSize;
    int         lastPageNo;

    TPageChain  inChain;
    TPageChain  outChain;

    bool        inChainNeedsHousekeeping;
    bool        outChainNeedsHousekeeping;

    inline void _printStats()
    {
        std::cout << "in: " << inChain.size() << "\tout:" << outChain.size() << std::endl;
    }
};


template <typename TValue, typename TDirection, typename TSpec = Async<> >
inline void
clear(FilePageTable<TValue, TDirection, TSpec> &pager)
{
    pager.lastPageNo = -1;
    pager.fileSize = 0;
    pager.inChainNeedsHousekeeping = false;
    pager.outChainNeedsHousekeeping = false;
}

// ============================================================================
// Functions
// ============================================================================

// Variant for File<TFileSpec> (fallback if MMap<> not available as on Windows).
template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_readBuffer(FilePageTable<TValue, TDirection, TSpec> &pager, TPageFrame &pf)
{
    // allocate page memory if not done already
    if (pf.begin == NULL)
        if (!allocPage(pf, pager.pageSize, pager.file))
            SEQAN_FAIL("ERROR: Could not allocate page of size %u.", (unsigned)pager.pageSize);

    // only read in read-mode
    if (IsSameType<TDirection, Output>::VALUE || (int)pf.pageNo > pager.lastPageNo)
    {
        resize(pf, 0);
        pf.status = READY;
        return true;
    }

    // shrink buffer size at the end of file
    if ((int)pf.pageNo < pager.lastPageNo)
        resize(pf, pager.pageSize);
    else
        resize(pf, pager.fileSize - pager.lastPageNo * pager.pageSize);

// TODO: move out
//    pager.setg(pf.begin, pf.begin, pf.end);

    // read next page into memory
    return readPage(pf, pager.file);
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_readBuffer(FilePageTable<TValue, TDirection, MMap<TConfig> > &pager, TPageFrame &pf)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

    _setCapacity(pf, pager.pageSize);
    TSize size = pager.pageSize;

    // Get the page size on the disk, adjusted for last page; required for setting itReadEnd.
    unsigned pageSizeOnDisk = size;

    // how to handle pages crossing/beyond the end of file
    if ((int)pf.pageNo >= pager.lastPageNo)
    {
        if (!IsSameType<TDirection, Input>::VALUE)
        {
            pageSizeOnDisk = pager.fileSize - pf.pageNo * pager.pageSize;

            // 1. buffer is writable:
            //
            // increase file size to next page boundary
            // and map the whole page
            pager.lastPageNo = pf.pageNo;
            pager.fileSize = (TSize)(pager.lastPageNo + 1) * (TSize)pager.pageSize;
            resize(pager.file, pager.fileSize);
        }
        else
        {
            // 2. buffer is read-only:
            //
            // adapt buffer to file size
            // map only the contents of the file
            if ((int)pf.pageNo == pager.lastPageNo)
            {
                pageSizeOnDisk = size;
                size = pager.fileSize - (TSize)pager.lastPageNo * (TSize)pager.pageSize;
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
        res = mapWritePage(pf, pager.file, size);
    else
        res = mapReadPage(pf, pager.file, size);

// TODO:move out
//    pager.setg(pf.begin, pf.begin, pf.begin + pageSizeOnDisk);

    return res;
}
#endif


template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_preprocessBuffer(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = READY;
    return true;
}

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_postprocessBuffer(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = POSTPROCESSED;
    return true;
}


template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_writeBuffer(FilePageTable<TValue, TDirection, TSpec> &pager, TPageFrame &pf)
{
    // only write in write-mode
    if (IsSameType<TDirection, Input>::VALUE)
    {
        pf.status = UNUSED;
        return true;
    }

    // shrink buffer size if it was not fully written (like the last buffer)
//TODO: move to streambuf
//    if (buffer.pptr() != buffer.epptr())
//        resize(pf, buffer.pptr() - pf.begin);

    return writePage(pf, pager.file);
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_writeBuffer(FilePageTable<TValue, TDirection, MMap<TConfig> > &pager, TPageFrame &pf)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

    TSize pageLen = length(pf);
    pf.status = UNUSED;
    unmapPage(pf, pager.file);

    if (!IsSameType<TDirection, Input>::VALUE && (pageLen < pager.pageSize))
    {
        // shrink the file again, it was enlarged earlier (in _readBuffer)
        pager.lastPageNo = pf.pageNo;
        pager.fileSize = (TSize)pf.pageNo * (TSize)pager.pageSize + pageLen;
        resize(pager.file, pager.fileSize);
    }
    return true;
}
#endif

// ----------------------------------------------------------------------------
// Function tryFetchBuffer()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
tryFetchBuffer(FilePageTable<TValue, TDirection, TSpec> &pager, TPageFrame &pf, bool &inProgress, bool doWait = false)
{
    inProgress = true;
    if (pf.status == UNUSED)
    {
        // buffer is idle, now read asynchronously
        pf.pageNo = pager._getNextFetchPageNo();
        if (!_readBuffer(pager, pf))
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
        if (!_preprocessBuffer(pager, pf, doWait))
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

template <typename TFilePageTable, typename TPageFrame>
inline bool
tryReleaseBuffer(TFilePageTable &pager, TPageFrame &pf, bool &inProgress, bool doWait = false)
{
    inProgress = true;
    #ifdef SEQAN_DEBUG_FILESTREAM
    std::cout << "Release:" << std::endl;
    #endif
    if (pf.status == READY)
    {
        pf.status = POSTPROCESSING;
        _postprocessBuffer(pager, pf, doWait);
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
        if (!_writeBuffer(pager, pf))
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

template <typename TFilePageTable, typename TPageFrame>
inline bool
_injectUnusedBuffer(TFilePageTable &pager, TPageFrame &frame)
{
    bool inProgress;
    frame.status = UNUSED;
    if (!tryFetchBuffer(pager, frame, inProgress))
        return false;

    #ifdef SEQAN_DEBUG_FILESTREAM
    std::cout << "ADDED INTO IN-CHAIN (" << pager.inChain.frames << "):" << std::endl;
    std::cout << frame << std::endl;
    #endif
    pager.inChain.pushBack(frame);
    pager.inChainNeedsHousekeeping |= inProgress;
    return true;
}

template <typename TValue, typename TDirection, typename TSpec>
inline bool
_houseKeeping(FilePageTable<TValue, TDirection, TSpec> &pager)
{
    typedef FilePageTable<TValue, TDirection, TSpec>    TFilePageTable;
    typedef typename TFilePageTable::TPageFrame         TPageFrame;

    bool inProgress;

    if (pager.outChainNeedsHousekeeping)
    {
        bool someAreInProgress = false;
        for (TPageFrame *p = pager.outChain.first; p != NULL;)
        {
            TPageFrame &page = *p;
            p = p->next;

            if (!tryReleaseBuffer(pager, page, inProgress))
                return false;

            if (!inProgress)
            {
                pager.outChain.erase(page);
                _injectUnusedBuffer(pager, page);
            }
            else
                someAreInProgress = true;
        }
        pager.outChainNeedsHousekeeping = someAreInProgress;
    }

    if (pager.inChainNeedsHousekeeping)
    {
        bool someAreInProgress = false;
        for (TPageFrame *p = pager.inChain.first; p != NULL; p = p->next)
        {
            if (!tryFetchBuffer(pager, *p, inProgress))
                return false;
            someAreInProgress |= inProgress;
        }
        pager.inChainNeedsHousekeeping = someAreInProgress;
    }
    return true;
}



// ----------------------------------------------------------------------------
// Class FileStreamBuffer
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec = Async<> >
struct FileStreamBuffer : public std::basic_streambuf<TValue>
{
    typedef std::basic_streambuf<TValue>                TBase;
    typedef std::char_traits<TValue>                    TTraits;

    typedef typename TBase::int_type                    TIntValue;
    typedef typename TBase::off_type                    TDifference;
    typedef typename TBase::pos_type                    TPosition;

    typedef FilePageTable<TValue, TDirection, TSpec>    TFilePageTable;
    typedef typename TFilePageTable::TFile              TFile;
    typedef typename TFilePageTable::TPageFrame         TPageFrame;

    typedef Buffer<TValue>                              TBuffer;
    typedef typename Size<TFile>::Type                  TSize;
    typedef	typename Iterator<TBuffer, Standard>::Type	TIterator;

    TFilePageTable  pager;
    TPageFrame      *framePtr;
    int             nextReadyPageNo;
    int             nextFetchPageNo;


    // Iterators to the current position in the page, the end of the page and the end of the page as it is on the disk.
    // itReadEnd is relevant if the file is opened in read+write mode on the last page.  itEnd can point behind the end
    // of the file on the disk whereas itReadEnd points to the end of the file on the disk.  itReadEnd will be
    // incremented if something is written on the last page.

    FileStreamBuffer()
    {
        clear(*this);
    }

    ~FileStreamBuffer()
    {
        close(*this);
    }
    
    inline int _getNextFetchPageNo()
    {
        return nextFetchPageNo++;
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
            resize(oldFrame, this->pptr() - this->pbase());

            // start to release the current buffer
            result &= tryReleaseBuffer(pager, oldFrame, inProgress);
            if (inProgress)
            {
                // writing is in progress
                // -> we put the frame into the out-chain
                #ifdef SEQAN_DEBUG_FILESTREAM
                std::cout << "ADDED INTO OUT-CHAIN (" << pager.outChain.frames << "):" << std::endl;
                std::cout << oldFrame << std::endl;
                #endif
                pager.outChain.pushBack(oldFrame);
                pager.outChainNeedsHousekeeping |= inProgress;
            }
            else
                // writing finished immediately
                // -> we reinject the now unused buffer into the in-chain
                _injectUnusedBuffer(pager, oldFrame);
        }

        // watch current progress
        // initiate next asynchronous task for each buffer
        // reinject written buffers for prefetching the next ones
        _houseKeeping(pager);

        // fetch new buffer
        if (pager.inChain.empty())
        {
            // we need more incoming buffers
            if (pager.inChain.size() + pager.outChain.size() < pager.numPages)
            {
                // we can afford to create new buffers
                _injectUnusedBuffer(pager, *(new TPageFrame()));
            }
            else
            {
                // all out-chain buffers are in progress (otherwise in-chain would not be empty)
                // so we take and wait for the first out-chain buffer
                TPageFrame &anyOutChainBuffer = pager.outChain.popFront();
                tryReleaseBuffer(pager, anyOutChainBuffer, inProgress, true);
                _injectUnusedBuffer(pager, anyOutChainBuffer);
            }
        }

        // we implicitly assume the front buffer
        // to be the next in sequence nextReadyPageNo, nextReadyPageNo+1, ...
        framePtr = &pager.inChain.popFront();
        SEQAN_ASSERT_EQ((int)framePtr->pageNo, nextReadyPageNo);
        ++nextReadyPageNo;

        // now fetch-and-wait for new buffer
        result &= tryFetchBuffer(pager, *framePtr, inProgress, true);

        // framePtr->begin+pageSize always points to the end of the page, whereas framePtr->end can point in the middle.
        this->setg(framePtr->begin, framePtr->begin, framePtr->end);
        this->setp(framePtr->begin, framePtr->begin + pager.pageSize);

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
        resize(oldFrame, this->pptr() - this->pbase());

        // start to release the current buffer
        bool dummy;
        bool result = tryReleaseBuffer(pager, oldFrame, dummy, true);
        freePage(oldFrame, pager.file);
        delete &oldFrame;

        while (!pager.outChain.empty())
        {
            result &= tryReleaseBuffer(pager, pager.outChain.front(), dummy, true);
            freePage(pager.outChain.front(), pager.file);
            delete &pager.outChain.popFront();
        }

        while (!pager.inChain.empty())
        {
            // cancel(inChain.front(), file);
            freePage(pager.inChain.front(), pager.file);
            delete &pager.inChain.popFront();
        }

        TSize _fileSize = pager.fileSize;
        clear(*this);
        if (preserveFileSize)
            pager.fileSize = _fileSize;

        return result;
    }


    virtual TIntValue
    overflow(TIntValue val)
    {
        if (SEQAN_UNLIKELY(!pager.file))
            return TTraits::eof();

        if (SEQAN_UNLIKELY(this->pptr() >= this->epptr() && !_advanceBuffer()))
        {
            this->setp(NULL, NULL);
            this->setg(NULL, NULL, NULL);
            return TTraits::eof();
        }

        if (!TTraits::eq_int_type(val, TTraits::eof()))
        {
            *this->pptr() = TTraits::to_char_type(val);
            this->pbump(1);
        }
        return TTraits::not_eof(val);
    }

    virtual TIntValue
    underflow()
    {
        if (SEQAN_UNLIKELY(!pager.file))
            return TTraits::eof();

        if (SEQAN_UNLIKELY(this->gptr() >= this->egptr()))
        {
            if (SEQAN_UNLIKELY(nextFetchPageNo > pager.lastPageNo || !_advanceBuffer()))
            {
                this->setp(NULL, NULL);
                this->setg(NULL, NULL, NULL);
                return TTraits::eof();
            }
        }

        return TTraits::to_int_type(*this->gptr());
    }

    TValue* eback() const   { return TBase::eback(); }
    TValue* gptr()  const   { return TBase::gptr();  }
    TValue* egptr() const   { return TBase::egptr(); }

    TValue* pbase() const   { return TBase::pbase(); }
    TValue* pptr()  const   { return TBase::pptr();  }
    TValue* epptr() const   { return TBase::epptr(); }

    void setg(TValue* new_eback, TValue* new_gptr, TValue* new_egptr)
    {
        TBase::setg(new_eback, new_gptr, new_egptr);
    }

    void setp(TValue* new_pbase, TValue* new_epptr)
    {
        TBase::setp(new_pbase, new_epptr);
    }

    TPosition _tell()
    {
        if (SEQAN_LIKELY(pager.file))
        {
            if (SEQAN_LIKELY(framePtr != NULL))
                return framePtr->pageNo * pager.pageSize + (gptr() - eback());
            else
                return 0;
        }
        else
        {
            return -1;
        }
    }

    TPosition _seek(TPosition pos)
    {
        if (SEQAN_UNLIKELY(pos < (TPosition)0))
            return -1;

        // Compute the page number and offset of target position.
        unsigned pageNo = pos / pager.pageSize;
        unsigned offset = pos % pager.pageSize;

        if (framePtr == NULL || framePtr->pageNo != pageNo)
        {
            if (framePtr == NULL && pos == 0)
                return 0;

            _stop(true);

            // Fetch the given buffer and wait for this.
            nextFetchPageNo = pageNo;
            nextReadyPageNo = pageNo;

            if (!_advanceBuffer())
            {
                this->setg(NULL, NULL, NULL);
                this->setp(NULL, NULL);
                return -1;
            }
        }

        // Seek to correct position in page.
        this->setg(framePtr->begin, framePtr->begin + offset, framePtr->end);
        this->setp(framePtr->begin, framePtr->begin + pager.pageSize);
        this->pbump(offset);

        return pos;
    }

    virtual TPosition
    seekpos(
        TPosition pos,
        std::ios::openmode which)
    {
        SEQAN_ASSERT_NEQ(which & IosOpenMode<TDirection>::VALUE, 0u);

        if (SEQAN_UNLIKELY(!pager.file))
            return -1;

        return _seek(pos);
    }

    virtual TPosition
    seekoff(
        TDifference off,
        std::ios::seekdir dir,
        std::ios::openmode which)
    {
        SEQAN_ASSERT_NEQ(which & IosOpenMode<TDirection>::VALUE, 0u);

        if (SEQAN_UNLIKELY(!pager.file))
            return -1;

        TPosition pos;

        if (dir == std::ios::beg)
            pos = off;
        else if (dir == std::ios::cur)
            pos = _tell() + off;
        else
            pos = pager.fileSize + off;

        return _seek(pos);
    }
};

template <typename TValue, typename TDirection, typename TSpec>
inline void
clear(FileStreamBuffer<TValue, TDirection, TSpec> &buffer)
{
    buffer.setp(NULL, NULL);
    buffer.setg(NULL, NULL, NULL);

    buffer.framePtr = NULL;
    buffer.nextReadyPageNo = 0;
    buffer.nextFetchPageNo = 0;
}




template <typename TValue, typename TDirection, typename TSpec>
class FileStream;


template <typename TValue, typename TSpec>
struct DefaultOpenMode<FileStream<TValue, Input, TSpec> >
{
    enum { VALUE = OPEN_RDONLY };
};

template <typename TValue, typename TSpec>
struct DefaultOpenMode<FileStream<TValue, Output, TSpec> >
{
    enum { VALUE = OPEN_WRONLY | OPEN_CREATE };
};


template <typename TValue, typename TSpec>
struct DefaultOpenMode<FileStream<TValue, Input, MMap<TSpec> > >
{
    enum { VALUE = OPEN_RDWR | OPEN_APPEND };
};

template <typename TValue, typename TSpec>
struct DefaultOpenMode<FileStream<TValue, Output, MMap<TSpec> > >
{
    enum { VALUE = OPEN_RDWR | OPEN_CREATE };
};




template <typename TValue, typename TDirection, typename TSpec>
class FileStream : public BasicStream<TValue, TDirection>::Type
{
public:
    typedef FileStreamBuffer<TValue, TDirection, TSpec>     TStreamBuffer;
    typedef typename BasicStream<TValue, TDirection>::Type  TBasicStream;

    TStreamBuffer buffer;

    FileStream():
        TBasicStream(&buffer)
    {}

    FileStream(const char *fileName, int openMode = DefaultOpenMode<FileStream>::VALUE):
        TBasicStream(&buffer)
    {
        open(*this, fileName, openMode);
    }

    ~FileStream()
    {
        close(*this);
    }
};

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

template <typename TValue, typename TDirection, typename TSpec>
struct Value<FileStream<TValue, TDirection, TSpec> >:
    Value<FileStreamBuffer<TValue, TDirection, TSpec> > {};

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

template <typename TValue, typename TDirection, typename TSpec>
struct Position<FileStream<TValue, TDirection, TSpec> >:
    Position<FileStreamBuffer<TValue, TDirection, TSpec> > {};

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

template <typename TValue, typename TDirection, typename TSpec>
struct Size<FileStream<TValue, TDirection, TSpec> >:
    Size<FileStreamBuffer<TValue, TDirection, TSpec> > {};

// ----------------------------------------------------------------------------
// Concepts
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
SEQAN_CONCEPT_IMPL((FileStream<TValue, Input, TSpec>), (InputStreamConcept));

template <typename TValue, typename TSpec>
SEQAN_CONCEPT_IMPL((FileStream<TValue, Output, TSpec>), (OutputStreamConcept));

template <typename TValue, typename TSpec>
SEQAN_CONCEPT_IMPL((FileStream<TValue, Bidirectional, TSpec>), (BidirectionalStreamConcept));


// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TFilename, typename TFlags>
inline bool
open(FilePageTable<TValue, TDirection, TSpec> & pager, TFilename const &filename, TFlags const &flags)
{
    clear(pager);
    bool result = open(pager.file, filename, flags);
    pager.fileSize = (result)? size(pager.file) / sizeof(TValue) : 0ul;
    if (pager.fileSize != 0)
        pager.lastPageNo = (pager.fileSize - 1) / pager.pageSize;
    else
        pager.lastPageNo = -1;
    return result;
}

template <typename TValue, typename TDirection, typename TSpec>
inline bool
open(FileStream<TValue, TDirection, TSpec> &stream, const char *fileName, int openMode = DefaultOpenMode<FileStream<TValue, TDirection, TSpec> >::VALUE)
{
    return open(stream.buffer.pager, fileName, openMode);
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
inline int
close(FileStreamBuffer<TValue, TDirection, TSpec> & buffer)
{
    if (buffer.pager.file)
    {
        buffer._stop();
        return close(buffer.pager.file);
    }
    return true;
}

template <typename TValue, typename TDirection, typename TSpec>
inline bool
close(FileStream<TValue, TDirection, TSpec> &stream)
{
    return close(stream.buffer);
}

} // namespace seqan

#endif  // #ifndef SEQAN_FILE_FILE_STREAM_H_
