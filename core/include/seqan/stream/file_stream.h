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
// Functions
// ============================================================================

template <typename TValue, typename TSpec, typename TSize>
inline void
clear(Buffer<TValue, TSpec> &me)
{
    free(me.begin);
    me.begin = me.end = NULL;
    _setCapacity(me, 0);
}

template <typename TValue, typename TSpec, typename TSize>
inline void
reserve(Buffer<TValue, TSpec> &me, TSize newCapacity)
{
    if (newCapacity < capacity(me))
        return;

    free(me.begin);
    me.begin = me.end = (TValue*) valloc(newCapacity * sizeof(TValue));

    if (me.begin == NULL)
    {
        _setCapacity(me, 0);
        throw std::bad_alloc();
    }

    _setCapacity(me, newCapacity);
}

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct IOException:
    public std::runtime_error
{
    IOException(const std::string& message):
        std::runtime_error(message)
    {}
};

// ----------------------------------------------------------------------------
// Class FilePage
// ----------------------------------------------------------------------------

enum PageCompletionStatus
{
    IS_READ = 1,
    IS_PREPROCESSED = 2,
    IS_POSTPROCESSED = 4,
    IS_WRITTEN = 8
};

template <typename TValue, typename TSpec>
struct FilePage
{
    typedef File<TSpec>                         TFile;
    typedef typename AsyncRequest<TFile>::Type  TAsyncFileRequest;
    typedef typename Position<TFile>::Type      TFilePos;
    typedef typename Size<TFile>::Type          TFileSize;

    Buffer<TValue>          rawBuffer;
    Buffer<TValue>          dataBuffer;
    TAsyncFileRequest       request;

    TFilePos                filePos;
    TFileSize               size;
    PageFrameStatus         status;
    PageCompletionStatus    completionStatus;
};

template <typename TValue, typename TConfig>
struct FilePage<TValue, MMap<TConfig> >
{
    typedef File<Async<> >                      TFile;
    typedef typename Position<TFile>::Type      TFilePos;
    typedef typename Size<TFile>::Type          TFileSize;

    Buffer<TValue>          rawBuffer;
    Buffer<TValue>          dataBuffer;

    TFilePos                filePos;
    TFileSize               size;
    PageFrameStatus         status;
    PageCompletionStatus    completionStatus;
};

// ============================================================================
// Functions
// ============================================================================


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FilePager
// ----------------------------------------------------------------------------

template <unsigned PAGE_SIZE = 4 * 1024 * 1024>
struct FixedPagingScheme
{
    enum { pageSize = PAGE_SIZE };
};

template <typename TValue, typename TDirection, typename TSpec, typename TPagingScheme>
struct Host<FilePageTable<TValue, TDirection, TSpec, TPagingScheme> >
{
    typedef File<TSpec> Type;
};

template <typename TValue, typename TDirection, typename TSpec, typename TPagingScheme>
struct Host<FilePageTable<TValue, TDirection, MMap<TSpec>, TPagingScheme> >
{
    typedef FileMapping<TSpec> Type;
};

template <typename TValue, typename TDirection, typename TSpec = Async<>, typename TPagingScheme = FixedPagingScheme<> >
struct FilePageTable
{
    typedef typename Host<FilePageTable>::Type          TFile;
    typedef typename Size<TFile>::Type                  TSize;

    typedef FilePage<TValue, TSpec>                     TPageFrame;
    typedef PageChain<TPageFrame>                       TPageChain;
    typedef Buffer<TValue>                              TBuffer;
    typedef	typename Iterator<TBuffer, Standard>::Type	TIterator;

    typedef short int                                   TPageId;
    typedef String<TPageId>                             TPageDir;

//    static const unsigned pageSize = 4*1024;      // 1M entries per page
    static const unsigned numPages = 2;                // double-buffering
    static const unsigned pageSize = 4*1024*1024;      // 1M entries per page

    TFile           file;
    TSize           fileSize;
    unsigned        filePages;

    TPageChain      inChain;
    TPageChain      outChain;

    bool            inChainNeedsHousekeeping;
    bool            outChainNeedsHousekeeping;

    TPagingScheme   scheme;

    inline void _printStats()
    {
        std::cout << "in: " << inChain.size() << "\tout:" << outChain.size() << std::endl;
    }
};


// ============================================================================
// Functions
// ============================================================================

template <typename TValue, typename TDirection, typename TSpec>
inline void
clear(FilePageTable<TValue, TDirection, TSpec> &pager)
{
    pager.filePages = 0;
    pager.fileSize = 0;
    pager.inChainNeedsHousekeeping = false;
    pager.outChainNeedsHousekeeping = false;
}

template <typename TValue, typename TDirection, typename TSpec, typename TPageNo, unsigned PAGE_SIZE>
inline Pair<__int64, unsigned>
_getPageOffsetAndLength(FilePageTable<TValue, TDirection, TSpec, FixedPagingScheme<PAGE_SIZE> > &pager, TPageNo pageNo)
{
    return Pair<__int64, unsigned>((__int64)pager.scheme.pageSize * (__int64)pageNo, pager.scheme.pageSize);
}



// Variant for POSIX file access
template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TOffset, typename TPageSize>
inline void
_readFilePage(FilePageTable<TValue, TDirection, TSpec> &pager, TPageFrame &pf, TOffset offset, TPageSize pageSize)
{
    typedef typename FilePageTable<TValue, TDirection, TSpec>::TSize TSize;

    // allocate required memory
    reserve(pf, pageSize);

    // do nothing in output-only mode or when there is nothing to read
    if (IsSameType<TDirection, Output>::VALUE || (TSize)offset >= pager.fileSize)
    {
        // no valid data read and we return immediately
        resize(pf, 0);
        pf.status = READY;
    }

    // shrink read buffer size at the end of file
    resize(pf, std::min((TSize)pageSize, pager.fileSize - (TSize)offset));

    // start asynchronous reading
	pf.filePos = offset;
    pf.status = READING;
    bool success = asyncReadAt(pager.file, pf.begin, length(pf), (TSize)offset, pf.request);

    // if an error occurred, throw an I/O exception
    if (!success)
        throw IOException((std::string)_pageFrameStatusString(pf) + " operation could not be initiated: \"" + strerror(errno) + '"');
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame, typename TOffset, typename TPageSize>
inline void
_readFilePage(FilePageTable<TValue, TDirection, MMap<TConfig> > &pager, TPageFrame &pf, TOffset offset, TPageSize pageSize)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

	TSize endOfs = (TSize)offset + (TSize)pageSize;
	TSize dataSize;

    // compute valid readable data
	if ((TSize)offset >= length(pager.file))
		dataSize = 0;	// no valid data to read from behind the end of file
	else
		dataSize = std::min((TSize)pageSize, length(pager.file) - (TSize)offset);
	
    // how to handle pages crossing/beyond the end of file
    if (endOfs > length(pager.file))
    {
        if (!IsSameType<TDirection, Input>::VALUE)
        {
            // increase file size to next page boundary and map the whole page
			resize(pager.file, endOfs);
        }
        else
        {
            // adapt page size to file size and map only the contents of the file
			pageSize = dataSize;
        }
    }

	pf.filePos = offset;
    pf.begin = mapFileSegment(pager.file, offset, pageSize);
	resize(pf, dataSize);
	pf.status = READY;
}
#endif


template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline void
_preprocessBuffer(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = READY;
}

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline void
_postprocessBuffer(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &pf, bool)
{
    // not used yet, could be usefull for external sorters/mappers
    pf.status = POSTPROCESSED;
}


template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline void
_writeBuffer(FilePageTable<TValue, TDirection, TSpec> &pager, TPageFrame &pf)
{
    // only write in write-mode
    if (IsSameType<TDirection, Input>::VALUE)
    {
        pf.status = UNUSED;
        return;
    }

    // start asynchronous writing
    pf.status = WRITING;
    bool success = asyncWriteAt(pager.file, pf.begin, length(pf), pf.offset, pf.request);

    // if an error occurred, throw an I/O exception
    if (!success)
        throw IOException((std::string)_pageFrameStatusString(pf) + " operation could not be initiated: \"" + strerror(errno) + '"');
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_writeBuffer(FilePageTable<TValue, TDirection, MMap<TConfig> > &pager, TPageFrame &pf)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

    pf.status = UNUSED;
    unmapFileSegment(pager.file, pf.begin, capacity(pf));

    if (!IsSameType<TDirection, Input>::VALUE)
    {
        pager.fileSize = std::max(pager.fileSize, (TSize)pf.filePos + (TSize)length(pf));

        // TODO: shrink file before closing
        //resize(pager.file, pager.fileSize);
    }
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
        _readBuffer(pager, pf);
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
        _preprocessBuffer(pager, pf, doWait);
    }

    if (pf.status == PREPROCESSING)
    {
        // we get here only if doWait==false and the
        // asynchronous postprocessing is in progress
        SEQAN_ASSERT_NOT(doWait);
        return;
    }

    if (pf.status == READY)
    {
        // all done
        inProgress = false;
        return;
    }

    SEQAN_FAIL("PageFrame has inconsistent state.");
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
        SEQAN_ASSERT_NEQ(which & IosOpenMode<TDirection>::VALUE, 0);

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
        SEQAN_ASSERT_NEQ(which & IosOpenMode<TDirection>::VALUE, 0);

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
