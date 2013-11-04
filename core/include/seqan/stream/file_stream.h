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

template <typename TValue, typename TSpec>
inline void
free(Buffer<TValue, TSpec> & me)
{
    ::free(me.begin);
}

template <typename TValue, typename TSpec>
inline void
clear(Buffer<TValue, TSpec> & me)
{
    me.begin = me.end = NULL;
    _setCapacity(me, 0);
}

template <typename TValue, typename TSpec, typename TSize>
inline void
reserve(Buffer<TValue, TSpec> & me, TSize newCapacity)
{
    if (newCapacity <= capacity(me))
        return;

    free(me);
    me.begin = me.end = (TValue *) valloc(newCapacity * sizeof(TValue));

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

struct IOException :
    public std::runtime_error
{
    IOException(const std::string & message) :
        std::runtime_error(message)
    {}
};

// ----------------------------------------------------------------------------
// Class FilePage
// ----------------------------------------------------------------------------

enum PageCompletionState
{
    IS_READ = 1,            // raw buffer is in sync with disk
    IS_PREPROCESSED = 2,    // data buffer is in sync (e.g. decompressed) with raw
    IS_POSTPROCESSED = 4,   //
    IS_WRITTEN = 8
};

template <typename TValue, typename TSpec>
struct FilePage
{
    typedef File<TSpec>                         TFile;
    typedef typename AsyncRequest<TFile>::Type  TAsyncFileRequest;
    typedef typename Position<TFile>::Type      TFilePos;
    typedef typename Size<TFile>::Type          TFileSize;

    FilePage                * next;
    unsigned                lockCount;

    Buffer<TValue>          raw;
    Buffer<TValue>          data;
    TAsyncFileRequest       request;

    TFilePos                filePos;
    TFileSize               size;
    PageFrameState          state;
    PageFrameState          targetState;
    PageCompletionState     completionState;
};

template <typename TValue, typename TConfig>
struct FilePage<TValue, MMap<TConfig> >
{
    typedef File<Async<> >                      TFile;
    typedef typename Position<TFile>::Type      TFilePos;
    typedef typename Size<TFile>::Type          TFileSize;

    FilePage                * next;
    unsigned                lockCount;

    Buffer<TValue>          raw;
    Buffer<TValue>          data;
    AsyncDummyRequest       request;

    TFilePos                filePos;
    TFileSize               size;
    PageFrameState          state;
    PageCompletionState     completionState;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TValue, typename TSpec>
inline void
clear(FilePage<TValue, TSpec> & me)
{
    free(me.raw);
    free(me.data);
    clear(me.raw);
    clear(me.data);
}

template <typename TValue, typename TConfig>
inline void
clear(FilePage<TValue, MMap<TConfig> > & me)
{
    free(me.data);
    clear(me.raw);
    clear(me.data);
}

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

    static const void * EMPTY;
    static const void * ON_DISK;

    String<void *> frameStart;
};

template <unsigned PAGE_SIZE>
const void * FixedPagingScheme<PAGE_SIZE>::EMPTY = NULL;

template <unsigned PAGE_SIZE>
const void * FixedPagingScheme<PAGE_SIZE>::ON_DISK = (void *)-1;

struct RandomPagingScheme
{
    static const void * EMPTY;
    static const void * ON_DISK;

    std::map<__uint64, void *> frameStart;
};

const void * RandomPagingScheme::EMPTY = NULL;
const void * RandomPagingScheme::ON_DISK = (void *)-1;



template <typename TValue, typename TDirection, typename TSpec = Async<>, typename TPagingScheme = FixedPagingScheme<> >
struct FilePageTable;



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


// incoming chain
// cached chain
// outgoing chain

template <typename TValue, typename TDirection, typename TSpec, typename TPagingScheme>
struct FilePageTable
{
    typedef typename Host<FilePageTable>::Type          TFile;
    typedef typename Size<TFile>::Type                  TSize;

    typedef FilePage<TValue, TSpec>                     TPageFrame;
    typedef PageChain<TPageFrame>                       TPageChain;
    typedef Buffer<TValue>                              TBuffer;
    typedef typename Iterator<TBuffer, Standard>::Type  TIterator;

    typedef short int                                   TPageId;
    typedef String<TPageId>                             TPageDir;

    static const unsigned numPages = 2;                // double-buffering

    TFile           file;
    TSize           fileSize;

    TPageChain      unused;
    TPageChain      ready;
    TPageChain      inProcess;

    TPagingScheme   table;

    inline void _printStats()
    {
        std::cout << "unused: " << unused.size() << "\tready:" << ready.size() << "\tinProcess:" << inProcess.size() << std::endl;
    }

};


// ============================================================================
// Functions
// ============================================================================

template <typename TValue, typename TDirection, typename TSpec>
inline void
clear(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    pager.fileSize = 0;
}

template <typename TValue, typename TDirection, typename TSpec, unsigned PAGE_SIZE, typename TPos>
inline Pair<__int64, unsigned>
_getPageOffsetAndLength(FilePageTable<TValue, TDirection, TSpec, FixedPagingScheme<PAGE_SIZE> > & pager, TPos pos)
{
    return Pair<__int64, unsigned>((__int64)pos & (__int64)(pager.table.pageSize - 1), pager.table.pageSize);
}

// Variant for POSIX file access
template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_readFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TPageFrame & pf)
{
    // allocate required memory
    reserve(pf.raw, pf.size);

    // do nothing in output-only mode or when there is nothing to read
    if (IsSameType<TDirection, Output>::VALUE || pf.filePos >= pager.fileSize)
    {
        // no valid data read and we return immediately
        resize(pf.raw, 0);
        return true;    // true = reading completed
    }

    // shrink read buffer size at the end of file
    resize(pf.raw, std::min(pf.size, pager.fileSize - pf.filePos));

    // start asynchronous reading
    bool success = asyncReadAt(pager.file, pf.begin, length(pf.raw), pf.filePos, pf.request);

    // if an error occurred, throw an I/O exception
    if (!success)
        throw IOException((std::string)_pageFrameStateString(pf) + " operation could not be initiated: \"" + strerror(errno) + '"');

    return false;   // false = reading in process
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_readFilePage(FilePageTable<TValue, TDirection, MMap<TConfig> > & pager, TPageFrame & pf)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

    TSize endOfs = pf.filePos + pf.size;
    TSize dataSize;

    // compute valid readable data
    if (pf.filePos >= length(pager.file))
        dataSize = 0;   // no valid data to read from behind the end of file
    else
        dataSize = std::min(pf.size, length(pager.file) - pf.filePos);

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
            pf.size = dataSize;
        }
    }

    pf.raw.begin = (TValue *) mapFileSegment(pager.file, pf.filePos, pf.size);
    _setCapacity(pf.raw, pf.size);
    resize(pf.raw, dataSize);
    return true;    // true = reading completed
}

#endif


template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_preprocessFilePage(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &, TBool const &)
{
    // not used yet, could be usefull for external sorters/mappers
    return true;
}

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_postprocessFilePage(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &, TBool const &)
{
    // not used yet, could be usefull for external sorters/mappers
    return true;
}

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_writeFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TPageFrame & pf)
{
    // only write in write-mode
    if (IsSameType<TDirection, Input>::VALUE)
    {
        pf.state = UNUSED;
        return true;    // true = writing completed
    }

    // start asynchronous writing
    pf.state = WRITING;
    bool success = asyncWriteAt(pager.file, pf.raw.begin, length(pf.raw), pf.filePos, pf.request);

    // if an error occurred, throw an I/O exception
    if (!success)
        throw IOException((std::string)_pageFrameStateString(pf) + " operation could not be initiated: \"" + strerror(errno) + '"');
    return false;   // false = writing in process
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_writeFilePage(FilePageTable<TValue, TDirection, MMap<TConfig> > & pager, TPageFrame & pf)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

    pf.state = UNUSED;
    unmapFileSegment(pager.file, pf.raw.begin, pf.size);

    if (!IsSameType<TDirection, Input>::VALUE)
    {
        pager.fileSize = std::max(pager.fileSize, (TSize)pf.filePos + (TSize)length(pf.raw));

        // TODO: shrink file before closing
        //resize(pager.file, pager.fileSize);
    }
    return true;    // true = writing completed
}

#endif

// ----------------------------------------------------------------------------
// Function _tryFetchFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_processFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TPageFrame & pf, TBool const & doWait)
{
    bool inProgress;
    while (pf.state != pf.targetState)
    {
        switch (pf.state)
        {
        case UNUSED:
            // page is ready, start reading
            pf.state = (_readFilePage(pager, pf)) ? READING_DONE : READING;
            break;

        case READING:
            inProgress = false;
            if (doWait)
                waitFor(pf.request);
            else
                waitFor(pf.request, 0, inProgress);

            // still in operation?
            if (inProgress)
                return false;

            // was reading, now preprocessing
            pf.state = READING_DONE;
            break;

        case READING_DONE:
            pf.state = (_preprocessFilePage(pager, pf, doWait)) ? PREPROCESSING_DONE : PREPROCESSING;
            break;

        case PREPROCESSING:
            // TODO: fill me

            pf.state = PREPROCESSING_DONE;
            break;

        case PREPROCESSING_DONE:
            // TODO: put me into ready chain
            pf.state = READY;
            break;

        case READY:
            // all done
            pf.state = (_postprocessFilePage(pager, pf, doWait)) ? POSTPROCESSING_DONE : POSTPROCESSING;
            break;

        case POSTPROCESSING:
            // TODO: fill me

            pf.state = PREPROCESSING_DONE;
            break;

        case PREPROCESSING_DONE:
            // page was postprocessed, now write asynchronously
            pf.state = (_writeFilePage(pager, pf)) ? WRITING_DONE : WRITING;
            break;

        case WRITING:
            inProgress = false;
            if (doWait)
                waitFor(pf.request);
            else
                waitFor(pf.request, 0, inProgress);

            // still in operation?
            if (inProgress)
                return false;

            pf.state = WRITING_DONE;
            break;

        case WRITING_DONE:
            // TODO: put me into ready chain
            pf.state = UNUSED;
            break;
        }
    }
    return true;
}

template <typename TFilePageTable, typename TPageFrame>
inline bool
_injectUnusedFilePage(TFilePageTable & pager, TPageFrame & frame)
{
    frame.state = UNUSED;
    pager.incomingNeedsHousekeeping |= !_tryFetchFilePage(pager, frame, False());
    pager.incoming.pushBack(frame);
    return true;
}

template <typename TValue, typename TDirection, typename TSpec>
inline bool
_houseKeeping(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    typedef FilePageTable<TValue, TDirection, TSpec>    TFilePageTable;
    typedef typename TFilePageTable::TPageFrame         TPageFrame;

    if (pager.inProcess.size() != 0)
    {
        for (TPageFrame * p = pager.inProcess.first; p != NULL; )
        {
            TPageFrame & page = *p;
            p = p->next;

            if (_processFilePage(pager, page, False()))
            {
                if (page.state == READY)
                {
                    pager.inProcess.erase(page);
                    pager.ready.pushBack(page);
                }
                else if (page.state == UNUSED)
                {
                    pager.inProcess.erase(page);
                    pager.unused.pushBack(page);
                }
            }
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function _newFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
inline FilePage<TValue, TSpec> *
_newFilePage(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    typedef FilePage<TValue, TSpec> TPage;
    if (unused.size())
    // fetch new page
    if (unused.size() < pager.numPages)
    {
        // we can afford to create new buffers
        return new TPage();
    }
    else
    {
        // all out-chain buffers are in progress (otherwise in-chain would not be empty)
        // so we take and wait for the first out-chain buffer
        TPage & anyOutChainFilePage = pager.outgoing.popFront();
        _tryReleaseFilePage(pager, anyOutChainFilePage, True());
        return &anyOutChainFilePage;
    }
}

// ----------------------------------------------------------------------------
// Function _lockFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TFilePos, typename TSize>
inline FilePage<TValue, TSpec> &
_lockFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TFilePos filePos, TSize size)
{
    typedef FilePage<TValue, TSpec> TPage;
    TPage * p;

    SEQAN_OMP_PRAGMA(critical(lockPageTable))
    {
        unsigned pageNo = filePos / pager.table.pageSize;
        if (length(pager.table.frameStart) <= pageNo)
            resize(pager.table.frameStart, pageNo + 1, pager.table.EMPTY);

        TPage * p = static_cast<TPage *>(pager.table.frameStart[pageNo]);

        if (p == pager.table.EMPTY || p == pager.table.ON_DISK)
        {
            p = _newFilePage(pager);
            p->filePos = filePos;
            p->size = size;
            pager.table[pageNo] = p;
        }
        ++p->count;
    }
    return *p;
}

// ----------------------------------------------------------------------------
// Function fetchFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TFilePos, typename TSize>
inline FilePage<TValue, TSpec> &
fetchFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TFilePos filePos, TSize size, bool readOnly = false)
{
    typedef FilePage<TValue, TSpec> TPage;

    TPage & page = _lockFilePage(pager, filePos, size);
    _tryFetchFilePage(pager, page, True());
}

template <typename TValue, typename TDirection, typename TSpec, typename TFilePage>
inline void
swapOutFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TFilePage & pf)
{
    SEQAN_OMP_PRAGMA(critical(swapOutFilePage))
    {
        if (pf.lockCount == 0)
        {
            // start to release the current buffer
            if (_tryReleaseFilePage(pager, pf, False()))
            {
                // writing finished immediately
                // -> we reinject the now unused buffer into the in-chain
                _injectUnusedFilePage(pager, pf);
            }
            else
            {
                // writing is in progress
                // -> we put the frame into the out-chain
                #ifdef SEQAN_DEBUG_FILESTREAM
                std::cout << "ADDED INTO OUT-CHAIN (" << pager.outgoing.frames << "):" << std::endl;
                std::cout << oldFrame << std::endl;
                #endif
                pager.outgoing.pushBack(pf);
                pager.outgoingNeedsHousekeeping = true;
            }
        }
    }
}

template <typename TValue, typename TDirection, typename TSpec, typename TFilePage>
inline void
releaseFilePage(FilePageTable<TValue, TDirection, TSpec> &, TFilePage & pf)
{
    atomicDec(pf.lockCount);
}

// ----------------------------------------------------------------------------
// Class FileStreamBuffer
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec = Async<> >
struct FileStreamBuffer :
    public std::basic_streambuf<TValue>
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
    typedef typename Iterator<TBuffer, Standard>::Type  TIterator;

    TFilePageTable  pager;
    TPageFrame      * framePtr;
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

    inline void _advanceFilePage()
    {
        // release old buffer
        if (framePtr != NULL)
        {
            // shrink current buffer size if it was not fully written (like the last buffer)
            TPageFrame & oldFrame = *framePtr;
            resize(oldFrame, this->pptr() - this->pbase());

            swapOutFilePage(pager, oldFrame);
        }

        // watch current progress
        // initiate next asynchronous task for each buffer
        // reinject written buffers for prefetching the next ones
        _houseKeeping(pager);

        // we implicitly assume the front buffer
        // to be the next in sequence nextReadyPageNo, nextReadyPageNo+1, ...
        framePtr = &pager.incoming.popFront();

        // now fetch-and-wait for new buffer
        _tryFetchFilePage(pager, *framePtr, True());

        // framePtr->begin+pageSize always points to the end of the page, whereas framePtr->end can point in the middle.
        this->setg(framePtr->begin, framePtr->begin, framePtr->end);
        this->setp(framePtr->begin, framePtr->begin + pager.pageSize);
    }

    // Flush everything to disk.  If preserveFileSize then the file size is preserved.
    inline void _stop(bool preserveFileSize = false)
    {
        // has anything been written?
        if (nextReadyPageNo == 0)
            return;

        // shrink current buffer size if it was not fully written (like the last buffer)
        TPageFrame & oldFrame = *framePtr;
        resize(oldFrame.raw, this->pptr() - this->pbase());

        // start to release the current buffer
        _tryReleaseFilePage(pager, oldFrame, True());
        clear(oldFrame.raw);
        delete &oldFrame;

        while (!pager.outgoing.empty())
        {
            _tryReleaseFilePage(pager, pager.outgoing.front(), True());
            clear(pager.outgoing.front());
            delete &pager.outgoing.popFront();
        }

        while (!pager.incoming.empty())
        {
            // cancel(incoming.front(), file);
            clear(pager.incoming.front());
            delete &pager.incoming.popFront();
        }

        TSize _fileSize = pager.fileSize;
        clear(*this);
        if (preserveFileSize)
            pager.fileSize = _fileSize;
    }

    virtual TIntValue
    overflow(TIntValue val)
    {
        if (SEQAN_UNLIKELY(!pager.file))
            return TTraits::eof();

        if (SEQAN_UNLIKELY(this->pptr() >= this->epptr() && !_advanceFilePage()))
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
            if (SEQAN_UNLIKELY(nextFetchPageNo > pager.lastPageNo || !_advanceFilePage()))
            {
                this->setp(NULL, NULL);
                this->setg(NULL, NULL, NULL);
                return TTraits::eof();
            }
        }

        return TTraits::to_int_type(*this->gptr());
    }

    TValue * eback() const   { return TBase::eback(); }
    TValue * gptr()  const   { return TBase::gptr();  }
    TValue * egptr() const   { return TBase::egptr(); }

    TValue * pbase() const   { return TBase::pbase(); }
    TValue * pptr()  const   { return TBase::pptr();  }
    TValue * epptr() const   { return TBase::epptr(); }

    void setg(TValue * new_eback, TValue * new_gptr, TValue * new_egptr)
    {
        TBase::setg(new_eback, new_gptr, new_egptr);
    }

    void setp(TValue * new_pbase, TValue * new_epptr)
    {
        TBase::setp(new_pbase, new_epptr);
    }

    TPosition _tell()
    {
        if (SEQAN_LIKELY(pager.file))
        {
            if (SEQAN_LIKELY(framePtr != NULL))
                return framePtr->filePos + (gptr() - eback());
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
        if (framePtr == NULL || pos < framePtr->filePos || framePtr->filePos + framePtr->size <= pos)
        {
            // Release former page.
            if (framePtr != NULL)
                releaseFilePage(pager, *framePtr);

            // Fetch new page.
            Pair<__int64, unsigned> ol = _getPageOffsetAndLength(pager, pos);
            framePtr = fetchFilePage(pager, ol.i1, ol.i2);
        }

        // Seek to correct position in page.
        this->setg(framePtr->data.begin, framePtr->data.begin + (pos - framePtr->filePos), framePtr->data.end);
        this->setp(framePtr->data.begin, framePtr->data.begin + capacity(*framePtr));
        this->pbump((pos - framePtr->filePos));

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
clear(FileStreamBuffer<TValue, TDirection, TSpec> & buffer)
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
class FileStream :
    public BasicStream<TValue, TDirection>::Type
{
public:
    typedef FileStreamBuffer<TValue, TDirection, TSpec>     TStreamBuffer;
    typedef typename BasicStream<TValue, TDirection>::Type  TBasicStream;

    TStreamBuffer buffer;

    FileStream() :
        TBasicStream(&buffer)
    {}

    FileStream(const char * fileName, int openMode = DefaultOpenMode<FileStream>::VALUE) :
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
    Value<FileStreamBuffer<TValue, TDirection, TSpec> >{};

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
    Position<FileStreamBuffer<TValue, TDirection, TSpec> >{};

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
    Size<FileStreamBuffer<TValue, TDirection, TSpec> >{};

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
open(FilePageTable<TValue, TDirection, TSpec> & pager, TFilename const & filename, TFlags const & flags)
{
    clear(pager);
    bool result = open(pager.file, filename, flags);
    pager.fileSize = (result) ? length(pager.file) / sizeof(TValue) : 0ul;
    return result;
}

template <typename TValue, typename TDirection, typename TSpec>
inline bool
open(FileStream<TValue, TDirection, TSpec> & stream, const char * fileName, int openMode = DefaultOpenMode<FileStream<TValue, TDirection, TSpec> >::VALUE)
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
close(FileStream<TValue, TDirection, TSpec> & stream)
{
    return close(stream.buffer);
}

} // namespace seqan

#endif  // #ifndef SEQAN_FILE_FILE_STREAM_H_
