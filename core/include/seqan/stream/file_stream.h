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

// TODO(weese):
// * skip reading of a page frame when ptr is EMPTY and not ONDISK
// * skip writing when page was not modified
// * cancel writing and reading if page is locked while writing


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
    typedef typename Size<Buffer<TValue, TSpec> >::Type TBufferSize;
    if ((TBufferSize)newCapacity <= capacity(me))
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
    typedef File<Async<> >                          TFile;
    typedef typename AsyncRequest<FilePage>::Type   TAsyncRequest;
    typedef typename Position<TFile>::Type          TFilePos;
    typedef typename Size<TFile>::Type              TFileSize;

    FilePage                * next;
    unsigned                lockCount;

    Buffer<TValue>          raw;
    Buffer<TValue>          data;
    TAsyncRequest           request;

    TFilePos                filePos;
    TFileSize               size;
    PageFrameState          state;
    PageFrameState          targetState;
    PageCompletionState     completionState;
};

template <typename TValue, typename TSpec>
struct AsyncRequest<FilePage<TValue, TSpec> >:
    AsyncRequest<File<TSpec> > {};

template <typename TValue, typename TConfig>
struct AsyncRequest<FilePage<TValue, MMap<TConfig> > >
{
    typedef AsyncDummyRequest Type;
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

    static void * EMPTY;
    static void * ON_DISK;

    String<void *> frameStart;
};

template <unsigned PAGE_SIZE>
void * FixedPagingScheme<PAGE_SIZE>::EMPTY = NULL;

template <unsigned PAGE_SIZE>
void * FixedPagingScheme<PAGE_SIZE>::ON_DISK = (void *)-1;

struct RandomPagingScheme
{
    static void * EMPTY;
    static void * ON_DISK;

    std::map<__uint64, void *> frameStart;
};

void * RandomPagingScheme::EMPTY = NULL;
void * RandomPagingScheme::ON_DISK = (void *)-1;



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

    ~FilePageTable()
    {
        flushAndFree(*this);
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

// ----------------------------------------------------------------------------
// Function _readFilePage()
// ----------------------------------------------------------------------------

// Variant for POSIX file access
template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_readFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TPageFrame & page)
{
    // allocate required memory
    reserve(page.raw, page.size);

    // do nothing in output-only mode or when there is nothing to read
    if (IsSameType<TDirection, Output>::VALUE || page.filePos >= pager.fileSize)
    {
        // no valid data read and we return immediately
        resize(page.raw, 0);
        return true;    // true = reading completed
    }

    // shrink read buffer size at the end of file
    resize(page.raw, std::min(page.size, pager.fileSize - page.filePos));

    // start asynchronous reading
    bool success = asyncReadAt(pager.file, begin(page.raw, Standard()), length(page.raw), page.filePos, page.request);

    // if an error occurred, throw an I/O exception
    if (!success)
        throw IOException((std::string)_pageFrameStatusString(page.state) + " operation could not be initiated: \"" + strerror(errno) + '"');

    return false;   // false = reading in process
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_readFilePage(FilePageTable<TValue, TDirection, MMap<TConfig> > & pager, TPageFrame & page)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

    TSize endOfs = page.filePos + page.size;
    TSize dataSize;

    // compute valid readable data
    if (page.filePos >= length(pager.file))
        dataSize = 0;   // no valid data to read from behind the end of file
    else
        dataSize = std::min(page.size, length(pager.file) - page.filePos);

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
            page.size = dataSize;
        }
    }

    page.raw.begin = (TValue *) mapFileSegment(pager.file, page.filePos, page.size);
    _setCapacity(page.raw, page.size);
    resize(page.raw, dataSize);
    return true;    // true = reading completed
}

#endif

// ----------------------------------------------------------------------------
// Function _preprocessFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_preprocessFilePage(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &page, TBool const &)
{
    // not used yet, could be usefull for external sorters/mappers
    page.data = page.raw;
    return true;
}

// ----------------------------------------------------------------------------
// Function _postprocessFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_postprocessFilePage(FilePageTable<TValue, TDirection, TSpec> &, TPageFrame &page, TBool const &)
{
    // not used yet, could be usefull for external sorters/mappers
    clear(page.data);
    return true;
}

// ----------------------------------------------------------------------------
// Function _writeFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame>
inline bool
_writeFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TPageFrame & page)
{
    // only write in write-mode
    if (IsSameType<TDirection, Input>::VALUE)
    {
        page.state = UNUSED;
        return true;    // true = writing completed
    }

    // start asynchronous writing
    page.state = WRITING;
    bool success = asyncWriteAt(pager.file, begin(page.raw, Standard()), length(page.raw), page.filePos, page.request);

    // if an error occurred, throw an I/O exception
    if (!success)
        throw IOException((std::string)_pageFrameStatusString(page.state) + " operation could not be initiated: \"" + strerror(errno) + '"');
    return false;   // false = writing in process
}

#ifndef PLATFORM_WINDOWS
template <typename TValue, typename TDirection, typename TConfig, typename TPageFrame>
inline bool
_writeFilePage(FilePageTable<TValue, TDirection, MMap<TConfig> > & pager, TPageFrame & page)
{
    typedef FilePageTable<TValue, TDirection, MMap<TConfig> >   TFilePageTable;
    typedef typename TFilePageTable::TFile                      TFile;
    typedef typename Size<TFile>::Type                          TSize;

    page.state = UNUSED;
    unmapFileSegment(pager.file, begin(page.raw, Standard()), page.size);

    if (!IsSameType<TDirection, Input>::VALUE)
    {
        pager.fileSize = std::max(pager.fileSize, (TSize)page.filePos + (TSize)length(page.raw));

        // TODO: shrink file before closing
        //resize(pager.file, pager.fileSize);
    }
    return true;    // true = writing completed
}

#endif

// ----------------------------------------------------------------------------
// Function _processFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_processFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TPageFrame & page, TBool const & doWait)
{
    bool inProgress;
    while (page.state != page.targetState)
    {
        switch (page.state)
        {
        case UNUSED:
            // page is ready, start reading
            page.state = (_readFilePage(pager, page)) ? READING_DONE : READING;
            break;

        case READING:
            inProgress = false;
            if (doWait)
                waitFor(page.request);
            else
                waitFor(page.request, 0, inProgress);

            // still in operation?
            if (inProgress)
                return false;

            // was reading, now preprocessing
            page.state = READING_DONE;
            break;

        case READING_DONE:
            page.state = (_preprocessFilePage(pager, page, doWait)) ? PREPROCESSING_DONE : PREPROCESSING;
            break;

        case PREPROCESSING:
            // TODO: fill me

            page.state = PREPROCESSING_DONE;
            break;

        case PREPROCESSING_DONE:
            // TODO: put me into ready chain
            page.state = READY;
            break;

        case READY:
            // all done
            page.state = (_postprocessFilePage(pager, page, doWait)) ? POSTPROCESSING_DONE : POSTPROCESSING;
            break;

        case POSTPROCESSING:
            // TODO: fill me

            page.state = PREPROCESSING_DONE;
            break;

        case POSTPROCESSING_DONE:
            // page was postprocessed, now write asynchronously
            page.state = (_writeFilePage(pager, page)) ? WRITING_DONE : WRITING;
            break;

        case WRITING:
            inProgress = false;
            if (doWait)
                waitFor(page.request);
            else
                waitFor(page.request, 0, inProgress);

            // still in operation?
            if (inProgress)
                return false;

            page.state = WRITING_DONE;
            break;

        case WRITING_DONE:
            // TODO: put me into ready chain
            page.state = UNUSED;
            break;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function _houseKeeping()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
inline bool
_houseKeeping(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    typedef FilePageTable<TValue, TDirection, TSpec>    TFilePageTable;
    typedef typename TFilePageTable::TPageFrame         TPageFrame;

    if (!empty(pager.inProcess))
    {
        for (TPageFrame * p = pager.inProcess.first; p != NULL; )
        {
            TPageFrame & page = *p;
            p = p->next;

            if (_processFilePage(pager, page, False()))
            {
                if (page.state == READY)
                {
                    erase(pager.inProcess, page);
                    pushBack(pager.ready, page);
                }
                else if (page.state == UNUSED)
                {
                    erase(pager.inProcess, page);
                    pushBack(pager.unused, page);
                }
            }
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function _swapOutFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
inline FilePage<TValue, TSpec> *
_swapOutFilePage(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    SEQAN_ASSERT_NOT(empty(pager.ready) && empty(pager.inProcess));

    typedef FilePage<TValue, TSpec> TPage;

    if (!empty(pager.ready))
        for (TPage * p = pager.ready.first; p != NULL; p = p->next)
            if (p->lockCount == 0)
            {
                p->targetState = UNUSED;
                _processFilePage(pager, *p, True());
                erase(pager.ready, *p);
                return p;
            }

    if (!empty(pager.inProcess))
        for (TPage * p = pager.inProcess.first; p != NULL; p = p->next)
            if (p->lockCount == 0)
            {
                p->targetState = UNUSED;
                _processFilePage(pager, *p, True());
                erase(pager.inProcess, *p);
                return p;
            }
    return NULL;
}

// ----------------------------------------------------------------------------
// Function _flush()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
inline void
flush(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    typedef FilePage<TValue, TSpec> TPage;

    while (!empty(pager.ready))
    {
        TPage &page = popFront(pager.ready);
        page.targetState = UNUSED;
        _processFilePage(pager, page, True());
        pushBack(pager.unused, page);
    }

    while (!empty(pager.inProcess))
    {
        TPage &page = popFront(pager.inProcess);
        page.targetState = UNUSED;
        _processFilePage(pager, page, True());
        pushBack(pager.unused, page);
    }
}

template <typename TValue, typename TDirection, typename TSpec, typename TFilePage>
inline void
_freeFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TFilePage & page)
{
    pager.table.frameStart[page.filePos / pager.table.pageSize] = pager.table.ON_DISK;
    clear(page);
}

template <typename TValue, typename TDirection, typename TSpec>
inline void
flushAndFree(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    flush(pager);
    while (!empty(pager.unused))
        _freeFilePage(pager, popFront(pager.unused));
}

// ----------------------------------------------------------------------------
// Function _newFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec>
inline FilePage<TValue, TSpec> *
_newFilePage(FilePageTable<TValue, TDirection, TSpec> & pager)
{
    typedef FilePage<TValue, TSpec> TPage;

    if (empty(pager.unused))
        _houseKeeping(pager);

    if (!empty(pager.unused))
        return &popFront(pager.unused);

    // fetch new page
    if (length(pager.ready) + length(pager.inProcess) < pager.numPages)
    {
        // we can afford to create new buffers
        return new TPage();
    }
    else
    {
        // all out-chain buffers are in progress (otherwise in-chain would not be empty)
        // so we take and wait for the first out-chain buffer
        if (!empty(pager.inProcess))
        {
            for (TPage *p = pager.inProcess.first; p != NULL; p = p->next)
            {
                if (p->targetState == UNUSED)
                {
                    _processFilePage(pager, *p, True());
                    erase(pager.inProcess, *p);
                    return p;
                }
            }
        }
        TPage *p = _swapOutFilePage(pager);
        if (p != NULL)
            return p;
        else
            return new TPage(); // we can't swap out anything, so we have to allocate more pages
    }
}

// ----------------------------------------------------------------------------
// Function _lockFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TFilePos, typename TSize>
inline FilePage<TValue, TSpec> &
_lockFilePage(FilePageTable<TValue, TDirection, TSpec> & pager, TFilePos filePos, TSize size, bool /*readOnly*/)
{
    typedef FilePage<TValue, TSpec> TPage;
    TPage * p;

    SEQAN_OMP_PRAGMA(critical(lockPageTable))
    {
        unsigned pageNo = filePos / pager.table.pageSize;
        if (length(pager.table.frameStart) <= pageNo)
            resize(pager.table.frameStart, pageNo + 1, pager.table.EMPTY);

        p = static_cast<TPage *>(pager.table.frameStart[pageNo]);

        if (p == pager.table.EMPTY || p == pager.table.ON_DISK)
        {
            p = _newFilePage(pager);
            p->filePos = filePos;
            p->size = size;
            pager.table.frameStart[pageNo] = p;
        }
        ++(p->lockCount);
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

    TPage & page = _lockFilePage(pager, filePos, size, readOnly);
    page.targetState = READY;
    _processFilePage(pager, page, True());
    pushBack(pager.ready, page);
    return page;
}

// ----------------------------------------------------------------------------
// Function releaseFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TFilePage>
inline void
releaseFilePage(FilePageTable<TValue, TDirection, TSpec> &pager, TFilePage & page)
{
    // we automatically flush the page to disk here,
    // to speed up the swap-out process
    if (atomicDec(page.lockCount) == 0)
    {
        page.targetState = UNUSED;
        erase(pager.ready, page);
        if (_processFilePage(pager, page, False()))
            pushBack(pager.unused, page);
        else
            pushBack(pager.inProcess, page);
    }
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

    typedef FilePageTable<TValue, TDirection, TSpec>    TFilePageTable;
    typedef typename TFilePageTable::TPageFrame         TPageFrame;
    typedef typename TFilePageTable::TFile              TFile;

    typedef typename Size<TFile>::Type                  TSize;
    typedef typename TBase::int_type                    TIntValue;
    typedef typename TBase::off_type                    TDifference;
    typedef typename TBase::pos_type                    TPosition;

    TFilePageTable  pager;
    TPageFrame      *readPage;
    TPageFrame      *writePage;

    TPosition       readPagePos;
    TPosition       writePagePos;

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

    inline void _stop()
    {
        this->setg(NULL, NULL, NULL);
        this->setp(NULL, NULL);

        if (readPage != NULL)
        {
            releaseFilePage(pager, *readPage);
            readPage = NULL;
        }
        if (writePage != NULL)
        {
            releaseFilePage(pager, *writePage);
            writePage = NULL;
        }
    }

    bool _nextReadPage()
    {
        if (readPage != NULL)
        {
            readPagePos += readPage->size;
            releaseFilePage(pager, *readPage);
        }

        if (pager.fileSize <= readPagePos)
            return false;

        Pair<__int64, unsigned> ol = _getPageOffsetAndLength(pager, readPagePos);
        readPage = &fetchFilePage(pager, ol.i1, ol.i2);
        this->setg(readPage->data.begin, readPage->data.begin, readPage->data.end);
        return true;
    }

    bool _nextWritePage()
    {
        if (writePage != NULL)
        {
            writePagePos += writePage->size;
            releaseFilePage(pager, *writePage);
        }

        Pair<__int64, unsigned> ol = _getPageOffsetAndLength(pager, writePagePos);
        writePage = &fetchFilePage(pager, ol.i1, ol.i2);
        this->setp(writePage->data.begin, writePage->data.begin + capacity(*writePage));
        return true;
    }

    virtual TIntValue
    overflow(TIntValue val)
    {
        if (SEQAN_UNLIKELY(!pager.file))
            return TTraits::eof();

        if (SEQAN_UNLIKELY(this->pptr() >= this->epptr() && !_nextWritePage()))
        {
            this->setp(NULL, NULL);
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

        if (SEQAN_UNLIKELY(this->gptr() >= this->egptr() && !_nextReadPage()))
        {
            this->setg(NULL, NULL, NULL);
            return TTraits::eof();
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

    TPosition _tell(std::ios::openmode which)
    {
        if (SEQAN_LIKELY(pager.file))
        {
            if (which == std::ios::in)
            {
                if (SEQAN_LIKELY(readPage != NULL))
                    return (TSize)readPagePos + (gptr() - eback());
                else
                    return 0;
            }
            else
            {
                if (SEQAN_LIKELY(writePage != NULL))
                    return (TSize)writePagePos + (pptr() - pbase());
                else
                    return 0;
            }
        }
        else
        {
            return -1;
        }
    }

    TPosition _seek(TPosition pos, Input)
    {
        if (readPage != NULL)
        {
            if (readPage->filePos <= pos && pos < readPage->filePos + readPage->size)
            {
                this->setg(readPage->data.begin, readPage->data.begin + (pos - readPage->filePos), readPage->data.end);
                return pos;
            }

            // Release former page.
            releaseFilePage(pager, *readPage);
            readPage = NULL;
        }

        if (SEQAN_UNLIKELY(pos < (TPosition)0))
            return -1;


        // Fetch new page.
        Pair<__int64, unsigned> ol = _getPageOffsetAndLength(pager, pos);
        readPage = &fetchFilePage(pager, ol.i1, ol.i2);
        readPagePos = readPage->filePos;

        // Seek to correct position in page.
        this->setg(readPage->data.begin, readPage->data.begin + (pos - readPage->filePos), readPage->data.end);

        return pos;
    }

    TPosition _seek(TPosition pos, Output)
    {
        if (writePage != NULL)
        {
            if (writePage->filePos <= pos && pos < writePage->filePos + writePage->size)
            {
                this->setg(writePage->data.begin, writePage->data.begin + (pos - writePage->filePos), writePage->data.end);
                return pos;
            }

            // Release former page.
            releaseFilePage(pager, *writePage);
            writePage = NULL;
        }

        if (SEQAN_UNLIKELY(pos < (TPosition)0))
            return -1;


        // Fetch new page.
        Pair<__int64, unsigned> ol = _getPageOffsetAndLength(pager, pos);
        writePage = &fetchFilePage(pager, ol.i1, ol.i2);
        writePagePos = writePage->filePos;

        // Seek to correct position in page.
        this->setp(writePage->data.begin, writePage->data.begin + capacity(*writePage));
        this->pbump(pos - writePage->filePos);

        return pos;
    }

    virtual TPosition
    seekpos(
        TPosition pos,
        std::ios::openmode which)
    {
        SEQAN_ASSERT_NEQ((int)(which & IosOpenMode<TDirection>::VALUE), 0);

        if (SEQAN_UNLIKELY(!pager.file))
            return -1;

        if (which == std::ios::in)
            return _seek(pos, Input());
        else
            return _seek(pos, Output());
    }

    virtual TPosition
    seekoff(
        TDifference off,
        std::ios::seekdir dir,
        std::ios::openmode which)
    {
        SEQAN_ASSERT_NEQ((int)(which & IosOpenMode<TDirection>::VALUE), 0);

        if (SEQAN_UNLIKELY(!pager.file))
            return -1;

        TPosition pos;

        if (dir == std::ios::beg)
            pos = off;
        else if (dir == std::ios::cur)
            pos = _tell(which) + off;
        else
            pos = pager.fileSize + off;

        if (which == std::ios::in)
            return _seek(pos, Input());
        else
            return _seek(pos, Output());
    }

};

template <typename TValue, typename TDirection, typename TSpec>
inline void
clear(FileStreamBuffer<TValue, TDirection, TSpec> & buffer)
{
    buffer.setp(NULL, NULL);
    buffer.setg(NULL, NULL, NULL);

    buffer.readPage = NULL;
    buffer.writePage = NULL;
    buffer.readPagePos = 0;
    buffer.writePagePos = 0;
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
