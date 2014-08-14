/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003   (original zlib stream)
Author: David Weese, dave.weese@gmail.com, 2014             (extension to parallel block-wise compression in bgzf format)
*/

#ifndef BGZFSTREAM_HPP
#define BGZFSTREAM_HPP

#include <vector>
#include <iostream>
#include <algorithm>
#include <zlib.h>
#include "zutil.h"

namespace seqan {

/// default gzip buffer size,
/// change this to suite your needs
const size_t default_buffer_size = 4096;

const unsigned BGZF_MAX_BLOCK_SIZE = 64 * 1024;
const unsigned BGZF_BLOCK_HEADER_LENGTH = 18;
const unsigned BGZF_BLOCK_FOOTER_LENGTH = 8;
const unsigned ZLIB_BLOCK_OVERHEAD = 5; // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)

// Reduce the maximal input size, such that the compressed data 
// always fits in one block even for level Z_NO_COMPRESSION.
const unsigned BGZF_BLOCK_SIZE = BGZF_MAX_BLOCK_SIZE - BGZF_BLOCK_HEADER_LENGTH - BGZF_BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;

/// Compression strategy, see bgzf doc.
enum EStrategy
{
	StrategyFiltered = 1,
	StrategyHuffmanOnly = 2,
	DefaultStrategy = 0
};


// new Event class

struct Event2
{
    enum { Infinite = LONG_MAX };

    pthread_cond_t  data_cond;
    Mutex &         mutex;

    explicit
    Event2(Mutex &mutex) :
        mutex(mutex)
    {
        int result = pthread_cond_init(&data_cond, NULL);
        SEQAN_ASSERT_EQ(result, 0);
    }

    ~Event2()
    {
        int result = pthread_cond_destroy(&data_cond);
        SEQAN_ASSERT_EQ(result, 0);
    }
};

inline void
waitFor(Event2 &event)
{
    int result = pthread_cond_wait(&event.data_cond, event.mutex.hMutex);
    SEQAN_ASSERT_EQ(result, 0);
}

inline void
waitFor(Event2 &event, long timeoutMilliSec, bool & inProgress)
{
    if (timeoutMilliSec != Event2::Infinite)
    {
        timespec ts;
        ts.tv_sec = timeoutMilliSec / 1000;
        ts.tv_nsec = (timeoutMilliSec % 1000) * 1000;
        int result = pthread_cond_timedwait(&event.data_cond, event.mutex.hMutex, &ts);
        inProgress = (result == ETIMEDOUT);
        SEQAN_ASSERT(result == 0 || inProgress);
    }
    else
    {
        inProgress = false;
        waitFor(event);
    }
}

inline void
signal(Event2 &event)
{
    int result = pthread_cond_broadcast(&event.data_cond);
    SEQAN_ASSERT_EQ(result, 0);
}



struct Suspendable_;
typedef Tag<Suspendable_> Suspendable;

template <typename TValue>
class ConcurrentQueue<TValue, Suspendable>
{
public:
    typedef typename Host<ConcurrentQueue>::Type    TString;
    typedef typename Size<TString>::Type            TSize;

    size_t      readerCount;
    size_t      writerCount;

    TString     data;
    TSize       occupied;
    TSize       nextIn;
    TSize       nextOut;

    Mutex mutex;
    Event2 more;
    Event2 less;

    ConcurrentQueue(TSize maxSize):
        readerCount(0),
        writerCount(0),
        occupied(0),
        nextIn(0),
        nextOut(0),
        mutex(false),
        more(mutex),
        less(more)
    {
        reserve(data, maxSize, Exact());
    }

    ~ConcurrentQueue()
    {
        TValue tmp;

        // wait for all pending readers to finish
        while (readerCount != 0u)
        {}

        while (popFront(tmp, *this))
        {}
    }
};

template <typename TValue>
inline void
unlockWriting(ConcurrentQueue<TValue, Suspendable> & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    if (--me.writerCount == 0u)
        signal(me.more);
}

template <typename TValue>
inline void
unlockReading(ConcurrentQueue<TValue, Suspendable> & me)
{
    ScopedLock<Mutex> lock(me.mutex);
    if (--me.readerCount == 0u)
        signal(me.less);
}

template <typename TValue>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, Suspendable> & me)
{
    typedef ConcurrentQueue<TValue, Suspendable>        TQueue;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    ScopedLock<Mutex> lock(me.mutex);
    TSize cap = capacity(me.data);

    while (me.occupied == 0u && me.writerCount > 0u)
        waitFor(me.more);

    if (me.occupied == 0u)
        return false;

    SEQAN_ASSERT_NEQ(me.occupied, 0u);

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + me.nextOut;
    std::swap(result, *it);
    valueDestruct(it);

    me.nextOut = (me.nextOut + 1) % cap;
    me.occupied--;

    /* now: either me.occupied > 0 and me.nextout is the index
       of the next occupied slot in the buffer, or
       me.occupied == 0 and me.nextout is the index of the next
       (empty) slot that will be filled by a producer (such as
       me.nextout == me.nextin) */

    signal(me.less);
    return true;
}

template <typename TValue, typename TValue2>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable> & me,
            TValue2 SEQAN_FORWARD_CARG val)
{
    typedef ConcurrentQueue<TValue, Suspendable>        TQueue;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;

    ScopedLock<Mutex> lock(me.mutex);
    TSize cap = capacity(me.data);

    while (me.occupied >= cap && me.readerCount > 0u)
        waitFor(me.less);

    if (me.occupied >= cap)
        return false;

    SEQAN_ASSERT_LT(me.occupied, cap);

    valueConstruct(begin(me.data, Standard()) + me.nextIn, val);
    me.nextIn = (me.nextIn + 1) % cap;
    me.occupied++;

    /* now: either me.occupied < BSIZE and me.nextin is the index
       of the next empty slot in the buffer, or
       me.occupied == BSIZE and me.nextin is the index of the
       next (occupied) slot that will be emptied by a consumer
       (such as me.nextin == me.nextout) */

    signal(me.more);
    return true;
}

template <typename TValue, typename TSize>
inline bool
waitForMinSize(ConcurrentQueue<TValue, Suspendable> & me,
               TSize minSize)
{
    ScopedLock<Mutex> lock(me.mutex);
    while (me.occupied < minSize && me.writerCount > 0u)
        waitFor(me.more);
    return me.occupied >= minSize;
}

/** \brief A stream decorator that takes raw input and zips it to a ostream.

The class wraps up the inflate method of the bgzf library 1.1.4 http://www.gzip.org/bgzf/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_streambuf : public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef typename Tr::char_type char_type;
	typedef typename Tr::int_type int_type;

    typedef ConcurrentQueue<size_t, Suspendable>  TJobQueue;

    struct Concatter
    {
        ostream_reference ostream;
        Mutex lock;
        unsigned nextWritableKey;
        bool stop;

        Concatter(ostream_reference ostream) :
            ostream(ostream),
            lock(false),
            nextWritableKey(0),
            stop(false)
        {}
    };

    // serialization of output
    Concatter concatter;

    struct CompressionJob
    {
        typedef std::vector<byte_type, byte_allocator_type> TOutputBuffer;
        typedef std::vector<char_type, char_allocator_type> TBuffer;

        TOutputBuffer   outputBuffer;
        TBuffer         buffer;
        size_t          size;
        unsigned        key;

        CompressionJob() :
            outputBuffer(BGZF_MAX_BLOCK_SIZE, 0),
            buffer(BGZF_BLOCK_SIZE / sizeof(char_type), 0),
            key(0)
        {}
    };

    // string of recycable jobs
    String<CompressionJob>  jobs;
    TJobQueue               jobQueue;
    TJobQueue               idleQueue;
    unsigned                currentJobKey;
    size_t                  currentJobId;
    bool                    currentJobAvail;


    struct CompressionThread
    {
        basic_bgzf_streambuf            *streamBuf;
        CompressionContext<BgzfFile>    compressionCtx;

        void operator()()
        {
            ScopedReadLock<TJobQueue> readLock(streamBuf->jobQueue);
            ScopedWriteLock<TJobQueue> writeLock(streamBuf->idleQueue);

            // wait for a new job to become available
            while (true)
            {
                size_t jobId;
                if (!popFront(jobId, streamBuf->jobQueue))
                    return;

                CompressionJob &job = streamBuf->jobs[jobId];

                // compress block with zlib
                size_t outputLen = _compressBlock(
                    &job.outputBuffer[0], capacity(job.outputBuffer),
                    &job.buffer[0], job.size, compressionCtx);

                while (true)
                {
                    ScopedLock<Mutex> scopedLock(streamBuf->concatter.lock);
                    if (streamBuf->concatter.nextWritableKey == job.key)
                    {
                        streamBuf->concatter.ostream.write(
                            (const char_type*) &(job.outputBuffer[0]),
                            outputLen);
                        if (!streamBuf->concatter.ostream.good())
                        {
                            streamBuf->concatter.stop = false;
                            return;
                        }
                        // move serializer to next packet and reset our job
                        streamBuf->concatter.nextWritableKey++;

                        appendValue(streamBuf->idleQueue, jobId);
                        break;
                    }
                }
            }
        }
    };

    // array of worker threads
    Thread<CompressionThread>   *threads;
    size_t                      numThreads;
    size_t                      numJobs;

    /** Construct a zip stream
     * More info on the following parameters can be found in the bgzf documentation.
     */
    basic_bgzf_streambuf(ostream_reference ostream_,
                         size_t numThreads = 4,
                         size_t numJobs = 16) :
		concatter(ostream_),
        jobQueue(numJobs),
        idleQueue(numJobs),
        numThreads(numThreads),
        numJobs(numJobs)
    {
        resize(jobs, numJobs, Exact());
        currentJobId = 0;
        currentJobKey = 0;

        lockWriting(jobQueue);
        lockReading(idleQueue);

        for (unsigned i = 0; i < numJobs; ++i)
        {
            bool success = appendValue(idleQueue, i);
            SEQAN_ASSERT(success);
        }

        threads = new Thread<CompressionThread>[numThreads];
        for (unsigned i = 0; i < numThreads; ++i)
        {
            threads[i].worker.streamBuf = this;
            run(threads[i]);
        }

        currentJobAvail = popFront(currentJobId, idleQueue);
        SEQAN_ASSERT(currentJobAvail);

        CompressionJob &job = jobs[currentJobId];
        job.key = currentJobKey++;
		this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
    }

    ~basic_bgzf_streambuf()
    {
        // the buffer is now (after addFooter()) and flush will append the empty EOF marker
        flush(true);

        unlockWriting(jobQueue);
        unlockReading(idleQueue);

        for (unsigned i = 0; i < numThreads; ++i)
            waitFor(threads[i]);
        delete[] threads;
    }

    bool compressBuffer(size_t size)
    {
        // submit current job
        if (currentJobAvail)
        {
            jobs[currentJobId].size = size;
            appendValue(jobQueue, currentJobId);
        }

        // recycle existing idle job
        if (!(currentJobAvail = popFront(currentJobId, idleQueue)))
            return false;

        CompressionJob &job = jobs[currentJobId];
        job.key = currentJobKey++;

        return !concatter.stop;
    }

    int_type overflow(int_type c)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if (c != EOF)
        {
            *this->pptr() = c;
            ++w;
        }
        if (compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
            return c;
        }
        else
        {
            return EOF;
        }
    }

	/** flushes the zip buffer and output buffer.

	This method should be called at the end of the compression. Calling flush multiple times, will lower the
	compression ratio.
	*/
	std::streamsize flush(bool flushEmptyBuffer = false)
    {
        int w = static_cast<int>(this->pptr() - this->pbase());
        if ((w != 0 || flushEmptyBuffer) && compressBuffer(w))
        {
            CompressionJob &job = jobs[currentJobId];
            this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
        }
        else
        {
            w = 0;
        }

        // wait for running compressor threads
        waitForMinSize(idleQueue, numJobs - 1);

		concatter.ostream.flush();
		return w;
    }

	int sync()
    {
		if (this->pptr() != this->pbase())
		{
			int c = overflow(EOF);
			if (c == EOF)
				return -1;
        }
        return 0;
    }

    void addFooter()
    {
        // we flush the filled buffer here, so that an empty (EOF) buffer is flushed in the d'tor
        if (this->pptr() != this->pbase())
            overflow(EOF);
    }

	/// returns a reference to the output stream
	ostream_reference get_ostream() const	{ return concatter.ostream; };
};

/** \brief A stream decorator that takes compressed input and unzips it to a istream.

The class wraps up the deflate method of the bgzf library 1.1.4 http://www.gzip.org/bgzf/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_unbgzf_streambuf : 
	public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
	typedef typename Tr::char_type char_type;
	typedef typename Tr::int_type int_type;

    typedef std::vector<char_type, char_allocator_type> TBuffer;
    typedef ConcurrentQueue<size_t, Suspendable>  TJobQueue;

    static const size_t MAX_PUTBACK = 4;

    struct Serializer
    {
        istream_reference istream;
        Mutex lock;
        unsigned waitForKey;
        unsigned nextKey;
        bool stop;

        Serializer(istream_reference istream) :
            istream(istream),
            lock(false),
            stop(false)
        {}
    };

    Serializer serializer;

    struct DecompressionJob
    {
        typedef std::vector<byte_type, byte_allocator_type> TInputBuffer;

        TInputBuffer    inputBuffer;
        TBuffer         buffer;
        size_t          size;

        Mutex           mutex;
        Event           readyEvent;
        bool            ready;

        DecompressionJob() :
            inputBuffer(BGZF_MAX_BLOCK_SIZE, 0),
            buffer(MAX_PUTBACK + BGZF_MAX_BLOCK_SIZE / sizeof(char_type), 0),
            ready(mutex)
        {}

        DecompressionJob(DecompressionJob const &) :
            inputBuffer(BGZF_MAX_BLOCK_SIZE, 0),
            buffer(MAX_PUTBACK + BGZF_MAX_BLOCK_SIZE / sizeof(char_type), 0),
            ready(mutex)
        {}
    };

    // string of recycable jobs
    String<DecompressionJob>    jobs;
    TJobQueue                   runningQueue;
    TJobQueue                   idleQueue;
    size_t                      currentJobId;
    bool                        currentJobAvail;

    struct DecompressionThread
    {
        basic_unbgzf_streambuf          *streamBuf;
        CompressionContext<BgzfFile>    compressionCtx;

        void operator()()
        {
            ScopedWriteLock<TJobQueue> readLock(streamBuf->idleQueue);
            ScopedReadLock<TJobQueue> writeLock(streamBuf->runningQueue);

            // wait for a new job to become available
            while (true)
            {
                size_t jobId;
                if (!popFront(jobId, streamBuf->idleQueue))
                    return;

                DecompressionJob &job = streamBuf->jobs[jobId];
                size_t tailLen;

                {
                    ScopedLock<Mutex> scopedLock(streamBuf->serializer.lock);

                    // read header
                    streamBuf->serializer->istream.read(
                        (char*)&(streamBuf->inputBuffer[0]),
                        BGZF_BLOCK_HEADER_LENGTH);

                    // check header
                    if (!streamBuf->serializer->istream.good() || !_bgzfCheckHeader(&(streamBuf->inputBuffer[0])))
                    {
                        streamBuf->serializer->stop = true;
                        return;
                    }

                    // extract length of compressed data
                    tailLen = _bgzfUnpack16(&(streamBuf->inputBuffer[16])) + 1u - BGZF_BLOCK_HEADER_LENGTH;

                    // read compressed data and tail
                    streamBuf->serializer->istream.read(
                        (char*)&(streamBuf->inputBuffer[BGZF_BLOCK_HEADER_LENGTH]),
                        tailLen);
                    
                    if (!appendValue(streamBuf->runningQueue, jobId))
                        return;
                }

                // decompress block
                streamBuf->size = _decompressBlock(
                    &job.buffer[MAX_PUTBACK], capacity(job.buffer),
                    &streamBuf->inputBuffer[0], BGZF_BLOCK_HEADER_LENGTH + tailLen, streamBuf->ctx);

                // move serializer to next packet and reset our job
                streamBuf->serializer.nextWritableKey++;

                appendValue(streamBuf->runningQueue, jobId);
                break;
            }
        }
    };
/*

    struct BgzfThreadContext
    {
        typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;

        struct BgzfDecompressor
        {
            BgzfThreadContext *threadCtx;

            void operator()()
            {
                while (true)
                {
                    size_t tailLen;
                    {
                        ScopedLock<Mutex> scopedLock(threadCtx->serializer->lock);

                        if (threadCtx->serializer->stop)
                            return;

                        if (threadCtx->serializer->waitForKey != threadCtx->waitForKey)
                            continue;

                        // read header
                        threadCtx->serializer->istream.read(
                            (char*)&(threadCtx->inputBuffer[0]),
                            BGZF_BLOCK_HEADER_LENGTH);

                        // check header
                        if (!threadCtx->serializer->istream.good() || !_bgzfCheckHeader(&(threadCtx->inputBuffer[0])))
                        {
                            threadCtx->serializer->stop = true;
                            return;
                        }

                        // extract length of compressed data
                        tailLen = _bgzfUnpack16(&(threadCtx->inputBuffer[16])) + 1u - BGZF_BLOCK_HEADER_LENGTH;

                        // read compressed data and tail
                        threadCtx->serializer->istream.read(
                            (char*)&(threadCtx->inputBuffer[BGZF_BLOCK_HEADER_LENGTH]),
                            tailLen);

                        if (!threadCtx->serializer->istream.good())
                        {
                            threadCtx->serializer->stop = true;
                            return;
                        }

                        threadCtx->serializer->waitForKey = threadCtx->waitForKey + 1;
                    }

                    // decompress block
                    threadCtx->size = _decompressBlock(
                        &threadCtx->buffer[MAX_PUTBACK], capacity(threadCtx->buffer),
                        &threadCtx->inputBuffer[0], BGZF_BLOCK_HEADER_LENGTH + tailLen, threadCtx->ctx);
                }
            }
        };

        byte_vector_type                inputBuffer;
        char_vector_type                buffer;
        int                             size;
        CompressionContext<BgzfFile>    ctx;
        Thread<BgzfDecompressor>        decompressor;

        Serializer *serializer;
        unsigned waitForKey;

        BgzfThreadContext():
            inputBuffer(BGZF_MAX_BLOCK_SIZE, 0),
            buffer(MAX_PUTBACK + BGZF_MAX_BLOCK_SIZE / sizeof(char_type), 0)
        {}
    };
*/


    // array of worker threads
    Thread<DecompressionThread> *threads;
    TBuffer                     putbackBuffer;
    size_t                      numThreads;
    size_t                      numJobs;

    /** Construct a unzip stream
    * More info on the following parameters can be found in the bgzf documentation.
    */
    basic_unbgzf_streambuf(istream_reference istream_,
                           size_t numThreads = 4,
                           size_t numJobs = 16) :
		serializer(istream_),
        runningQueue(numJobs),
        idleQueue(numJobs),
        numThreads(numThreads),
        numJobs(numJobs),
        putbackBuffer(MAX_PUTBACK)
    {
        resize(jobs, numJobs, Exact());
        currentJobId = 0;

        lockReading(runningQueue);
        lockWriting(idleQueue);

        for (unsigned i = 0; i < numJobs; ++i)
        {
            bool success = appendValue(idleQueue, i);
            SEQAN_ASSERT(success);
        }

        threads = new Thread<DecompressionThread>[numThreads];
        for (unsigned i = 0; i < threads; ++i)
        {
            threads[i].worker.streamBuf = this;
            run(threads[i]);
        }

        currentJobAvail = popFront(currentJobId, runningQueue);
        SEQAN_ASSERT(currentJobAvail);

        // wait for first job
        DecompressionJob &job = jobs[currentJobId];
        {
            ScopedLock<Mutex> lock(job.mutex);
            if (!job.ready)
                waitFor(job.readyEvent);
        }
		this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
    }

	~basic_unbgzf_streambuf()
    {
        unlockWriting(jobQueue);
        unlockReading(idleQueue);

        for (unsigned i = 0; i < numThreads; ++i)
            waitFor(threads[i]);
        delete[] threads;
    }

    int_type underflow()
    {
//
//        // submit current job
//        if (currentJobAvail)
//        {
//            jobs[currentJobId].size = size;
//            appendValue(jobQueue, currentJobId);
//        }
//
//        // recycle existing idle job
//        if (!(currentJobAvail = popFront(currentJobId, idleQueue)))
//            return false;
//
//        CompressionJob &job = jobs[currentJobId];
//        job.key = currentJobKey++;
//
//        return !concatter.stop;


        if (this->gptr() && this->gptr() < this->egptr())
            return *this->gptr();

        size_t putback = this->gptr() - this->eback();
        if (putback > MAX_PUTBACK)
            putback = MAX_PUTBACK;

        // save at most MAX_PUTBACK characters from previous page to putback buffer
        if (putback != 0)
            std::copy(
                this->gptr() - putback,
                this->gptr(),
                &putbackBuffer[0]);

        appendValue(idleQueue, currentJobId);

        do
        {
//            currentJobAvail = popFront(currentJobId, idleQueue);
//            SEQAN_ASSERT(currentJobAvail);
//
//            // wait for first job
//            DecompressionJob &job = jobs[currentJobId];
//            {
//                ScopedLock<Mutex> lock(job.mutex);
//                if (!job.ready)
//                    waitFor(job.readyEvent);
//            }
//            this->setp(&(job.buffer[0]), &(job.buffer[job.buffer.size() - 1]));
//
            ctx = &threadCtx[nextThread];

            do {
                threadCtx[nextRunThread].size = -1;
                threadCtx[nextRunThread].waitForKey = serializer.nextKey++;
                run(threadCtx[nextRunThread].decompressor);
                nextRunThread = (nextRunThread + 1) % threads;
            } while (nextRunThread != nextThread);

            nextThread = (nextThread + 1) % threads;

            // restore putback buffer
            if (putback != 0)
                std::copy(
                    &putbackBuffer[0],
                    &putbackBuffer[putback],
                    &ctx->buffer[MAX_PUTBACK - putback]);

            waitFor(ctx->decompressor);

            if (ctx->size == -1) // EOF
               return EOF;

            // reset buffer pointers
            this->setg( 
                  &ctx->buffer[MAX_PUTBACK - putback],      // beginning of putback area
                  &ctx->buffer[MAX_PUTBACK],                // read position
                  &ctx->buffer[MAX_PUTBACK + ctx->size]);   // end of buffer
        } while (ctx->size == 0);

        // return next character
        return *this->gptr();
    }

	/// returns the compressed input istream
	istream_reference get_istream()	{ return serializer.istream;};
};

/*! \brief Base class for zip ostreams

Contains a basic_bgzf_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
	typedef basic_bgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > bgzf_streambuf_type;

    /** Construct a zip stream
     * More info on the following parameters can be found in the bgzf documentation.
     */
	basic_bgzf_ostreambase(ostream_reference ostream_)
		: m_buf(ostream_)
	{
		this->init(&m_buf );
	};
	
	/// returns the underlying zip ostream object
	bgzf_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the bgzf error state
	int get_zerr() const					{	return m_buf.get_err();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the compressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the uncompressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	bgzf_streambuf_type m_buf;
};

/*! \brief Base class for unzip istreams

Contains a basic_unbgzf_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
	typedef basic_unbgzf_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > unbgzf_streambuf_type;

	basic_bgzf_istreambase(istream_reference ostream_)
		: m_buf(ostream_)
	{
		this->init(&m_buf );
	};
	
	/// returns the underlying unzip istream object
	unbgzf_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the bgzf error state
	int get_zerr() const					{	return m_buf.get_zerr();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the uncompressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the compressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	unbgzf_streambuf_type m_buf;
};

/*! \brief A zipper ostream

This class is a ostream decorator that behaves 'almost' like any other ostream.

At construction, it takes any ostream that shall be used to output of the compressed data.

When finished, you need to call the special method zflush or call the destructor 
to flush all the intermidiate streams.

Example:
\code
// creating the target zip string, could be a fstream
ostringstream ostringstream_;
// creating the zip layer
bgzf_ostream zipper(ostringstream_);

	
// writing data	
zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// zip ostream needs special flushing...
zipper.zflush();
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_bgzf_ostream : 
	public basic_bgzf_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_ostream<Elem,Tr>
{
public:
	typedef basic_bgzf_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_ostreambase_type;
	typedef std::basic_ostream<Elem,Tr> ostream_type;
    typedef ostream_type& ostream_reference;

	/** Constructs a zipper ostream decorator
	 *
	 * \param ostream_ ostream where the compressed output is written

	 When is_gbgzf_ is true, a gzip header and footer is automatically added.
	 */
	basic_bgzf_ostream(ostream_reference ostream_)
	: 
		bgzf_ostreambase_type(ostream_),
		ostream_type(this->rdbuf())
	{}

	/// flush inner buffer and zipper buffer
	basic_bgzf_ostream<Elem,Tr>& zflush()	
	{	
		this->flush(); this->rdbuf()->flush(); return *this;
	};

    ~basic_bgzf_ostream()
    {
        this->rdbuf()->addFooter();
    }


private:
    static void put_long(ostream_reference out_, unsigned long x_);
};

/*! \brief A zipper istream

This class is a istream decorator that behaves 'almost' like any other ostream.

At construction, it takes any istream that shall be used to input of the compressed data.

Simlpe example:
\code
// create a stream on zip string
istringstream istringstream_( ostringstream_.str());
// create unzipper istream
bgzf_istream unzipper( istringstream_);

// read and unzip
unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bgzf_istream : 
	public basic_bgzf_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_istream<Elem,Tr>
{
public:
	typedef basic_bgzf_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bgzf_istreambase_type;
	typedef std::basic_istream<Elem,Tr> istream_type;
    typedef istream_type& istream_reference;
	typedef char byte_type;

	/** Construct a unzipper stream
	 *
	 * \param istream_ input buffer
	 */
	basic_bgzf_istream(istream_reference istream_)
	  : 
		bgzf_istreambase_type(istream_),
		istream_type(this->rdbuf()),
		m_is_gzip(false),
		m_gbgzf_data_size(0)
	{};

	/// returns true if it is a gzip file
	bool is_gzip() const				{	return m_is_gzip;};
	/// return data size check
	bool check_data_size() const		{	return this->get_out_size() == m_gbgzf_data_size;};

	/// return the data size in the file 
	long get_gbgzf_data_size() const		{	return m_gbgzf_data_size;};
protected:
    static void read_long(istream_reference in_, unsigned long& x_);

	int check_header();
	bool m_is_gzip;
	unsigned long m_gbgzf_data_size;
};

/// A typedef for basic_bgzf_ostream<char>
typedef basic_bgzf_ostream<char> bgzf_ostream;
/// A typedef for basic_bgzf_ostream<wchar_t>
typedef basic_bgzf_ostream<wchar_t> bgzf_wostream;
/// A typedef for basic_bgzf_istream<char>
typedef basic_bgzf_istream<char> bgzf_istream;
/// A typedef for basic_bgzf_istream<wchart>
typedef basic_bgzf_istream<wchar_t> bgzf_wistream;

}  // namespace seqan

#include "bgzfstream_impl.h"

#endif

