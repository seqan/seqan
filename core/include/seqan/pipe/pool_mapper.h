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

#ifndef SEQAN_HEADER_POOL_MAPPER_H
#define SEQAN_HEADER_POOL_MAPPER_H

namespace SEQAN_NAMESPACE_MAIN
{

    // external synchronous permutation mapping

/**
.Spec.MapperConfigSize:
..cat:Pipelining
..general:Spec.MapperSpec
..summary:Configuration of Mapper.
..signature:MapperConfigSize<TMap, TSize, TFile>
..param.TMap:The destination function (see STL's $unary_function$).
...remarks:This functions maps a stream element to its destined position. The result type of this unary function should convertible to $TSize$.
...note:The destination function must be bijective.
..param.TSize:The Mapper's size type.
..param.TFile:The underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..see:Spec.MapperConfig
..include:seqan/pipe.h
*/

    template < typename TMap,
		       typename TSize,
		       typename TFile = File<> >						// default file type
    struct MapperConfigSize {
        typedef TMap        Map;
		typedef TSize		SizeType;
        typedef TFile       File;
    };

/**
.Spec.MapperConfig:
..cat:Pipelining
..general:Spec.MapperSpec
..summary:Configuration of Mapper.
..signature:MapperConfig<TMap, TFile>
..param.TMap:The destination function (see STL's $unary_function$).
...remarks:This functions maps a stream element to its destined position. The result type of this unary function should convertible to $TSize$.
...note:The destination function must be bijective.
..param.TFile:The underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:Using this configuration spec., the Mapper's size type is $Size<TFile>::Type$. To use a custom size type @Spec.MapperConfigSize@ should be used.
..see:Spec.MapperConfigSize
..include:seqan/pipe.h
*/

    template < typename TMap,
		       typename TFile = File<> >
    struct MapperConfig {
        typedef TMap						Map;
		typedef typename Size<TFile>::Type	SizeType;
        typedef TFile						File;
    };

/**
.Spec.MapperSpec:
..cat:Pipelining
..general:Class.Pool
..summary:Permutes all elements using a custom destination function.
..signature:Pool<TValue, MapperSpec<TConfig> >
..param.TValue:The value type, that is the type of the stream elements.
..param.TConfig:Configuration Spec. Defines destination function, size type, and file type.
...type:Spec.MapperConfig
...type:Spec.MapperConfigSize
..remarks:The Pool's input/output type is $TValue$ and the size type is determined by the $TConfig$.
..include:seqan/pipe.h
*/

    template < typename TConfig >
    struct MapperSpec {
        typedef TConfig Config;
    };

    template < typename TValue,
			   typename TConfig >
    struct HandlerArgs< Pool< TValue, MapperSpec<TConfig> > >
    {
        typedef typename TConfig::Map Type;
    };

    struct ReadMapperSpec_;
	typedef Tag<ReadMapperSpec_> ReadMapperSpec;

/*
    // mapping phase 2
    template < typename TValue,
			   typename TConfig >
    inline SimpleBuffer< TValue > & processBuffer(
        PageFrame< TValue, typename TConfig::File, Dynamic<> > &buf,
        BufferHandler< Pool< TValue, MapperSpec<TConfig> >, ReadFileSpec > &_me)
    {
		typedef BufferHandler< Pool< TValue, MapperSpec<TConfig> >, ReadMapperSpec > Handler;
        Handler *me = static_cast<Handler*>(&_me);
        
        typename Size< Pool< TValue, MapperSpec<TConfig> > >::Type offset = buf.pageNo;
	    offset *= (unsigned)pageSize(buf);

		typename TConfig::Map M = me->pool.handlerArgs;
        for(TValue *cur = buf.begin; cur != buf.end; ++cur) {
            #ifdef SEQAN_DEBUG
                if (!(M(*cur) >= offset && M(*cur) < offset + pageSize(buf))) {
					::std::cerr << "Mapper assertion failed: " << ::std::hex << M(*cur);
					::std::cerr << " not in [" << offset << "," << (offset + pageSize(buf)) << ") at " << (cur - buf.begin);
                    ::std::cerr << " element is " << ::std::dec << *cur << ::std::endl;
                }
            #endif
            SEQAN_ASSERT(M(*cur) >= offset && M(*cur) < offset + pageSize(buf));
            me->mapBuffer[M(*cur) - offset] = *cur;
        }
		resize(me->mapBuffer, size(buf));
        return me->mapBuffer;
    }
*/

    // mapping phase 2 (in-place)
    template < typename TValue,
			   typename TConfig >
    inline SimpleBuffer< TValue > & processBuffer(
        PageFrame< TValue, typename TConfig::File, Dynamic<> > &buf,
        BufferHandler< Pool< TValue, MapperSpec<TConfig> >, ReadFileSpec > &me)
    {
        typedef typename Size< Pool< TValue, MapperSpec<TConfig> > >::Type TSize;

		TSize undefinedPos = 0;
		TSize offset = buf.pageNo;
		TSize dstPos;
	    offset *= (unsigned)pageSize(buf);

		typename TConfig::Map M = me.pool.handlerArgs;
		bool partiallyFilled = me.pool._partiallyFilled;

		if (partiallyFilled)
			undefinedPos = M(me.pool.undefinedValue);

		for(TValue *cur = buf.begin; cur != buf.end; ++cur) {

			dstPos = M(*cur);

			if (partiallyFilled && dstPos == undefinedPos) 
				continue;							// don't move undefined values

			#ifdef SEQAN_DEBUG
                if (!(dstPos >= offset && dstPos < offset + (TSize)pageSize(buf))) {
					::std::cerr << "Mapper assertion failed: " << ::std::hex << dstPos;
					::std::cerr << " not in [" << offset << "," << (offset + pageSize(buf)) << ") at " << (cur - buf.begin);
                    ::std::cerr << " element is " << ::std::dec << *cur << ::std::endl;
                }
			#endif
            SEQAN_ASSERT(dstPos >= offset && dstPos < offset + (TSize)pageSize(buf));

            TValue *I = buf.begin + (dstPos - offset);
            if (I != cur) {
				TValue tmp;
				TValue *refPrev = cur, *refNext = &tmp;

				do {
					*refNext = *I;
					*I = *refPrev;

					dstPos = M(*refNext);

					if (partiallyFilled && dstPos == undefinedPos)
						I = cur;					// move the undefined value to an arbitrary free position (*cur)
					else 
					{
						#ifdef SEQAN_DEBUG
							if (!(dstPos >= offset && dstPos < offset + (TSize)pageSize(buf))) {
								::std::cerr << "Mapper assertion failed: " << ::std::hex << dstPos;
								::std::cerr << " not in [" << offset << "," << (offset + pageSize(buf)) << ") at " << (refNext - buf.begin);
								::std::cerr << " element is " << ::std::dec << *refNext << ::std::endl;
							}
							TValue *oldI = I;
						#endif
						SEQAN_ASSERT(dstPos >= offset && dstPos < offset + (TSize)pageSize(buf));

						I = buf.begin + (dstPos - offset);

						#ifdef SEQAN_DEBUG
							if (!partiallyFilled && I < cur) {
								::std::cerr << "Mapper assertion failed: I=" << ::std::hex << I;
								::std::cerr << " < cur=" << cur << ::std::dec << ::std::endl;
								break;
							}
							if (I == oldI) {
								::std::cerr << "Mapper assertion failed: I=" << ::std::hex << I;
								::std::cerr << " in endless loop" << ::std::dec << ::std::endl;
								break;
							}
						#endif
					}

					TValue *_swap = refNext;
					refNext = refPrev;
					refPrev = _swap;
				} while (I != cur);

				if (refPrev != cur)
					*cur = *refPrev;
            }
        }
        return buf;
    }

/*
    // inherits buffered file reader and uses a map buffer
    template < typename TPool >
    struct BufferHandler< TPool, ReadMapperSpec >: public BufferHandler< TPool, ReadFileSpec > 
    {
        typedef BufferHandler< TPool, ReadFileSpec >    Base;
        typedef typename Base::Type                     Type;
        typedef typename Base::Buffer					Buffer;

        Buffer	mapBuffer;

        BufferHandler(TPool &_pool):
            Base(_pool)
        {
            allocPage(mapBuffer, _pool.pageSize, *this);
        }

        BufferHandler(TPool &_pool, unsigned _requestedBufferSize, unsigned _readAheadBuffers = 1):
            Base(_pool, _requestedBufferSize, _readAheadBuffers)
        {
            allocPage(mapBuffer, alignSize(_min((_pool.size(), _requestedBufferSize), _pool.pageSize), *this);
        }

        ~BufferHandler()
        {
            freePage(mapBuffer, *this);
        }
    };
*/

    //////////////////////////////////////////////////////////////////////////////
	// generic adapter for buffered memory writers
    struct MapperMemAdapter;

	template < typename TBufferHandler >
	struct Handler< TBufferHandler, MapperMemAdapter >
    {
		typedef typename TBufferHandler::Pool	Pool;
        typedef typename TBufferHandler::Type   Type;
        typedef typename TBufferHandler::Buffer	Buffer;

		Pool			&pool;
        TBufferHandler  handler;
        Buffer			buffer;

        template < typename TPool >
        Handler(TPool &_pool):
			pool(_pool),
            handler(_pool) { }

        inline void _initializeBuffer() {
            if (empty(buffer)) return;
            arrayFill(seqan::begin(buffer, Standard()), seqan::end(buffer, Standard()), pool.undefinedValue);
        }
        
        inline bool begin() {
            buffer = handler.first();
            _initializeBuffer();
            return buffer.begin != NULL;
        }

        inline void push(Type const & Val_) {
            buffer[pool.handlerArgs(Val_)] = Val_;
        }

        inline void end() {
            handler.end();
        }

		inline bool eof() { return false; }
		inline void process() {}
    };





    //////////////////////////////////////////////////////////////////////////////
	// cache bucket based synchronous write handler
    struct MapperSyncWriter;

	template < typename TPool >
	struct Handler< TPool, MapperSyncWriter >
    {
        typedef typename TPool::Type		Type;
        typedef typename TPool::Buffer		Buffer;
        typedef PageBucket<Type>            TPageBucket;
        typedef ::std::vector<TPageBucket>	Cache;

        TPool   &pool;
        Buffer	bucketBuffer;
        Cache   cache;

        Handler(TPool &_pool):
            pool(_pool) { }

        ~Handler() {
            cancel();
        }
        
		struct insertBucket : public ::std::unary_function<TPageBucket,void> {
			Handler &me;
			insertBucket(Handler &_me): me(_me) {}

			inline void operator() (TPageBucket const &cb) const {
                me.cache.push_back(cb);
			}
		};

        bool begin() {
    		cache.reserve(pool.pages);
            return equiDistantDistribution(
                bucketBuffer, pool.bucketBufferSize, *this,
                pool._size, pool.pageSize,
                insertBucket(*this));
        }

        inline void push(Type const &item) {
			unsigned pageNo = pool.handlerArgs(item) / pool.pageSize;
            #ifdef SEQAN_DEBUG
                if (!(pageNo < cache.size())) {
					::std::cerr << "Mapper push assertion failed: " << pageNo << " >= " << cache.size();
                    ::std::cerr << " element is " << item << ::std::endl;
                }
			#endif
    		SEQAN_ASSERT_LT(pageNo, cache.size());
			TPageBucket &cb = cache[pageNo];

			*cb.cur = item;
			if (++cb.cur == cb.end)
				writeBucket(cb, pageNo, pool.pageSize, pool.file);
        }

        inline void end() {
            // flush all cache buckets to disk and compact cache
			pool._partiallyFilled = false;
			unsigned pageNo = 0;
			unsigned curOfs, endOfs;
			for(typename Cache::iterator cb = cache.begin(); cb != cache.end(); ++cb, ++pageNo) {
				curOfs = cb->pageOfs + (cb->cur - cb->begin);
				endOfs = pool.dataSize(pageNo);
				if (curOfs < endOfs) {
					// this page is partially filled and needs to be filled up with undefined entries
					for(unsigned i = curOfs; i < endOfs; ++i) {
						*cb.cur = pool.undefinedValue;
						if (++cb.cur == cb.end)
							writeBucket(cb, pageNo, pool.pageSize, pool.file);
					}
					pool._partiallyFilled = true;
				}
				writeBucket(*cb, pageNo, pool.pageSize, pool.file);
			}
			flush(pool.file);
            cancel();
        }

        inline void cancel() {
            cache.clear();
			cache.reserve(0);
            freePage(bucketBuffer, *this);
        }
        
		inline bool eof() { return false; }
        inline void process() {}
	};


    //////////////////////////////////////////////////////////////////////////////
	// cache bucket based synchronous write handler
    struct MapperAsyncWriter;

	template < typename TPool >
	struct Handler< TPool, MapperAsyncWriter >
    {
        typedef typename TPool::Type                Type;
        typedef typename TPool::File                File;

        typedef SimpleBuffer<Type>					Buffer;
        typedef PageFrame<Type, File, Dynamic<> >   TPageFrame;
        typedef PageBucket<Type>                    TPageBucket;
        typedef ::std::vector<TPageBucket>	        Cache;
        typedef PageChain<TPageFrame>	            TPageChain;

        TPool       &pool;
        Buffer		bucketBuffer;
        TPageChain	chain;
        Buffer		writeCache;
        Cache       cache;
        unsigned    clusterSize;

        Handler(TPool &_pool):
            pool(_pool),
            chain(_pool.writeBackBuckets) {}

        Handler(TPool &_pool, unsigned _writeBackBuckets):
            pool(_pool),
            chain(_writeBackBuckets) {}

        ~Handler() {
            cancel();
        }
        
		struct insertBucket : public ::std::unary_function<TPageBucket,void> {
			Handler &me;
			insertBucket(Handler &_me): me(_me) {}

			inline void operator() (TPageBucket const &cb) const {
                me.cache.push_back(cb);
			}
		};

        bool begin() {
    		cache.reserve(pool.pages());
            clusterSize = equiDistantAlignedDistribution(
                bucketBuffer, sectorSize(pool.file), pool.bucketBufferSize, pool.file,
                pool._size, pool.pageSize,
                insertBucket(*this));

			if (clusterSize == 0) {
				clusterSize = UINT_MAX;
				#ifdef SEQAN_DEBUG
					::std::cerr << "mapper switched to synchronous mode" << ::std::endl;
				#endif
				return equiDistantDistribution(
					bucketBuffer, pool.bucketBufferSize, pool.file,
					pool._size, pool.pageSize,
					insertBucket(*this));
			}

			#ifdef SEQAN_VERBOSE
				::std::cerr << "async mapper clustersize " << clusterSize << ::std::endl;
			#endif
            allocPage(writeCache, chain.maxFrames * clusterSize, pool.file);

            // distribute write back buffers
            Type *cur = writeCache.begin;
            TPageFrame *p = chain.first;
            while (p) {
                p->begin = cur; cur += clusterSize;
                p->end = cur;
                setPageSize(*p, pool.pageSize);
                p = p->next;
            }
            return true;
        }

        inline void push(Type const &item) {
			unsigned pageNo = pool.handlerArgs(item) / pool.pageSize;
            #ifdef SEQAN_DEBUG
                if (!(pageNo < cache.size())) {
					::std::cerr << "Mapper push assertion failed: " << pageNo << " >= " << cache.size();
                    ::std::cerr << " element is " << item << ::std::endl;
					::std::cerr << ::std::hex << pool.handlerArgs(item) << " / " << pool.pageSize;
					::std::cerr << " = " << pageNo << ::std::dec << ::std::endl;
                }
			#endif
            SEQAN_ASSERT_LT(pageNo, cache.size());
			TPageBucket &cb = cache[pageNo];

			*cb.cur = item;
			if (++cb.cur == cb.end)
				_writeBucket(cb, pageNo);
        }

        inline void end() {
            // flush all cache buckets to disk and compact cache
			pool._partiallyFilled = false;
			unsigned pageNo = 0;
			unsigned curOfs, endOfs;
			for(typename Cache::iterator cb = cache.begin(); cb != cache.end(); ++cb, ++pageNo) {
				curOfs = cb->pageOfs + (cb->cur - cb->begin);
				endOfs = pool.dataSize(pageNo);
				if (curOfs < endOfs) {
					// this page is partially filled and needs to be filled up with undefined entries
					for(unsigned i = curOfs; i < endOfs; ++i) {
						*cb->cur = pool.undefinedValue;
						if (++cb->cur == cb->end)
							_writeBucket(*cb, pageNo);
					}
					pool._partiallyFilled = true;
				} else
					_writeBucket(*cb, pageNo);
			}
            chain.waitForAll();
			flush(pool.file);
            cancel();
        }

        inline void cancel() {
            chain.cancelAll(pool.file);
            cache.clear();
			cache.reserve(0);
            freePage(writeCache, pool.file);
            freePage(bucketBuffer, pool.file);
        }
        
		inline bool eof() { return false; }
        inline void process() {}

    protected:

        inline void _swap(TPageFrame &pf, TPageBucket &pb) {
            Type *tmp = pf.begin;
            pf.begin = pb.begin;
            pf.end   = pb.cur;

            pb.begin = tmp;
            pb.cur   = tmp;
            pb.end   = tmp + clusterSize;
        }

        bool _writeBucket(TPageBucket &cb, unsigned pageNo) {
            if ((unsigned)(cb.cur - cb.begin) != clusterSize)
                return writeBucket(cb, pageNo, pool.pageSize, pool.file);

            TPageFrame *pf = chain.getReadyPage();
            _swap(*pf, cb);
            pf->pageNo = pageNo;
            return writeBucket(*pf, cb.pageOfs, pool.file);
        }
	};


	//////////////////////////////////////////////////////////////////////////////
	// character and buffer based handler definitions
	template < typename TValue,
			   typename TConfig >
    struct BufReadHandler< Pool< TValue, MapperSpec<TConfig> > >
    {
        typedef BufferHandler< Bundle2<
			BufferHandler< Pool< TValue, MapperSpec<TConfig> >, MemorySpec >,
//			BufferHandler< Pool< TValue, MapperSpec<TConfig> >, ReadMapperSpec >
			BufferHandler< Pool< TValue, MapperSpec<TConfig> >, ReadFileSpec >
		>, MultiplexSpec > Type;
    };

    template < typename TValue,
			   typename TConfig >
    struct WriteHandler< Pool< TValue, MapperSpec<TConfig> > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler < Pool< TValue, MapperSpec<TConfig> >, MemorySpec >, MapperMemAdapter >,
			Handler< Pool< TValue, MapperSpec<TConfig> >, MapperSyncWriter >
		>, MultiplexSpec > Type;
	};


	// use async MapperHandler for all async files classes

    template < typename TValue,
   	           typename TMap,
		       typename TSize,
			   typename TSpec >
    struct WriteHandler< Pool< TValue, MapperSpec< MapperConfigSize< TMap, TSize, File<Async<TSpec> > > > > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler	< Pool< TValue, MapperSpec<MapperConfigSize< TMap, TSize, File<Async<TSpec> > > > >, MemorySpec >, MapperMemAdapter >,
			Handler< Pool< TValue, MapperSpec<MapperConfigSize< TMap, TSize, File<Async<TSpec> > > > >, MapperAsyncWriter >
		>, MultiplexSpec > Type;
	};

    template < typename TValue,
   	           typename TMap,
			   typename TSpec >
    struct WriteHandler< Pool< TValue, MapperSpec< MapperConfig< TMap, File<Async<TSpec> > > > > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler < Pool< TValue, MapperSpec<MapperConfig< TMap, File<Async<TSpec> > > > >, MemorySpec >, MapperMemAdapter >,
			Handler< Pool< TValue, MapperSpec<MapperConfig< TMap, File<Async<TSpec> > > > >, MapperAsyncWriter >
		>, MultiplexSpec > Type;
	};


	// async file arrays

    template < typename TValue,
   	           typename TMap,
		       typename TSize,
			   __int64  FileSize_,
			   typename TSpec >
    struct WriteHandler< Pool< TValue, MapperSpec< MapperConfigSize< TMap, TSize, File< Chained<FileSize_, File< Async<TSpec> > > > > > > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler	< Pool< TValue, MapperSpec<MapperConfigSize< TMap, TSize, File< Chained<FileSize_, File< Async<TSpec> > > > > > >, MemorySpec >, MapperMemAdapter >,
			Handler< Pool< TValue, MapperSpec<MapperConfigSize< TMap, TSize, File< Chained<FileSize_, File< Async<TSpec> > > > > > >, MapperAsyncWriter >
		>, MultiplexSpec > Type;
	};

    template < typename TValue,
   	           typename TMap,
			   __int64  FileSize_,
			   typename TSpec >
    struct WriteHandler< Pool< TValue, MapperSpec< MapperConfig< TMap, File< Chained< FileSize_, File< Async<TSpec> > > > > > > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler < Pool< TValue, MapperSpec<MapperConfig< TMap, File< Chained<FileSize_, File< Async<TSpec> > > > > > >, MemorySpec >, MapperMemAdapter >,
			Handler< Pool< TValue, MapperSpec<MapperConfig< TMap, File< Chained<FileSize_, File< Async<TSpec> > > > > > >, MapperAsyncWriter >
		>, MultiplexSpec > Type;
	};

    template < typename TValue,
   	           typename TMap,
		       typename TSize,
			   unsigned FileCount_,
			   typename TSpec >
    struct WriteHandler< Pool< TValue, MapperSpec< MapperConfigSize< TMap, TSize, File< Striped<FileCount_, File< Async<TSpec> > > > > > > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler	< Pool< TValue, MapperSpec<MapperConfigSize< TMap, TSize, File< Striped<FileCount_, File< Async<TSpec> > > > > > >, MemorySpec >, MapperMemAdapter >,
			Handler< Pool< TValue, MapperSpec<MapperConfigSize< TMap, TSize, File< Striped<FileCount_, File< Async<TSpec> > > > > > >, MapperAsyncWriter >
		>, MultiplexSpec > Type;
	};

    template < typename TValue,
   	           typename TMap,
			   unsigned FileCount_,
			   typename TSpec >
    struct WriteHandler< Pool< TValue, MapperSpec< MapperConfig< TMap, File< Striped<FileCount_, File< Async<TSpec> > > > > > > >
    {
        typedef Handler< Bundle2<
			Handler< BufferHandler < Pool< TValue, MapperSpec<MapperConfig< TMap, File< Striped<FileCount_, File< Async<TSpec> > > > > > >, MemorySpec >, MapperMemAdapter >,
			Handler< Pool< TValue, MapperSpec<MapperConfig< TMap, File< Striped<FileCount_, File< Async<TSpec> > > > > > >, MapperAsyncWriter >
		>, MultiplexSpec > Type;
	};


/*
    // WARNING:
    // Using template based polymorphy with oop inheritance causes
    // funny side effects, i.e. the VS compiler can't find the right
    // global template function/structure implicitly. And instead
	// of redefining every Pool structure for Mapper/Sorter/Pool
	// we now only use Pool and rename it to Pool


    template < typename TValue,
               typename TMap,
			   typename TFile = File<> >
    struct Mapper: public Pool< TValue, MapperSpec< MapperConfig<TMap, TFile> > >
    {
        typedef Pool< TValue, MapperSpec< MapperConfig<TMap, TFile> > > Base;

        typedef TValue	                    Type;
        typedef TFile                       File;
        typedef typename Size<TFile>::Type  SizeType;

        Mapper(TMap const & M)
        {
            Base::handlerArgs = M;
        }

        template < typename TInput, typename TSpec >
        Mapper(Pipe<TInput, TSpec> &src, TMap const & M):
            Base(src)
        {
            Base::handlerArgs = M;
        }

    };
*/
}

#endif
