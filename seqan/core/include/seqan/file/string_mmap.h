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

#ifndef SEQAN_HEADER_STRING_MMAP_H
#define SEQAN_HEADER_STRING_MMAP_H


/* IOREV
 * _tested_
 * _windows_
 *
 *
 * tested in library/demos/howto/efficiently_import_sequences.cpp and stellar
 *
 * relation to file_format_mmap.h unclear
 *
 * relation to string_external unclear, what benifit does string_mmap provide?
 *
 * contains windows-specific code
 */


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

    template < typename TFile_ = File<>,				// default file type
               typename TSize_ = size_t >				// size type
    struct MMapConfig {
//IOREV _nodoc_ doc says using MMap<ExternalConfig> is correct, whats this for?
        typedef TFile_ TFile;
        typedef TSize_ TSize;
    };

    template < typename TConfig = MMapConfig<> >
    struct MMap {};
//IOREV
	
	
	//////////////////////////////////////////////////////////////////////////////
    // Memory Mapped String
    //////////////////////////////////////////////////////////////////////////////

/**
.Spec.MMap String:
..cat:Strings
..general:Class.String
..summary:String that is stored in external memory. Uses memory mapping.
..signature:String<TValue, MMap<> >
..signature:String<TValue, MMap<TConfig> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.TConfig:A structure to configure the memory mapped string.
...type:Tag.ExternalConfig
...type:Tag.ExternalConfigLarge
...type:Tag.ExternalConfigSize
...default:@Tag.ExternalConfig@
..remarks:The MMap String enables to access sequences larger than the available physical memory (RAM) by using
external memory (e.g. Hard disk, Network storage, ...) mapped into memory.
The size of the string is limited by external memory and the logical address space (4GB on a 32bit OS).
See the @Memfunc.ExtString#String.constructor@ for more details.
..remarks:This String also supports fast appending and removing of values at the end (see @Spec.Block String@, @Function.appendValue@)
..include:seqan/file.h
*/

/**
.Enum.MMapAdviseScheme
..cat:Sequences
..cat:Input / Output
..summary:Enum with mmap advise values.
..value.MMAP_NORMAL:There is no advise on the given address range.
..value.MMAP_RANDOM:The address range will be accessed with a random access memory pattern.
..value.MMAP_SEQUENTIAL:The address range will be accessed sequentially.
..value.MMAP_WILLNEED:The address range in the advise will be needed in the future.
..value.MMAP_DONTNEED:The address range in the advise will not be needed any more.
..see:Function.mmapAdvise
..include:seqan/file.h
 */

#ifdef PLATFORM_WINDOWS

		enum MMapAdviseScheme {
			MMAP_NORMAL = 0,
			MMAP_RANDOM = 0,
			MMAP_SEQUENTIAL = 0,
			MMAP_WILLNEED = 0,
			MMAP_DONTNEED = 0
		};

#else

		enum MMapAdviseScheme {
			MMAP_NORMAL = POSIX_MADV_NORMAL,
			MMAP_RANDOM = POSIX_MADV_RANDOM,
			MMAP_SEQUENTIAL = POSIX_MADV_SEQUENTIAL,
			MMAP_WILLNEED = POSIX_MADV_WILLNEED,
			MMAP_DONTNEED = POSIX_MADV_DONTNEED
		};

#endif

    template < typename TValue,
               typename TConfig >
	class String<TValue, MMap<TConfig> >
	{
//IOREV
	public:

        typedef typename TConfig::TFile		TFile;
        typedef typename TConfig::TSize		TSize;

		TValue				*data_begin;
		TValue				*data_end;
		TSize				data_capacity;

		TFile				file;
		int					_openMode;
        bool                _temporary, _ownFile;

#ifdef PLATFORM_WINDOWS
        HANDLE				handle;
#endif
		MMapAdviseScheme	scheme;

        // TODO(holtgrew): Explicit?
		String(TSize size = 0):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL),
			scheme(MMAP_NORMAL)
        {
            _temporary = true;
            _ownFile = false;

			resize(*this, size);
        }

        // TODO(holtgrew): Explicit?
		String(TFile &_file):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL),
			scheme(MMAP_NORMAL)
        {
			open(*this, _file);
        }

		explicit String(const char *fileName, int openMode = DefaultOpenMode<TFile>::VALUE):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL),
			scheme(MMAP_NORMAL)
        {
			open(*this, fileName, openMode);
        }

		template <typename TSource>
		String & operator =(TSource const & source)
		{
	SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}
		String & operator =(String const & source)
		{
	SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}

		~String() 
		{
			close(*this);
		}

//____________________________________________________________________________

		template <typename TPos>
		inline typename Reference<String>::Type
		operator [] (TPos pos)
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<String const>::Type 
		operator [] (TPos pos) const
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

//____________________________________________________________________________

        inline operator bool() 
   	{
            return file;
        }

//____________________________________________________________________________

};

   template < typename TValue, typename TConfig >
	inline typename Iterator<String<TValue, MMap<TConfig> > , Standard>::Type
	begin(String<TValue, MMap<TConfig> > & me,
		Standard)
	{
//IOREV
SEQAN_CHECKPOINT
		return me.data_begin;
	}
   template < typename TValue, typename TConfig >
	inline typename Iterator<String<TValue, MMap<TConfig> >  const, Standard>::Type
	begin(String<TValue, MMap<TConfig> > const & me,
		Standard)
	{
//IOREV
SEQAN_CHECKPOINT
		return me.data_begin;
	}

//____________________________________________________________________________

   template < typename TValue, typename TConfig >
   inline typename Iterator<String<TValue, MMap<TConfig> > , Standard>::Type
	end(String<TValue, MMap<TConfig> > & me,
		Standard)
	{
//IOREV
SEQAN_CHECKPOINT
		return me.data_end;
	}
   template < typename TValue, typename TConfig >
	inline typename Iterator<String<TValue, MMap<TConfig> >  const, Standard>::Type
	end(String<TValue, MMap<TConfig> > const & me,
		Standard)
	{
//IOREV
SEQAN_CHECKPOINT
		return me.data_end;
	}

//____________________________________________________________________________

   template < typename TValue, typename TConfig >
	inline typename Size<String<TValue, MMap<TConfig> > >::Type
	capacity(String<TValue, MMap<TConfig> > & me) 
	{
//IOREV
SEQAN_CHECKPOINT
		return me.data_capacity;
	}

   template < typename TValue, typename TConfig >
	inline typename Size<String<TValue, MMap<TConfig> > >::Type
	capacity(String<TValue, MMap<TConfig> > const & me) 
	{
//IOREV
SEQAN_CHECKPOINT
		return me.data_capacity;
	}

//____________________________________________________________________________

   template < typename TValue, typename TConfig >
	inline void 
	_setLength(
		String<TValue, MMap<TConfig> > & me, 
		size_t new_length)
	{
//IOREV
SEQAN_CHECKPOINT
		me.data_end = me.data_begin + new_length;
	}

//____________________________________________________________________________

	template < typename TValue, typename TConfig >
   inline void 
	_setCapacity(
		String<TValue, MMap<TConfig> > & me, 
		size_t new_capacity)
	{
//IOREV
SEQAN_CHECKPOINT
		me.data_capacity = new_capacity;
	}


    //////////////////////////////////////////////////////////////////////////////
    // meta-function interface

    template < typename TValue, typename TConfig >
    struct Size< String<TValue, MMap<TConfig> > >
    {
//IOREV
        typedef typename TConfig::TSize Type;
    };

    template < typename TValue, typename TConfig >
    struct Difference< String<TValue, MMap<TConfig> > >
    {
//IOREV
		typedef typename MakeSigned_<typename TConfig::TSize>::Type Type;
    };
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	struct DefaultOverflowExplicit<String<TValue, MMap<TConfig> > >
	{
//IOREV
		typedef Generous Type;
	};

    template < typename TValue, typename TConfig >
	struct DefaultOverflowImplicit<String<TValue, MMap<TConfig> > >
	{
//IOREV
		typedef Generous Type;
	};
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	struct IsContiguous< String<TValue, MMap<TConfig> > >
	{
//IOREV
		typedef True Type;
		enum { VALUE = true };
	};

    template < typename TValue, typename TConfig >
	struct AllowsFastRandomAccess< String<TValue, MMap<TConfig> > >
	{
//IOREV
		typedef False Type;
		enum { VALUE = false };
	};


	//////////////////////////////////////////////////////////////////////////////
    // global interface

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	inline void 
	waitForAll(String<TValue, MMap<TConfig> > &)
	{
//IOREV _nodoc_
	}
	
#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES MMapStringDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    template < typename TValue, typename TConfig, typename TScheme, typename TBeginPos, typename TEndPos >
	inline int
	mmapAdvise(String<TValue, MMap<TConfig> > &/*me*/, TScheme const & /*scheme*/, TBeginPos const &, TEndPos const &)
	{
//IOREV  _nodoc_ _windows_
		return 0;
	}
		
    template < typename TValue, typename TConfig >
    inline void 
	flush(String<TValue, MMap<TConfig> > &) 
	{
//IOREV _nodoc_ _windows_
    }

	// cancel all transactions
    template < typename TValue, typename TConfig >
	inline void 
	cancel(String<TValue, MMap<TConfig> > &)
	{
//IOREV _windows_
	}

	// flush and free all allocated pages
    template < typename TValue, typename TConfig >
	inline int
	flushAndFree(String<TValue, MMap<TConfig> > &)
	{
//IOREV _nodoc_ _windows_
		return 0;
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    _map(String<TValue, MMap<TConfig> > &me, size_t new_capacity) 
	{
//IOREV _windows_
		if (new_capacity > 0) 
		{
			resize(me.file, new_capacity * sizeof(TValue));
			DWORD prot = 0;
			DWORD access = 0;
			if ((me._openMode & OPEN_MASK) == OPEN_RDONLY) 
			{
				prot = PAGE_READONLY;
				access = FILE_MAP_READ;
			} else {
				prot = PAGE_READWRITE;
				access = FILE_MAP_ALL_ACCESS;
			}
            LARGE_INTEGER largeSize;
			largeSize.QuadPart = new_capacity;
			largeSize.QuadPart *= sizeof(TValue);

			me.handle = CreateFileMapping(me.file.handle, &MMapStringDefaultAttributes, prot, largeSize.HighPart, largeSize.LowPart, NULL);
			if (me.handle == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "CreateFileMapping failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}

			void *addr = MapViewOfFile(me.handle, access, 0, 0, 0);	
			if (addr == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "MapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}
				
			me.data_begin = (TValue *) addr;
			_setLength(me, new_capacity);
			_setCapacity(me, new_capacity);
		}
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    _unmap(String<TValue, MMap<TConfig> > &me) 
	{
//IOREV _windows_
		bool result = true;
		if (me.data_begin) 
		{
			if (UnmapViewOfFile(me.data_begin) == 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "UnmapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				result = false;
			}
			
			if (CloseHandle(me.handle) == 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "CloseHandle failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				result = false;
			}
			
			resize(me.file, length(me) * sizeof(TValue));
			me.data_begin = NULL;
			me.data_end = NULL;
			me.data_capacity = 0;
		}
		return result;
	}

	template < typename TValue, typename TConfig, typename TCapSize >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TCapSize new_capacity) 
	{
//IOREV _windows_
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;

		if (me.data_begin) 
		{
			if (new_capacity > 0) 
			{
				TSize seq_length = length(me);
				
				DWORD prot = 0;
				DWORD access = 0;
				if ((me._openMode & OPEN_MASK) == OPEN_RDONLY) 
				{
					prot = PAGE_READONLY;
					access = FILE_MAP_READ;
				} else {
					prot = PAGE_READWRITE;
					access = FILE_MAP_ALL_ACCESS;
				}
				LARGE_INTEGER largeSize;
				largeSize.QuadPart = new_capacity;
				largeSize.QuadPart *= sizeof(TValue);

				bool result = true;
				result &= (UnmapViewOfFile(me.data_begin) != 0);
				result &= (CloseHandle(me.handle) != 0);

				HANDLE handle = CreateFileMapping(me.file.handle, &MMapStringDefaultAttributes, prot, largeSize.HighPart, largeSize.LowPart, NULL);
				if (handle == NULL)
				{
				#ifdef SEQAN_DEBUG
					::std::cerr << "CreateFileMapping failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
				#endif
					return false;
				}

				void *addr = MapViewOfFile(handle, access, 0, 0, 0);	
				if (addr == NULL)
				{
				#ifdef SEQAN_DEBUG
					::std::cerr << "MapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
				#endif
					return false;
				}

				if (capacity(me) > new_capacity)
					resize(me.file, new_capacity * sizeof(TValue));

				me.handle = handle;
				me.data_begin = (TValue*) addr;
				_setLength(me, seq_length);
				_setCapacity(me, new_capacity);
				return true;
			} else
				return _unmap(me);
		} else
			return _allocateStorage(me, new_capacity) != NULL;
	}

#else

///.Function.flush.param.string.type:Spec.MMap String
///.Function.flush.class:Spec.MMap String

/**
.Function.mmapAdvise
..class:Spec.MMap String
..cat:Sequences
..summary:Call advise function for memory mapped files.
..signature:mmapAdvise(mmapString, scheme, beginPos, endPos)
..param.mmapString:The @Spec.MMap String@ that contains the location in the advise call.
...type:Spec.MMap String
..param.scheme:The memory access scheme to use.
...type:Enum.MMapAdviseScheme
..param.beginPos:Begin position in the string for the advise call.
..param.endPos:End position in the string for the advise call.
..returns:$int$, return code 0 on success.
..see:Enum.MMapAdviseScheme
..include:seqan/file.h
 */
    template < typename TValue, typename TConfig, typename TScheme, typename TBeginPos, typename TEndPos >
	inline int
	mmapAdvise(String<TValue, MMap<TConfig> > &me, TScheme const & scheme, TBeginPos const & beginPos, TEndPos const & endPos)
	{
//IOREV _nodoc_
		me.scheme = scheme;
//		posix_fadvise(me.file.handle, beginPos * sizeof(TValue), (endPos - beginPos) * sizeof(TValue), scheme);
		if (scheme == MMAP_DONTNEED)
			msync(me.data_begin + beginPos, (endPos - beginPos) * sizeof(TValue), MS_INVALIDATE);
		return posix_madvise(me.data_begin + beginPos, (endPos - beginPos) * sizeof(TValue), scheme);
	}
		
    template < typename TValue, typename TConfig >
    inline void 
	flush(String<TValue, MMap<TConfig> > &me) 
	{
//IOREV
		msync(me.data_begin, length(me) * sizeof(TValue), MS_SYNC);
    }

	// cancel all transactions
    template < typename TValue, typename TConfig >
	inline void 
	cancel(String<TValue, MMap<TConfig> > &me)
	{
//IOREV
		msync(me.data_begin, capacity(me) * sizeof(TValue), MS_INVALIDATE);
	}

	// flush and free all allocated pages
    template < typename TValue, typename TConfig >
	inline int
	flushAndFree(String<TValue, MMap<TConfig> > &me)
	{
//IOREV _nodoc_
		return posix_madvise(me.data_begin, capacity(me) * sizeof(TValue), MADV_DONTNEED);
	}
	
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    _map(String<TValue, MMap<TConfig> > &me, size_t new_capacity) 
	{
//IOREV
		if (new_capacity > 0) 
		{
			_ensureFileIsOpen(me);
			resize(me.file, new_capacity * sizeof(TValue));
			int prot = 0;
			if (me._openMode & OPEN_RDONLY) prot |= PROT_READ;
			if (me._openMode & OPEN_WRONLY) prot |= PROT_WRITE;
			void *addr = mmap(NULL, new_capacity * sizeof(TValue), prot, MAP_SHARED, me.file.handle, 0);
			mmapAdvise(me, me.scheme);
			
			if (addr == MAP_FAILED)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "mmap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
			#endif
				return false;
			}
				
			me.data_begin = (TValue *) addr;
			_setLength(me, new_capacity);
			_setCapacity(me, new_capacity);
		}
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    _unmap(String<TValue, MMap<TConfig> > &me) 
	{
//IOREV
		if (me.data_begin) 
		{
			int error = munmap(me.data_begin, capacity(me) * sizeof(TValue));
			if (error != 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "munmap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
			#endif
				return false;
			}
			
			resize(me.file, length(me) * sizeof(TValue));
			me.data_begin = NULL;
			me.data_end = NULL;
			me.data_capacity = 0;
		}
		return true;
	}

	template < typename TValue, typename TConfig, typename TCapSize >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TCapSize new_capacity) 
	{
//IOREV
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;

		if (me.data_begin)
		{
			if (new_capacity > 0) 
			{
				TSize seq_length = length(me);
				
				_ensureFileIsOpen(me);
				if (capacity(me) < new_capacity)
					resize(me.file, new_capacity * sizeof(TValue));

#ifdef MREMAP_MAYMOVE
				void *addr = mremap(me.data_begin, capacity(me) * sizeof(TValue), new_capacity * sizeof(TValue), MREMAP_MAYMOVE);
				mmapAdvise(me, me.scheme);
#else
				// for BSD systems without mremap(..) like Mac OS X ...
				int prot = 0;
				if (me._openMode & OPEN_RDONLY) prot |= PROT_READ;
				if (me._openMode & OPEN_WRONLY) prot |= PROT_WRITE;
	//			void *addr = mmap(me.data_begin, new_capacity * sizeof(TValue), prot, MAP_SHARED | MAP_FIXED, me.file.handle, 0);
				munmap(me.data_begin, capacity(me) * sizeof(TValue));
				void *addr = mmap(NULL, new_capacity * sizeof(TValue), prot, MAP_SHARED, me.file.handle, 0);
				mmapAdvise(me, me.scheme);
#endif

				if (addr == MAP_FAILED) 
				{
				#ifdef SEQAN_DEBUG
					::std::cerr << "mremap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
				#endif
					return false;
				}

				if (capacity(me) > new_capacity)
					resize(me.file, new_capacity * sizeof(TValue));

				me.data_begin = (TValue*) addr;
				_setLength(me, seq_length);
				_setCapacity(me, new_capacity);
				return true;
			} else
				return _unmap(me);
		} else
			return _map(me, new_capacity);
	}

#endif

	template < typename TValue, typename TConfig >
    inline void 
    clear(String<TValue, MMap<TConfig> > &me) 
	{
//IOREV
		cancel(me);
		_unmap(me);
		resize(me.file, 0);
	}
//____________________________________________________________________________

    template < typename TValue, typename TConfig, typename TScheme >
	inline int
	mmapAdvise(String<TValue, MMap<TConfig> > &me, TScheme scheme)
	{
//IOREV _nodoc_
		return mmapAdvise(me, scheme, 0, capacity(me));
	}
		
//____________________________________________________________________________

	template < typename TValue, typename TConfig, typename TSize >
    inline typename Value<String<TValue, MMap<TConfig> > >::Type * 
    _allocateStorage(String<TValue, MMap<TConfig> > &me, TSize new_capacity) 
	{
//IOREV
		TSize size = _computeSizeForCapacity(me, new_capacity);
		_map(me, size);
		return NULL;
	}

	template < typename TValue, typename TConfig, typename TSize >
    inline typename Value<String<TValue, MMap<TConfig> > >::Type * 
    _reallocateStorage(
		String<TValue, MMap<TConfig> > &me,
		TSize new_capacity) 
	{
//IOREV
		TSize size = _computeSizeForCapacity(me, new_capacity);
		_remap(me, size);
		return NULL;
	}

	template < typename TValue, typename TConfig, typename TSize >
    inline void
    _deallocateStorage(String<TValue, MMap<TConfig> > &/*me*/, TValue * /*ptr*/, TSize /*capacity*/)
	{
//IOREV
	}
//____________________________________________________________________________
///.Function.open.param.string.type:Spec.MMap String
///.Function.open.class:Spec.MMap String

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, const char *fileName, int openMode) 
	{
//IOREV
		close(me);
		me._temporary = false;
				
		if ((me._ownFile = open(me.file, fileName, openMode))) 
		{
			me._openMode = openMode;
			return _map(me, (size_t)(size(me.file) / sizeof(TValue)));
		}

		return false;
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, const char *fileName) 
	{
//IOREV
		typedef typename String<TValue, MMap<TConfig> >::TFile	TFile;
		return open(me, fileName, DefaultOpenMode<TFile>::VALUE);
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, typename TConfig::TFile file) 
	{
//IOREV
		close(me);
		me.file = file;
        me._temporary = false;
        me._ownFile = false;

		if (me.file) 
		{
			me._openMode = OPEN_RDWR;
			return _map(me, size(me.file) / sizeof(TValue));
		}

		return false;
    }

///.Function.openTemp.param.string.type:Spec.MMap String
///.Function.openTemp.class:Spec.MMap String

	template < typename TValue, typename TConfig >
    inline bool 
    openTemp(String<TValue, MMap<TConfig> > &me) 
	{
//IOREV
		close(me);
        me._temporary = true;
		me._openMode = OPEN_RDWR;
		
		return me._ownFile = openTemp(me.file);
    }
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
	inline void _ensureFileIsOpen(String<TValue, MMap<TConfig> > &me) 
	{
//IOREV
		if (!me.file)
		{
			me._temporary = true;
			me._openMode = OPEN_RDWR;

			if (!(me._ownFile = openTemp(me.file)))
				::std::cerr << "Memory Mapped String couldn't open temporary file" << ::std::endl;
		}
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &/*me*/, const char * /*fileName*/, int /*openMode*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &/*me*/, const char * /*fileName*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &/*me*/, typename TConfig::TFile /*file*/) {
//IOREV _nodoc_ shouldn't we flush here? in case of abnormal termination...
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}
//____________________________________________________________________________
///.Function.flush.close.string.type:Spec.MMap String
///.Function.flush.class:Spec.MMap String

	template < typename TValue, typename TConfig >
    inline bool 
    close(String<TValue, MMap<TConfig> > &me) 
	{
//IOREV shouldn't we cancel/flush in any case, even if not temp?
		if (me.file) 
		{
			// close associated file
			if (me._temporary) 
			{
				me._temporary = false;
				cancel(me);
			}

			_unmap(me);

			if (me._ownFile) 
			{
				me._ownFile = false;
				return close(me.file);
			}
		}
		return true;
    }

//////////////////////////////////////////////////////////////////////////////



} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
