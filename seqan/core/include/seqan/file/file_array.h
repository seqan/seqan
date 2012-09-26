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

#ifndef SEQAN_HEADER_FILE_ARRAY_H
#define SEQAN_HEADER_FILE_ARRAY_H

#include <sstream>
#include <iomanip>

/* IOREV
 * _nottested_
 * _doc_
 * 
 * 
 * Basically contains stuff fro striped and Chained file access
 * (a little) documentation for both is in file_base
 * Chained-tag has no test-case but is used in different demos
 * Striped-tag has no test-case and AFAICT is not used in any app
 * functions used are those of file_cstyle.h not cstream.h
 * not clear why certain things are here and other in file_base
 * 
 * I am unsure whether Striped works at all, there seems to be no open
 * functions or other functions for acutally assigning multiple files
 * to the Striped virtual file
 */

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	//template < __int64 FILE_SIZE = 2*1024*1024*1024-1, typename TFile = File<> >
	//struct Chained;

	//template < unsigned FileCount_ = 2, typename TFile = File<> >
	//struct Striped;


    template < __int64 FILE_SIZE, typename TFile >
    struct Size< File< Chained<FILE_SIZE, TFile> > >
    {
//IOREV 
        typedef __int64 Type;
    };

    template < __int64 FILE_SIZE, typename TFile >
    struct Position< File< Chained<FILE_SIZE, TFile> > >
    {
//IOREV 
        typedef __int64 Type;
    };

    template < __int64 FILE_SIZE, typename TFile >
    struct Difference< File< Chained<FILE_SIZE, TFile> > >
    {
//IOREV 
        typedef __int64 Type;
    };

    template < __int64 FILE_SIZE, typename TFile >
    struct AsyncRequest< File< Chained<FILE_SIZE, TFile> > >
    {
//IOREV 
		typedef typename AsyncRequest<TFile>::Type Type;
    };


    template < unsigned FileCount_, typename TFile >
    struct Size< File< Striped<FileCount_, TFile> > >
    {
//IOREV
        typedef __int64 Type;
    };

    template < unsigned FileCount_, typename TFile >
    struct Position< File< Striped<FileCount_, TFile> > >
    {
//IOREV
        typedef __int64 Type;
    };

    template < unsigned FileCount_, typename TFile >
    struct Difference< File< Striped<FileCount_, TFile> > >
    {
//IOREV
        typedef __int64 Type;
    };

    template < unsigned FileCount_, typename TFile >
    struct AsyncRequest< File< Striped<FileCount_, TFile> > >
    {
//IOREV
		typedef typename AsyncRequest<TFile>::Type Type;
    };


	template < unsigned FileCount_, typename TFile >
	class File< Striped<FileCount_, TFile> >: public Tuple< TFile, FileCount_ > {
//IOREV _stub_ Is this a stub? see header of file_array.h
		File(void * /*dummy = NULL*/) {}	// to be compatible with the FILE*(NULL) constructor
		operator bool() const { return (*this)[0]; }
	};

    template < __int64 FILE_SIZE, typename TFile >
	class File< Chained<FILE_SIZE, TFile> >: public String< TFile > {
//IOREV
		typedef String< TFile > Base;

		::std::string	baseName;
		int				openMode;
		__int64			fileSize;
		bool			temporary;

		File(void * /*dummy = NULL*/) :	// to be compatible with the FILE*(NULL) constructor
			fileSize(0),
			_realign(false) {}

	private:
		
		bool _realign;
	
		template < typename TSize, typename TValue >
		inline void _alignFloor(TSize _size, TValue const *) {
			__int64 alignment = sizeof(TValue) * sectorSize(TFile());
			fileSize = (_size / alignment) * alignment;
		}

		template < typename TSize, typename TValue >
		inline void _alignCeil(TSize _size, TValue const *) {
			__int64 alignment = sizeof(TValue) * sectorSize(TFile());
			fileSize = ((_size + alignment - 1) / alignment) * alignment;
		}

	public:
	
		inline ::std::string getFileName(int i) const { 
			::std::stringstream strm;
			strm << baseName << '.' << ::std::setfill('0') << ::std::setw(3) << i;
			return strm.str();
		}

		inline operator bool() const { 
			return (*this)[0]; 
		}

		inline unsigned fileCount() const {
			return length(*(Base*)this);
		}

		inline TFile& getFile(int fileNo) {
			unsigned _oldFileCount = fileCount();
			if (fileNo > 0 && static_cast<unsigned>(fileNo) >= _oldFileCount) {
				resize(*(Base*)this, fileNo + 1);
				for(unsigned i = _oldFileCount; i <= static_cast<unsigned>(fileNo); ++i)  // Cast OK since fileNo > 0 checked above.
					if (temporary)
						openTemp((*this)[i], openMode);
					else
						open((*this)[i], getFileName(i).c_str(), openMode);
			}
			return (*this)[fileNo];
		}

		inline void tryOpen() {
			unsigned fileCount = 0;
			while (fileExists(getFileName(fileCount).c_str())) ++fileCount;
			if (fileCount) {
				fileSize = size(getFile(0));
				_realign = (fileCount == 1);
				getFile(fileCount - 1);
			} 
		}

		// fileSize dependent functions

		template < typename TValue >
		inline void adjustFileSize(TValue const *dummy) {
			if (_realign) {
				_alignCeil(fileSize, dummy);
				_realign = false;
				if (fileSize < FILE_SIZE)
					fileSize = 0;
			}
			if (!fileSize) 
				_alignFloor(FILE_SIZE, dummy);
		}

		template < typename TPos, typename TOffset, typename TValue >
		inline TFile& getFileAndOffset(TPos offset, TOffset &fileOffset, TValue const *dummy) {
			adjustFileSize(dummy);
			offset *= sizeof(TValue);
			fileOffset = (offset % fileSize) / sizeof(TValue);
			return getFile(offset / fileSize);
		}

		template < typename TOffset, typename TValue >
		inline __int64 restAt(TOffset fileOffset, TValue const *dummy) {
			adjustFileSize(dummy);
			__int64 restBytes = fileSize;
			restBytes -= fileOffset * sizeof(TValue);
			return restBytes / sizeof(TValue);
		}

		inline void resizeArray(__int64 _newSize) {
			if (fileSize) {
				unsigned _oldFileCount = fileCount();
				unsigned _newFileCount = enclosingBlocks(_newSize, fileSize);
				for(unsigned i = _newFileCount; i < _oldFileCount; ++i) {
					close((*this)[i]);
					if (!temporary) fileUnlink(getFileName(i).c_str());
				}
				resize(*(Base*)this, _newFileCount);
				if (_newFileCount) {
					typename Size<TFile>::Type lastFileSize = _newSize % fileSize;
					if (fileSize) resize((*this)[_newFileCount - 1], lastFileSize);
				}
			}
		}

        inline void clearInternals() {
			clear(*(Base*)this);
            fileSize = 0;
            _realign = false;
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // generic open/close interface
    template < typename TFileArray >
    inline bool _openTempFArray(TFileArray &me, int openMode) {
//IOREV
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= openTemp(me[i], openMode);
		return result;
    }

    template < typename TFileArray >
    inline bool _openTempFArray(TFileArray &me) {
//IOREV
		return _openTempFArray(me, DefaultOpenTempMode<TFileArray>::VALUE);
	}

    template < typename TFileArray >
    inline bool _reopenFArray(TFileArray &me, int openMode) {
//IOREV
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= reopen(me[i], openMode);
		return result;
    }

    template < typename TFileArray >
    inline bool _closeFArray(TFileArray &me) {
//IOREV
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			if (me[i]) result &= close(me[i]);
		return result;
    }

    template < typename TFileArray >
    inline unsigned _sectorSizeFArray(TFileArray &me, int /*openMode*/) {
//IOREV
		return sectorSize(me[0]);
    }

    template < typename TFileArray >
    inline typename Size<TFileArray>::Type
	_sizeFArray(TFileArray &me) {
//IOREV
        typename Size<TFileArray>::Type sum = 0;
		for(int i = 0; i < length(me); ++i)
			sum += size(me[i]);
		return sum;
    }

    template < typename TFileArray >
    inline bool _flushFArray(TFileArray &me) {
//IOREV
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= flush(me[i]);
		return result;
    }

    template < typename TFileArray, typename TRequest >
    inline bool _cancelFArray(TFileArray &me, TRequest &request) {
//IOREV
		bool result = true;
		for(int i = 0; i < length(me); ++i)
			result &= cancel(me[i], &request);
		return result;
    }


    //////////////////////////////////////////////////////////////////////////////
    // standard file array wrappers

    template < __int64 FILE_SIZE, typename TFile >
	inline unsigned length(File< Chained<FILE_SIZE, TFile> > const &me) {
//IOREV
		return me.fileCount();
	}

    template < unsigned FileCount_, typename TFile >
	inline unsigned length(File< Striped<FileCount_, TFile> > const &/*me*/) {
//IOREV
		return FileCount_;
	}

    template < __int64 FILE_SIZE, typename TFile >
	inline bool open(File< Chained<FILE_SIZE, TFile> > &me, const char *fileName, int openMode) {
//IOREV
		me.baseName = fileName;
		me.openMode = openMode;
		me.temporary = false;
		me.tryOpen();
		return true;
	}

    template < __int64 FILE_SIZE, typename TFile >
	inline bool openTemp(File< Chained<FILE_SIZE, TFile> > &me, int openMode) {
//IOREV
		me.openMode = openMode;
		me.temporary = true;
		return true;
	}

    template < unsigned FileCount_, typename TFile >
	inline bool openTemp(File< Striped<FileCount_, TFile> > &me, int openMode) {
//IOREV I dont think this does what intended
		return _openTempFArray(me, openMode);
	}

    template < __int64 FILE_SIZE, typename TFile >
	inline bool close(File< Chained<FILE_SIZE, TFile> > &me) {
//IOREV
        _closeFArray(me);
        me.clearInternals();
        return true;
    }

    template < unsigned FileCount_, typename TFile >
	inline bool close(File< Striped<FileCount_, TFile> > &me) {	return _closeFArray(me); }
//IOREV

    template < __int64 FILE_SIZE, typename TFile >
	__int64 size(File< Chained<FILE_SIZE, TFile> > &me) {
//IOREV
		return _sizeFArray(me);
	}

    template < unsigned FileCount_, typename TFile >
	__int64 size(File< Striped<FileCount_, TFile> > &me) {
//IOREV
		return _sizeFArray(me);
	}

    template < __int64 FILE_SIZE, typename TFile, typename TSize >
    inline void resize(File< Chained<FILE_SIZE, TFile> > &me, TSize new_length) {
//IOREV
		me.resizeArray(new_length);
    }

    template < __int64 FILE_SIZE, typename TFile, typename TValue, typename TSize >
	inline void allocate(File< Chained<FILE_SIZE, TFile> > const &me, TValue* &data, TSize count) {
//IOREV
		allocate(me[0], data, count);
	}

    template < __int64 FILE_SIZE, typename TFile, typename TValue, typename TSize >
	inline void deallocate(File< Chained<FILE_SIZE, TFile> > const &me, TValue* &data, TSize count) {
//IOREV
		deallocate(me[0], data, count);
	}

    template < unsigned FileCount_, typename TFile, typename TValue, typename TSize >
	inline void allocate(File< Striped<FileCount_, TFile> > const &me, TValue* &data, TSize count) {
//IOREV
		allocate(me[0], data, count);
	}

    template < unsigned FileCount_, typename TFile, typename TValue, typename TSize >
	inline void deallocate(File< Striped<FileCount_, TFile> > const &me, TValue* &data, TSize count) {
//IOREV
		deallocate(me[0], data, count);
	}


    //////////////////////////////////////////////////////////////////////////////
    // read/write wrappers

    template < __int64 FILE_SIZE, typename TFile, typename TValue, typename TSize, typename TOffset >
    inline bool readAt(File< Chained<FILE_SIZE, TFile> > &me, TValue *memPtr, TSize count, TOffset offset) {
//IOREV
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (!readAt(file, memPtr, xmitSize, fileOfs)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }
    
    template < __int64 FILE_SIZE, typename TFile, typename TValue, typename TSize, typename TOffset >
    inline bool writeAt(File< Chained<FILE_SIZE, TFile> > &me, TValue const *memPtr, TSize count, TOffset offset) {
//IOREV
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (!writeAt(file, memPtr, xmitSize, fileOfs)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }

    template < __int64 FILE_SIZE, typename TFile, typename TValue, typename TSize, typename TOffset, typename TRequest >
    inline bool asyncReadAt(File< Chained<FILE_SIZE, TFile> > &me, TValue *memPtr, TSize count, TOffset offset, TRequest &req) {
//IOREV
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (count != xmitSize) {
				if (!readAt(file, memPtr, xmitSize, fileOfs)) return false;
			} else
				if (!asyncReadAt(file, memPtr, xmitSize, fileOfs, req)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }
    
    template < __int64 FILE_SIZE, typename TFile, typename TValue, typename TSize, typename TOffset, typename TRequest  >
    inline bool asyncWriteAt(File< Chained<FILE_SIZE, TFile> > &me, TValue const *memPtr, TSize count, TOffset offset, TRequest &req) {
//IOREV
		TOffset fileOfs = 0;
		while (count) {
			TFile &file = me.getFileAndOffset(offset, fileOfs, memPtr);
			TSize xmitSize = _min(me.restAt(fileOfs, memPtr), (__int64)count);
			if (count != xmitSize) {
				if (!writeAt(file, memPtr, xmitSize, fileOfs)) return false;
			} else
				if (!asyncWriteAt(file, memPtr, xmitSize, fileOfs, req)) return false;
			count -= xmitSize;
			offset += xmitSize;
			memPtr += xmitSize;
		}
		return true;
    }

}

#endif
