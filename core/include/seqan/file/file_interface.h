// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_FILE_INTERFACE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_FILE_INTERFACE_H_

namespace seqan {

/**
.Spec.Sync:
..cat:Files
..general:Class.File
..summary:File structure supporting synchronous input/output access.
..signature:File<Sync<> >
..remarks:This class suports pseudo-asynchronous access methods, i.e. the methods to initiate a I/O request return after request completion.
..include:seqan/file.h
*/

	template <typename TSpec = void>
    struct Sync;
//IOREV

/**
.Spec.Async:
..cat:Files
..general:Class.File
..summary:File structure supporting synchronous and asynchronous input/output access.
..signature:File<Async<> >
..include:seqan/file.h
*/

	template <typename TSpec = void>
    struct Async;
//IOREV


/**
.Class.File:
..cat:Input/Output
..summary:Represents a file.
..signature:File<TSpec>
..param.TSpec:The specializing type.
...default:$Async<>$, see @Spec.Async@.
..include:seqan/file.h
*/

	template <typename TSpec = Async<> >
    class File;
//IOREV

/**
.Spec.Chained:
..cat:Files
..general:Class.File
..summary:Splits a large file into a chain of smaller files.
..signature:File<Chained<FileSize, TFile> >
..param.FileSize:The maximal split file size in byte.
...default:2^31-1 (~2GB)
..param.TFile:Underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:This file type uses a chain of $TFile$ files, whose file sizes are at most $FileSize$ bytes.
Chained Files should be used for file systems or $TFile$ types that don't support large files (e.g. FAT32, C-style FILE*).
..remarks:The chain can be used as if it were one contiguous file.
..include:seqan/file.h
*/

	// chained file's default filesize is 2gb-1byte (fat16 filesize limitation)
	template < __int64 FileSize_ = ~(((__int64)1) << 63), typename TFile = File<> >
	struct Chained;
//IOREV

/**
.Spec.Striped:
..cat:Files
..general:Class.File
..summary:Stripes a file across multiple files.
..signature:File<Chained<FileCount, TFile> >
..param.FileCount:The number of files used for striping.
...default:2
..param.TFile:Underlying @Class.File@ type.
...default:$File<>$, see @Class.File@.
..remarks:This file type uses a software striping without redundance (see RAID0) to accelerate I/O access when using more than one disks.
..remarks:Striped files should only be used in @Class.Pool@s or external Strings as they only support block operations and no random accesses.
..include:seqan/file.h
*/

	template < unsigned FileCount_ = 2, typename TFile = File<> >
	struct Striped;
//IOREV not known if working, see file_array.h

    enum FileOpenMode {
        OPEN_RDONLY     = 1,
        OPEN_WRONLY     = 2,
        OPEN_RDWR       = 3,
        OPEN_MASK       = 3,
        OPEN_CREATE     = 4,
        OPEN_APPEND     = 8,
        OPEN_ASYNC      = 16,
		OPEN_TEMPORARY	= 32,
		OPEN_QUIET		= 128
    }; //IOREV is it intended that two labels share the same value? What is OPEN_MASK anyway?

	template <typename T>
	struct DefaultOpenMode {
//IOREV
		enum { VALUE = OPEN_RDWR | OPEN_CREATE | OPEN_APPEND };
	};

	template <typename T>
	struct DefaultOpenTempMode {
//IOREV
		enum { VALUE = OPEN_RDWR | OPEN_CREATE };
	};

    enum FileSeekMode {
        SEEK_BEGIN   = 0,
        SEEK_CURRENT = 1
#ifndef SEEK_END
      , SEEK_END     = 2
#endif
    }; //IOREV why not use constants SEEK_SET, SEEK_CUR, SEEK_END from cstdio?

    //////////////////////////////////////////////////////////////////////////////
    // result type of asynch. functions
    // you have to call release(AsyncRequest<T>) after a finished *event based* transfer
	struct AsyncDummyRequest {};
//IOREV

/**
.Class.AsyncRequest:
..cat:Input/Output
..summary:Associated with an asynchronous I/O request.
..signature:AsyncRequest<TFile>
..param.TFile:A File type.
..remarks:This structure is used to identify asynchronous requests after their initiation.
..include:seqan/file.h
*/

    template < typename T >
    struct AsyncRequest
    {
//IOREV _stub_ this seems not to be implemented at all, most functions are commented
        typedef AsyncDummyRequest Type;
    };
}  // namespace seqan;

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_FILE_INTERFACE_H_
