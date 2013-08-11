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
// Defines basic tags and flags for file access.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_FILE_INTERFACE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_FILE_INTERFACE_H_

namespace seqan {

/*!
 * @class SyncFile
 * @headerfile <seqan/file.h>
 * @brief File structure supporting synchronous input/output access.
 *
 * @signature template <[typename TSpec]>
 *            class File<Sync<TSpec> >;
 *
 * @tparam TSpec Further specializing type.  Default: <tt>void</tt>.
 *
 * @section Remarks
 *
 * This class supports pseudo-asynchronous access methods, i.e. the method to initiate an I/O request blocks until the
 * request completion.
 */

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

/*!
 * @class AsyncFile
 * @headerfile <seqan/file.h>
 * @brief File structure supporting asynchronous input/output access.
 *
 * @signature template <[typename TSpec]>
 *            class File<Async<TSpec> >;
 *
 * @tparam TSpec Further specializing type.  Default: <tt>void</tt>.
 */

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

/*!
 * @class File
 * @headerfile <seqan/file.h>
 * @brief Represents a file.
 *
 * @signature template <[typename TSpec]>
 *            class File;
 *
 * @tparam TSpec Specializing type.  Default: <tt>Async&lt;&gt;</tt>.
 */

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

/*!
 * @enum FileOpenMode
 * @headerfile <seqan/file.h>
 * @brief Flags to select th eopen mode of a @link File @endlink or external string.
 *
 * @section Remarks
 *
 * These flags can be comined via the <tt>|</tt> operator (bitwise OR).  The defualt open mode is <tt>OPEN_RDWR |
 * OPEN_CREATE | OPEN_APPEND</tt>.
 *
 * If you omit the <tt>OPEN_APPEND</tt> flag in write mode, the file will be truncated to size 0 when opened.
 *
 * @var FileOpenMode OPEN_RDONLY
 * @brief Open in read-only mode.
 *
 * @var FileOpenMode OPEN_WRONLY
 * @brief Open in write-only mode.
 *
 * @var FileOpenMode OPEN_RDWR
 * @brief Open for reading and writing.
 *
 * @var FileOpenMode OPEN_CREATE
 * @brief Create the file if it does not yet exists.
 *
 * @var FileOpenMode OPEN_APPEND
 * @brief Keep the existing data.  If this flag is not given then the file is cleared in write mode.
 *
 * @var FileOpenMode OPEN_QUIET
 * @brief Don't print any warning message if the file could not be opened.
 *
 * @var FileOpenMode OPEN_MASK
 * @brief (Internal) Bitmask to extract the read/write open mode.
 *
 * @section Examples
 *
 * Code example to test for read-only mode:
 *
 * @code{.cpp}
 * if (openMode & OPEN_MASK == OPEN_READ)
 *     // do something if opened in read-only mode
 * @endlink
 *
 * @var FileOpenMode OPEN_ASYNC
 * @brief (Internal) Open the file for asynchronous file access.  For asynchronous file access, use the @link AsyncFile
 *                   @endlink.
 *
 * @var FileOpenMode OPEN_TEMPORARY
 * @brief (Internal) Open automatically delete the file after close.  Use the <tt>openTemp</tt> methods to open
 *        temporary files.
 */

/**
.Enum.FileOpenMode
..cat:Input/Output
..summary:Flags to select the open mode of a @Class.File@ or external string.
..value.OPEN_RDONLY:Open for only reading.
..value.OPEN_WRONLY:Open for only writing.
..value.OPEN_RDWR:Open for reading and writing.
..value.OPEN_CREATE:Create a file if it not yet exists.
..value.OPEN_APPEND:Keep the existing data. If this flag is not given, the file is cleared in write mode.
..value.OPEN_QUIET:Don't print any warning message if the file could not be opened.
..value.OPEN_MASK:(Internal) Bitmask to extract the read/write open mode.
..example.text:Code example to test for read-only mode.
..example.code:
if (openMode & OPEN_MASK == OPEN_READ)
    // do something if opened in read-only mode
..value.OPEN_ASYNC:(Internal) Open the file for asynchronous file access. For asynchronous file access, use @Spec.Async@.
..value.OPEN_TEMPORARY:(Internal) Automatically delete the file after close. Use @Function.openTemp@ to open temporary files.
..remarks:These flags can be combined via the $|$ operator. The default open mode is $OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$.
..remarks:If you omit the $OPEN_APPEND$ flag in write mode, the file will be cleared when opened.
..include:seqan/seq_io.h
*/
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

/*!
 * @class AsyncRequest
 * @headerfile <seqan/file.h>
 * @brief Associated with an asynchronous I/O request.
 *
 * @signature template <typenam TFile>
 *            struct AsyncRequest;
 *
 * @tparam TFile The file type.
 *
 * @section Remarks
 *
 * This structure is used to identify asynchronous requests after their initiation.
 */

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
