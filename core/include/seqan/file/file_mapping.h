// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Provides a platform independent access to memory mapping of files.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_FILE_MAPPING_H_
#define SEQAN_CORE_INCLUDE_SEQAN_FILE_MAPPING_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @enum FileMappingMode
 * @headerfile <seqan/file.h>
 * @brief Flags to determine the mapping mode of mapFileSegment.
 *
 * @var FileMappingMode MAP_RDONLY = 1;
 * @brief Map the segment in read-only mode.
 *
 * @var FileMappingMode MAP_WRONLY = 2;
 * @brief Map the segment in write-only mode.
 *
 * @var FileMappingMode MAP_RDWR = 3;
 * @brief Map the segment for reading and writing.
 *
 * @var FileMappingMode MAP_COPYONWRITE = 4;
 * @brief Write accesses are not written back to file and not shared among different mappings.
 *
 * @section Remarks
 *
 * The mapping mode must be compatible to the open mode of a @link FileMapping @endlink, e.g. <tt>MAP_RDWR</tt> is not
 * allowed if the file mapping was opened with <tt>OPEN_RDONLY</tt>.
 */

/**
.Enum.FileMappingMode
..cat:Sequences
..cat:Input/Output
..summary:Flags to define the mapping mode of @Function.mapFileSegment@.
..value.MAP_COPYONWRITE:Write accesses are not written back to file and not shared among different mappings.
..value.MAP_RDONLY:Map the segment in read-only mode.
..value.MAP_RDWR:Map the segment for reading and writing.
..value.MAP_WRONLY:Map the segment in write-only mode.
..remarks:The mapping mode must be compatible to the open mode of a @Class.FileMapping@,
e.g. $MAP_RDWR$ is not allowed if the file mapping was opened with $OPEN_RDONLY$.
..include:seqan/file.h
*/

enum FileMappingMode {
    MAP_RDONLY = 1,
    MAP_WRONLY = 2,
    MAP_RDWR = 3,
    MAP_COPYONWRITE = 4
};

/*!
 * @enum FileMappingAdvise
 * @headerfile <seqan/file.h>
 * @brief Enum with MMAP advise values.
 *
 * @var FileMappingAdvise MAP_NORMAL;
 * @brief There is no advise on the given address range.
 *
 * @var FileMappingAdvise MAP_RANDOM;
 * @brief The address range will be accessed with random access memory pattern.
 *
 * @var FileMappingAdvise MAP_SEQUENTIAL;
 * @brief The address range will be accessed sequentially.
 *
 * @var FileMappingAdvise MAP_WILLNEED;
 * @brief The address range in the advise will be needed in the future.
 *
 * @var FileMappingAdvise MAP_DONTNEED;
 * @brief The address range in the advise will not be needed any more.
 */

/**
.Enum.FileMappingAdvise
..cat:Sequences
..cat:Input/Output
..summary:Enum with mmap advise values.
..value.MAP_NORMAL:There is no advise on the given address range.
..value.MAP_RANDOM:The address range will be accessed with a random access memory pattern.
..value.MAP_SEQUENTIAL:The address range will be accessed sequentially.
..value.MAP_WILLNEED:The address range in the advise will be needed in the future.
..value.MAP_DONTNEED:The address range in the advise will not be needed any more.
..include:seqan/file.h
*/

#ifdef PLATFORM_WINDOWS

enum FileMappingAdvise {
    MAP_NORMAL = 0,
    MAP_RANDOM = 0,
    MAP_SEQUENTIAL = 0,
    MAP_WILLNEED = 0,
    MAP_DONTNEED = 0
};

#else

enum FileMappingAdvise {
    MAP_NORMAL = POSIX_MADV_NORMAL,
    MAP_RANDOM = POSIX_MADV_RANDOM,
    MAP_SEQUENTIAL = POSIX_MADV_SEQUENTIAL,
    MAP_WILLNEED = POSIX_MADV_WILLNEED,
    MAP_DONTNEED = POSIX_MADV_DONTNEED
};

#endif

/*!
 * @class FileMapping
 * @headerfile <seqan/file.h>
 * @brief A structure to memory-map a file.
 *
 * @signature template <typename TSpec>
 *            stuct FileMapping;
 *
 * @tparam TSpec The specializing type.  Default: <tt>void</tt>.
 *
 * @section Remarks
 *
 * This structure represents both a file and its memory mapping.
 */

/**
.Class.FileMapping:
..cat:File
..summary:A structure to memory-map a file.
..signature:FileMapping<TSpec>
..param.TSpec:The specializing type.
...default:$void$
..remarks:This structure represents both a file and its memory mapping.
..include:seqan/file.h
*/

#ifdef PLATFORM_WINDOWS
static SECURITY_ATTRIBUTES FileMappingDefaultAttributes =
{
    sizeof(SECURITY_ATTRIBUTES),
    NULL,
    true
};
#endif

template <typename TSpec = void>
struct FileMapping
{
    // -----------------------------------------------------------------------
    // Typedefs
    // -----------------------------------------------------------------------

    typedef File<Async<> >              TFile;
    typedef typename Size<TFile>::Type  TFileSize;

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------

#ifdef PLATFORM_WINDOWS
    HANDLE      handle;
#endif

    TFileSize   fileSize;
    TFile       file;
    int         openMode;
    bool        ownFile;
    bool        temporary;


    FileMapping()
    {
        _initialize(*this);
    }

//____________________________________________________________________________

    inline operator bool()
    {
        return file;
    }
//____________________________________________________________________________
};

template <typename TSpec>
struct Size<FileMapping<TSpec> >:
    public Size<typename FileMapping<TSpec>::TFile> {};


template <typename TSpec>
inline void
_initialize(FileMapping<TSpec> &mapping)
{
#ifdef PLATFORM_WINDOWS
    mapping.handle = NULL;
#endif
    mapping.fileSize = 0;
    mapping.openMode = OPEN_RDWR;
    mapping.ownFile = false;
    mapping.temporary = true;
}


template <typename TSpec, typename TSize>
inline bool
_mapFile(FileMapping<TSpec> &mapping, TSize mappingSize)
{
    ignoreUnusedVariableWarning(mapping);
    ignoreUnusedVariableWarning(mappingSize);

    bool result = true;
#ifdef PLATFORM_WINDOWS
    if (mappingSize == 0)
    {
        mapping.handle = NULL;
        return true;
    }

    DWORD prot = ((mapping.openMode & OPEN_MASK) == OPEN_RDONLY) ? PAGE_READONLY : PAGE_READWRITE;
    LARGE_INTEGER largeSize;
    largeSize.QuadPart = mappingSize;   // 0 = map the whole file

    mapping.handle = CreateFileMapping(
        mapping.file.handle,            // _In_     HANDLE hFile,
        &FileMappingDefaultAttributes,  // _In_opt_ LPSECURITY_ATTRIBUTES lpAttributes,
        prot,                           // _In_     DWORD flProtect,
        largeSize.HighPart,             // _In_     DWORD dwMaximumSizeHigh,
        largeSize.LowPart,              // _In_     DWORD dwMaximumSizeLow,
        NULL                            // _In_opt_ LPCTSTR lpName
    );
    result &= (mapping.handle != NULL);

    if (mapping.handle == NULL)
    {
        LPVOID lpMsgBuf;
        FormatMessage(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
            NULL,
            GetLastError(),
            MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
            (LPTSTR) &lpMsgBuf,
            0,
            NULL);
    
        SEQAN_FAIL("CreateFileMapping failed in resize: \"%s\"", lpMsgBuf /*strerror(GetLastError())*/);
        LocalFree(lpMsgBuf);
    }
#endif
    return result;
}

template <typename TSpec>
inline bool
_unmapFile(FileMapping<TSpec> &mapping)
{
    ignoreUnusedVariableWarning(mapping);

    bool result = true;
#ifdef PLATFORM_WINDOWS
    if (mapping.handle != NULL)
    {
        result &= (CloseHandle(mapping.handle) != 0);
        if (!result)
            SEQAN_FAIL("CloseHandle failed in unmap: \"%s\"", strerror(errno));
        mapping.handle = NULL;
    }
#endif
    return result;
}

/*!
 * @fn FileMapping#open
 * @brief Open a file to be mapped into memory.
 *
 * @signature bool open(fileMapping, fileName[, openMode]);
 *
 * @param[in,out] fileMapping The FileMapping to open.
 * @param[in]     fileName    The path to the fie.
 * @param[in]     openMode    The mode to open the file in, flags from @link FileOpenMode @endlink to combine
 *                            using OR.  Write-only mode is not supported, use <tt>OPEN_RDWR</tt> if you need
 *                            write access.  If you omit the <tt>OPEN_APPEND</tt> flag in write mode, the file
 *                            will be cleared when opened.  Default: <tt>OPEN_RDWR | OPEN_CREATe | OPEN_APPEND</tt>.
 *
 * @return bool <tt>true</tt> if the opening was successful, <tt>false</tt> otherwise.
 */

/**
.Function.FileMapping#open:
..class:Class.FileMapping
..summary:Open a file to be mapped into memory.
..cat:Input/Output
..signature:open(mapping, fileName[, openMode])
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened. See @Enum.FileOpenMode@ for more details.
...remarks:Write-only mode is not supported, use OPEN_RDWR if you need write access.
If you omit the $OPEN_APPEND$ flag in write mode, the file will be cleared when opened.
...type:Enum.FileOpenMode
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

template <typename TSpec>
inline bool
open(FileMapping<TSpec> &mapping, const char *filename, int openMode)
{
    _initialize(mapping);
    bool result = open(mapping.file, filename, openMode);
    mapping.openMode = openMode;
    mapping.ownFile = true;
    mapping.temporary = false;
    mapping.fileSize = (result)? size(mapping.file) : 0ul;
    result &= _mapFile(mapping, mapping.fileSize);
    return result;
}

template <typename TSpec, typename TFile>
inline bool
open(FileMapping<TSpec> &mapping, TFile const &file)
{
    _initialize(mapping);
    mapping.file = file;
    mapping.openMode = OPEN_RDWR;
    mapping.ownFile = false;
    mapping.temporary = false;
    if (mapping.file)
    {
        mapping.fileSize = size(mapping.file);
        return _mapFile(mapping, mapping.fileSize);
    }
    return false;
}

/*!
 * @fn FileMapping#openTemp
 * @brief Open a temporary file to be mapped into memory.
 *
 * @signature bool openTemp(fileMapping);
 *
 * @param[in,out] fileMapping The FileMapping to open using a temporary file.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> otherwise.
 */

/**
.Function.FileMapping#openTemp:
..class:Class.FileMapping
..summary:Open a temporary file to be mapped into memory.
..cat:Input/Output
..signature:openTemp(mapping)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

template <typename TSpec>
inline bool
openTemp(FileMapping<TSpec> &mapping)
{
    _initialize(mapping);
    bool result = openTemp(mapping.file);
    mapping.openMode = OPEN_RDWR;
    mapping.ownFile = true;
    mapping.temporary = true;
    return result;
}

/*!
 * @fn FileMapping#close
 * @brief Close a file and its memory mapping.
 * 
 * @signature bool close(mapping);
 *
 * @param[in,out] mapping The FileMapping to close
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> otherwise.
 */

/**
.Function.FileMapping#close:
..class:Class.FileMapping
..summary:Close a file and its memory mapping.
..cat:Input/Output
..signature:close(mapping)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

template <typename TSpec>
inline bool
close(FileMapping<TSpec> &mapping)
{
    bool result = _unmapFile(mapping);
    if (mapping.ownFile)
        result &= close(mapping.file);
    _initialize(mapping);
    return result;
}

/*!
 * @fn FileMapping#closeAndResize
 * @brief Close a memory mapping and resize and close the underlying file.
 *
 * @signature bool closeAndResize(mapping, newFileSize);
 *
 * @param[in,out] mapping     The FileMapping to close.
 * @param[in]     newFileSize The size the file should have after closing.
 *
 * @return bool <tt>true</tt> indicating success, <tt>false</tt> failure.
 */

/**
.Function.FileMapping#closeAndResize:
..class:Class.FileMapping
..summary:Close a memory mapping and resize and close the underlying file.
..cat:Input/Output
..signature:closeAndResize(mapping, newFileSize)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.newFileSize:The new file size.
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

template <typename TSpec, typename TSize>
inline bool
closeAndResize(FileMapping<TSpec> &mapping, TSize newFileSize)
{
    bool result = _unmapFile(mapping);
    if (newFileSize != mapping.fileSize)
        resize(mapping.file, newFileSize);
    if (mapping.ownFile)
        result &= close(mapping.file);
    _initialize(mapping);
    return result;
}

/*!
 * @fn FileMapping#length
 * @brief Return the file size of a memory mapping.
 * 
 * @signature TSize length(mapping);
 *
 * @param[in] mapping The FileMapping to return the length for.
 *
 * @return TSize The file size.  TSize is the size type of the file, as returned by the Size metafunction.
 */

/**
.Function.FileMapping#length:
..class:Class.FileMapping
..summary:Return the file size of a memory mapping.
..cat:Input/Output
..signature:length(mapping)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..returns:The size of the underlying file.
..include:seqan/file.h
*/

template <typename TSpec>
inline typename Size<FileMapping<TSpec> >::Type
length(FileMapping<TSpec> &mapping)
{
    return mapping.fileSize;
}

template <typename TSpec>
inline typename Size<FileMapping<TSpec> >::Type
length(FileMapping<TSpec> const &mapping)
{
    return mapping.fileSize;
}

/*!
 * @fn FileMapping#resize
 * @brief Resize the underlying file.
 *
 * @signature bool resize(mapping, newFileSize);
 *
 * @param[in,out] mapping     The FileMapping to resize.
 * @param[in]     newFileSize The new file size to set.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> otherwise.
 *
 * @section Remarks
 *
 * On Windows, all existing file mappings must be unmapped via @link FileMapping#unmapFileSegment @endlink before
 * claling this function.
 */

/**
.Function.FileMapping#resize:
..class:Class.FileMapping
..summary:Resize the underlying file.
..cat:Input/Output
..signature:resize(mapping, newFileSize)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.newFileSize:The new file size.
..returns:A $bool$ which is $true$ on success.
..remarks:Under Windows, all existing file mappings must be unmapped via @Function.unmapFileSegment@ prior calling this function.
..include:seqan/file.h
*/

template <typename TSpec, typename TSize>
inline bool
resize(FileMapping<TSpec> &mapping, TSize newFileSize)
{
    typedef typename Size<FileMapping<TSpec> >::Type TFileSize;
    if ((TFileSize)newFileSize == mapping.fileSize)
        return true;

    bool result = _unmapFile(mapping);
    resize(mapping.file, newFileSize);
    mapping.fileSize = newFileSize;
    result &= _mapFile(mapping, newFileSize);
    return result;
}

/*!
 * @fn FileMapping#flushFileSegment
 * @brief Wait for all outstanding transactions of a memory-mapped file segment.
 *
 * @signature bool flushFileSegment(fileMapping, addr, beginPos, size);
 *
 * @param[in,out] fileMapping A FileMapping object.
 * @param[in]     addr        A pointer to the beginning of memory-mapped segment in memory (returned by a prior
 *                            call of @link FileMapping#mapFileSegment @endlink.
 * @param[in]     beginPos    The absolute start address of the segment in bytes.
 * @param[in]     size        The segment length in bytes.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 *
 * @section Remarks
 *
 * This function has no effect under Windows.  On all other platforms it calls <tt>msync</tt>.  This function is only
 * needed to synchronize file accesses in non-shared-memory environments.
 */

/**
.Function.flushFileSegment:
..class:Class.FileMapping
..summary:Wait for all outstanding transactions of a memory-mapped file segment.
..cat:Input/Output
..signature:flushFileSegment(mapping, addr, fileOfs, size)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.addr:A pointer to the beginning of the memory-mapped segment in memory (returned by a prior call of @Function.mapFileSegment@).
..param.fileOfs:The absolute start address of the segment in bytes.
..param.size:The segment length in bytes.
..returns:A $bool$ which is $true$ on success.
..remarks:This function has no effect under Windows. On all other platforms it calls $msync$.
This function is only needed to synchronize file accesses in non-shared-memory environments.
..include:seqan/file.h
*/

template <typename TSpec, typename TPos, typename TSize>
inline bool
flushFileSegment(FileMapping<TSpec> &, void *addr, TPos beginPos, TSize size)
{
#ifdef PLATFORM_WINDOWS
    ignoreUnusedVariableWarning(addr);
    ignoreUnusedVariableWarning(beginPos);
    ignoreUnusedVariableWarning(size);
    return true;
#else
    return (msync(static_cast<char*>(addr) + beginPos, size, MS_SYNC) == 0);
#endif
}

/*!
 * @fn FileMapping#cancelFileSegment
 * @brief Cancel all outstanding transactions of a memory-mapped file segment.
 *
 * @signature bool cancelFileSegment(fileMapping, addr, beginPos, size);
 *
 * @param[in,out] fileMapping The FileMapping object.
 * @param[in]     addr        A pointer to the beginning of the memory-mapped segment in memory (returned by a prior
 *                            call of FileMapping#mapFileSegment).
 * @param[in]     beginPos    the absolute start address of the segment in bytes.
 * @param[in]     size        The segment length in bytes.
 *
 * @returns bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

/**
.Function.cancelFileSegment:
..class:Class.FileMapping
..summary:Cancel all outstanding transactions of a memory-mapped file segment.
..cat:Input/Output
..signature:cancelFileSegment(mapping, addr, fileOfs, size)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.addr:A pointer to the beginning of the memory-mapped segment in memory (returned by a prior call of @Function.mapFileSegment@).
..param.fileOfs:The absolute start address of the segment in bytes.
..param.size:The segment length in bytes.
..returns:A $bool$ which is $true$ on success.
..remarks:This function has no effect under Windows. On all other platforms it calls $msync$.
..include:seqan/file.h
*/

template <typename TSpec, typename TPos, typename TSize>
inline bool
cancelFileSegment(FileMapping<TSpec> &, void *addr, TPos fileOfs, TSize size)
{
#ifdef PLATFORM_WINDOWS
    ignoreUnusedVariableWarning(addr);
    ignoreUnusedVariableWarning(fileOfs);
    ignoreUnusedVariableWarning(size);
    return true;
#else
    return (msync(addr, size + fileOfs, MS_INVALIDATE) == 0);
#endif
}

/*!
 * @fn FileMapping#adviseFileSegment
 * @brief Give advise about use of a memory-mapped file segment.
 *
 * @signature bool adviseFileSegment(fileMapping, advise, addr, fileOfs, size);
 *
 * @param[in,out] fileMapping The FileMapping object to adsive segment in.
 * @param[in]     advise      The advise type.  Type: @link FileMappingAdvise @endlink.
 * @param[in]     addr        A pointer t othe beginning of the memory-mapped segment in memory (returned by a
 *                            prior call of @link FileMapping#mapFileSegment @endlink).
 * @param[in]     fileOfs     The absolute start address of the segment in bytes.
 * @param[in]     size        The segment length in bytes.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 *
 * @section Remarks
 *
 * This function has no effect on Windows.  On all other platforms it calls <tt>posix_madvise</tt>.
 */

/**
.Function.adviseFileSegment:
..class:Class.FileMapping
..summary:Give advice about use of a memory-mapped file segment.
..cat:Input/Output
..signature:adviseFileSegment(mapping, advise, addr, fileOfs, size)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.advise:Advise flags
...type:Enum.FileMappingAdvise
..param.addr:A pointer to the beginning of the memory-mapped segment in memory (returned by a prior call of @Function.mapFileSegment@).
..param.fileOfs:The absolute start address of the segment in bytes.
..param.size:The segment length in bytes.
..returns:A $bool$ which is $true$ on success.
..remarks:This function has no effect under Windows. On all other platforms it calls $posix_madvise$.
..include:seqan/file.h
*/
template <typename TSpec, typename TPos, typename TSize>
inline bool
adviseFileSegment(FileMapping<TSpec> &, FileMappingAdvise advise, void *addr, TPos fileOfs, TSize size)
{
#ifdef PLATFORM_WINDOWS
    ignoreUnusedVariableWarning(advise);
    ignoreUnusedVariableWarning(addr);
    ignoreUnusedVariableWarning(fileOfs);
    ignoreUnusedVariableWarning(size);
    return true;
#else
//		posix_fadvise(mapping.file.handle, beginPos, size, advise);
    return (posix_madvise(static_cast<char*>(addr) + fileOfs, size, advise) == 0);
#endif
}

/*!
 * @fn FileMapping#mapFileSegment
 * @brief Map a segment of a file into memory.
 *
 * @signature TPtr mapFileSegment(fileMapping, fileOfs[, size[, mode]]);
 *
 * @param[in,out] fileMapping A FileMapping object.
 * @param[in]     fileOfs     The absolute start address of the segment in bytes.
 * @param[in]     size        The segment length in bytes.
 * @param[in]     mode        The mapping access mode.  Default rread/write open mode of the underlying file.
 *                            Type: @link FileMappingMode @endlink.
 *
 * @return TPtr A pointer to the beginning of the memory-mapped segment in memory or <tt>NULL</tt> on error.  TPtr is
 *              <tt>void *</tt>.
 */

/**
.Function.mapFileSegment:
..class:Class.FileMapping
..summary:Map a segment of a file into memory.
..cat:Input/Output
..signature:mapFileSegment(mapping, fileOfs[, size[, mode]])
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.fileOfs:The absolute start address of the segment in bytes.
..param.size:The segment length in bytes.
...default:The rest of the file.
..param.mode:The mapping access mode.
...type:Enum.FileMappingMode
...default:The read/write open mode of the underlying file.
..returns:A pointer to the beginning of the memory-mapped segment in memory or NULL on error. 
...type:nolink:$void *$
..include:seqan/file.h
*/

template <typename TSpec, typename TPos, typename TSize, typename TFileMappingMode>
inline void *
mapFileSegment(FileMapping<TSpec> &mapping, TPos fileOfs, TSize size, TFileMappingMode mode)
{
    SEQAN_ASSERT_EQ((int)OPEN_RDONLY, (int)MAP_RDONLY);
    SEQAN_ASSERT_EQ((int)OPEN_WRONLY, (int)MAP_WRONLY);

    void *addr;
    if (size == 0)
        return NULL;
    mode = (FileMappingMode)(mode & (mapping.openMode & OPEN_MASK));

#ifdef PLATFORM_WINDOWS

    DWORD access = ((mode & OPEN_MASK) == OPEN_RDONLY) ? FILE_MAP_READ : FILE_MAP_ALL_ACCESS;
    LARGE_INTEGER largeOfs;
    largeOfs.QuadPart = fileOfs;

    addr = MapViewOfFile(
        mapping.handle,                 // _In_     HANDLE hFileMappingObject,
        access,                         // _In_     DWORD dwDesiredAccess,
        largeOfs.HighPart,              // _In_     DWORD dwFileOffsetHigh,
        largeOfs.LowPart,               // _In_     DWORD dwFileOffsetLow,
        size);                          // _In_     SIZE_T dwNumberOfBytesToMap (0 = map the whole file)
#else
    int prot = 0;

    if (mode & MAP_RDONLY) prot |= PROT_READ;
    if (mode & MAP_WRONLY) prot |= PROT_WRITE;

    int flags = 0;
    if ((mode & MAP_COPYONWRITE) != 0)
        flags |= MAP_PRIVATE;
    else
        flags |= MAP_SHARED;

    addr = mmap(
        NULL,
        size,
        prot,
        flags,
        mapping.file.handle,
        fileOfs);
//    std::cerr << "mmap(0,"<<size<<','<<prot<<','<<flags<<','<<mapping.file.handle<<','<<fileOfs<<")="<<addr<<std::endl;

    if (addr == MAP_FAILED)
        addr = NULL;
#endif
    if (addr == NULL)
    {
        SEQAN_FAIL("mapFileSegment(%i,%i,%i) failed (filesize=%i): \"%s\"", fileOfs, size, mode, seqan::size(mapping.file), strerror(errno));
    }
    return addr;
}

template <typename TSpec, typename TPos, typename TSize>
inline void *
mapFileSegment(FileMapping<TSpec> &mapping, TPos fileOfs, TSize size)
{
    return mapFileSegment(mapping, fileOfs, size, MAP_RDWR);
}

template <typename TSpec, typename TPos, typename TSize>
inline void *
mapFileSegment(FileMapping<TSpec> &mapping, TPos fileOfs = 0)
{
    return mapFileSegment(mapping, fileOfs, length(mapping) - fileOfs);
}

/*!
 * @fn FileMapping#unmapFileSegment
 * @brief Unmap a memory-mapped file segment.
 *
 * @signature bool unmappedFileSegment(fileMapping, addr, size);
 *
 * @param[in,out] fileMapping A FileMapping object.
 * @param[in]     addr        The pointer to the beginning of the memory-mapped segment in memory.  Type: <tt>void*</tt>.
 * @param[in]     size        The segment length in bytes.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

/**
.Function.unmapFileSegment:
..class:Class.FileMapping
..summary:Unmap a memory-mapped file segment.
..cat:Input/Output
..signature:unmapFileSegment(mapping, addr, size)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.addr:A pointer to the beginning of the memory-mapped segment in memory.
..param.size:The segment length in bytes.
..returns:A $bool$ which is $true$ on success.
..include:seqan/file.h
*/

template <typename TSpec, typename TSize>
inline bool
unmapFileSegment(FileMapping<TSpec> &, void *addr, TSize size)
{
    bool result;
#ifdef PLATFORM_WINDOWS
    ignoreUnusedVariableWarning(size);
    result = (UnmapViewOfFile(addr) != 0);
#else
    result = (munmap(addr, size) == 0);
//    std::cerr << "munmap("<<addr<<','<<size<<")="<<result<<std::endl;

#endif
    if (!result)
        SEQAN_FAIL("unmapFileSegment(%x,%i) failed: \"%s\"", (unsigned long)addr, size, strerror(errno));
    return result;
}

/*!
 * @fn FileMapping#remapFileSegment
 * @brief Change the size of a memory-mapped file segment.
 *
 * @signature TPtr remapFileSegment(fileMapping, oldAddr, oldFileOfs, oldSize, newSize);
 *
 * @param[in,out] fileMapping The FileMapping object.
 * @param[in]     oldAddr     The address returned by @link FileMapping#mapFileSegment @endlink.
 * @param[in]     oldFileOfs  The <tt>fileOfs</tt> parameter used in @link FileMapping#mapFileSegment @endlink.
 * @param[in]     oldSize     The <tt>size</tt> parameter used in @link FileMapping#mapFileSegment @endlink.
 * @param[in]     newSize     the new segment length in bytes.
 *
 * @return TPtr A pointer to the beginning of the memory-mapped segment in memory or <tt>NULL</tt> on error.  Type:
 *              <tt>void*</tt>.
 */

/**
.Function.remapFileSegment:
..class:Class.FileMapping
..summary:Change the size of a memory-mapped file segment.
..cat:Input/Output
..signature:remapFileSegment(mapping, oldAddr, oldFileOfs, oldSize, newSize)
..param.mapping:A file mapping object.
...type:Class.FileMapping
..param.oldAddr:The address returned by @Function.mapFileSegment@.
..param.oldFileOfs:The fileOfs parameter used in @Function.mapFileSegment@.
..param.oldSize:The size parameter used in @Function.mapFileSegment@.
..param.newSize:The new segment length in bytes.
..returns:A pointer to the beginning of the memory-mapped segment in memory or NULL on error. 
...type:nolink:$void *$
..include:seqan/file.h
*/

template <typename TSpec, typename TPos, typename TSize>
inline void *
remapFileSegment(FileMapping<TSpec> &mapping, void *oldAddr, TPos oldFileOfs, TSize oldSize, TSize newSize)
{
    void *addr;
#if !defined(PLATFORM_WINDOWS) && defined(MREMAP_MAYMOVE)
    ignoreUnusedVariableWarning(mapping);
    ignoreUnusedVariableWarning(oldFileOfs);
    addr = mremap(oldAddr, oldSize, newSize, MREMAP_MAYMOVE);
//    std::cerr << "mremap("<<oldAddr<<','<<oldSize<<','<<newSize<<','<<MREMAP_MAYMOVE<<")="<<addr<<std::endl;
    if (addr == MAP_FAILED)
        addr = NULL;
#else
    // for BSD systems without mremap(..) like Mac OS X ...
    unmapFileSegment(mapping, oldAddr, oldSize);
    addr = mapFileSegment(mapping, oldFileOfs, newSize);
#endif
    if (addr == NULL)
        SEQAN_FAIL("remapFileSegment failed: \"%s\"", strerror(errno));
    return addr;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_FILE_MAPPING_H_
