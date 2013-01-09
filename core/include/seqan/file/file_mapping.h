
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
// Author: David Weese <david.weese@fu-berlin.de>
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

/**
.Enum.FileMappingAdviseScheme
..cat:Sequences
..cat:Input / Output
..summary:Enum with mmap advise values.
..value.MAP_NORMAL:There is no advise on the given address range.
..value.MAP_RANDOM:The address range will be accessed with a random access memory pattern.
..value.MAP_SEQUENTIAL:The address range will be accessed sequentially.
..value.MAP_WILLNEED:The address range in the advise will be needed in the future.
..value.MAP_DONTNEED:The address range in the advise will not be needed any more.
..see:Function.adviseFileSegment
..include:seqan/file.h
 */

enum FileMappingMode {
    MAP_RDONLY = 1,
    MAP_WRONLY = 2,
    MAP_RDWR = 3,
    MAP_COPYONWRITE = 4
};

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


/**
.Class.FileMapping:
..cat:File
..summary:A structure used by a @Class.FilePager@ to represent a file and its memory mapping.
..FileMapping<TFile>
..param.TFile:The file type.
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
        SEQAN_FAIL("CreateFileMapping failed in resize: \"%s\"", strerror(errno));
#endif
    return result;
}

template <typename TSpec>
inline bool
_mapFile(FileMapping<TSpec> &mapping)
{
    return _mapFile(mapping, 0);
}

template <typename TSpec>
inline bool
_unmapFile(FileMapping<TSpec> &mapping)
{
    ignoreUnusedVariableWarning(mapping);

    bool result = true;
#ifdef PLATFORM_WINDOWS
    result &= (CloseHandle(mapping.handle) != 0);
    if (!result)
        SEQAN_FAIL("CloseHandle failed in unmap: \"%s\"", strerror(errno));
#endif
    return result;
}

template <typename TSpec, typename TFilename, typename TOpenMode>
inline bool
open(FileMapping<TSpec> &mapping, TFilename const &filename, TOpenMode const &openMode)
{
    _initialize(mapping);
    bool result = open(mapping.file, filename, openMode);
    mapping.openMode = openMode;
    mapping.ownFile = true;
    mapping.temporary = false;
    mapping.fileSize = (result)? size(mapping.file) : 0ul;
    result &= _mapFile(mapping);
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
        return _mapFile(mapping);
    }
    return false;
}

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

template <typename TSpec, typename TSize>
inline bool
closeAndResize(FileMapping<TSpec> &mapping, TSize newFileSize)
{
    bool result = _unmapFile(mapping);
    resize(mapping.file, newFileSize);
    if (mapping.ownFile)
        result &= close(mapping.file);
    _initialize(mapping);
    return result;
}


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

template <typename TSpec, typename TSize>
inline bool
resize(FileMapping<TSpec> &mapping, TSize newFileSize)
{
    bool result = _unmapFile(mapping);
    resize(mapping.file, newFileSize);
    mapping.fileSize = newFileSize;
    result &= _mapFile(mapping, newFileSize);
    return result;
}

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

// cancel all transactions
template <typename TSpec, typename TPos, typename TSize>
inline bool
cancelFileSegment(FileMapping<TSpec> &, void *addr, TPos beginPos, TSize size)
{
#ifdef PLATFORM_WINDOWS
    ignoreUnusedVariableWarning(addr);
    ignoreUnusedVariableWarning(beginPos);
    ignoreUnusedVariableWarning(size);
    return true;
#else
    return (msync(addr, size + beginPos, MS_INVALIDATE) == 0);
#endif
}

template <typename TSpec, typename TPos, typename TSize>
inline bool
adviseFileSegment(FileMapping<TSpec> &, FileMappingAdvise advise, void *addr, TPos beginPos, TSize size)
{
#ifdef PLATFORM_WINDOWS
    ignoreUnusedVariableWarning(advise);
    ignoreUnusedVariableWarning(addr);
    ignoreUnusedVariableWarning(beginPos);
    ignoreUnusedVariableWarning(size);
    return true;
#else
//		posix_fadvise(mapping.file.handle, beginPos, size, advise);
    return (posix_madvise(static_cast<char*>(addr) + beginPos, size, advise) == 0);
#endif
}

template <typename TSpec, typename TPos, typename TSize>
inline void *
mapFileSegment(FileMapping<TSpec> &mapping, TPos fileOfs, TSize size, FileMappingMode mode)
{
    SEQAN_ASSERT_EQ(OPEN_RDONLY, MAP_RDONLY);
    SEQAN_ASSERT_EQ(OPEN_WRONLY, MAP_WRONLY);

    void *addr;
    if (size == 0)
        return NULL;
    mode = (FileMappingMode)(mode & (mapping.openMode & OPEN_RDWR));

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

    if (mode & OPEN_RDONLY) prot |= PROT_READ;
    if (mode & OPEN_WRONLY) prot |= PROT_WRITE;

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
        SEQAN_FAIL("mapFileSegment failed: \"%s\"", strerror(errno));
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
        SEQAN_FAIL("unmapFileSegment failed: \"%s\"", strerror(errno));
    return result;
}

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
