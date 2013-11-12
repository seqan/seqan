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
// BGZF File.
// ==========================================================================

#ifndef SEQAN_FILE_FILE_BGZF_H_
#define SEQAN_FILE_FILE_BGZF_H_

#include <zlib.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BgzfFile_;
typedef Tag<BgzfFile_> BgzfFile;

template <typename T>
struct MagicHeader<BgzfFile, T>
{
    static unsigned char const VALUE[3];
};

template <typename T>
unsigned char const MagicHeader<BgzfFile, T>::VALUE[3] = { 0x1f, 0x8b, 0x08 };  // gzip's magic number


template <typename TSpec = MMap<> >
struct Bgzf {};


template <typename TValue, typename TDirection, typename TSpec>
struct Host<FilePageTable<TValue, TDirection, Bgzf<TSpec> > >:
    public Host<FilePageTable<TValue, TDirection, TSpec> > {};


struct RandomPagingScheme
{
    typedef std::map<__uint64, void *> TMap;

    enum { pageSize = 64 * 1024 };

    static void * EMPTY;
    static void * ON_DISK;

    TMap frameStart;
};

void * RandomPagingScheme::EMPTY = NULL;
void * RandomPagingScheme::ON_DISK = (void *)-1;

template <typename TValue, typename TDirection, typename TSpec>
struct PagingScheme<FilePageTable<TValue, TDirection, Bgzf<TSpec> > >
{
    typedef RandomPagingScheme Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _getPageOffsetAndLength()
// ----------------------------------------------------------------------------

template <typename TPos>
inline Pair<__int64, unsigned>
_getPageOffsetAndLength(RandomPagingScheme & table, TPos pos)
{
    // as we don't know the actual compressed size of a BGZF page here,
    // we have to read the whole 64K page
    return Pair<__int64, unsigned>(pos, table.pageSize);
}

// ----------------------------------------------------------------------------
// Function _getFrameStart()
// ----------------------------------------------------------------------------

template <typename TFilePos, typename TSize>
inline void *
_getFrameStart(RandomPagingScheme const &table, TFilePos filePos, TSize)
{
    typedef RandomPagingScheme::TMap TMap;
    TMap::const_iterator iter = table.frameStart.find(filePos);

    if (SEQAN_LIKELY(iter != table.frameStart.cend()))
        return iter->second;
    else
        return table.EMPTY;
}

// ----------------------------------------------------------------------------
// Function _setFrameStart()
// ----------------------------------------------------------------------------

template <typename TFilePos, typename TSize>
inline void
_setFrameStart(RandomPagingScheme &table, TFilePos filePos, TSize, void * frameStart)
{
    table.frameStart[filePos] = frameStart;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfUnpackInt16()
// ----------------------------------------------------------------------------

inline unsigned short
_bgzfUnpackInt16(unsigned char const * buffer)
{
    return (unsigned short)buffer[0] | ((unsigned short)buffer[1] << 8);
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfPackInt16()
// ----------------------------------------------------------------------------

inline
void
_bgzfPackInt16(unsigned char * buffer, unsigned short value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfPackInt32()
// ----------------------------------------------------------------------------

inline void
_bgzfPackInt32(unsigned char * buffer, unsigned value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfCheckHeader()
// ----------------------------------------------------------------------------

inline bool
_bgzfCheckHeader(unsigned char const * header)
{
    const char FLG_FEXTRA = 4;
    const char BGZF_ID1 = 'B';
    const char BGZF_ID2 = 'C';
    const char BGZF_LEN = 2;
    const char BGZF_XLEN = 6;  // BGZF_LEN+4

    return (header[0] == MagicHeader<BgzfFile>::VALUE[0] &&
            header[1] == MagicHeader<BgzfFile>::VALUE[1] &&
            header[2] == MagicHeader<BgzfFile>::VALUE[2] &&
            (header[3] & FLG_FEXTRA) != 0 &&
            _bgzfUnpackInt16(header + 10) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            _bgzfUnpackInt16(header + 14) == BGZF_LEN);
}

// ----------------------------------------------------------------------------
// Function _preprocessFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_preprocessFilePage(FilePageTable<TValue, TDirection, Bgzf<TSpec> > &, TPageFrame &page, TBool const &)
{
    const unsigned MAX_BLOCK_SIZE = 64 * 1024;
    const unsigned BLOCK_HEADER_LENGTH = 18;
    const int GZIP_WINDOW_BITS = -15;  // no zlib header

    Bytef *header = static_cast<Bytef *>(static_cast<void *>(begin(page.raw, Standard())));

    if (length(page.raw) < BLOCK_HEADER_LENGTH)
        throw IOException("BGZF block header too short.");

    if (!_bgzfCheckHeader(header))
        throw IOException("Invalid BGZF block header.");

    // Make sure there is enough space in the buffer for decompressed data.
    resize(page.raw, _bgzfUnpackInt16(header + 16) + 1);
    reserve(page.data, MAX_BLOCK_SIZE);

    z_stream zs;
	int status;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = header + BLOCK_HEADER_LENGTH;
    zs.avail_in = length(page.raw) - (BLOCK_HEADER_LENGTH - 2);
    zs.next_out = static_cast<Bytef *>(static_cast<void *>(begin(page.data, Standard())));
    zs.avail_out = capacity(page.data);

    status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if (status != Z_OK)
        throw IOException("BGZF inflateInit2() failed.");

    status = inflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END)
    {
        inflateEnd(&zs);
        throw IOException("BGZF inflate() failed.");
    }

    status = inflateEnd(&zs);
    if (status != Z_OK)
        throw IOException("BGZF inflateEnd() failed.");

    resize(page.data, zs.total_out);
    return true;
}

template <typename TValue, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_preprocessFilePage(FilePageTable<TValue, Output, Bgzf<TSpec> > &, TPageFrame &page, TBool const &)
{
    const unsigned MAX_BLOCK_SIZE = 64 * 1024;
    const unsigned BLOCK_HEADER_LENGTH = 18;
    const unsigned BLOCK_FOOTER_LENGTH = 8;
    const unsigned ZLIB_BLOCK_OVERHEAD = 5; // 5 bytes block overhead (see 3.2.4. at http://www.gzip.org/zlib/rfc-deflate.html)

    // Reduce the maximal input size, such that the compressed data always fits in one block even for level Z_NO_COMPRESSION.
    reserve(page.data, MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD);
    return true;
}

// ----------------------------------------------------------------------------
// Function _postprocessFilePage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDirection, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_postprocessFilePage(FilePageTable<TValue, TDirection, Bgzf<TSpec> > &, TPageFrame &page, TBool const &)
{
    static int compressionLevels[] = { Z_DEFAULT_COMPRESSION, Z_BEST_COMPRESSION, Z_BEST_SPEED, Z_NO_COMPRESSION };

    const unsigned MAX_BLOCK_SIZE = 64 * 1024;
    const unsigned BLOCK_HEADER_LENGTH = 18;
    const unsigned BLOCK_FOOTER_LENGTH = 8;

    const int GZIP_WINDOW_BITS = -15; // no zlib header
    const int Z_DEFAULT_MEM_LEVEL = 8;

    // Make sure there is enough space in the buffer for compressed data.
    reserve(page.raw, MAX_BLOCK_SIZE);
    resize(page.raw, 0);

    Bytef *header = static_cast<Bytef *>(static_cast<void *>(begin(page.raw, Standard())));

    // Init gzip header
    header[0] = MagicHeader<BgzfFile>::VALUE[0];
    header[1] = MagicHeader<BgzfFile>::VALUE[1];
    header[2] = MagicHeader<BgzfFile>::VALUE[2];
    header[3] = 4;      // FLG_FEXTRA;
    header[4] = 0;      // mtime
    header[5] = 0;
    header[6] = 0;
    header[7] = 0;
    header[8] = 0;
    header[9] = -1;     // OS_UNKNOWN;
    header[10] = 6;     // BGZF_LEN+4
    header[11] = 0;
    header[12] = 'B';
    header[13] = 'C';
    header[14] = 2;     // BGZF_LEN
    header[15] = 0;
    header[16] = 0;     // Placeholder for block length.
    header[17] = 0;

    // Loop to retry for blocks that do not compress enough.
    for (unsigned clIdx = 0; clIdx < sizeof(compressionLevels) / sizeof(int); ++clIdx)
    {
        z_stream zs;
        zs.zalloc = NULL;
        zs.zfree = NULL;
        zs.next_in = static_cast<Bytef *>(static_cast<void *>(begin(page.data, Standard())));
        zs.avail_in = length(page.data);
        zs.next_out = header + BLOCK_HEADER_LENGTH;
        zs.avail_out = capacity(page.raw) - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

        int status = deflateInit2(&zs, compressionLevels[clIdx], Z_DEFLATED,
                                  GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
        if (status != Z_OK)
            throw IOException("BGZF deflateInit2() failed.");

        status = deflate(&zs, Z_FINISH);
        bool rawDataTooBig = (status != Z_STREAM_END);

        status = deflateEnd(&zs);
        if (status != Z_OK)
            throw IOException("BGZF deflateEnd() failed.");

        if (!rawDataTooBig)
        {
            resize(page.raw, zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH);
            break;
        }
        // If we are here, compression was too low. Try another compression level.
    }
    if (empty(page.raw))
    {
        // Should never happen.
        throw IOException("Deflation failed. Compressed BGZF data is bigger than uncompressed data.");
    }

    // Set compressed length into buffer, compute CRC and write CRC into buffer.
    uLong crc = crc32(0L, NULL, 0L);
    crc = crc32(crc, static_cast<Bytef const *>(static_cast<void *>(begin(page.data, Standard()))), length(page.data));

    _bgzfPackInt16(header + 16, length(page.raw) - 1);
    _bgzfPackInt32(header + length(page.raw) - 8, crc);
    _bgzfPackInt32(header + length(page.raw) - 4, length(page.data));

    return true;
}

template <typename TValue, typename TSpec, typename TPageFrame, typename TBool>
inline bool
_postprocessFilePage(FilePageTable<TValue, Input, Bgzf<TSpec> > &, TPageFrame &, TBool const &)
{
    // Do nothing.
    return true;
}

} // namespace seqan

#endif  // #ifndef SEQAN_FILE_FILE_BGZF_H_
