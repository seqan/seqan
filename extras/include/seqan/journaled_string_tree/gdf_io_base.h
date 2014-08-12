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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Basic io data structures and methods for gdf format.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_BASE_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_BASE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(rmaerker): Rename the tag.
struct Gdf_;
typedef Tag<Gdf_> Gdf;

// ----------------------------------------------------------------------------
// Class GdfIO
// ----------------------------------------------------------------------------

struct GdfIO
{
    enum FILE_VERSION
    {
        FILE_VERSION_BIG = 0,
        FILE_VERSION_LITTLE = 1
    };

    enum ERROR
    {
        PARSE_OK = 0,
        UNSUPPORTED_PREFIX_ERROR = 1,
        UNSUPPORTED_FILE_VERSION_ERROR = 2,
        UNSUPPORTED_FILE_FORMAT_ERROR = 3,
        UNSUPPORTED_REFERNCE_INFORMATION_ERROR = 4
    };

    static const CharString FILE_VERSION_VALUE_PREFIX;// = "GDFv";
    static const CharString FILE_VERSION_VALUE_SEPARATOR;// = ".";
    static const CharString FILE_VERSION_KEY; //"FileVersion";
    static const CharString FILE_ENDIANNESS_KEY;  //"Endianness";
    static const CharString FILE_ENDIANNESS_LITTLE;  //"littleEndian";
    static const CharString FILE_ENDIANNESS_BIG;  //"bigEndian";
    static const CharString FILE_SNP_COMPRESSION_KEY;
    static const CharString FILE_SNP_COMPRESSION_2BIT;
    static const CharString FILE_SNP_COMPRESSION_GENERIC;
    static const CharString FILE_BLOCKSIZE_KEY;
    static const CharString REFERENCE_ID_KEY; //"Reference";
    static const CharString REFERENCE_FILE_KEY; //"ReferenceLocation";
    static const CharString REFERENCE_HASH_KEY; //"ReferenceHash";
    static const CharString HEADER_PREFIX;// = "##";
    static const CharString SEQ_NAMES_PREFIX;// = "!!";
    static const CharString SEQ_NAMES_SEPARATOR;// = ",";
    static const CharString KEY_VALUE_SEPARATOR;// = "=";
};

const CharString GdfIO::FILE_VERSION_VALUE_PREFIX = "GDFv";
const CharString GdfIO::FILE_VERSION_VALUE_SEPARATOR = ".";
const CharString GdfIO::FILE_VERSION_KEY = "file version";
const CharString GdfIO::FILE_ENDIANNESS_KEY = "endianness";
const CharString GdfIO::FILE_ENDIANNESS_LITTLE = "little endian";
const CharString GdfIO::FILE_ENDIANNESS_BIG = "big endian";
const CharString GdfIO::FILE_SNP_COMPRESSION_KEY = "compression alphabet";
const CharString GdfIO::FILE_SNP_COMPRESSION_2BIT = "2bit";
const CharString GdfIO::FILE_SNP_COMPRESSION_GENERIC = "generic";
const CharString GdfIO::FILE_BLOCKSIZE_KEY = "block size";
const CharString GdfIO::REFERENCE_ID_KEY = "reference ID";
const CharString GdfIO::REFERENCE_FILE_KEY = "reference filename";
const CharString GdfIO::REFERENCE_HASH_KEY = "reference hash";
const CharString GdfIO::HEADER_PREFIX = "##";
const CharString GdfIO::KEY_VALUE_SEPARATOR = "=";
const CharString GdfIO::SEQ_NAMES_PREFIX = "!!";
const CharString GdfIO::SEQ_NAMES_SEPARATOR = ",";

// ----------------------------------------------------------------------------
// Struct SystemByteOrder
// ----------------------------------------------------------------------------

// Checks the byte order of the running system.
struct SystemByteOrder
{
   static bool IS_LITTLE_ENDIAN()
   {
       union  {
           __int16 i;
           char c[2];
       } tmp;
       tmp.i = 0x006C;
       return tmp.c[0] == 'l';
   }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Struct GdfFileBlockSize_
// ----------------------------------------------------------------------------

// Number of variants before the block ends.
// Which also means we write the coverage into the blocks.

template <typename T>
struct GdfFileBlockSize_
{
    enum
    {
        VALUE = 1000000
    };
};

template <typename TValue>
struct MsbIndex
{
    enum MSB_INDEX
    {
        VALUE = sizeof(TValue) << 3
    };
};

//template <typename T>
//const T SetHighBit<T>::VALUE = static_cast<T>(1) << (sizeof(T) * 8 - 1);

template <typename TValue>
struct SetSnpBits
{
    static const TValue VALUE;
};

template <typename T>
const T SetSnpBits<T>::VALUE = static_cast<T>(1) << (sizeof(T) * 8 - 3);

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _arrayMoveForwarReverse() using pointers
// ----------------------------------------------------------------------------

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayMoveForwardReverseDefault(TSource1 source_begin,
                                TSource2 source_end,
                                TTarget target_begin)
{
    while (source_end != source_begin)
    {
        move(*target_begin, *(--source_end));
        ++target_begin;
    }
}

template<typename TValue>
inline void
_arrayMoveForwardReversePointer(TValue * source_begin,
                                TValue * source_end,
                                TValue * target_begin,
                                True)
{
    ::std::reverse(source_begin, source_end);
    ::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}

template<typename TValue>
inline void
_arrayMoveForwardReversePointer(TValue * source_begin,
                                TValue * source_end,
                                TValue * target_begin,
                                False)
{
    _arrayMoveForwardReverseDefault(source_begin, source_end, target_begin);
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void
arrayMoveForwardReverse(TSource1 source_begin,
                        TSource2 source_end,
                        TTarget target_begin)
{
    _arrayMoveForwardReverseDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void
arrayMoveForwardReverse(TValue * source_begin,
                        TValue * source_end,
                        TValue * target_begin)
{
    _arrayMoveForwardReversePointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_BASE_H_
