// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

struct Gdf_;
typedef Tag<Gdf_> Gdf;

struct GdfIOException : public Exception
{
    std::string message;

    const char* what() const noexcept
    {
        return message.c_str();
    }

    GdfIOException(std::string const & errMessage)
    {
        message = "Gdf_IO_Exception: (" + errMessage + ")";
    }
};

struct GdfIOWrongReferenceException : public GdfIOException
{

    GdfIOWrongReferenceException() : GdfIOException("Wrong reference")
    {}
};

// ----------------------------------------------------------------------------
// Class GdfIO
// ----------------------------------------------------------------------------

struct GdfIO
{
    // TODO(rmaerker): Replace by exceptions.
    enum ERROR
    {
        PARSE_OK = 0,
        UNSUPPORTED_PREFIX_ERROR = 1,
        UNSUPPORTED_FILE_VERSION_ERROR = 2,
        UNSUPPORTED_FILE_FORMAT_ERROR = 3,
        UNSUPPORTED_REFERNCE_INFORMATION_ERROR = 4
    };

    enum CompressionMode
    {
        COMPRESSION_MODE_2_BIT_SNP_COMPRESSION,
        COMPRESSION_MODE_NO_SNP_COMPRESSION
    };

    enum ReferenceMode
    {
        REFERENCE_MODE_WRITE_ENABLED,  // Writes the reference under the given filename.
        REFERENCE_MODE_WRITE_DISABLED  // Dosen't write the reference.
    };

//    static const TKeyType FILE_VERSION_VALUE_PREFIX;// = "GDFv";
//    static const TKeyType FILE_VERSION_VALUE_SEPARATOR;// = ".";
//    static const TKeyType FILE_VERSION_KEY;// = "file version";
//    static const TKeyType FILE_ENDIANNESS_KEY;// = "endianness";
//    static const TKeyType FILE_ENDIANNESS_LITTLE;// = "little endian";
//    static const TKeyType FILE_ENDIANNESS_BIG;// = "big endian";
//    static const TKeyType FILE_SNP_COMPRESSION_KEY;// = "compression alphabet";
//    static const TKeyType FILE_SNP_COMPRESSION_2BIT;// = "2bit";
//    static const TKeyType FILE_SNP_COMPRESSION_GENERIC;// = "generic";
//    static const TKeyType FILE_BLOCKSIZE_KEY;// = "block size";
//    static const TKeyType REFERENCE_ID_KEY;// = "reference ID";
//    static const TKeyType REFERENCE_FILE_KEY;// = "reference filename";
//    static const TKeyType REFERENCE_HASH_KEY;// = "reference hash";
//    static const TKeyType HEADER_PREFIX;// = "##";
//    static const TKeyType SEQ_NAMES_PREFIX;// = "!!";
//    static const TKeyType SEQ_NAMES_SEPARATOR;// = ",";
//    static const TKeyType KEY_VALUE_SEPARATOR;// = "=";
};

#define GDF_IO_FILE_VERSION_VALUE_PREFIX "GDFv"
#define GDF_IO_FILE_VERSION_VALUE_SEPARATOR "."
#define GDF_IO_FILE_VERSION_KEY "file version"
#define GDF_IO_FILE_ENDIANNESS_KEY "endianness"
#define GDF_IO_FILE_ENDIANNESS_LITTLE "little endian"
#define GDF_IO_FILE_ENDIANNESS_BIG "big endian"
#define GDF_IO_FILE_SNP_COMPRESSION_KEY "compression alphabet"
#define GDF_IO_FILE_SNP_COMPRESSION_2BIT "2bit"
#define GDF_IO_FILE_SNP_COMPRESSION_GENERIC "generic"
#define GDF_IO_FILE_BLOCKSIZE_KEY "block size"
#define GDF_IO_REFERENCE_ID_KEY "reference ID"
#define GDF_IO_REFERENCE_FILE_KEY "reference filename"
#define GDF_IO_REFERENCE_HASH_KEY "reference hash"
#define GDF_IO_HEADER_PREFIX "##"
#define GDF_IO_KEY_VALUE_SEPARATOR "="
#define GDF_IO_SEQ_NAMES_PREFIX "!!"
#define GDF_IO_SEQ_NAMES_SEPARATOR ","
#define GDF_IO_FILE_VERSION_BIG 0
#define GDF_IO_FILE_VERSION_LITTLE 1

//const GdfIO::TKeyType GdfIO::FILE_VERSION_VALUE_PREFIX = "GDFv";
//const GdfIO::TKeyType GdfIO::FILE_VERSION_VALUE_SEPARATOR = ".";
//const GdfIO::TKeyType GdfIO::FILE_VERSION_KEY = "file version";
//const GdfIO::TKeyType GdfIO::FILE_ENDIANNESS_KEY = "endianness";
//const GdfIO::TKeyType GdfIO::FILE_ENDIANNESS_LITTLE = "little endian";
//const GdfIO::TKeyType GdfIO::FILE_ENDIANNESS_BIG = "big endian";
//const GdfIO::TKeyType GdfIO::FILE_SNP_COMPRESSION_KEY = "compression alphabet";
//const GdfIO::TKeyType GdfIO::FILE_SNP_COMPRESSION_2BIT = "2bit";
//const GdfIO::TKeyType GdfIO::FILE_SNP_COMPRESSION_GENERIC = "generic";
//const GdfIO::TKeyType GdfIO::FILE_BLOCKSIZE_KEY = "block size";
//const GdfIO::TKeyType GdfIO::REFERENCE_ID_KEY = "reference ID";
//const GdfIO::TKeyType GdfIO::REFERENCE_FILE_KEY = "reference filename";
//const GdfIO::TKeyType GdfIO::REFERENCE_HASH_KEY = "reference hash";
//const GdfIO::TKeyType GdfIO::HEADER_PREFIX = "##";
//const GdfIO::TKeyType GdfIO::KEY_VALUE_SEPARATOR = "=";
//const GdfIO::TKeyType GdfIO::SEQ_NAMES_PREFIX = "!!";
//const GdfIO::TKeyType GdfIO::SEQ_NAMES_SEPARATOR = ",";

// ----------------------------------------------------------------------------
// Struct SystemByteOrder
// ----------------------------------------------------------------------------

// TODO(rmaerker): No compile time check possible with C'98 & C'11
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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function computeHash()
// ----------------------------------------------------------------------------

#ifdef __SSE4_2__
template <typename TReference>
inline unsigned int computeReferenceCrc(TReference const & ref)
{
    typedef typename Iterator<TReference const, Standard>::Type TIterator;

    unsigned int crc;
    for (TIterator it = begin(ref, Standard()); it != end(ref, Standard()); ++it)
        crc = _mm_crc32_u8(crc, static_cast<unsigned char>(*it));
    return crc;
}
#else  //__SSE4_2_
template <typename TReference>
inline unsigned int computeReferenceCrc(TReference const & /*ref*/)
{
    return 0;
}
#endif

#ifdef __SSE4_2__
template <typename TReference, typename TCrc>
inline bool checkReferenceCrc(TReference const & ref, TCrc crc)
{
    return computeReferenceCrc(ref) == crc;
}
#else  //__SSE4_2_
template <typename TReference, typename TCrc>
inline bool checkReferenceCrc(TReference const & /*ref*/, TCrc /*crc*/)
{
    return true;  // Always return true.
}
#endif  //__SSE4_2_

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


// ----------------------------------------------------------------------------
// Function readFromBinary()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Replace by standard write method when io-module is updated in develop.
template <typename TStream, typename TValue>
inline void readFromBinary(TStream & stream, TValue & val)
{
    streamWriteBlock(stream, reinterpret_cast<const char *>(&val), sizeof(val));
}

// ----------------------------------------------------------------------------
// Function writeToBinary()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Replace by standard write method when io-module is updated in develop.
template <typename TStream, typename TValue>
inline void writeBinary(TStream & stream, TValue const & val)
{
    streamWriteBlock(stream, reinterpret_cast<const char *>(&val), sizeof(val));
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_BASE_H_
