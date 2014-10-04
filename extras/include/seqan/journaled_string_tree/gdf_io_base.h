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

/*!
 * @defgroup GdfIO GDF I/O
 * @brief Functionality for GDF I/O.
 */

// ----------------------------------------------------------------------------
// Tag Gdf
// ----------------------------------------------------------------------------

/*!
 * @tag GdfIO#Gdf
 * @brief The tag used for GDF I/O functionaliy.
 * @headerfile <seqan/journaled_string_tree.h>
 * 
 * @signature struct Gdf_;
 *            typedef Tag<Gdf_> Gdf;
 */

struct Gdf_;
typedef Tag<Gdf_> Gdf;

// ----------------------------------------------------------------------------
// Class GdfIOException
// ----------------------------------------------------------------------------

/*!
 * @class GdfIOException
 * @extends Exception
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Exception thrown by errors during GDF I/O operations.
 *
 * @signature struct GdfIOException : public Exception;
 */

/*!
 * @fn GdfIOException#what
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Outputs the error message.
 *
 * @signature const char* e.what();
 * @parm e The exception.
 *
 * @return "const char*" The error message.
 */

struct GdfIOException
{
    /*!
     * @var std::string GdfIOException::message
     * @brief The error message.
     */
    std::string message;

    /*!
    * @fn GdfIOException::what
    * @brief Outputs the associated error message.
    *
    * @signature const char* GdfIOException::what();
    */
    const char* what()
    {
        return message.c_str();
    }

    /*!
    * @fn GdfIOException::GdfIOException
    * @brief The constructor.
    *
    * @signature GdfIOException()
    * @signature GdfIOException(message)
    * @param message The error message to be printed. Must be of type <a href="http://www.cplusplus.com/reference/string/string/?kw=string">std::string</a>.
    */
    GdfIOException()
    {}

    template <typename TString>
    GdfIOException(TString const & errMessage)
    {
        std::string tmp(errMessage);
        message = "Gdf_IO_Exception: (" + tmp + ")";
    }
};

struct GdfIOWrongReferenceException : public GdfIOException
{

    GdfIOWrongReferenceException() : GdfIOException("Wrong reference!")
    {}

    template <typename TString, typename TCrc>
    GdfIOWrongReferenceException(TString const & isId, TString const & shouldId, TCrc isCrc, TCrc shouldCrc) :
        GdfIOException()
    {
        std::stringstream errMessage;
        errMessage << "The id of the reference is \'" << isId << "\' but should be \'" << shouldId <<
                      "\' and the crc is \'" << isCrc << "\' but should be \'" << shouldCrc << "\'!";
        message = GdfIOException(errMessage.str()).message;
    }
};

// ----------------------------------------------------------------------------
// Class GdfIO
// ----------------------------------------------------------------------------

 /*!
 * @enum GdfIOMode::CompressionMode
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief SNP compression modes.
 *
 * @val GdfIOMode::CompressionMode GdfIOMode::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION
 * @brief Enables 2 bit compression of SNPs only if the SNP value fits into 2 bit.
 *
 * @val GdfIOMode::CompressionMode GdfIOMode::COMPRESSION_MODE_NO_SNP_COMPRESSION
 * @brief Disables snp compression.
 */

/*!
 * @enum GdfIOMode::CoverageCompression
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Coverage compression modes.
 *
 * @val GdfIOMode::CoverageCompression GdfIOMode::COVERAGE_COMPRESSION_1_BYTE_PER_VALUE
 * @brief Uses 1 byte sized host value if the length of the coverage is less than max. value of <tt>__uint8</tt>.
 *
 * @val GdfIOMode::CoverageCompression GdfIOMode::COVERAGE_COMPRESSION_2_BYTE_PER_VALUE
 * @brief Uses 2 byte sized host value if the length of the coverage is less than max. value of <tt>__uint16</tt>.
 *
 * @val GdfIOMode::CoverageCompression GdfIOMode::COVERAGE_COMPRESSION_4_BYTE_PER_VALUE
 * @brief Uses 4 byte sized host value if the length of the coverage is less than max. value of <tt>__uint32</tt>.
 *
 * @val GdfIOMode::CoverageCompression GdfIOMode::COVERAGE_COMPRESSION_8_BYTE_PER_VALUE
 * @brief Uses 8 byte sized host value if the length of the coverage is less than max. value of <tt>__uint64</tt>.
 */

/*!
 * @enum GdfIOMode::SaveReferenceMode
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Mode to enable or disable write of the reference sequence.
 *
 * @val GdfIOMode::SaveReferenceMode GdfIOMode::SAVE_REFERENCE_MODE_ENABLED
 * @brief Explicitly saves the reference at given filename.
 *
 * @val GdfIOMode::SaveReferenceMode GdfIOMode::SAVE_REFERENCE_MODE_DISABLED
 * @brief Reference sequence is not saved explicitly. Note only use this if you are sure, that the proper reference 
 *        sequence is already stored somewhere.
 */

struct GdfIOMode
{
    enum CompressionMode
    {
        COMPRESSION_MODE_2_BIT_SNP_COMPRESSION,
        COMPRESSION_MODE_NO_SNP_COMPRESSION
    };

    enum CoverageCompression
    {
        COVERAGE_COMPRESSION_1_BYTE_PER_VALUE = 1,
        COVERAGE_COMPRESSION_2_BYTE_PER_VALUE = 2,
        COVERAGE_COMPRESSION_4_BYTE_PER_VALUE = 4,
        COVERAGE_COMPRESSION_8_BYTE_PER_VALUE = 8
    };

    enum SaveReferenceMode
    {
        SAVE_REFERENCE_MODE_ENABLED,  // Writes the reference under the given filename.
        SAVE_REFERENCE_MODE_DISABLED  // Dosen't write the reference.
    };
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
#define GDF_IO_FILE_COVERAGE_COMPRESSION "coverage compression"
#define GDF_IO_REFERENCE_ID_KEY "reference ID"
#define GDF_IO_REFERENCE_FILE_KEY "reference filename"
#define GDF_IO_REFERENCE_HASH_KEY "reference hash"
#define GDF_IO_HEADER_PREFIX "##"
#define GDF_IO_KEY_VALUE_SEPARATOR "="
#define GDF_IO_SEQ_NAMES_PREFIX "!!"
#define GDF_IO_SEQ_NAMES_SEPARATOR ","
#define GDF_IO_FILE_VERSION_BIG 0
#define GDF_IO_FILE_VERSION_LITTLE 1

// ----------------------------------------------------------------------------
// Tag GdfIO2BitSnpCompression
// ----------------------------------------------------------------------------

struct GdfIO2BitSnpCompression_;
typedef Tag<GdfIO2BitSnpCompression_> GdfIO2BitSnpCompression;

struct GdfIOGenericSnpCompression_;
typedef Tag<GdfIOGenericSnpCompression_> GdfIOGenericSnpCompression;

// ----------------------------------------------------------------------------
// Tag GdfIOCoverageCompression
// ----------------------------------------------------------------------------

template <typename TType>
struct GdfIOCoverageCompression{};

// ----------------------------------------------------------------------------
// Class BitCompressedInDel_
// ----------------------------------------------------------------------------

struct BitCompressedInDel_
{
    __uint32 isDel : 1;
    __uint32 isSV : 1;
    __uint32 value : BitsPerValue<__uint32>::VALUE - 2;


    BitCompressedInDel_() : isDel(0), isSV(0), value(0)
    {}

    template <typename TValue>
    BitCompressedInDel_(bool _isDel, bool _isSV, TValue _val) : isDel(_isDel), isSV(_isSV), value(_val)
    {}

    inline __uint32 toWord()
    {
        return isDel << (BitsPerValue<__uint32>::VALUE - 1) | isSV << (BitsPerValue<__uint32>::VALUE - 2) | value;
    }

    inline void fromWord(__uint32 word)
    {
        isDel = isBitSet(word, BitsPerValue<__uint32>::VALUE - 1);
        isSV = isBitSet(word, BitsPerValue<__uint32>::VALUE - 2);
        value = word & (~static_cast<__uint32>(0) >> 2);
    }
};

// ----------------------------------------------------------------------------
// Class BitCompressedDeltaPos_                                       [Generic]
// ----------------------------------------------------------------------------

template <typename TCompressionMode>
struct BitCompressedDeltaPos_
{
    __uint32 isSnp : 1;
    __uint32 pos : BitsPerValue<__uint32>::VALUE - 1;

    BitCompressedDeltaPos_() : isSnp(0), pos(0)
    {}

    template <typename TPos>
    BitCompressedDeltaPos_(TPos _pos) : isSnp(0), pos(_pos)
    {}

    inline __uint32 toWord()
    {
        return isSnp << (BitsPerValue<__uint32>::VALUE - 1) | pos;
    }

    inline void fromWord(__uint32 word)
    {
        isSnp = isBitSet(word, BitsPerValue<__uint32>::VALUE - 1);
        pos = word & (~static_cast<__uint32>(0) >> 1);
    }
};

// ----------------------------------------------------------------------------
// Class BitCompressedDeltaPos_                                          [2Bit]
// ----------------------------------------------------------------------------

template <>
struct BitCompressedDeltaPos_<GdfIO2BitSnpCompression>
{
    __uint32 isSnp : 1;
    __uint32 snp   : 2;
    __uint32 pos   : BitsPerValue<__uint32>::VALUE - 3;

    BitCompressedDeltaPos_() : isSnp(0), snp(0), pos(0)
    {}

    template <typename TPos>
    BitCompressedDeltaPos_(TPos _pos) : isSnp(0), snp(0), pos(_pos)
    {}

    inline __uint32 toWord()
    {
        return isSnp << (BitsPerValue<__uint32>::VALUE - 1) | snp << (BitsPerValue<__uint32>::VALUE - 3) | pos;
    }

    inline void fromWord(__uint32 word)
    {
        isSnp = isBitSet(word, BitsPerValue<__uint32>::VALUE - 1);
        snp = (word >> (BitsPerValue<__uint32>::VALUE - 3)) & 3;
        pos = word & (~static_cast<__uint32>(0) >> 3);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ValueSize
// ----------------------------------------------------------------------------

template <typename TMode>
struct ValueSize<BitCompressedDeltaPos_<TMode> >
{
    static const __uint32 VALUE = static_cast<__uint32>(~0) >> 1;
};

template <>
struct ValueSize<BitCompressedDeltaPos_<GdfIO2BitSnpCompression> >
{
    static const __uint32 VALUE = static_cast<__uint32>(~0) >> 3;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _copyToBuffer()                                           [Generic]
// ----------------------------------------------------------------------------

template <typename TTargetValue, typename TSrcValue, typename TSize>
inline TSize _copyToBuffer(TTargetValue * target, TSrcValue const * source, TSize copySize)
{
    for (unsigned i = 0; i < copySize; ++i)
        target[i] = convert<TTargetValue>(source[i]);
    return copySize;
}

// ----------------------------------------------------------------------------
// Function _copyToBuffer()                                      [ConstUInt<1>]
// ----------------------------------------------------------------------------

inline unsigned _copyToBuffer(char * target, const char * source, ConstUInt<1> /*value size*/)
{
    *target = *source;
    return 1;
}

// ----------------------------------------------------------------------------
// Function _copyToBuffer()                                      [ConstUInt<2>]
// ----------------------------------------------------------------------------

inline unsigned _copyToBuffer(char * target, const char * source, ConstUInt<2> /*value size*/)
{
    _copyToBuffer(target, source, ConstUInt<1>());
    _copyToBuffer(target + 1, source + 1, ConstUInt<1>());
    return 2;
}

// ----------------------------------------------------------------------------
// Function _copyToBuffer()                                      [ConstUInt<4>]
// ----------------------------------------------------------------------------

inline unsigned _copyToBuffer(char * target, const char * source, ConstUInt<4> /*value size*/)
{
    _copyToBuffer(target, source, ConstUInt<2>());
    _copyToBuffer(target + 2, source + 2, ConstUInt<2>());
    return 4;
}

// ----------------------------------------------------------------------------
// Function _copyToBuffer()                                      [ConstUInt<8>]
// ----------------------------------------------------------------------------

inline unsigned _copyToBuffer(char * target, const char * source, ConstUInt<8> /*value size*/)
{
    _copyToBuffer(target, source, ConstUInt<4>());
    _copyToBuffer(target + 4, source + 4, ConstUInt<4>());
    return 8;
}

// ----------------------------------------------------------------------------
// Function _copyToBuffer()                                   [ConstUInt<SIZE>]
// ----------------------------------------------------------------------------

template <unsigned SIZE>
inline unsigned _copyToBuffer(char * target, const char * source, ConstUInt<SIZE> /*value size*/)
{
    return _copyToBuffer(target, source, SIZE);
}

// ----------------------------------------------------------------------------
// Function _copyFromBuffer()
// ----------------------------------------------------------------------------

template <typename TTargetValue, typename TSize>
inline TSize _copyFromBuffer(TTargetValue & target, const char * source, TSize copySize)
{
    for (unsigned i = 0; i < copySize; ++i)
        target[i] = convert<char>(source[i]);
    return copySize;
}

// ----------------------------------------------------------------------------
// Function _toBinary()
// ----------------------------------------------------------------------------

template <typename TBuffer, typename TValue, typename TFromByteOrder, typename TToByteOrder>
inline unsigned _toBinary(TBuffer & buffer, TValue val, TFromByteOrder /*tag*/, TToByteOrder /*tag*/)
{
    val = endianSwap(val, TFromByteOrder(), TToByteOrder());
    return _copyToBuffer(buffer, reinterpret_cast<const char *>(&val), ConstUInt<sizeof(TValue)>());
}

template <typename TBuffer, typename TValue>
inline unsigned _toBinary(TBuffer & buffer, TValue val)
{
    return _toBinary(buffer, val, HostByteOrder(), HostByteOrder());
}

// ----------------------------------------------------------------------------
// Function _fromBinary()
// ----------------------------------------------------------------------------

template <typename TValue, typename TBuffer, typename TFromByteOrder, typename TToByteOrder>
inline unsigned _fromBinary(TValue & val, TBuffer const & buffer, TFromByteOrder /*tag*/, TToByteOrder /*tag*/)
{
    _copyToBuffer(reinterpret_cast<char *>(&val), buffer, ConstUInt<sizeof(TValue)>());
    val = endianSwap(val, TFromByteOrder(), TToByteOrder());
    return sizeof(TValue);
}

template <typename TValue, typename TBuffer>
inline unsigned _fromBinary(TValue val, TBuffer const & buffer)
{
    return _fromBinary(val, buffer, HostByteOrder(), HostByteOrder());
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_BASE_H_
