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
// Implements the write methods to write out the journal strings as
// binary version.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_JOURNALED_STRING_TREE_DELTA_MAP_IO_WRITE_H_
#define EXTRAS_INCLUDE_JOURNALED_STRING_TREE_DELTA_MAP_IO_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

static const char DIFF_FILE_SEPARATOR = '\t';

#define GDF_IO_WRITE_KEY_VALUE(stream,key,...) stream << GDF_IO_HEADER_PREFIX << key << GDF_IO_KEY_VALUE_SEPARATOR << (__VA_ARGS__) << '\n'

// ----------------------------------------------------------------------------
// Class ByteCounterHelper_
// ----------------------------------------------------------------------------

template <typename TMapIterator, typename TBitCompression>
struct ByteCounterHelper_
{
    unsigned byteCount;
    TMapIterator it;

    ByteCounterHelper_() : byteCount(0), it()
    {}

    template <typename TTag>
    inline bool operator()(TTag /*tag*/)
    {
        if (isDeltaType(deltaType(it), TTag()))
        {
            byteCount += _countBytes(it, TTag(), TBitCompression());
            return true;
        }
        return false;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _countBytes()                                        [DeltaTypeSnp]
// ----------------------------------------------------------------------------

template <typename TIterator>
inline unsigned _countBytes(TIterator /*it*/, DeltaTypeSnp /*deltaType*/, GdfIO2BitSnpCompression /*tag*/)
{
    return 4;  // DeltaPosition and Snp value are compressed in 4 bytes.
}

template <typename TIterator>
inline unsigned _countBytes(TIterator /*it*/, DeltaTypeSnp /*deltaType*/, GdfIOGenericSnpCompression /*tag*/)
{
    return 5;  // DeltaPosition + 1 byte for the SNP.
}

// ----------------------------------------------------------------------------
// Function _countBytes()                                        [DeltaTypeDel]
// ----------------------------------------------------------------------------

template <typename TIterator, typename TTag>
inline unsigned _countBytes(TIterator /*it*/, DeltaTypeDel /*deltaType*/, TTag /*tag*/)
{
    return 8;  // 4 Bytes DeltaPosition and 4 Bytes for the deletion size.
}

// ----------------------------------------------------------------------------
// Function _countBytes()                                        [DeltaTypeIns]
// ----------------------------------------------------------------------------

template <typename TIterator, typename TTag>
inline unsigned _countBytes(TIterator const & it, DeltaTypeIns /*deltaType*/, TTag /*tag*/)
{
    return 8 + length(deltaValue(it, DeltaTypeIns()));  // 4 Bytes DeltaPosition and 4 Bytes for the insertion size + number of stored values.
}

// ----------------------------------------------------------------------------
// Function _countBytes()                                         [DeltaTypeSV]
// ----------------------------------------------------------------------------

template <typename TIterator, typename TTag>
inline unsigned _countBytes(TIterator const & it, DeltaTypeSV /*deltaType*/, TTag /*tag*/)
{
    return 12 + length(deltaValue(it, DeltaTypeSV()).i2);  // 4 Bytes DeltaPosition + 4 Bytes for the deletion + 4 Bytes for the insertion + number of stored values.
}

// ----------------------------------------------------------------------------
// Function _encodeBitVector()
// ----------------------------------------------------------------------------

template <typename TTarget, typename THostSpec>
inline void
_encodeBitVector(TTarget & target,
                 String<bool, Packed<THostSpec> > const & bv)
{
    typedef String<bool, Packed<THostSpec> >  TPackedString;
    typedef typename Position<TPackedString>::Type TPosition;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Iterator<TPackedHost const, Standard>::Type TConstPackedHostIterator;
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename TTraits::THostValue THostValue;
    typedef typename THostValue::TBitVector TBitVector;
    typedef typename Value<TTarget>::Type TEncodedValue;
    typedef Pair<TEncodedValue, TEncodedValue> TTargetValue;

    if (empty(host(bv)))  // Simply return unmodified target.
        return;

    TConstPackedHostIterator itBlockBegin = begin(host(bv), Standard()) + 1;
    TConstPackedHostIterator it = itBlockBegin;
    TConstPackedHostIterator itEnd = end(host(bv), Standard()) - 1;

    TPosition removedBits, zeroCount, oneCount = 0;
    TTargetValue val(0,0);
    bool lastBitSet = false;

    while (true)
    {
        // Check for zero entries.
        if (lastBitSet)
            for (; it != itEnd && testAllZeros(it->i); ++it)
            {}
        else
            for (; it != itEnd && testAllZeros(it->i); ++it)
            {}

        val.i1 += ((it - itBlockBegin) * BitsPerValue<TBitVector>::VALUE);
        removedBits = 0;
        TBitVector tmp = (it != itEnd) ? it->i : it->i &
                         (~static_cast<TBitVector>(0) <<  (BitsPerValue<TBitVector>::VALUE -
                                                          (length(bv) % BitsPerValue<TBitVector>::VALUE)));
        while (true)  // Process the current word.
        {
            if (testAllZeros(tmp))  // Check if no bit is set.
            {
                val.i1 += BitsPerValue<TBitVector>::VALUE - removedBits;
                lastBitSet = false;
                break;
            }

            zeroCount = BitsPerValue<TBitVector>::VALUE - 1 - bitScanReverse(tmp);
            val.i1 += zeroCount;

            // All bits flipped and leading zeros are shifted out of the word.
            tmp = ~(tmp << zeroCount);
            removedBits += zeroCount;

            if (testAllZeros(tmp))  // Check if no bit is set -> Hence all bits are set, since the word is negated.
            {  // Add number of ones to current block.
                val.i2 += BitsPerValue<TBitVector>::VALUE;
                lastBitSet = true;
                break;
            }
            // Get the block of ones.
            oneCount = BitsPerValue<TBitVector>::VALUE - 1 - bitScanReverse(tmp);
            if (removedBits + oneCount == BitsPerValue<TBitVector>::VALUE)  // Stop here since everything is set to ones.
            {
                val.i2 += oneCount;
                lastBitSet = true;
                break;
            }

            val.i2 += oneCount;
            SEQAN_ASSERT_LT(val.i1, MaxValue<TEncodedValue>::VALUE);
            SEQAN_ASSERT_LT(val.i2, MaxValue<TEncodedValue>::VALUE);
            appendValue(target, val.i1);
            appendValue(target, val.i2);
            val = TTargetValue(0,0);
            tmp = (~tmp) << oneCount;
            removedBits += oneCount;
        }

        ++it;
        itBlockBegin = it;

        if (it == end(host(bv), Standard()))
        {
            if (lastBitSet)
            {
                SEQAN_ASSERT_LT(val.i1, MaxValue<TEncodedValue>::VALUE);
                SEQAN_ASSERT_LT(val.i2, MaxValue<TEncodedValue>::VALUE);
                appendValue(target, val.i1);
                appendValue(target, val.i2);
            }
            break;
        }
        // Check condition of previous value.
        if (lastBitSet && !isBitSet(it->i, BitsPerValue<TBitVector>::VALUE - 1))
        {
            SEQAN_ASSERT_LT(val.i1, MaxValue<TEncodedValue>::VALUE);
            SEQAN_ASSERT_LT(val.i2, MaxValue<TEncodedValue>::VALUE);
            appendValue(target, val.i1);
            appendValue(target, val.i2);
            val = TTargetValue(0,0);
        }
    }
    appendValue(target, MaxValue<TEncodedValue>::VALUE);  // Add stop signal
}

// ----------------------------------------------------------------------------
// Function _encodeAndComputeByteCount()
// ----------------------------------------------------------------------------

template <typename TStringSet, typename TIterator, typename TSnpCompression>
inline unsigned int
_encodeAndComputeByteCount(TStringSet & set,
                           TIterator const & itBlockBegin,
                           TIterator const & itBlockEnd,
                           TSnpCompression /*tag*/)
{
    typedef typename Value<TStringSet>::Type TSegment;
    typedef typename Host<TSegment>::Type TEncodedString;
    typedef typename Value<TEncodedString>::Type TEncodeValue;

    ByteCounterHelper_<TIterator, TSnpCompression> counterHelper;
    counterHelper.it = itBlockBegin;

    TEncodedString target;
    for (; counterHelper.it != itBlockEnd; ++counterHelper.it)
    {
        clear(target);
        tagApply(counterHelper, DeltaTypes());
        _encodeBitVector(target, deltaCoverage(counterHelper.it));
        appendValue(set, target);
    }
    return counterHelper.byteCount + length(concat(set)) * sizeof(TEncodeValue);
}

// ----------------------------------------------------------------------------
// Function _writeGdfFileConfigInfo()
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue>
inline void
_writeGdfFileConfigInfo(TStream & stream, GdfFileConfiguration<TValue> const & config)
{
    // Write the file version information.
    std::stringstream fileVersion;
    fileVersion << GDF_IO_FILE_VERSION_VALUE_PREFIX << GDF_IO_FILE_VERSION_BIG << GDF_IO_FILE_VERSION_VALUE_SEPARATOR << GDF_IO_FILE_VERSION_LITTLE;
    GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_VERSION_KEY, fileVersion.str());

    // Write system's endianness.
    if (IsLittleEndian::VALUE)
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_ENDIANNESS_KEY, GDF_IO_FILE_ENDIANNESS_LITTLE);
    else
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_ENDIANNESS_KEY, GDF_IO_FILE_ENDIANNESS_BIG);
    // Write Compression mode.
    if (config.compressionMode == GdfIOMode::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION)
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_SNP_COMPRESSION_KEY, GDF_IO_FILE_SNP_COMPRESSION_2BIT);
    else
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_SNP_COMPRESSION_KEY, GDF_IO_FILE_SNP_COMPRESSION_GENERIC);

    // Write Coverage Compression Mode
    GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_COVERAGE_COMPRESSION, config.coverageCompression);
}

// ----------------------------------------------------------------------------
// Function _writeGdfHeaderReferenceInfo()
// ----------------------------------------------------------------------------

template <typename TStream, typename TConfig>
inline void
_writeGdfHeaderReferenceInfo(TStream & stream, GdfHeader const & jseqHeader, TConfig const & config)
{
    // Write the reference id.
      GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_REFERENCE_ID_KEY, jseqHeader.referenceId);
    // Write the reference file.
      GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_REFERENCE_FILE_KEY, jseqHeader.referenceFilename);
      // Write the reference hash.
      GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_REFERENCE_HASH_KEY, config.refHash);
}

// ----------------------------------------------------------------------------
// Function _writeBitVector()
// ----------------------------------------------------------------------------

// @deprecated.
template <typename TStream, typename TSpec>
inline void
_writeBitVector(TStream & stream, String<bool, Packed<TSpec> > const & bitVec)
{
    typedef String<bool, Packed<TSpec> > const TConstBitVec;
    typedef typename Host<TConstBitVec>::Type THost;
    typedef typename Value<THost>::Type THostValue;
    typedef typename THostValue::TBitVector TBitVector;
    typedef typename Iterator<THost>::Type THostIterator;

    THostIterator it = begin(host(bitVec));
    THostIterator itEnd = end(host(bitVec));

    for (; it != itEnd; ++it)
        streamWriteBlock(stream, reinterpret_cast<const char *>(&it->i), sizeof(TBitVector));
}

// ----------------------------------------------------------------------------
// Function _writeGdfHeader()
// ----------------------------------------------------------------------------

template <typename TStream, typename TConfig>
inline void _writeGdfHeader(TStream & stream, GdfHeader const & header, TConfig const & config)
{
    // Write the file information.
    _writeGdfFileConfigInfo(stream, config);
    // Write reference information.
    _writeGdfHeaderReferenceInfo(stream, header, config);

    // Write additional data.
    stream << GDF_IO_SEQ_NAMES_PREFIX;
    for (unsigned i = 0; i < length(header.nameStore); ++i)
        stream << header.nameStore[i] << GDF_IO_SEQ_NAMES_SEPARATOR;
    stream << '\n';
}

// ----------------------------------------------------------------------------
// Function _bufferSnp()
// ----------------------------------------------------------------------------

// Writes SNP in separat value after delta pos.
template <typename TBuffer, typename TValue>
inline void
_bufferSnp(TBuffer & buffer, BitCompressedDeltaPos_<GdfIOGenericSnpCompression> deltaPos, TValue val)
{
    buffer += _toBinary(buffer, deltaPos.toWord(), HostByteOrder(), BigEndian());
    *buffer = val;
    ++buffer;
}

// ----------------------------------------------------------------------------
// Function _bufferSnp()                                                  [Dna]
// ----------------------------------------------------------------------------

// Encodes SNP in delta position assuming, that 28 bits are sufficient to store the delta to the previous variant.
template <typename TBuffer, typename TValue>
inline void
_bufferSnp(TBuffer & buffer, BitCompressedDeltaPos_<GdfIO2BitSnpCompression> deltaPos, TValue val)
{
    deltaPos.snp = val;
    buffer += _toBinary(buffer, deltaPos.toWord(), HostByteOrder(), BigEndian());
}

// ----------------------------------------------------------------------------
// Function _bufferDataBlock()
// ----------------------------------------------------------------------------

template <typename TBuffer, typename TDeltaMapIter, typename TValue, typename TAlphabet, typename TSpec,
          typename TCompressionMode>
inline void _bufferDataBlock(TBuffer & buffer,
                             TDeltaMapIter & it,
                             TDeltaMapIter & itEnd,
                             DeltaMap<TValue, TAlphabet, TSpec> const & /*deltaMap*/,
                             TCompressionMode /*tag*/)
{
    typedef BitCompressedDeltaPos_<TCompressionMode> TDeltaPos;

    __uint32 lastRefPos = deltaPosition(it);
    // Write the reference offset of the current block.
    buffer += _toBinary(buffer, lastRefPos);

    TDeltaMapIter itDelta = it;
    while(itDelta != itEnd)
    {
        SEQAN_ASSERT_LEQ(deltaPosition(itDelta) - lastRefPos, +ValueSize<TDeltaPos>::VALUE);
        // Store the reference position.
        TDeltaPos deltaPos(deltaPosition(itDelta) - lastRefPos);

        // Write SNP Data.
        if (deltaType(itDelta) == DELTA_TYPE_SNP)
        {
            deltaPos.isSnp = 1;
            _bufferSnp(buffer, deltaPos, deltaValue(itDelta, DeltaTypeSnp()));
        }
        else  // Write Indel
        {
            buffer += _toBinary(buffer, deltaPos.toWord(), HostByteOrder(), BigEndian());

            if (deltaType(itDelta) == DELTA_TYPE_DEL)  // Write deletion.
            {
                buffer += _toBinary(buffer, BitCompressedInDel_(1, 0, deltaValue(itDelta, DeltaTypeDel())).toWord(),
                                    HostByteOrder(), BigEndian());
            }
            else
            {
                // Handle Indel.
                if (deltaType(itDelta) == DELTA_TYPE_SV)
                {
                    buffer += _toBinary(buffer,
                                        BitCompressedInDel_(0, 1, deltaValue(itDelta, DeltaTypeSV()).i1).toWord(),
                                        HostByteOrder(), BigEndian());
                    buffer += _toBinary(buffer, static_cast<__uint32>(length(deltaValue(itDelta, DeltaTypeSV()).i2)),
                                        HostByteOrder(), BigEndian());
                    buffer += _copyToBuffer(buffer, &deltaValue(itDelta, DeltaTypeSV()).i2[0],
                                            length(deltaValue(itDelta, DeltaTypeSV()).i2));
                }
                else  // Handle Insertion.
                {
                    SEQAN_ASSERT(deltaType(itDelta) == DELTA_TYPE_INS);
                    buffer += _toBinary(buffer,
                                        BitCompressedInDel_(0, 0, length(deltaValue(itDelta, DeltaTypeIns()))).toWord(),
                                        HostByteOrder(), BigEndian());
                    buffer += _copyToBuffer(buffer, &deltaValue(itDelta, DeltaTypeIns())[0],
                                            length(deltaValue(itDelta, DeltaTypeIns())));
                }
            }
        }
        lastRefPos = deltaPosition(itDelta);
        ++itDelta;
    }
}

// ----------------------------------------------------------------------------
// Function _bufferCoverages()
// ----------------------------------------------------------------------------

template <typename TBuffer, typename TStringSet>
inline void _bufferCoverages(TBuffer & buffer, TStringSet const & set)
{
    typedef typename Concatenator<TStringSet>::Type TConcatenate;
    typedef typename Iterator<TConcatenate const, Standard>::Type TConcatIterator;

    for (TConcatIterator it = begin(concat(set), Standard()); it != end(concat(set), Standard()); ++it)
        buffer += _toBinary(buffer, *it);
}

// ----------------------------------------------------------------------------
// Function _writeGdfData()
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TAlphabet, typename TSpec, typename TConfig,
          typename TSnpCompression, typename TCoverageEncodingType>
inline void _writeGdfData(TStream & stream,
                          DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
                          GdfFileConfiguration<TConfig> const & config,
                          TSnpCompression /*tag*/,
                          GdfIOCoverageCompression<TCoverageEncodingType> /*tag*/)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Size<TDeltaMap>::Type TSize;
    typedef typename Iterator<TDeltaMap const, Standard>::Type TIterator;

    typedef String<TCoverageEncodingType> TEncodedString;
    typedef StringSet<TEncodedString, Owner<ConcatDirect<> > > TEncodedSet;

    TSize maxNumOfNodes = length(deltaMap);
    __uint32 numOfBlocks = (maxNumOfNodes + config.blockSize - 1) / config.blockSize;

    // Write the block containing the number of blocks to read.
    char tmp[sizeof(numOfBlocks)];
    _toBinary(tmp, numOfBlocks);
    streamWriteBlock(stream, tmp, sizeof(numOfBlocks));

    TEncodedSet set;

    for (__uint32 i = 0; i < numOfBlocks; ++i)  // For each block:
    {
        TIterator it = begin(deltaMap, Standard()) + (config.blockSize * i);
        TIterator itEnd = _min(end(deltaMap, Standard()), it + config.blockSize);
        clear(set);

        // Encode the bit vectors and compute buffer size.
        unsigned int bufferSize = _encodeAndComputeByteCount(set, it, itEnd, TSnpCompression()) + (sizeof(__uint32) * 3);
        unsigned int chunkCoverage = length(concat(set)) * sizeof(TCoverageEncodingType);
        unsigned int chunkDelta = bufferSize - chunkCoverage - (sizeof(__uint32) * 2);
        char * buffer = new char [bufferSize];
        buffer += _toBinary(buffer, chunkDelta);  // Write the buffer size of delta block including the reference position.
        _bufferDataBlock(buffer, it, itEnd, deltaMap, TSnpCompression());  // Buffer the variants.
        buffer += _toBinary(buffer, chunkCoverage);  // Write the buffer size of the coverage block.
        _bufferCoverages(buffer, set);  // Buffer the encoded bit vectors.

        buffer -= bufferSize;  // Set back to beginning of buffer;
        streamWriteBlock(stream, buffer, bufferSize);  // Write buffer to stream.
        delete[] buffer;
    }
}

// ----------------------------------------------------------------------------
// Function _writeGdfData()
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TAlphabet, typename TSpec, typename TConfig,
          typename TSnpCompressionMode>
inline void _writeGdfData(TStream & stream,
                          DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
                          GdfFileConfiguration<TConfig> const & config,
                          TSnpCompressionMode /*tag*/)
{
    switch (config.coverageCompression)
    {
        case GdfIOMode::COVERAGE_COMPRESSION_1_BYTE_PER_VALUE:
            _writeGdfData(stream, deltaMap, config, TSnpCompressionMode(), GdfIOCoverageCompression<__uint8>()); break;
        case GdfIOMode::COVERAGE_COMPRESSION_2_BYTE_PER_VALUE:
            _writeGdfData(stream, deltaMap, config, TSnpCompressionMode(), GdfIOCoverageCompression<__uint16>()); break;
        case GdfIOMode::COVERAGE_COMPRESSION_4_BYTE_PER_VALUE:
            _writeGdfData(stream, deltaMap, config, TSnpCompressionMode(), GdfIOCoverageCompression<__uint32>()); break;
        default:
            _writeGdfData(stream, deltaMap, config, TSnpCompressionMode(), GdfIOCoverageCompression<__uint64>()); break;
    }
}

template <typename TStream, typename TValue, typename TAlphabet, typename TSpec, typename TConfig>
inline void _writeGdfData(TStream & stream,
                          DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
                          GdfFileConfiguration<TConfig> const & config)
{
    if (config.compressionMode == GdfIOMode::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION)
        _writeGdfData(stream, deltaMap, config, GdfIO2BitSnpCompression());
    else
        _writeGdfData(stream, deltaMap, config, GdfIOGenericSnpCompression());
}

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------

/*!
 * @fn GdfIO#write
 * @brief Writes the specified @link DeltaMap @endlink to disk in the GDF file format.
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @signature write(stream, deltaMap, header, config, tag);
 *
 * @param[out] stream   The opened file stream to write the delta map to.
 * @param[in]  deltaMap The delta map to write to file.
 * @param[in]  header   The header containing the header information. Of type @link GdfHeader @endlink.
 * @param[in]  config   The config object used to invoke file dependent configurations. Of type @link GdfFileConfiguration @endlink.
 * @param[in]  tag      The tag to determine the correct file format. Must be of type @link GdfIO#Gdf @endlink.
 *
 * @throw GdfIOException The exception type thrown in case of an unhandled runtime error.
 *
 * @see GdfIO#read
 */

template <typename TStream, typename TValue, typename TAlphabet, typename TConfig>
inline void
write(TStream & stream,
      DeltaMap<TValue, TAlphabet> const & deltaMap,
      GdfHeader const & gdfHeader,
      TConfig const & config,
      Gdf const & /*tag*/)
{
    SEQAN_TRY
    {
        _writeGdfHeader(stream, gdfHeader, config);
        _writeGdfData(stream, deltaMap, config);
    }
    SEQAN_CATCH(Exception e)
    {
        SEQAN_THROW(GdfIOException(e.what()));
    }
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_JOURNALED_STRING_TREE_DELTA_MAP_IO_WRITE_H_
