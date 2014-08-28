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
// Implements the read methods to load journal sequence format.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_READ_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_READ_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TConfig, typename TStream>
inline void
_readGdfFileConfigInfo(GdfFileConfiguration<TConfig> & config,
                       RecordReader<TStream, SinglePass<> > & reader)
{
    CharString buffer;

    // Skip the file prefix.
    skipNChars(reader, 2);
    // Read file version.
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_FILE_VERSION_KEY));
    if (!startsWith(buffer, GDF_IO_FILE_VERSION_KEY))
        SEQAN_THROW(GdfIOException("Unsupported file version key!"));

    // Skip = sign
    skipNChars(reader, 1);
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_FILE_VERSION_VALUE_PREFIX));
    if (!startsWith(buffer, GDF_IO_FILE_VERSION_VALUE_PREFIX))
        SEQAN_THROW(GdfIOException("Unsupported file version prefix!"));
    // Read big version number.
    clear(buffer);
    int fileVersionBig;
    clear(buffer);
    readNChars(buffer, reader, 1);
    lexicalCast2(fileVersionBig, buffer);

    if (fileVersionBig != GDF_IO_FILE_VERSION_BIG)
        SEQAN_THROW(GdfIOException("Unsupported file version!"));
    // Read version number separator.
    skipNChars(reader, 1);
    // Read little version number.
    clear(buffer);
    int fileVersionLittle;
    clear(buffer);
    readNChars(buffer, reader, 1);
    lexicalCast2(fileVersionLittle, buffer);

    if (fileVersionLittle != GDF_IO_FILE_VERSION_LITTLE)
        SEQAN_THROW(GdfIOException("Unsupported file version!"));

    // Skip until next line.
    skipLine(reader);

    // Read byte order.
    skipNChars(reader, 2);  // skip prefix.
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_FILE_ENDIANNESS_KEY));
    if (buffer != GDF_IO_FILE_ENDIANNESS_KEY)
        SEQAN_THROW(GdfIOException("Unsupported endianness key!"));
    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);
    config.isLittleEndian = buffer == GDF_IO_FILE_ENDIANNESS_LITTLE;

    //Read snp compression status.
    clear(buffer);
    skipNChars(reader, 2);  // skip prefix.
    readNChars(buffer, reader, length(GDF_IO_FILE_SNP_COMPRESSION_KEY));
    if (buffer != GDF_IO_FILE_SNP_COMPRESSION_KEY)
        SEQAN_THROW(GdfIOException("Unsupported snp compression key!"));

    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);
    if (buffer == GDF_IO_FILE_SNP_COMPRESSION_2BIT)
        config.compressionMode = GdfIO::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION;
    else
    {
        SEQAN_ASSERT_EQ(buffer, GDF_IO_FILE_SNP_COMPRESSION_GENERIC);
        config.compressionMode = GdfIO::COMPRESSION_MODE_NO_SNP_COMPRESSION;
    }

    //Read coverage compression status.
    clear(buffer);
    skipNChars(reader, 2);  // skip prefix.
    readNChars(buffer, reader, length(GDF_IO_FILE_COVERAGE_COMPRESSION));
    if (buffer != GDF_IO_FILE_COVERAGE_COMPRESSION)
        SEQAN_THROW(GdfIOException("Unsupported coverage compression key!"));

    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);
    unsigned val;
    lexicalCast2(val, buffer);  // Store the coverage compression value.
    switch(val)
    {
        case 1: config.coverageCompression = GdfIO::COVERAGE_COMPRESSION_1_BYTE_PER_VALUE; break;
        case 2: config.coverageCompression = GdfIO::COVERAGE_COMPRESSION_2_BYTE_PER_VALUE; break;
        case 4: config.coverageCompression = GdfIO::COVERAGE_COMPRESSION_4_BYTE_PER_VALUE; break;
        case 8: config.coverageCompression = GdfIO::COVERAGE_COMPRESSION_8_BYTE_PER_VALUE; break;
        default: SEQAN_THROW(GdfIOException("Unkown Coverage Compression Value")); break;
    }
}

template <typename TConfig, typename TStream>
inline void
_readGdfHeaderRefInfo(GdfHeader & gdfHeader,
                      GdfFileConfiguration<TConfig> & config,
                      RecordReader<TStream, SinglePass<> >  & reader)
{
    CharString buffer;
    skipNChars(reader, 2);  // skip prefix.
    // Read refernece ID.
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_REFERENCE_ID_KEY));
    if (!startsWith(buffer, GDF_IO_REFERENCE_ID_KEY))
        SEQAN_THROW(GdfIOException("Unsupported reference id key!"));

    // Skip = sign
    skipNChars(reader, 1);
    clear(gdfHeader.referenceId);
    readLine(gdfHeader.referenceId, reader);  // Read value.

    skipNChars(reader, 2);  // skip prefix.
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_REFERENCE_FILE_KEY));
    if (!startsWith(buffer, GDF_IO_REFERENCE_FILE_KEY))
        SEQAN_THROW(GdfIOException("Unsupported reference filename key!"));

    // Skip = sign
    skipNChars(reader, 1);
    clear(gdfHeader.referenceFilename);
    readLine(gdfHeader.referenceFilename, reader);  // Read value.

    skipNChars(reader, 2);  // skip prefix.
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_REFERENCE_HASH_KEY));
    if (!startsWith(buffer, GDF_IO_REFERENCE_HASH_KEY))
        SEQAN_THROW(GdfIOException("Unsupported reference hash key!"));

    // Skip = sign
    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);  // Read value.
    lexicalCast2(config.refHash, buffer);
}

// ----------------------------------------------------------------------------
// Function _readSeqNames()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
_readSeqNames(GdfHeader & gdfHeader, RecordReader<TStream, SinglePass<> >  & reader)
{
    CharString buffer;
    // Read the sequence names.
    SEQAN_ASSERT_EQ(value(reader), '!');
    skipNChars(reader, 1);
    SEQAN_ASSERT_EQ(value(reader), '!');

    while (true)
    {
        goNext(reader);
        if (value(reader) == GDF_IO_SEQ_NAMES_SEPARATOR[0])
        {
            appendValue(gdfHeader.nameStore, buffer);
            clear(buffer);
            skipNChars(reader, 1);
        }
        if (value(reader) == '\n')
        {
            skipNChars(reader, 1);
            break;
        }
        append(buffer, value(reader));
    }
}

template <typename TConfig, typename TStream>
inline void
readHeader(GdfHeader & gdfHeader,
           GdfFileConfiguration<TConfig> & config,
           RecordReader<TStream, SinglePass<> > & reader,
           Gdf const & /*tag*/)
{
    _readGdfFileConfigInfo(config, reader);
    _readGdfHeaderRefInfo(gdfHeader, config, reader);

    _readSeqNames(gdfHeader, reader);
}

// ----------------------------------------------------------------------------
// Function _addDeltaWithoutCoverage()
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TValue, typename TSpec, typename TDeltaPos, typename TDeltaValue, typename TTag>
inline void _addDeltaWithoutCoverage(DeltaMap<TRefPos, TValue, TSpec> & deltaMap,
                                     TDeltaPos deltaPos,
                                     TDeltaValue const & val,
                                     TTag const & /*deltaType*/)
{
    typedef DeltaMap<TRefPos, TValue, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename Member<TDeltaMap, DeltaMapStoreMember>::Type TDeltaStore;
    typedef typename Size<TDeltaStore>::Type TStoreSize;
    typedef typename DeltaRecord<TEntry>::Type TDeltaRecord;

    TStoreSize storePos = addDeltaValue(deltaMap._deltaStore, val, TTag());
    TEntry tmp;
    tmp.deltaPosition = deltaPos;
    tmp.deltaRecord = TDeltaRecord(selectDeltaType(TTag()), storePos);
    resize(tmp.deltaCoverage, getCoverageSize(deltaMap), false, Exact());
    appendValue(deltaMap._entries, tmp);
}

// ----------------------------------------------------------------------------
// Function _readDeltaCoverage()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TStream>
inline int
_readDeltaCoverage(String<bool, Packed<TSpec> > & bitVector,
                   RecordReader<TStream, SinglePass<> > & reader)
{
    typedef String<bool, Packed<TSpec> > TPackedString;
    typedef typename Host<TPackedString>::Type TPackedStringHost;
    typedef typename Iterator<TPackedStringHost>::Type THostIterator;

    typedef typename PackedHostValue_<TPackedString>::Type TPackedHostValue;
    typedef typename TPackedHostValue::TBitVector TBitVector;

    SEQAN_ASSERT(!empty(host(bitVector)));

    THostIterator it = begin(host(bitVector));
    THostIterator itEnd = end(host(bitVector));
    CharString buffer;
    for (; it != itEnd; ++it)
    {
        clear(buffer);
        readNChars(buffer, reader, sizeof(TBitVector));
        it->i = *reinterpret_cast<TBitVector*>(&buffer[0]);  // TODO(rmaerker): Try move function here?
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function _decodeBitVector()
// ----------------------------------------------------------------------------

template <typename THostSpec, typename TSource>
inline void
_decodeBitVector(String<bool, Packed<THostSpec> > & target,
                 TSource const & source)
{
    typedef String<bool, Packed<THostSpec> >  TPackedString;
    typedef typename Position<TPackedString>::Type TPosition;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Iterator<TPackedHost, Standard>::Type TConstPackedHostIterator;
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename TTraits::THostValue THostValue;
    typedef typename THostValue::TBitVector TBitVector;

    typedef typename Iterator<TSource const, Standard>::Type TSourceIterator;

    SEQAN_ASSERT_NOT(empty(host(target)));

    TSourceIterator srcIt = begin(source, Standard());
    TSourceIterator srcItEnd = end(source, Standard());

    TConstPackedHostIterator it = begin(host(target), Standard()) + 1;

    bool beginWithBit = false;
    TPosition blockSize = *srcIt;  // First value corresponds to leading zero count.

    while (srcIt != srcItEnd)
    {
        // Skip all completely filled blocks of either '0's or '1's.
        TPosition skipFullWords = blockSize / BitsPerValue<TBitVector>::VALUE;
        blockSize -= BitsPerValue<TBitVector>::VALUE * skipFullWords;
        if (beginWithBit)
            for (unsigned i = 0; i < skipFullWords; ++i, ++it)
                it->i = ~static_cast<TBitVector>(0);
        else
            it += skipFullWords;

        TPosition remainingBits = BitsPerValue<TBitVector>::VALUE;
        while (true)  // Process the current word.
        {
            if (remainingBits <= blockSize)  // Stop processing current word.
            {
                blockSize -= remainingBits;
                beginWithBit = (blockSize == 0) ? true : false;  // Handle case when remainingBits == blockSize.
                break;
            }

            it->i |= (~static_cast<TBitVector>(0)) >> (BitsPerValue<TBitVector>::VALUE - remainingBits) + blockSize;  // Skip remaining zeros at the beginning of the current word.
            remainingBits -= blockSize;
            blockSize = *(++srcIt);  // Extract number of ones.

            if (remainingBits <= blockSize)
            {
                blockSize -= remainingBits;
                beginWithBit = (blockSize == 0) ? false : true;  // Handle case when remainingBits == blockSize.
                break;
            }

            remainingBits -= blockSize;
            it->i &= (~static_cast<TBitVector>(0)) << remainingBits;
            if (++srcIt == srcItEnd)  // Check special end condition.
                break;
            blockSize = *srcIt;  // Extract number of zeros.
        }
        ++it;  // Move to next word in bit vector.
        if (!remainingBits)
            if(++srcIt != srcItEnd)  // Check special end condition.
                blockSize = *srcIt;
    }
}

// ----------------------------------------------------------------------------
// Function _readCoverageBlock()
// ----------------------------------------------------------------------------

template <typename TDeltaMapIterator, typename TStream, typename TDecodeType, typename TInputOrder>
inline void
_readCoverageBlock(TDeltaMapIterator & it,
                   TDeltaMapIterator & itEnd,
                   RecordReader<TStream, SinglePass<> > & reader,
                   GdfIOCoverageCompression<TDecodeType> /*tag*/,
                   TInputOrder /*tag*/)
{
    // Read the byte size of the coverage block.
    CharString tmp;
    readNChars(tmp, reader, sizeof(__uint32));  // Read the byte size of the current block.
    __uint32 blockSize = endianSwap(*reinterpret_cast<__uint32*>(&tmp[0]), TInputOrder(), HostByteOrder());
    clear(tmp);
    resize(tmp, blockSize, Exact());
    readNChars(tmp, reader, blockSize);  // Change to read directly from the stream.

    // Process the delta buffer.
    const char * buffer = &tmp[0];
    __uint32 deltaRef = 0;
    buffer += _fromBinary(deltaRef, buffer, TInputOrder(), HostByteOrder());

    while (it != itEnd)  // Process the current block.
    {
        String<TDecodeType> coverageBuffer;

        TDecodeType tmpVal = 0;
        while (tmpVal != MaxValue<TDecodeType>::VALUE)
        {
            buffer += _fromBinary(tmpVal, buffer, TInputOrder(), HostByteOrder());
            appendValue(coverageBuffer, tmpVal);
        }
        _decodeBitVector(deltaCoverage(it), coverageBuffer);
        ++it;
    }
    SEQAN_ASSERT_EQ(buffer - &tmp[0], blockSize);
}

// ----------------------------------------------------------------------------
// Function _readSnp()
// ----------------------------------------------------------------------------

// Writes SNP in separat value after delta pos.
template <typename TAlphabet, typename TBuffer>
inline void
_readSnp(TAlphabet & snp,
         BitCompressedDeltaPos_<GdfIOGenericSnpCompression> & /*deltaPos*/,
         TBuffer & buffer)
{
    buffer += _copyToBuffer(&snp, buffer, sizeof(TAlphabet));
    // Read the offset for the encoded delta.
//    setBitTo(buffer[0], BitsPerValue<__uint8>::VALUE -1, false);  // Reset the MSB to extract delta.
//    TPosition tmpPos = *reinterpret_cast<__uint32*>(&buffer[0]);
//    endianSwap(tmpPos, BigEndian(), HostByteOrder());
//    deltaPos += tmpPos;
//    clear(buffer);
//    readNChars(buffer, reader, sizeof(TAlphabet));
    // Read the value for the snp.
//    snp = *reinterpret_cast<TAlphabet*>(&buffer[0]);
//    return sizeof(__uint32) + sizeof(TAlphabet);
}

// ----------------------------------------------------------------------------
// Function _readSnp()                                                    [Dna]
// ----------------------------------------------------------------------------

// Encodes SNP in delta position assuming, that 28 bits are sufficient to store the delta to the previous variant.
template <typename TAlphabet, typename TBuffer>
inline void
_readSnp(TAlphabet & snp,
         BitCompressedDeltaPos_<GdfIO2BitSnpCompression> & deltaPos,
         TBuffer /*buffer*/)
{
    snp = convert<TAlphabet>(deltaPos.snp);
//    CharString buffer;
//    readNChars(buffer, reader, sizeof(__uint32));
//    snp = static_cast<TAlphabet>((buffer[0] >> 5) & 3);  // Extract the SNP;
//    buffer[0] &= static_cast<__uint8>(~0) >> 3;  // Disable leading three bits.
//    TPosition tmp = *reinterpret_cast<__uint32*>(&buffer[0]);
//    endianSwap(tmp, BigEndian(), HostByteOrder());
//    deltaPos += tmp;
//    return sizeof(__uint32);
}

template <typename TValue, typename TAlphabet, typename TStream, typename TGdfSnpCompression, typename TInputOrder>
inline void _readGdfBlock(DeltaMap<TValue, TAlphabet> & deltaMap,
                          RecordReader<TStream, SinglePass<> > & reader,
                          TGdfSnpCompression /*isDnaCompressed*/,
                          TInputOrder /*tag*/)
{
    typedef DeltaMap<TValue, TAlphabet>                        TDeltaMap;
    typedef typename Size<TDeltaMap>::Type                     TSize;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnp;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDel;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSV>::Type  TIndel;
    typedef typename Iterator<TDeltaMap, Standard>::Type       TMapIter;
    typedef BitCompressedDeltaPos_<TGdfSnpCompression> TDeltaPos;

    // Read the byte size of the delta block.
    CharString tmp;
    readNChars(tmp, reader, sizeof(__uint32));  // Read the byte size of the current block.
    __uint32 blockSize = endianSwap(*reinterpret_cast<__uint32*>(&tmp[0]), TInputOrder(), HostByteOrder());
    clear(tmp);
    resize(tmp, blockSize, Exact());
    readNChars(tmp, reader, blockSize);  // Change to read directly from the stream.

    // Process the delta buffer.
    const char * buffer = &tmp[0];
    __uint32 deltaRef = 0;
    buffer += _fromBinary(deltaRef, buffer, TInputOrder(), HostByteOrder());
//    __uint32 blockRef = endianSwap(*reinterpret_cast<__uint32*>(&buffer[0]), TInputOrder(), HostByteOrder());
//    clear(buffer);
//    readNChars(buffer, reader, sizeof(__uint32));
//    __uint32 blockSize = *reinterpret_cast<__uint32*>(&buffer[0]);
//    endianSwap(blockSize, TInputOrder(), HostByteOrder());

    // TODO(rmaerker): Rewrite to first load the block.
    while (buffer - &tmp[0] < blockSize)
    {
        TDeltaPos deltaPos;
        __uint32 tmp = 0;
        buffer += _fromBinary(tmp, buffer, BigEndian(), HostByteOrder());
        deltaPos.fromWord(tmp);
        deltaRef += deltaPos.pos;

        if (deltaPos.isSnp)  // Dealing with a SNP.
        {
            TSnp snp;
            _readSnp(snp, deltaPos, buffer);  // Read alphabet dependent snp value.
            _addDeltaWithoutCoverage(deltaMap, deltaRef, snp, DeltaTypeSnp());
        }
        else  // Is an InDel.
        {
            buffer += _fromBinary(tmp, buffer, BigEndian(), HostByteOrder());
            BitCompressedInDel_ indel;
            indel.fromWord(tmp);

            if(indel.isDel)  // Read deletion info.
            {
//                __uint32 delSize;
//                buffer += _fromBinary(delSize, buffer, BigEndian(), HostByteOrder());
//                readNChars(buffer, reader, sizeof(__uint32));
                  //                = *reinterpret_cast<__uint32*>(&buffer[0]);
//                endianSwap(delSize, BigEndian(), HostByteOrder());
//                setBitTo(delSize, BitsPerValue<__uint32>::VALUE - 1, false);
                _addDeltaWithoutCoverage(deltaMap, deltaRef, static_cast<TDel>(indel.value), DeltaTypeDel());
//                blockSize -= length(buffer);
            }
            else
            {
                __uint32 insSize = indel.value;
                if (indel.isSV)
                    buffer += _fromBinary(insSize, buffer, BigEndian(), HostByteOrder());

                TIns insBuffer;
                resize(insBuffer, insSize, Exact());
                buffer += _copyToBuffer(&insBuffer[0], buffer, insSize);

                if (indel.isSV)  // RecordSV
                    _addDeltaWithoutCoverage(deltaMap, deltaRef, TIndel(indel.value, insBuffer), DeltaTypeSV());
                else  // Record insertion.
                    _addDeltaWithoutCoverage(deltaMap, deltaRef, insBuffer, DeltaTypeIns());

//                bool isIndel = false;
//                readNChars(buffer, reader, sizeof(__uint32));
//                __uint32 insSize = *reinterpret_cast<__uint32*>(&buffer[0]);
//                endianSwap(insSize, BigEndian(), HostByteOrder());
//                blockSize -= length(buffer);
//                if (isBitSet(insSize, BitsPerValue<__uint32>::VALUE - 2))
//                {  // Handle indel.
//                    setBitTo(insSize, BitsPerValue<__uint32>::VALUE - 2, false);
//                    isIndel = true;
//                }
//                clear(buffer);
//                readNChars(buffer, reader, insSize);
//                TIns insSegment = buffer;
//                blockSize -= length(buffer);
//                if (!isIndel)
//                {
//                    _addDeltaWithoutCoverage(deltaMap, deltaRef, insSegment, DeltaTypeIns());
//                }
//                else
//                {
//                    clear(buffer);
//                    readNChars(buffer, reader, sizeof(__uint32));
//                    __uint32 delSize = *reinterpret_cast<__uint32*>(&buffer[0]);
//                    endianSwap(delSize, BigEndian(), HostByteOrder());
//                    _addDeltaWithoutCoverage(deltaMap, deltaRef, TIndel(delSize, insSegment), DeltaTypeSV());
//                    blockSize -= length(buffer);
//                }
            }
        }
    }
    buffer -= blockSize;  // Set buffer pointer to begin.
    delete[] buffer;
}

// ----------------------------------------------------------------------------
// Function _readGdfData()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TStream, typename TSnpCompressionTag,
          typename TCoverageEncodingType, typename TInputEndian>
inline void
_readGdfData(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
             RecordReader<TStream, SinglePass<> > & reader,
             TSnpCompressionTag /*tag*/,
             GdfIOCoverageCompression<TCoverageEncodingType> /*tag*/,
             TInputEndian /*tag*/)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TIterator;
    // We need to store the number of blocks in the
    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));

    __int32 numOfBlocks = *reinterpret_cast<__uint32*>(&buffer[0]);
    endianSwap(numOfBlocks, TInputEndian(), HostByteOrder());

    unsigned numDeltas = 0;
    while (numOfBlocks != 0)
    {
        _readGdfBlock(deltaMap, reader, TSnpCompressionTag(), TInputEndian());
        TIterator itBlockBegin = begin(deltaMap, Standard()) + numDeltas;
        TIterator itBlockEnd = end(deltaMap, Standard());
        _readCoverageBlock(itBlockBegin, itBlockEnd, reader, GdfIOCoverageCompression<TCoverageEncodingType>(), TInputEndian());
        numDeltas = length(deltaMap);
        --numOfBlocks;
    }
}

// ----------------------------------------------------------------------------
// Function _readGdfData()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TStream, typename TConfig,
          typename TSnpCompressionTag, typename TCoverageEncodingType>
inline void
_readGdfData(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
             RecordReader<TStream, SinglePass<> > & reader,
             GdfFileConfiguration<TConfig> const & config,
             TSnpCompressionTag /*tag*/,
             GdfIOCoverageCompression<TCoverageEncodingType> /*tag*/)
{
    if (config.isLittleEndian)
        _readGdfData(deltaMap, reader, TSnpCompressionTag(), GdfIOCoverageCompression<TCoverageEncodingType>(), LittleEndian());
    else
        _readGdfData(deltaMap, reader, TSnpCompressionTag(), GdfIOCoverageCompression<TCoverageEncodingType>(), BigEndian());
}

// ----------------------------------------------------------------------------
// Function _readGdfData()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TStream, typename TConfig,
          typename TSnpCompressionTag>
inline void
_readGdfData(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
             RecordReader<TStream, SinglePass<> > & reader,
             GdfFileConfiguration<TConfig> const & config,
             TSnpCompressionTag /*tag*/)
{
    switch(config.coverageCompression)
    {
        case GdfIO::COVERAGE_COMPRESSION_1_BYTE_PER_VALUE:
            _readGdfData(deltaMap, reader, config, TSnpCompressionTag(), GdfIOCoverageCompression<__uint8>()); break;
        case GdfIO::COVERAGE_COMPRESSION_2_BYTE_PER_VALUE:
            _readGdfData(deltaMap, reader, config, TSnpCompressionTag(), GdfIOCoverageCompression<__uint16>()); break;
        case GdfIO::COVERAGE_COMPRESSION_4_BYTE_PER_VALUE:
            _readGdfData(deltaMap, reader, config, TSnpCompressionTag(), GdfIOCoverageCompression<__uint32>()); break;
        default:
            _readGdfData(deltaMap, reader, config, TSnpCompressionTag(), GdfIOCoverageCompression<__uint64>()); break;
    }
}

// ----------------------------------------------------------------------------
// Function read()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TConfig, typename TStream>
inline void
read(DeltaMap<TValue, TAlphabet> & deltaMap,
     GdfHeader & gdfHeader,
     GdfFileConfiguration<TConfig> & config,
     TStream  & stream,
     Gdf const & /*tag*/)
{
    SEQAN_TRY
    {
        RecordReader<TStream, SinglePass<> > reader(stream);
        readHeader(gdfHeader, config, reader, Gdf());
        setCoverageSize(deltaMap, length(gdfHeader.nameStore));

        if (config.compressionMode == GdfIO::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION)
            _readGdfData(deltaMap, reader, config, GdfIO2BitSnpCompression());
        else
            _readGdfData(deltaMap, reader, config, GdfIOGenericSnpCompression());
    }
    SEQAN_CATCH(Exception e)
    {
        SEQAN_THROW(GdfIOException(e.what()));
    }
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_READ_H_
