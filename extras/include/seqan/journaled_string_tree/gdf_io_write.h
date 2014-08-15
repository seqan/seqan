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

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#define GDF_IO_WRITE_KEY_VALUE(stream,key,...) stream << GDF_IO_HEADER_PREFIX << key << GDF_IO_KEY_VALUE_SEPARATOR << (__VA_ARGS__) << '\n'

template <typename TStream, typename TValue>
inline void
_writeGdfHeaderFileInfo(TStream & stream, GdfHeader const & jseqHeader, GdfFileConfiguration<TValue> const & config)
{
//    streamWriteBlock(stream, toCString(GDF_IO_HEADER_PREFIX), length(GDF_IO_HEADER_PREFIX));
//    streamWriteBlock(stream, toCString(GDF_IO_FILE_VERSION_KEY), length(GDF_IO_FILE_VERSION_KEY));
//    streamWriteBlock(stream, toCString(GDF_IO_KEY_VALUE_SEPARATOR), length(GDF_IO_KEY_VALUE_SEPARATOR));
//    streamWriteBlock(stream, toCString(GDF_IO_FILE_VERSION_VALUE_PREFIX), length(GDF_IO_FILE_VERSION_VALUE_PREFIX));
//    std::stringstream versionNumber;
//    versionNumber << GDF_IO_FILE_VERSION_BIG << GDF_IO_FILE_VERSION_VALUE_SEPARATOR << GDF_IO_FILE_VERSION_LITTLE;
//    streamWriteBlock(stream, versionNumber.str().c_str(), versionNumber.str().size());
//    streamWriteChar(stream, '\n');

    // Write the file version information.
    GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_VERSION_KEY,
                           GDF_IO_FILE_VERSION_VALUE_PREFIX << GDF_IO_FILE_VERSION_BIG << GDF_IO_FILE_VERSION_VALUE_SEPARATOR << GDF_IO_FILE_VERSION_LITTLE);

    // Write system's endianness.
    if (IsLittleEndian::VALUE)
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_ENDIANNESS_KEY, GDF_IO_FILE_ENDIANNESS_LITTLE);
    else
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_ENDIANNESS_KEY, GDF_IO_FILE_ENDIANNESS_BIG);
//    stream << GDF_IO_HEADER_PREFIX << GDF_IO_FILE_VERSION_KEY << GDF_IO_KEY_VALUE_SEPARATOR <<
//              GDF_IO_FILE_VERSION_VALUE_PREFIX << GDF_IO_FILE_VERSION_BIG << GDF_IO_FILE_VERSION_VALUE_SEPARATOR <<
//              GDF_IO_FILE_VERSION_LITTLE << '\n';
    // Write the endianess.  -> We always write in Network Byte Order (BgEndian)
//    streamWriteBlock(stream, toCString(GDF_IO_HEADER_PREFIX), length(GDF_IO_HEADER_PREFIX));
//    streamWriteBlock(stream, toCString(GDF_IO_FILE_ENDIANNESS_KEY), length(GDF_IO_FILE_ENDIANNESS_KEY));
//    streamWriteBlock(stream, toCString(GDF_IO_KEY_VALUE_SEPARATOR), length(GDF_IO_KEY_VALUE_SEPARATOR));
//    if (SystemByteOrder::IS_LITTLE_ENDIAN())
//        streamWriteBlock(stream, toCString(GDF_IO_FILE_ENDIANNESS_LITTLE), length(GDF_IO_FILE_ENDIANNESS_LITTLE));
//    else
//        streamWriteBlock(stream, toCString(GDF_IO_FILE_ENDIANNESS_BIG), length(GDF_IO_FILE_ENDIANNESS_BIG));
//    streamWriteChar(stream, '\n');

    // Write the block information.
//    streamWriteBlock(stream, toCString(GDF_IO_HEADER_PREFIX), length(GDF_IO_HEADER_PREFIX));
//    streamWriteBlock(stream, toCString(GDF_IO_FILE_BLOCKSIZE_KEY), length(GDF_IO_FILE_BLOCKSIZE_KEY));
//    streamWriteBlock(stream, toCString(GDF_IO_KEY_VALUE_SEPARATOR), length(GDF_IO_KEY_VALUE_SEPARATOR));
//    std::stringstream blockSize;
//    blockSize << jseqHeader._fileInfos._blockSize;
//    streamWriteBlock(stream, blockSize.str().c_str(), blockSize.str().size());
//    streamWriteChar(stream, '\n');

    // Write Compression mode.

//    streamWriteBlock(stream, toCString(GDF_IO_HEADER_PREFIX), length(GDF_IO_HEADER_PREFIX));
//    streamWriteBlock(stream, toCString(GDF_IO_FILE_SNP_COMPRESSION_KEY), length(GDF_IO_FILE_SNP_COMPRESSION_KEY));
//    streamWriteBlock(stream, toCString(GDF_IO_KEY_VALUE_SEPARATOR), length(GDF_IO_KEY_VALUE_SEPARATOR));
    if (config.compressionMode == GdfIO::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION)
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_SNP_COMPRESSION_KEY, GDF_IO_FILE_SNP_COMPRESSION_2BIT);
    else
        GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_FILE_SNP_COMPRESSION_KEY, GDF_IO_FILE_SNP_COMPRESSION_GERNERIC);
}

template <typename TStream, typename TConfig>
inline void
_writeGdfHeaderReferenceInfo(TStream & stream, GdfHeader const & jseqHeader, TConfig const & config)
{
    // Write the reference id.
      GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_REFERENCE_ID_KEY, jseqHeader.referenceId);
//    streamWriteBlock(stream, toCString(GDF_IO_HEADER_PREFIX), length(GDF_IO_HEADER_PREFIX));
//    streamWriteBlock(stream, toCString(GDF_IO_REFERENCE_ID_KEY), length(GDF_IO_REFERENCE_ID_KEY));
//    streamWriteBlock(stream, toCString(GDF_IO_KEY_VALUE_SEPARATOR), length(GDF_IO_KEY_VALUE_SEPARATOR));
//    streamWriteBlock(stream, toCString(jseqHeader._refInfos._refId), length(jseqHeader._refInfos._refId));
//    streamWriteChar(stream, '\n');
    // Write the reference file.
      GDF_IO_WRITE_KEY_VALUE(stream, GDF_IO_REFERENCE_FILE_KEY, jseqHeader.referenceFilename);
//      streamWriteBlock(stream, toCString(GDF_IO_HEADER_PREFIX), length(GDF_IO_HEADER_PREFIX));
//      streamWriteBlock(stream, toCString(GDF_IO_REFERENCE_FILE_KEY), length(GDF_IO_REFERENCE_FILE_KEY));
//      streamWriteBlock(stream, toCString(GDF_IO_KEY_VALUE_SEPARATOR), length(GDF_IO_KEY_VALUE_SEPARATOR));
//      streamWriteBlock(stream, toCString(jseqHeader._refInfos._refFile), length(jseqHeader._refInfos._refFile));
//      streamWriteChar(stream, '\n');
      // Write the reference hash.
      GDF_IO_WRITE_KEY_VALUE(stram, GDF_IO_REFERENCE_HASH_KEY, config.refHash);
//      streamWriteBlock(stream, toCString(GDF_IO_HEADER_PREFIX), length(GDF_IO_HEADER_PREFIX));
//      streamWriteBlock(stream, toCString(GDF_IO_REFERENCE_HASH_KEY), length(GDF_IO_REFERENCE_HASH_KEY));
//      streamWriteBlock(stream, toCString(GDF_IO_KEY_VALUE_SEPARATOR), length(GDF_IO_KEY_VALUE_SEPARATOR));
//      std::stringstream referenceHash;
//      referenceHash << jseqHeader._refInfos._refHash;
//      streamWriteBlock(stream, referenceHash.str().c_str(), referenceHash.str().size());
//      streamWriteChar(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function
// ----------------------------------------------------------------------------

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

template <typename TStream, typename TConfig>
inline void _writeGdfHeader(TStream & stream, GdfHeader const & header, TConfig const & config)
{
    // Write the file information.
    _writeGdfHeaderFileInfo(stream, header, config);
    // Write reference information.
    _writeGdfHeaderReferenceInfo(stream, header, config);

    // Write additional data.
    stream << GDF_IO_SEQ_NAMES_PREFIX;
    for (unsigned i = 0; i < length(header.nameStore); ++i)
        stream << header.nameStore)[i] << GDF_IO_SEQ_NAMES_SEPARATOR;
//        streamWriteBlock(stream, toCString((header.nameStore)[i]), length((header._nameStorePtr)[i]));
//        streamWriteBlock(stream, toCString(GDF_IO_SEQ_NAMES_SEPARATOR), length(GDF_IO_SEQ_NAMES_SEPARATOR));
    stream << '\n';
}

// ----------------------------------------------------------------------------
// Function _writeSnp()
// ----------------------------------------------------------------------------

// Writes SNP in separat value after delta pos.
template <typename TTarget, typename TPosition, typename TAlphabet>
inline void
_writeSnp(TTarget & blockBuffer, TPosition deltaPos, TAlphabet snp)
{
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -1));  // The last bit should not be set.

    setBit(deltaPos, BitsPerValue<__uint32>::VALUE -1);  // Set bit to indicate SNP
    register __uint32 offset = sizeof(__uint32) + sizeof(TAlphabet);
    resize(blockBuffer, length(blockBuffer) + offset);
    endianSwap(deltaPos, HostByteOrder(), BigEndian());  // swaps endianness if HostByteOrder is little endian.
    char* refBuffer = reinterpret_cast<char*>(&deltaPos);
    arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - offset);
    char* snpBuffer = reinterpret_cast<char*>(&snp);
    arrayMoveForward(snpBuffer, snpBuffer + sizeof(TAlphabet), end(blockBuffer, Standard()) - sizeof(TAlphabet));
}

// ----------------------------------------------------------------------------
// Function _writeSnp()                                                   [Dna]
// ----------------------------------------------------------------------------

// Encodes SNP in delta position assuming, that 28 bits are sufficient to store the delta to the previous variant.
template <typename TTarget, typename TPosition>
inline void
_writeSnp(TTarget & blockBuffer, TPosition deltaPos, Dna snp)
{
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -1));  // Bit at index 31 not set.
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -2));  // Bit at index 30 not set.
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -3));  // Bit at index 29 not set.

    setBit(deltaPos, BitsPerValue<__uint32>::VALUE -1);
    deltaPos |= static_cast<__uint32>(snp) << ((sizeof(__uint32) << 3) - 3);
    resize(blockBuffer, length(blockBuffer) + sizeof(__uint32));
    endianSwap(deltaPos, HostByteOrder(), BigEndian());
    char* refBuffer = reinterpret_cast<char*>(&deltaPos);
    arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - sizeof(__uint32));
}

// ----------------------------------------------------------------------------
// Function _writeDataBlock()
// ----------------------------------------------------------------------------

template <typename TStream, typename TDeltaMapIter, typename TValue, typename TAlphabet>
inline int _writeDataBlock(TStream & stream,
                           TDeltaMapIter & it,
                           TDeltaMapIter & itEnd,
                           DeltaMap<TValue, TAlphabet> const & /*deltaMap*/)
{
    typedef DeltaMap<TValue, TAlphabet> TDeltaMap;

    CharString blockBuffer;
    typename Iterator<CharString>::Type blockBuffIt;
    __uint32 lastRefPos = deltaPosition(it);
    // Write the reference offset of the current block.
    writeBinary(stream, lastRefPos);
    TDeltaMapIter itDelta = it;
    while(itDelta != itEnd)
    {
        // Store the reference position.
        __uint32 deltaPos = deltaPosition(itDelta) - lastRefPos;

        // Write SNP Data.
        if (deltaType(itDelta) == DELTA_TYPE_SNP)
        {
            _writeSnp(blockBuffer, deltaPos, deltaValue(itDelta, DeltaTypeSnp()));
        }
        else  // Write Indel
        {
            resize(blockBuffer, length(blockBuffer) + sizeof(__uint32));
            blockBuffIt = end(blockBuffer) - sizeof(__uint32);
            endianSwap(deltaPos, HostByteOrder(), BigEndian());
            char * refBuffer = reinterpret_cast<char *>(&deltaPos);
            arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);

            if (deltaType(itDelta) == DELTA_TYPE_DEL)  // Write deletion.
            {
                __uint32 del = static_cast<__uint32>(deltaValue(itDelta, DeltaTypeDel()));
                setBit(del, BitsPerValue<__uint32>::VALUE - 1);
                resize(blockBuffer, length(blockBuffer) + sizeof(del));
                blockBuffIt = end(blockBuffer) - sizeof(del);
                endianSwap(del, HostByteOrder(), BigEndian());
                const char * delBuffer = reinterpret_cast<const char *>(&del);
                arrayMoveForward(delBuffer, delBuffer + sizeof(del), blockBuffIt);
            }
            else
            {
                typedef typename DeltaValue<TDeltaMap const, DeltaTypeIns>::Type TIns;

                // Handle Indel.
                if (deltaType(itDelta) == DELTA_TYPE_SV)
                {
                    TIns ins = deltaValue(itDelta, DeltaTypeSV()).i2;
                    __uint32 insLength = length(ins);
                    setBit(insLength, BitsPerValue<__uint32>::VALUE - 2);
                    resize(blockBuffer, length(blockBuffer) + sizeof(insLength) + sizeof(del) + length(ins));
                    blockBuffIt = end(blockBuffer) - sizeof(insLength) - sizeof(del) - length(ins);
                    // Write size of insertion.
                    endianSwap(insLength, HostByteOrder(), BigEndian());
                    const char * insBuffer = reinterpret_cast<const char *>(&insLength);
                    arrayMoveForward(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
                    // Write inserted characters.
                    blockBuffIt = end(blockBuffer) - length(ins) - sizeof(del);
                    arrayMoveForward(begin(ins, Standard()), end(ins, Standard()), blockBuffIt);
                    // Write size of deletion
                    __uint32 del = static_cast<__uint32>(deltaValue(itDelta, DeltaTypeSV()).i1);
                    endianSwap(del, HostByteOrder(), BigEndian());
                    const char * delBuffer = reinterpret_cast<const char *>(&del);
                    blockBuffIt = end(blockBuffer) - sizeof(del);
                    arrayMoveForward(delBuffer, delBuffer + sizeof(del), blockBuffIt);
                }
                else  // Handle Insertion.
                {
                    SEQAN_ASSERT(deltaType(itDelta) == DELTA_TYPE_INS);
                    TIns ins = deltaValue(itDelta, DeltaTypeIns());
                    __uint32 insLength = length(ins);
                    SEQAN_ASSERT_NOT(isBitSet(insLength, BitsPerValue<__uint32>::VALUE - 2));
                    endianSwap(insLength, HostByteOrder(), BigEndian());
                    const char * insBuffer = reinterpret_cast<const char *>(&insLength);
                    resize(blockBuffer, length(blockBuffer) + sizeof(insLength) + length(ins));
                    blockBuffIt = end(blockBuffer) - sizeof(insLength) - length(ins);
                    arrayMoveForward(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
                    blockBuffIt = end(blockBuffer) - length(ins);
                    arrayMoveForward(begin(ins, Standard()), end(ins, Standard()), blockBuffIt);
                }
            }
        }
        lastRefPos = deltaPosition(itDelta);
        ++itDelta;
    }

    // Write the block Data.
    unsigned blockLength = length(blockBuffer);
    writeBinary(stream, blockLength);
    streamWriteBlock(stream, &blockBuffer[0], blockLength);

    // Write the coverage of the block
    for (; it != itEnd; ++it)
        _writeBitVector(stream, deltaCoverage(it));

    return 0;
}

// Write the io context.
template <typename TStream, typename TDeltaStore, typename TDeltaCoverageStore, typename TValue>
inline int _writeGdfData(TStream & stream,
                          DeltaMap<TDeltaStore, TDeltaCoverageStore> const & deltaMap,
                          GdfFileConfiguration<TValue> const & /*config*/)
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore> TDeltaMap;
    typedef typename Size<TDeltaMap>::Type TSize;
    typedef typename Iterator<TDeltaMap const, Standard>::Type TIterator;
    typedef GdfFileConfiguration<TValue> TConfig;

    TSize maxNumOfNodes = length(deltaMap);
    TSize numOfBlocks = (maxNumOfNodes + TConfig::BLOCK_SIZE - 1) / TConfig::BLOCK_SIZE;

    // Write the block containing the number of blocks to read.
    writeBinary(stream, numOfBlocks);
    for (TSize i = 0; i < numOfBlocks; ++i)  // For each block:
    {
        TIterator it = begin(deltaMap, Standard()) + (blockSize * i);
        TIterator itEnd = _min(end(deltaMap, Standard()), it + blockSize);
        _writeDataBlock(stream, it, itEnd, deltaMap);
    }
    return 0;
}

template <typename TStream, typename TValue, typename TAlphabet, typename TConfig>
inline void
write(TStream & stream,
      DeltaMap<TValue, TAlphabet> const & deltaMap,
      GdfHeader const & gdfHeader,
      TConfig const & config,
      Gdf const & /*tag*/)
{
    _writeGdfHeader(stream, gdfHeader, config);
    _writeGdfData(stream, deltaMap, config);
}

///*!
// * @fn journalSet#write
// * @deprecated
// */
//template <typename TStream, typename TJournalSequence>
//inline int
//write(TStream & stream,
//      StringSet<TJournalSequence, Owner<JournaledSet> > const & journalSet,
//      GdfHeader const & jseqHeader,
//      Gdf const & /*tag*/)
//{
//    typedef typename Value<TJournalSequence>::Type TAlphabet;
//    typedef typename Position<TJournalSequence>::Type TPosition;
//
//    if (empty(host(journalSet)))
//        return -1;
//    DeltaMap<TPosition, TAlphabet> deltaMap;
//    adaptTo(deltaMap, journalSet);
//    write(stream, deltaMap, jseqHeader, Gdf());
//    return 0;
//}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_JOURNALED_STRING_TREE_DELTA_MAP_IO_WRITE_H_
