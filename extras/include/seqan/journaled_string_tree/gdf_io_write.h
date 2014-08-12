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

template <typename TStream, typename TConfig>
inline int
_writeGdfHeaderFileInfo(TStream & stream, GdfHeader<TConfig> const & jseqHeader)
{
    // Write the file version information.
    streamWriteBlock(stream, toCString(GdfIO::HEADER_PREFIX), length(GdfIO::HEADER_PREFIX));
    streamWriteBlock(stream, toCString(GdfIO::FILE_VERSION_KEY), length(GdfIO::FILE_VERSION_KEY));
    streamWriteBlock(stream, toCString(GdfIO::KEY_VALUE_SEPARATOR), length(GdfIO::KEY_VALUE_SEPARATOR));
    streamWriteBlock(stream, toCString(GdfIO::FILE_VERSION_VALUE_PREFIX), length(GdfIO::FILE_VERSION_VALUE_PREFIX));
    std::stringstream versionNumber;
    versionNumber << GdfIO::FILE_VERSION_BIG << GdfIO::FILE_VERSION_VALUE_SEPARATOR << GdfIO::FILE_VERSION_LITTLE;
    streamWriteBlock(stream, versionNumber.str().c_str(), versionNumber.str().size());
    streamWriteChar(stream, '\n');
    // Write the endianess.
    streamWriteBlock(stream, toCString(GdfIO::HEADER_PREFIX), length(GdfIO::HEADER_PREFIX));
    streamWriteBlock(stream, toCString(GdfIO::FILE_ENDIANNESS_KEY), length(GdfIO::FILE_ENDIANNESS_KEY));
    streamWriteBlock(stream, toCString(GdfIO::KEY_VALUE_SEPARATOR), length(GdfIO::KEY_VALUE_SEPARATOR));
    if (SystemByteOrder::IS_LITTLE_ENDIAN())
        streamWriteBlock(stream, toCString(GdfIO::FILE_ENDIANNESS_LITTLE), length(GdfIO::FILE_ENDIANNESS_LITTLE));
    else
        streamWriteBlock(stream, toCString(GdfIO::FILE_ENDIANNESS_BIG), length(GdfIO::FILE_ENDIANNESS_BIG));
    streamWriteChar(stream, '\n');
    // Write the block information.
    streamWriteBlock(stream, toCString(GdfIO::HEADER_PREFIX), length(GdfIO::HEADER_PREFIX));
    streamWriteBlock(stream, toCString(GdfIO::FILE_BLOCKSIZE_KEY), length(GdfIO::FILE_BLOCKSIZE_KEY));
    streamWriteBlock(stream, toCString(GdfIO::KEY_VALUE_SEPARATOR), length(GdfIO::KEY_VALUE_SEPARATOR));
    std::stringstream blockSize;
    blockSize << jseqHeader._fileInfos._blockSize;
    streamWriteBlock(stream, blockSize.str().c_str(), blockSize.str().size());
    streamWriteChar(stream, '\n');
    streamWriteBlock(stream, toCString(GdfIO::HEADER_PREFIX), length(GdfIO::HEADER_PREFIX));
    streamWriteBlock(stream, toCString(GdfIO::FILE_SNP_COMPRESSION_KEY), length(GdfIO::FILE_SNP_COMPRESSION_KEY));
    streamWriteBlock(stream, toCString(GdfIO::KEY_VALUE_SEPARATOR), length(GdfIO::KEY_VALUE_SEPARATOR));
    if (jseqHeader._fileInfos._snpCompression)
        streamWriteBlock(stream, toCString(GdfIO::FILE_SNP_COMPRESSION_2BIT), length(GdfIO::FILE_SNP_COMPRESSION_2BIT));
    else
        streamWriteBlock(stream, toCString(GdfIO::FILE_SNP_COMPRESSION_GENERIC), length(GdfIO::FILE_SNP_COMPRESSION_GENERIC));
    streamWriteChar(stream, '\n');

    return 0;
}

template <typename TStream, typename TConfig>
inline int
_writeGdfHeaderReferenceInfo(TStream & stream, GdfHeader<TConfig> const & jseqHeader)
{
    // Write the reference id.
    streamWriteBlock(stream, toCString(GdfIO::HEADER_PREFIX), length(GdfIO::HEADER_PREFIX));
    streamWriteBlock(stream, toCString(GdfIO::REFERENCE_ID_KEY), length(GdfIO::REFERENCE_ID_KEY));
    streamWriteBlock(stream, toCString(GdfIO::KEY_VALUE_SEPARATOR), length(GdfIO::KEY_VALUE_SEPARATOR));
    streamWriteBlock(stream, toCString(jseqHeader._refInfos._refId), length(jseqHeader._refInfos._refId));
    streamWriteChar(stream, '\n');
    // Write the reference file.
    streamWriteBlock(stream, toCString(GdfIO::HEADER_PREFIX), length(GdfIO::HEADER_PREFIX));
    streamWriteBlock(stream, toCString(GdfIO::REFERENCE_FILE_KEY), length(GdfIO::REFERENCE_FILE_KEY));
    streamWriteBlock(stream, toCString(GdfIO::KEY_VALUE_SEPARATOR), length(GdfIO::KEY_VALUE_SEPARATOR));
    streamWriteBlock(stream, toCString(jseqHeader._refInfos._refFile), length(jseqHeader._refInfos._refFile));
    streamWriteChar(stream, '\n');
    // Write the reference hash.
    streamWriteBlock(stream, toCString(GdfIO::HEADER_PREFIX), length(GdfIO::HEADER_PREFIX));
    streamWriteBlock(stream, toCString(GdfIO::REFERENCE_HASH_KEY), length(GdfIO::REFERENCE_HASH_KEY));
    streamWriteBlock(stream, toCString(GdfIO::KEY_VALUE_SEPARATOR), length(GdfIO::KEY_VALUE_SEPARATOR));
    std::stringstream referenceHash;
    referenceHash << jseqHeader._refInfos._refHash;
    streamWriteBlock(stream, referenceHash.str().c_str(), referenceHash.str().size());
    streamWriteChar(stream, '\n');

    return 0;
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
inline int _writeGdfHeader(TStream & stream, GdfHeader<TConfig> const & header)
{
    // Write the file information.
    _writeGdfHeaderFileInfo(stream, header);
    // Write reference information.
    _writeGdfHeaderReferenceInfo(stream, header);

    // Write additional data.
    streamWriteBlock(stream, toCString(GdfIO::SEQ_NAMES_PREFIX), length(GdfIO::SEQ_NAMES_PREFIX));
    for (unsigned i = 0; i < length(*header._nameStorePtr); ++i)
    {
        streamWriteBlock(stream, toCString((*header._nameStorePtr)[i]), length((*header._nameStorePtr)[i]));
        streamWriteBlock(stream, toCString(GdfIO::SEQ_NAMES_SEPARATOR), length(GdfIO::SEQ_NAMES_SEPARATOR));
    }
    streamWriteChar(stream, '\n');
    return 0;
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
    char* refBuffer = reinterpret_cast<char*>(&deltaPos);
    resize(blockBuffer, length(blockBuffer) + offset);
    if (SystemByteOrder::IS_LITTLE_ENDIAN())
        arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - offset);
    else
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
    char* refBuffer = reinterpret_cast<char*>(&deltaPos);
    if (SystemByteOrder::IS_LITTLE_ENDIAN())
        arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), end(blockBuffer, Standard()) - sizeof(__uint32));
    else
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
    __uint32 lastRefPos = *it;
    // Write the reference offset of the current block.
    streamWriteBlock(stream, reinterpret_cast<char*>(&lastRefPos), sizeof(lastRefPos));

    TDeltaMapIter itDelta = it;

    while(itDelta != itEnd)
    {
        // Store the reference position.
        __uint32 deltaPos = *itDelta - lastRefPos;

        // Write SNP Data.
//        TMappedDelta deltaInfo = mappedDelta(deltaMap, itDelta - itBegin);
        if (deltaType(itDelta) == DeltaType::DELTA_TYPE_SNP)
        {
//            resize(blockBuffer, length(blockBuffer) + sizeof(__uint32));
//            blockBuffIt = end(blockBuffer) - sizeof(__uint32);
            _writeSnp(blockBuffer, deltaPos, deltaSnp(itDelta));
//            setBit(deltaPos, BitsPerValue<__uint32>::VALUE -1);
//            deltaPos |= static_cast<__uint32>(deltaSnp(deltaMap, deltaPosition(deltaInfo))) << (sizeof(__uint32) * 8 - 3);
//            char* refBuffer = reinterpret_cast<char*>(&deltaPos);
//
//            if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);
//            else
//                arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);
        }
        else  // Write Indel
        {
            char * refBuffer = reinterpret_cast<char *>(&deltaPos);
            resize(blockBuffer, length(blockBuffer) + sizeof(__uint32));
            blockBuffIt = end(blockBuffer) - sizeof(__uint32);
            if (SystemByteOrder::IS_LITTLE_ENDIAN())
                arrayMoveForwardReverse(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);
            else
                arrayMoveForward(refBuffer, refBuffer + sizeof(__uint32), blockBuffIt);

            if (deltaType(itDelta) == DeltaType::DELTA_TYPE_DEL)  // Write deletion.
            {
                __uint32 del = static_cast<__uint32>(deltaDel(itDelta));
                setBit(del, BitsPerValue<__uint32>::VALUE - 1);
                const char * delBuffer = reinterpret_cast<const char *>(&del);
                resize(blockBuffer, length(blockBuffer) + sizeof(del));
                blockBuffIt = end(blockBuffer) - sizeof(del);
                if (SystemByteOrder::IS_LITTLE_ENDIAN())
                    arrayMoveForwardReverse(delBuffer, delBuffer + sizeof(del), blockBuffIt);
                else
                    arrayMoveForward(delBuffer, delBuffer + sizeof(del), blockBuffIt);
            }
            else
            {
                typedef typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_INS>::Type TIns;

                // Handle Indel.
                if (deltaType(itDelta) == DeltaType::DELTA_TYPE_INDEL)
                {
                    TIns ins = deltaIndel(itDelta).i2;
                    __uint32 insLength = length(ins);
                    setBit(insLength, BitsPerValue<__uint32>::VALUE - 2);
                    __uint32 del = static_cast<__uint32>(deltaIndel(itDelta).i1);
                    const char * delBuffer = reinterpret_cast<const char *>(&del);
                    const char * insBuffer = reinterpret_cast<const char *>(&insLength);
                    resize(blockBuffer, length(blockBuffer) + sizeof(insLength) + sizeof(del) + length(ins));
                    blockBuffIt = end(blockBuffer) - sizeof(insLength) - sizeof(del) - length(ins);
                    // Write size of insertion.
                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
                        arrayMoveForwardReverse(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
                    else
                        arrayMoveForward(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
                    // Write inserted characters.
                    blockBuffIt = end(blockBuffer) - length(ins) - sizeof(del);
                    arrayMoveForward(begin(ins, Standard()), end(ins, Standard()), blockBuffIt);
                    // Write size of deletion
                    blockBuffIt = end(blockBuffer) - sizeof(del);
                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
                        arrayMoveForwardReverse(delBuffer, delBuffer + sizeof(del), blockBuffIt);
                    else
                        arrayMoveForward(delBuffer, delBuffer + sizeof(del), blockBuffIt);
                }
                else  // Handle Insertion.
                {
                    SEQAN_ASSERT(deltaType(itDelta) == DeltaType::DELTA_TYPE_INS);
                    TIns ins = deltaIns(itDelta);
                    __uint32 insLength = length(ins);
                    SEQAN_ASSERT_NOT(isBitSet(insLength, BitsPerValue<__uint32>::VALUE - 2));

                    const char * insBuffer = reinterpret_cast<const char *>(&insLength);
                    resize(blockBuffer, length(blockBuffer) + sizeof(insLength) + length(ins));
                    blockBuffIt = end(blockBuffer) - sizeof(insLength) - length(ins);
                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
                        arrayMoveForwardReverse(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
                    else
                        arrayMoveForward(insBuffer, insBuffer + sizeof(insLength), blockBuffIt);
                    blockBuffIt = end(blockBuffer) - length(ins);
                    arrayMoveForward(begin(ins, Standard()), end(ins, Standard()), blockBuffIt);
                }
            }
        }
        lastRefPos = *itDelta;
        ++itDelta;
    }

    // Write the block Data.
    unsigned blockLength = length(blockBuffer);
    streamWriteBlock(stream, reinterpret_cast<char*>(&blockLength), sizeof(blockLength));
    streamWriteBlock(stream, &blockBuffer[0], blockLength);

    // Write the coverage of the block
    for (; it != itEnd; ++it)
        _writeBitVector(stream, deltaCoverage(it));

    return 0;
}

// Write the io context.
template <typename TStream, typename TDeltaStore, typename TDeltaCoverageStore, typename TSize>
inline int _writeGdfData(TStream & stream,
                          DeltaMap<TDeltaStore, TDeltaCoverageStore> const & deltaMap,
                          TSize const & blockSize)
{
    typedef DeltaMap<TDeltaStore, TDeltaCoverageStore> TDeltaMap;
    typedef typename Iterator<TDeltaMap const, Standard>::Type TIterator;

    unsigned maxNumOfNodes = length(deltaMap);
    unsigned numOfBlocks = std::ceil(static_cast<double>(maxNumOfNodes)/static_cast<double>(blockSize));

    // Write the block containing the number of blocks to read.
    streamWriteBlock(stream, reinterpret_cast<char*>(&numOfBlocks), sizeof(numOfBlocks));

    for (unsigned i = 0; i < numOfBlocks; ++i)  // For each block:
    {
        TIterator it = begin(deltaMap, Standard()) + (blockSize * i);
        TIterator itEnd = _min(end(deltaMap, Standard()), it + blockSize);
        _writeDataBlock(stream, it, itEnd, deltaMap);
    }
    return 0;
}

template <typename TStream, typename TValue, typename TAlphabet, typename TConfig>
inline int
write(TStream & stream,
      DeltaMap<TValue, TAlphabet> const & deltaMap,
      GdfHeader<TConfig> const & gdfHeader,
      Gdf const & /*tag*/)
{
    int res = _writeGdfHeader(stream, gdfHeader);
    if (res != 0)
        return res;
    return _writeGdfData(stream, deltaMap, gdfHeader._fileInfos._blockSize);
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
