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
_readGdfHeaderFileInfo(GdfHeader & gdfHeader,
                       GdfFileConfiguration<TConfig> & config,
                       RecordReader<TStream, SinglePass<> > & reader)
{
    CharString buffer;

    // Skip the file prefix.
    skipNChars(reader, 2);
    // Read file version.
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_FILE_VERSION_KEY));
//        streamReadBlock(&buffer[0], stream, length(GDF_IO_FILE_VERSION_KEY));
    if (!startsWith(buffer, GDF_IO_FILE_VERSION_KEY))
        SEQAN_THROW(GdfIOException("Unsupported file version key!"));

    // Skip = sign
    skipNChars(reader, 1);
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_FILE_VERSION_VALUE_PREFIX));
    if (!startsWith(buffer, GDF_IO_FILE_VERSION_VALUE_PREFIX))
        SEQAN_THROW(GdfIOException("Unsupported file version prefix!"));
    // Read big version number.
//        readD(buffer[0], stream);
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

    // Read block size.
//    clear(buffer);
//    skipNChars(reader, 2);  // skip prefix.
//    readNChars(buffer, reader, length(GDF_IO_FILE_BLOCKSIZE_KEY));
//    if (buffer != GDF_IO_FILE_BLOCKSIZE_KEY)
//        return -1;      // TODO(rmaerker): Either use error codes or use exceptions.
//    skipNChars(reader, 1);
//    clear(buffer);
//    readLine(buffer, reader);
//    lexicalCast2(gdfHeader._fileInfos._blockSize, buffer);

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
    return 0;
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
    clear(gdfHeader._refInfos._refId);
    readLine(gdfHeader._refInfos._refId, reader);  // Read value.

    skipNChars(reader, 2);  // skip prefix.
    clear(buffer);
    readNChars(buffer, reader, length(GDF_IO_REFERENCE_FILE_KEY));
    if (!startsWith(buffer, GDF_IO_REFERENCE_FILE_KEY))
        SEQAN_THROW(GdfIOException("Unsupported reference filename key!"));

    // Skip = sign
    skipNChars(reader, 1);
    clear(gdfHeader._refInfos._refFile);
    readLine(gdfHeader._refInfos._refFile, reader);  // Read value.

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
inline int
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
    return GdfIO::PARSE_OK;
}

template <typename TConfig, typename TStream>
inline void
readHeader(GdfHeader & gdfHeader,
           GdfFileConfiguration<TConfig> & config,
           RecordReader<TStream, SinglePass<> > & reader,
           Gdf const & /*tag*/)
{
    _readGdfHeaderFileInfo(gdfHeader, config, reader);
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
// Function _readSnp()
// ----------------------------------------------------------------------------

// Writes SNP in separat value after delta pos.
template <typename TAlphabet, typename TPosition, typename TReader>
inline unsigned
_readSnp(TAlphabet & snp,
         TPosition & deltaPos,
         TReader & reader,
         False /*isDnaCompressed*/)
{
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -1));  // The last bit should not be set.

    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));
    // Read the offset for the encoded delta.
    setBitTo(buffer[0], BitsPerValue<__uint8>::VALUE -1, false);  // Reset the MSB to extract delta.
    TPosition tmpPos = *reinterpret_cast<__uint32*>(&buffer[0]);
    endianSwap(tmpPos, BigEndian(), HostByteOrder());
    deltaPos += tmpPos;
    clear(buffer);
    readNChars(buffer, reader, sizeof(TAlphabet));
    // Read the value for the snp.
//    if (gdfHeader._fileInfos._byteOrder != SystemByteOrder::IS_LITTLE_ENDIAN())  // Reverse buffer if byte order between source and target os doesn't match.
//        reverse(buffer);
    snp = *reinterpret_cast<TAlphabet*>(&buffer[0]);
    return sizeof(__uint32) + sizeof(TAlphabet);
}

// ----------------------------------------------------------------------------
// Function _readSnp()                                                    [Dna]
// ----------------------------------------------------------------------------

// Encodes SNP in delta position assuming, that 28 bits are sufficient to store the delta to the previous variant.
template <typename TAlphabet, typename TPosition, typename TReader>
inline unsigned
_readSnp(TAlphabet & snp,
         TPosition & deltaPos,
         TReader & reader,
         True /*isDnaCompressed*/)
{
    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));
    snp = static_cast<TAlphabet>((buffer[0] >> 5) & 3);  // Extract the SNP;
    buffer[0] &= static_cast<__uint8>(~0) >> 3;  // Disable leading three bits.
    if (SystemByteOrder::IS_LITTLE_ENDIAN())  // Reverse buffer to access correct position.
        reverse(buffer);
    TPosition tmp = *reinterpret_cast<__uint32*>(&buffer[0]);
    endianSwap(tmp, BigEndian(), HostByteOrder());
    deltaPos += tmp;
    return sizeof(__uint32);
}

template <typename TValue, typename TAlphabet, typename TStream, typename TConfig, typename TBoolFlag,
          typename TByteOrder>
inline void _readGdfBlock(DeltaMap<TValue, TAlphabet> & deltaMap,
                           RecordReader<TStream, SinglePass<> > & reader,
                           TBoolFlag isDnaCompressed,
                           TByteOrder byteOrder)
{
    typedef DeltaMap<TValue, TAlphabet> TDeltaMap;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnp;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDel;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSV>::Type TIndel;

    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIter;

//    typedef typename GetDeltaCoverageStore_<TDeltaMap>::Type TDeltaCoverageStore;
//    typedef typename Iterator<TDeltaCoverageStore, Standard>::Type TCoverageStoreIter;

    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));
    __uint32 blockRef = *reinterpret_cast<__uint32*>(&buffer[0]);
    endianSwap(blockRef, byteOrder, HostByteOrder());

    clear(buffer);
    readNChars(buffer, reader, sizeof(__uint32));
    __uint32 blockSize = *reinterpret_cast<__uint32*>(&buffer[0]);
    endianSwap(blockSize, byteOrder, HostByteOrder());

    // TODO(rmaerker): Rewrite to first load the block.
    __uint32 deltaRef = blockRef;
    unsigned oldMapLength = length(deltaMap);
//    unsigned counter = 0;
    while (blockSize != 0)
    {
        if (isBitSet(value(reader), BitsPerValue<char>::VALUE - 1))  // Dealing with a SNP.
        {
            TSnp snp;
            blockSize -= _readSnp(snp, deltaRef, reader, isDnaCompressed, byteOrder);  // Read alphabet dependent snp value.
//            _insert(deltaMap, deltaRef, length(deltaMap), snp);  // Record the snp.
            _addDeltaWithoutCoverage(deltaMap, deltaRef, snp, DeltaTypeSnp());
        }
        else  // Is an insertion or deletion.
        {
            clear(buffer);
            readNChars(buffer, reader, sizeof(__uint32));
//            if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                reverse(buffer);
            __uint32 tmp = *reinterpret_cast<__uint32*>(&buffer[0]);
            endianSwap(tmp, BigEndian(), HostByteOrder());
            deltaRef += tmp;
            blockSize -= length(buffer);
            clear(buffer);
            if(isBitSet(value(reader), BitsPerValue<char>::VALUE - 1))  // Read deletion info.
            {
                readNChars(buffer, reader, sizeof(__uint32));
//                if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                    reverse(buffer);
                __uint32 delSize = *reinterpret_cast<__uint32*>(&buffer[0]);
                endianSwap(delSize, BigEndian(), HostByteOrder());
                setBitTo(delSize, BitsPerValue<__uint32>::VALUE - 1, false);
//                _insert(deltaMap, deltaRef, length(deltaMap), static_cast<TDel>(delSize));  // Record the deletion.
                _addDeltaWithoutCoverage(deltaMap, deltaRef, static_cast<TDel>(delSize), DeltaTypeDel());
                blockSize -= length(buffer);
            }
            else
            {
                bool isIndel = false;
                readNChars(buffer, reader, sizeof(__uint32));
//                if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                    reverse(buffer);
                __uint32 insSize = *reinterpret_cast<__uint32*>(&buffer[0]);
                endianSwap(insSize, BigEndian(), HostByteOrder());
                blockSize -= length(buffer);
                if (isBitSet(insSize, BitsPerValue<__uint32>::VALUE - 2))
                {  // Handle indel.
                    setBitTo(insSize, BitsPerValue<__uint32>::VALUE - 2, false);
                    isIndel = true;
                }
                clear(buffer);
                readNChars(buffer, reader, insSize);
                TIns insSegment = buffer;  // TODO(rmaerker): Do we have to copy the data first?
                blockSize -= length(buffer);
                if (!isIndel)
                {
                    _addDeltaWithoutCoverage(deltaMap, deltaRef, insSegment, DeltaTypeIns());
//                    _insert(deltaMap, deltaRef, length(deltaMap), insSegment);
                }
                else
                {
                    clear(buffer);
                    readNChars(buffer, reader, sizeof(__uint32));
//                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
//                        reverse(buffer);
                    __uint32 delSize = *reinterpret_cast<__uint32*>(&buffer[0]);
                    endianSwap(delSize, BigEndian(), HostByteOrder());
                    _addDeltaWithoutCoverage(deltaMap, deltaRef, TIndel(delSize, insSegment), DeltaTypeSV());
//                    _insert(deltaMap, deltaRef, length(deltaMap), TIndel(delSize, insSegment));
                    blockSize -= length(buffer);
                }
            }
        }
    }

    // Afterwards we read the coverages and store them as well.
//    unsigned oldLength = length(deltaMap._deltaCoverageStore);
//    resize(deltaMap._deltaCoverageStore, length(deltaMap._deltaCoverageStore) + counter, Exact());
//    TCoverageStoreIter it = begin(deltaMap._deltaCoverageStore, Standard()) + oldLength;
//    TCoverageStoreIter itEnd = end(deltaMap._deltaCoverageStore, Standard());

    TMapIter it = begin(deltaMap, Standard()) + oldMapLength;
    TMapIter itEnd = end(deltaMap, Standard());

    for (; it != itEnd; ++it)
    {
        resize(deltaCoverage(it), getCoverageSize(deltaMap), Exact());
        _readDeltaCoverage(deltaCoverage(it), reader);
    }
}

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TStream, typename TConfig, typename TBoolFlag,
          typename TByteOrder>
inline void
_readGdfData(DeltaMap<TDeltaStore, TDeltaCoverageStore> & varStore,
             RecordReader<TStream, SinglePass<> > & reader,
             GdfHeader<TConfig> & gdfHeader,
             TBoolFlag is2BitCompressed,
             TByteOrder /*byteOrder*/)
{

//    unsigned lastRefPos = 0;
    // We need to store the number of blocks in the
    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));

    __int32 numOfBlocks = *reinterpret_cast<__uint32*>(&buffer[0]);
    endianSwap(numOfBlocks, TByteOrder(), HostByteOrder());

    while (numOfBlocks != 0)
    {
        _readGdfBlock(varStore, reader, gdfHeader, is2BitCompressed);
        --numOfBlocks;
    }
}

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TStream, typename TConfig, typename TBoolFlag>
inline void
_readGdfData(DeltaMap<TDeltaStore, TDeltaCoverageStore> & varStore,
             RecordReader<TStream, SinglePass<> > & reader,
             GdfHeader & gdfHeader,
             GdfFileConfiguration<TConfig> const & config,
             TBoolFlag is2BitCompressed)
{
    if (config.isLittleEndian)
        _readGdfData(varStore, reader, gdfHeader, is2BitCompressed, LittleEndian());
    else
        _readGdfData(varStore, reader, gdfHeader, is2BitCompressed, BigEndian());
}

template <typename TValue, typename TAlphabet, typename TConfig, typename TStream>
inline void
read(DeltaMap<TValue, TAlphabet> & deltaMap,
     GdfHeader & gdfHeader,
     GdfFileConfiguration<TConfig> & config,
     TStream  & stream,
     Gdf const & /*tag*/)
{
    RecordReader<TStream, SinglePass<> > reader(stream);
    readHeader(gdfHeader, reader, config, Gdf());
    setCoverageSize(deltaMap, length(gdfHeader.nameStore));

    if (config.ompressionMode == GdfIO::COMPRESSION_MODE_2_BIT_SNP_COMPRESSION)
        _readGdfData(deltaMap, reader, gdfHeader, config, True());
    else
        _readGdfData(deltaMap, reader, gdfHeader, config, False());
}

//template <typename TJournalSeq, typename TStream>
//inline int
//read(StringSet<TJournalSeq, Owner<JournaledSet> > & journalSet,
//     GdfHeader & jseqHeader,
//     RecordReader<TStream, SinglePass<> > & reader,
//     Gdf const & /*tag*/)
//{
//    typedef typename Value<TJournalSeq>::Type TAlphabet;
//    typedef typename Position<TJournalSeq>::Type TPosition;
//
//    if (empty(host(journalSet)))
//        return -1;
//
//    DeltaMap<TPosition, TAlphabet> deltaMap;
//    read(deltaMap, jseqHeader, reader, Gdf());
//    adaptTo(journalSet, deltaMap, 0, length(deltaMap), Serial());
//    return 0;
//}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_IO_READ_H_
