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
inline int
_readGdfHeaderFileInfo(GdfHeader<TConfig> & gdfHeader, RecordReader<TStream, SinglePass<> >  & reader)
{
    CharString buffer;

    // Skip the file prefix.
    skipNChars(reader, 2);
    // Read file version.
    clear(buffer);
    readNChars(buffer, reader, length(GdfIO::FILE_VERSION_KEY));
//        streamReadBlock(&buffer[0], stream, length(GdfIO::FILE_VERSION_KEY));
    if (!startsWith(buffer, GdfIO::FILE_VERSION_KEY))
        return GdfIO::UNSUPPORTED_FILE_VERSION_ERROR;      // TODO(rmaerker): Either use error codes or use exceptions.

    // Skip = sign
    skipNChars(reader, 1);
    clear(buffer);
    readNChars(buffer, reader, length(GdfIO::FILE_VERSION_VALUE_PREFIX));
    if (!startsWith(buffer, GdfIO::FILE_VERSION_VALUE_PREFIX))
        return GdfIO::UNSUPPORTED_FILE_VERSION_ERROR;      // TODO(rmaerker): Either use error codes or use exceptions.
    // Read big version number.
//        readD(buffer[0], stream);
    clear(buffer);
    int fileVersionBig;
    clear(buffer);
    readNChars(buffer, reader, 1);
    lexicalCast2(fileVersionBig, buffer);

    if (fileVersionBig != GdfIO::FILE_VERSION_BIG)
        return GdfIO::UNSUPPORTED_FILE_VERSION_ERROR;
    // Read version number separator.
    skipNChars(reader, 1);
    // Read little version number.
    clear(buffer);
    int fileVersionLittle;
    clear(buffer);
    readNChars(buffer, reader, 1);
    lexicalCast2(fileVersionLittle, buffer);

    if (fileVersionLittle != GdfIO::FILE_VERSION_LITTLE)
        return GdfIO::UNSUPPORTED_FILE_VERSION_ERROR;

    // Skip until next line.
    skipLine(reader);

    // Read byte order.
    skipNChars(reader, 2);  // skip prefix.
    clear(buffer);
    readNChars(buffer, reader, length(GdfIO::FILE_ENDIANNESS_KEY));
    if (buffer != GdfIO::FILE_ENDIANNESS_KEY)
        return -1;      // TODO(rmaerker): Either use error codes or use exceptions.
    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);
    gdfHeader._fileInfos._byteOrder = buffer == GdfIO::FILE_ENDIANNESS_LITTLE;

    // Read block size.
    clear(buffer);
    skipNChars(reader, 2);  // skip prefix.
    readNChars(buffer, reader, length(GdfIO::FILE_BLOCKSIZE_KEY));
    if (buffer != GdfIO::FILE_BLOCKSIZE_KEY)
        return -1;      // TODO(rmaerker): Either use error codes or use exceptions.
    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);
    lexicalCast2(gdfHeader._fileInfos._blockSize, buffer);

    //Read snp compression status.
    clear(buffer);
    skipNChars(reader, 2);  // skip prefix.
    readNChars(buffer, reader, length(GdfIO::FILE_SNP_COMPRESSION_KEY));
    if (buffer != GdfIO::FILE_SNP_COMPRESSION_KEY)
        return -1;      // TODO(rmaerker): Either use error codes or use exceptions.
    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);
    if (buffer == GdfIO::FILE_SNP_COMPRESSION_2BIT)
        gdfHeader._fileInfos._snpCompression = true;
    else
    {
        SEQAN_ASSERT_EQ(buffer, GdfIO::FILE_SNP_COMPRESSION_GENERIC);
        gdfHeader._fileInfos._snpCompression = false;
    }
    return 0;
}

template <typename TConfig, typename TStream>
inline int
_readGdfHeaderRefInfo(GdfHeader<TConfig> & gdfHeader, RecordReader<TStream, SinglePass<> >  & reader)
{
    CharString buffer;
    skipNChars(reader, 2);  // skip prefix.
    // Read refernece ID.
    clear(buffer);
    readNChars(buffer, reader, length(GdfIO::REFERENCE_ID_KEY));
    if (!startsWith(buffer, GdfIO::REFERENCE_ID_KEY))
        return GdfIO::UNSUPPORTED_REFERNCE_INFORMATION_ERROR;      // TODO(rmaerker): Either use error codes or use exceptions.

    // Skip = sign
    skipNChars(reader, 1);
    clear(gdfHeader._refInfos._refId);
    readLine(gdfHeader._refInfos._refId, reader);  // Read value.

    skipNChars(reader, 2);  // skip prefix.
    clear(buffer);
    readNChars(buffer, reader, length(GdfIO::REFERENCE_FILE_KEY));
    if (!startsWith(buffer, GdfIO::REFERENCE_FILE_KEY))
        return GdfIO::UNSUPPORTED_REFERNCE_INFORMATION_ERROR;      // TODO(rmaerker): Either use error codes or use exceptions.

    // Skip = sign
    skipNChars(reader, 1);
    clear(gdfHeader._refInfos._refFile);
    readLine(gdfHeader._refInfos._refFile, reader);  // Read value.

    skipNChars(reader, 2);  // skip prefix.
    clear(buffer);
    readNChars(buffer, reader, length(GdfIO::REFERENCE_HASH_KEY));
    if (!startsWith(buffer, GdfIO::REFERENCE_HASH_KEY))
        return GdfIO::UNSUPPORTED_REFERNCE_INFORMATION_ERROR;      // TODO(rmaerker): Either use error codes or use exceptions.

    // Skip = sign
    skipNChars(reader, 1);
    clear(buffer);
    readLine(buffer, reader);  // Read value.
    lexicalCast2(gdfHeader._refInfos._refHash, buffer);
    return 0;
}

// ----------------------------------------------------------------------------
// Function _readSeqNames()
// ----------------------------------------------------------------------------

template <typename TConfig, typename TStream>
inline int
_readSeqNames(GdfHeader<TConfig> & gdfHeader, RecordReader<TStream, SinglePass<> >  & reader)
{
    CharString buffer;
    // Read the sequence names.
    SEQAN_ASSERT_EQ(value(reader), '!');
    skipNChars(reader, 1);
    SEQAN_ASSERT_EQ(value(reader), '!');

    while (true)
    {
        goNext(reader);
        if (value(reader) == GdfIO::SEQ_NAMES_SEPARATOR[0])
        {
            appendValue(*gdfHeader._nameStorePtr, buffer);
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
inline int
readHeader(GdfHeader<TConfig> & gdfHeader,
           RecordReader<TStream, SinglePass<> > & reader,
           Gdf const & /*tag*/)
{
    _readGdfHeaderFileInfo(gdfHeader, reader);
    _readGdfHeaderRefInfo(gdfHeader, reader);

    // Read variable header informations.
//    while (value(reader) == '#')
//    {
//        skipNChars(reader, 2);
//        GdfHeaderRecord record;
//        clear(record._key);
//        clear(record._value);
//        readUntilChar(record._key, reader, '=');
//        skipNChars(reader, 1);
//        readLine(record._value, reader);
//        appendValue(jseqHeader._headerRecords, record);
//    }

    return _readSeqNames(gdfHeader, reader);
}

// ----------------------------------------------------------------------------
// Function _readSnp()
// ----------------------------------------------------------------------------

// Writes SNP in separat value after delta pos.
template <typename TAlphabet, typename TPosition, typename TReader, typename TConfig>
inline unsigned
_readSnp(TAlphabet & snp,
         TPosition & deltaPos,
         TReader & reader,
         GdfHeader<TConfig> const & gdfHeader,
         False /*isDnaCompressed*/)
{
    SEQAN_ASSERT_NOT(isBitSet(deltaPos, BitsPerValue<__uint32>::VALUE -1));  // The last bit should not be set.

    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));
    // Read the offset for the encoded delta.
    setBitTo(buffer[0], BitsPerValue<__uint8>::VALUE -1, false);  // Reset the MSB to extract delta.
    if (SystemByteOrder::IS_LITTLE_ENDIAN())  // Reverse buffer if system is little endian.
        reverse(buffer);
    deltaPos += *reinterpret_cast<__uint32*>(&buffer[0]);  // Translate into integer.
    clear(buffer);
    readNChars(buffer, reader, sizeof(TAlphabet));
    // Read the value for the snp.
    if (gdfHeader._fileInfos._byteOrder != SystemByteOrder::IS_LITTLE_ENDIAN())  // Reverse buffer if byte order between source and target os doesn't match.
        reverse(buffer);
    snp = *reinterpret_cast<TAlphabet*>(&buffer[0]);
    return sizeof(__uint32) + sizeof(TAlphabet);
}

// ----------------------------------------------------------------------------
// Function _readSnp()                                                    [Dna]
// ----------------------------------------------------------------------------

// Encodes SNP in delta position assuming, that 28 bits are sufficient to store the delta to the previous variant.
template <typename TAlphabet, typename TPosition, typename TReader, typename TConfig>
inline unsigned
_readSnp(TAlphabet & snp,
         TPosition & deltaPos,
         TReader & reader,
         GdfHeader<TConfig> const & /*gdfHeader*/,
         True /*isDnaCompressed*/)
{
    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));
    snp = static_cast<TAlphabet>((buffer[0] >> 5) & 3);  // Extract the SNP;
    buffer[0] &= static_cast<__uint8>(~0) >> 3;  // Disable leading three bits.
    if (SystemByteOrder::IS_LITTLE_ENDIAN())  // Reverse buffer to access correct position.
        reverse(buffer);
    deltaPos += *reinterpret_cast<__uint32*>(&buffer[0]);  // Translate into integer.
    return sizeof(__uint32);
}

template <typename TValue, typename TAlphabet, typename TStream, typename TConfig, typename TBoolFlag>
inline void _readGdfBlock(DeltaMap<TValue, TAlphabet> & deltaMap,
                           RecordReader<TStream, SinglePass<> > & reader,
                           GdfHeader<TConfig> & gdfHeader,
                           TBoolFlag isDnaCompressed)
{
    typedef DeltaMap<TValue, TAlphabet> TDeltaMap;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type TSnp;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDel;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type TIndel;

    typedef typename GetDeltaCoverageStore_<TDeltaMap>::Type TDeltaCoverageStore;
    typedef typename Iterator<TDeltaCoverageStore, Standard>::Type TCoverageStoreIter;

    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));
    if (gdfHeader._fileInfos._byteOrder != SystemByteOrder::IS_LITTLE_ENDIAN())
        reverse(buffer);
    __uint32 blockRef = *reinterpret_cast<__uint32*>(&buffer[0]);

    clear(buffer);
    readNChars(buffer, reader, sizeof(__uint32));
    if (gdfHeader._fileInfos._byteOrder != SystemByteOrder::IS_LITTLE_ENDIAN())
        reverse(buffer);
    __uint32 blockSize = *reinterpret_cast<__uint32*>(&buffer[0]);

    __uint32 deltaRef = blockRef;
    unsigned counter = 0;
    for (; blockSize != 0; ++counter)
    {
        if (isBitSet(value(reader), BitsPerValue<char>::VALUE - 1))  // Dealing with a SNP.
        {
            TSnp snp;
            blockSize -= _readSnp(snp, deltaRef, reader, gdfHeader, isDnaCompressed);  // Read alphabet dependent snp value.
            _insert(deltaMap, deltaRef, length(deltaMap), snp);  // Record the snp.
        }
        else  // Is an insertion or deletion.
        {
            clear(buffer);
            readNChars(buffer, reader, sizeof(__uint32));
            if (SystemByteOrder::IS_LITTLE_ENDIAN())
                reverse(buffer);
            deltaRef += *reinterpret_cast<__uint32*>(&buffer[0]);  // Translate into integer.
            blockSize -= length(buffer);
            clear(buffer);
            if(isBitSet(value(reader), BitsPerValue<char>::VALUE - 1))  // Read deletion info.
            {
                readNChars(buffer, reader, sizeof(__uint32));
                if (SystemByteOrder::IS_LITTLE_ENDIAN())
                    reverse(buffer);
                __uint32 delSize = *reinterpret_cast<__uint32*>(&buffer[0]);
                setBitTo(delSize, BitsPerValue<__uint32>::VALUE - 1, false);
                _insert(deltaMap, deltaRef, length(deltaMap), static_cast<TDel>(delSize));  // Record the deletion.
                blockSize -= length(buffer);
            }
            else
            {
                bool isIndel = false;
                readNChars(buffer, reader, sizeof(__uint32));
                if (SystemByteOrder::IS_LITTLE_ENDIAN())
                    reverse(buffer);
                __uint32 insSize = *reinterpret_cast<__uint32*>(&buffer[0]);
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
                    _insert(deltaMap, deltaRef, length(deltaMap), insSegment);
                else
                {
                    clear(buffer);
                    readNChars(buffer, reader, sizeof(__uint32));
                    if (SystemByteOrder::IS_LITTLE_ENDIAN())
                        reverse(buffer);
                    __uint32 delSize = *reinterpret_cast<__uint32*>(&buffer[0]);
                    _insert(deltaMap, deltaRef, length(deltaMap), TIndel(delSize, insSegment));
                    blockSize -= length(buffer);
                }
            }
        }
    }

    unsigned oldLength = length(deltaMap._deltaCoverageStore);
    resize(deltaMap._deltaCoverageStore, length(deltaMap._deltaCoverageStore) + counter, Exact());
    TCoverageStoreIter it = begin(deltaMap._deltaCoverageStore, Standard()) + oldLength;
    TCoverageStoreIter itEnd = end(deltaMap._deltaCoverageStore, Standard());

    for (; it != itEnd; ++it)
    {
        resize(*it, coverageSize(deltaMap), Exact());
        _readDeltaCoverage(*it, reader);
    }
}

template <typename TDeltaStore, typename TDeltaCoverageStore, typename TStream, typename TConfig, typename TBoolFlag>
inline void
_readGdfData(DeltaMap<TDeltaStore, TDeltaCoverageStore> & varStore,
             RecordReader<TStream, SinglePass<> > & reader,
             GdfHeader<TConfig> & gdfHeader,
             TBoolFlag isDnaCompressed)
{

//    unsigned lastRefPos = 0;
    // We need to store the number of blocks in the
    CharString buffer;
    readNChars(buffer, reader, sizeof(__uint32));

    // Reverse array if endianness doesn't match.
    if (gdfHeader._fileInfos._byteOrder != SystemByteOrder::IS_LITTLE_ENDIAN())
        reverse(buffer);
    unsigned numOfBlocks = *reinterpret_cast<__uint32*>(&buffer[0]);

    while (numOfBlocks != 0)
    {
        _readGdfBlock(varStore, reader, gdfHeader, isDnaCompressed);
        --numOfBlocks;
    }
}

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

template <typename TValue, typename TAlphabet, typename TConfig, typename TStream>
inline int
read(DeltaMap<TValue, TAlphabet> & deltaMap,
     GdfHeader<TConfig> & gdfHeader,
     TStream  & stream,
     Gdf const & /*tag*/)
{
    RecordReader<TStream, SinglePass<> > reader(stream);
    readHeader(gdfHeader, reader, Gdf());
    setCoverageSize(deltaMap, length(*gdfHeader._nameStorePtr));

    if (gdfHeader._fileInfos._snpCompression)
        _readGdfData(deltaMap, reader, gdfHeader, True());
    else
        _readGdfData(deltaMap, reader, gdfHeader, False());
    return 0;
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
