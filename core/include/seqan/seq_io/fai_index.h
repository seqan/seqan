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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Code for loading and writing FAI FASTA index files.
// ==========================================================================

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#ifndef SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_
#define SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FaiIndexEntry_
// ----------------------------------------------------------------------------

// Stores one entry in the FAI file.

class FaiIndexEntry_
{
public:
    // Name of the sequences.
    CharString name;
    // Length of the sequence.
    __uint64 sequenceLength;
    // Offset in the file.
    __uint64 offset;
    // Number of sequence characters per line.
    unsigned lineLength;
    // Number of bytes per line, including newline character(s).
    unsigned overallLineLength;

    FaiIndexEntry_() :
        sequenceLength(0), offset(0), lineLength(0), overallLineLength(0)
    {}
};

inline void
clear(FaiIndexEntry_ &entry)
{
    clear(entry.name);
    entry.sequenceLength = 0;
    entry.offset = 0;
    entry.lineLength = 0;
    entry.overallLineLength = 0;
}

// ----------------------------------------------------------------------------
// Class FaiIndex
// ----------------------------------------------------------------------------

/*!
 * @class FaiIndex
 * @headerfile <seqan/seq_io.h>
 * @brief Data structure for access to FAI indices.
 *
 * @signature class FaiIndex;
 *
 * FAI indices allow the rast random access to sequences or parts of sequences in a FASTA file.  Originally, they were
 * introduced in the <a href="http://samtools.sourceforge.net/samtools.shtml">samtools</a> program.
 *
 * Also see the <a href="http://trac.seqan.de/wiki/Tutorial/IndexedFastaIO">Indexed FASTA I/O Tutorial</a>.
 *
 * @section Example
 *
 * The following example demonstrates the usage of the FaiIndex class.
 *
 * @include demos/seq_io/fai_index_example.cpp
 *
 * The output is as follows:
 *
 * @include demos/seq_io/fai_index_example.cpp.stdout
 *
 * @fn FaiIndex::FaiIndex
 * @brief Constructor.
 *
 * @signature FaiIndex::FaiIndex();
 */

/**
.Class.FaiIndex
..cat:Input/Output
..signature:FaiIndex
..summary:Data type for storing FAI indices.
..wiki:Tutorial/IndexedFastaIO|Tutorial: Indexed FASTA I/O
..example.text:The following example demonstrate the usage of the FAIIndex class.
..example.file:demos/seq_io/fai_index_example.cpp
..include:seqan/seq_io.h

.Memvar.FaiIndex#FaiIndex
..class:Class.FaiIndex
..signature:FaiIndex()
..summary:The @Class.FaiIndex@ class only provides the default constructor.
*/

class FaiIndex
{
public:
    // TODO(holtgrew): The members are not publically documented for now. Keep it this way?

    // The name of the FASTA file.
    CharString fastaFilename;
    // The name of the FAI file.
    CharString faiFilename;

    // The index entries.
    String<FaiIndexEntry_> indexEntryStore;
    // A store for the sequence names.
    StringSet<CharString> seqNameStore;
    // A cache for fast access to the sequence name store.
    NameStoreCache<StringSet<CharString> > seqNameStoreCache;

    // We use this memory mapped string (opened read-only) to read from the file.
    String<char, MMap<> > mmapString;

    FaiIndex() :
        seqNameStoreCache(seqNameStore)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#clear
 * @brief Reset a FaiIndex object to the state after default construction.
 *
 * @signature void clear(faiIndex);
 *
 * @param[in,out] faiIndex The FaiIndex to clear.
 */

/**
.Function.FaiIndex#clear
..cat:Input/Output
..class:Class.FaiIndex
..signature:clear(faiIndex)
..param.faiIndex:The @Class.FaiIndex@ to reset.
...type:Class.FaiIndex
..summary:Reset a @Class.FaiIndex@ object to the state after default construction.
..include:seqan/seq_io.h
*/

inline void clear(FaiIndex & index)
{
    clear(index.fastaFilename);
    clear(index.faiFilename);
    clear(index.indexEntryStore);
    clear(index.seqNameStore);
    clear(index.seqNameStoreCache);
}

// ----------------------------------------------------------------------------
// Function getIdByName()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#getIdByName
 * @brief Return reference ID (numeric index in the file) of a sequence in a FAI file.
 *
 * @signature bool getIdByName(rID, faiIndex, name);
 *
 * @param[in]  faiIndex The FaiIndex to query.
 * @param[in]  name     The name of the sequence to look the id up for.  Type: @link ContainerConcept @endlink.
 * @param[out] rID      The id of the sequence is written here.
 *
 * @return bool true if a sequence with the given name is known in the index.
 */

/**
.Function.FaiIndex#getIdByName
..cat:Input/Output
..class:Class.FaiIndex
..signature:getIdByName(faiIndex, name, rID)
..summary:Return id (numeric index in the file) of a sequence in a FAI file.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.name:The name of the sequence to get the id for.
...type:Shortcut.CharString
..param.rID:The id of the sequence is written here.
...type:nolink:$unsigned$
..return:$bool$ indicating whether a sequence with this name is known in the index.
..include:seqan/seq_io.h
*/

template <typename TName, typename TId>
inline bool getIdByName(TId & rID, FaiIndex const & index, TName const & name)
{
    return getIdByName(rID, index.seqNameStoreCache, name);
}

// ----------------------------------------------------------------------------
// Function sequenceLength()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#sequenceLength
 * @brief Return length of the sequence with the given id in the FaiIndex.
 *
 * @signature __uint64 sequenceLength(faiIndex, rID);
 *
 * @param[in] faiIndex The FaiIndex to query.
 * @param[in] rID    The id of the sequence to get the length of.
 *
 * @return __uint64 The length of the sequence with index rID in faiIndex.
 */

/**
.Function.FaiIndex#sequenceLength
..cat:Input/Output
..class:Class.FaiIndex
..signature:__uint64 sequenceLength(faiIndex, rID)
..summary:Return length of the sequence with the given id in the @Class.FaiIndex@.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.rID:The id of the sequence to get the length of.
...type:nolink:$unsigned$
..return:$__uint64$ with the length of the sequence.
..include:seqan/seq_io.h
*/

template <typename TSeqId>
inline __uint64 sequenceLength(FaiIndex const & index, TSeqId rID)
{
    return index.indexEntryStore[rID].sequenceLength;
}

// TODO(holtgrew): Wrapper and template only here because sequenceLength in string_set_base.h is weird.

template <typename TSeqId>
inline __uint64 sequenceLength(FaiIndex & index, TSeqId rID)
{
    return index.indexEntryStore[rID].sequenceLength;
}

// ----------------------------------------------------------------------------
// Function sequenceName()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#sequenceName
 * @brief Return the name of the sequence with the given id in the FaiIndex.
 *
 * @signature CharString sequenceName(faiIndex, rID);
 *
 * @param[in] faiIndex The FaiIndex to query.
 * @param[in] rID    The index of the sequence.
 *
 * @return CharString The name of the sequence with the given id.
 */

/**
.Function.FaiIndex#sequenceName
..cat:Input/Output
..class:Class.FaiIndex
..signature:CharString sequenceName(faiIndex, rID)
..summary:Return the name of the sequence with the given id in the @Class.FaiIndex@.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.rID:The id of the sequence to get the name of.
...type:Shortcut.CharString
..return:@Shortcut.CharString@ with the name of the sequence.
..include:seqan/seq_io.h
*/

inline CharString const & sequenceName(FaiIndex const & index, unsigned rID)
{
    return index.indexEntryStore[rID].name;
}

// ----------------------------------------------------------------------------
// Function numSeqs()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#numSeqs
 * @brief Return the number of sequences known to a FaiIndex.
 *
 * @signature __uint64 numSeqs(faiIndex);
 *
 * @param[in] faiIndex The FaiIndex to query.
 *
 * @return __uint64 The number of sequences in the index.
 */

/**
.Function.FaiIndex#numSeqs
..cat:Input/Output
..class:Class.FaiIndex
..signature:numSeqs(faiIndex)
..summary:Return number of sequences known to an @Class.FaiIndex@.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..return:$__uint64$ with the number of sequences.
..include:seqan/seq_io.h
*/

inline __uint64 numSeqs(FaiIndex const & index)
{
    return length(index.indexEntryStore);
}

// ----------------------------------------------------------------------------
// Function readRegion()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#readRegion
 * @brief Read a region through an FaiIndex.
 *
 * @signature int readRegion(str, faiIndex, rID, beginPos, endPos);
 * @signature int readRegion(str, faiIndex, region);
 *
 * @param[out] str      The @link String @endlink to read the sequence into.
 * @param[in]  faiIndex The FaiIndex to read from.
 * @param[in]  rID    The id of the sequence to read (Type: <tt>unsigned).
 * @param[in]  beginPos The begin position of the region to read (Type: <tt>unsigned).
 * @param[in]  endPos   The end position of the region to read  (Type: <tt>unsigned).
 * @param[in]  region   The @link GenomicRegion @endlink to read.
 *
 * @return int 0 on success, non-0 on errors.
 */

/**
.Function.FaiIndex#readRegion
..cat:Input/Output
..class:Class.FaiIndex
..signature:readRegion(str, faiIndex, rID, beginPos, endPos)
..signature:readRegion(str, faiIndex, region);
..summary:Load the infix of a sequence from a @Class.FaiIndex@.
..description:The given region is loaded from the indexed FASTA file.
..param.str:The sequence infix is written into this string.
...type:Class.String
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.rID:The index of the sequence in the file.
...type:nolink:$unsigned$
..param.beginPos:The begin position of the infix to write to $str$.
...type:nolink:$unsigned$
..param.endPos:The end position of the infix to write to $str$.
...type:nolink:$unsigned$
..param.region:The @Class.GenomicRegion@ to read.
...type:Class.GenomicRegion
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

template <typename TValue, typename TSpec, typename TSeqId, typename TBeginPos, typename TEndPos>
inline void readRegion(String<TValue, TSpec> & str,
                       FaiIndex const & index,
                       TSeqId rID,
                       TBeginPos beginPos,
                       TEndPos endPos)
{
    typedef String<char, MMap<> > const TFastaFile;

    FaiIndexEntry_ const & entry = index.indexEntryStore[rID];

    // Limit region to the infix, make sure that beginPos < endPos, compute character to read.
    TEndPos seqLen = entry.sequenceLength;;
    beginPos = std::min((TEndPos)beginPos, seqLen);
    endPos = std::min(std::max((TEndPos)beginPos, endPos), seqLen);
    TEndPos toRead = endPos - beginPos;

    // seek to start position
    DirectionIterator<TFastaFile, Input>::Type reader = directionIterator(index.mmapString, Input());
    setPosition(
        reader,
        index.indexEntryStore[rID].offset +
        (beginPos / entry.lineLength) * entry.overallLineLength);

    // set up countdowns
    CountDownFunctor<NotFunctor<IsWhitespace> > countDownIgnore(beginPos % entry.lineLength);
    CountDownFunctor<NotFunctor<IsWhitespace> > countDownData(toRead);
    IsWhitespace ignWhiteSpace;

    // read characters
    clear(str);
    skipUntil(reader, countDownIgnore);
    readUntil(str, reader, countDownData, ignWhiteSpace);
    if (!countDownData)
        SEQAN_THROW(UnexpectedEnd());
//
//
//    typedef typename Iterator<String<char, MMap<> > const, Standard>::Type TSourceIter;
//    typedef typename Iterator<String<TValue, TSpec>, Standard>::Type TTargetIter;
//    TSourceIter itSource = begin(index.mmapString, Standard());
//    __uint64 offset = index.indexEntryStore[rID].offset;
//    // First, compute offset of the completely filled lines.
//    unsigned numLines = beginPos / index.indexEntryStore[rID].lineLength;
//    unsigned numBytes = numLines * index.indexEntryStore[rID].overallLineLength;
//    // Then, compute overall offset by adding remaining bytes, too.
//    numBytes += beginPos % index.indexEntryStore[rID].lineLength;
//    offset += numBytes;
//    // Advance iterator in MMap file.
//    itSource += offset;
//
//    // Copy out the characters from FASTA file and convert via iterator assignment to target string's type.
//    resize(str, toRead, TValue());
//    TTargetIter itTarget = begin(str, Standard());
//    for (TEndPos i = 0; i < toRead; )
//    {
//        if (isspace(*itSource))
//        {
//            ++itSource;
//            continue;  // Skip spaces.
//        }
//        *itTarget = *itSource;
//        ++itTarget;
//        ++itSource;
//        ++i;
//    }
//
//    return true;
}

template <typename TValue, typename TSpec>
inline bool readRegion(String<TValue, TSpec> & str,
                       FaiIndex const & index,
                       GenomicRegion const & region)
{
    unsigned rID = region.rID;
    if (rID == GenomicRegion::INVALID_ID)
        if (!getIdByName(rID, index, region.seqName))
            return false;  // Sequence with this name could not be found.

    unsigned beginPos = (region.beginPos != GenomicRegion::INVALID_POS)? region.beginPos : 0;
    unsigned endPos = (region.endPos != GenomicRegion::INVALID_POS)? region.endPos : sequenceLength(index, rID);
    readRegion(str, index, rID, beginPos, endPos);
    return true;
}

// ----------------------------------------------------------------------------
// Function readSequence()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#readSequence
 * @brief Load a whole sequence from a FaiIndex.
 *
 * @signature int readSequence(str, faiIndex, rID);
 *
 * @param[out] str      The @link String @endlink to read into.
 * @param[in]  faiIndex The FaiIndex to read from.
 * @param[in]  seqID    The index of the sequence in the file.
 *
 * @return int 0 on success, non-0 on errors.
 */

/**
.Function.FaiIndex#readSequence
..cat:Input/Output
..class:Class.FaiIndex
..signature:readSequence(str, faiIndex, rID)
..summary:Load a whole sequence from an @Class.FaiIndex@.
..param.str:The sequence is written into this string.
...type:Class.String
..param.faiIndex:The @Class.FaiIndex@ to use.
...type:Class.FaiIndex
..param.rID:The index of the sequence in the file.
...type:nolink:$unsigned$
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

template <typename TValue, typename TSpec>
inline void readSequence(String<TValue, TSpec> & str, FaiIndex const & index, unsigned rID)
{
    readRegion(str, index, rID, 0u, sequenceLength(index, rID));
}

// ----------------------------------------------------------------------------
// Function read()
// ----------------------------------------------------------------------------

/*!
 * @fn FaiIndex#read
 * @brief Read a FAI index from file.
 *
 * @signature int read(faiIndex, fastaFileName[, faiFileName]);
 *
 * @param[out] faiIndex      The FaiIndex to read into.
 * @param[in]  fastaFileName Path to the FASTA file to read.  Type: <tt>char const *</tt>.
 * @param[in]  faiFileName   Path to the FAI file to read.  Type: <tt>char const *</tt>.  Defaults to
 *                           <tt>"${fastaFileName}.fai"</tt>.
 *
 * @return int 0 on success, non-0 on errors.
 */

/**
.Function.FaiIndex#read
..cat:Input/Output
..class:Class.FaiIndex
..signature:read(faiIndex, fastaFileName[, faiFileName])
..summary:Read a FAI index.
..param.faiIndex:The @Class.FaiIndex@ object to read the file into.
...type:Class.FaiIndex
..param.fastaFileName:The name of the FASTA file to read.
...type:nolink:$char const *$
..param.faiFileName:The name of the FAI file to read.
...default:$fastaFileName + ".fai"$.
...type:nolink:$char const *$
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

template <typename TFwdIterator>
inline void
readRecord(FaiIndexEntry_ & entry, TFwdIterator & reader, CharString & buffer)
{
    clear(entry);

    // Read REF_NAME.
    readUntil(entry.name, reader, IsTab());

    skipOne(reader, IsTab());   // Must have been on tab

    // Read SEQ_LENGTH.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(entry.sequenceLength, buffer);

    skipOne(reader, IsTab());   // Must have been on tab

    // Read OFFSET.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(entry.offset, buffer);

    skipOne(reader, IsTab());   // Must have been on tab

    // Read LINE_LENGTH.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(entry.lineLength, buffer);

    skipOne(reader, IsTab());   // Must have been on tab

    // Read OVERALL_LINE_LENGTH.
    clear(buffer);
    readUntil(buffer, reader, OrFunctor<IsTab, IsNewline>());
    lexicalCast(entry.overallLineLength, buffer);

    skipLine(reader);           // Skip over line ending.
}

inline bool open(FaiIndex & index, char const * fastaFilename, char const * faiFilename)
{
    clear(index);  // Also clears filename, thus backup above and restore below.
    index.fastaFilename = fastaFilename;
    index.faiFilename = faiFilename;

    if (!open(index.mmapString, toCString(fastaFilename), OPEN_RDONLY))
        return false;  // Could not open file.

    // Open file.
    std::ifstream faiStream(toCString(index.faiFilename));
    if (!faiStream.good())
        return false;  // problem opening
    DirectionIterator<std::ifstream, Input>::Type reader = directionIterator(faiStream, Input());

    // Read FAI file.
    CharString buffer;
    FaiIndexEntry_ entry;
    while (!atEnd(reader))
    {
        readRecord(entry, reader, buffer);
        appendValue(index.seqNameStore, entry.name);
        appendValue(index.indexEntryStore, entry);
    }

    // Recreate name store cache.
    refresh(index.seqNameStoreCache);
    return true;
}

inline bool open(FaiIndex & index, char const * fastaFilename)
{
    std::string faiFilename = fastaFilename;
    faiFilename += ".fai";
    return open(index, fastaFilename, toCString(faiFilename));
}

// ---------------------------------------------------------------------------
// Function save()
// ---------------------------------------------------------------------------

/*!
 * @fn FaiIndex#save
 * @brief Write out an FaiIndex object.
 *
 * @signature int save(faiIndex[, faiFileName]);
 *
 * @param[in] faiIndex    The FaiIndex to write out.
 * @param[in] faiFileName The name of the FAI file to write to.  This parameter is optional only if the FAI index knows
 *                        the FAI file name from a previous @link FaiIndex#build @endlink call.  By default, the FAI
 *                        file name from the previous call to @link FaiIndex#build @endlink is used.  Type: <tt>char
 *                        const *</tt>.
 *
 * @return int 0 on success, 1 on errors.
 */

/**
.Function.FaiIndex#save
..cat:Input/Output
..class:Class.FaiIndex
..signature:write(faiIndex[, faiFileName])
..summary:Write out an @Class.FaiIndex@ object.
..param.faiIndex:The @Class.FaiIndex@ object to write out.
...type:Class.FaiIndex
..param.faiFileName:The name of the FAI file to write to.
...remarks:This is optional only if the FAI index knows the FAI file name from a previous @Function.FaiIndex#build@ call.
...default:The FAI file name from the previous call to @Function.FaiIndex#build@, if any.
...type:nolink:$char const *$
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

inline bool save(FaiIndex const & index, char const * faiFilename)
{
    // Open index files.
    std::ofstream file(faiFilename);
    if (!file.good())
        return false;

    file.exceptions(std::ios_base::badbit);
    for (unsigned i = 0; i < length(index.indexEntryStore); ++i)
    {
        FaiIndexEntry_ const & entry = index.indexEntryStore[i];
        file << entry.name << '\t' << entry.sequenceLength << '\t' << entry.offset << '\t'
             << entry.lineLength << '\t' << entry.overallLineLength << '\n';
    }
    return true;
}

inline bool save(FaiIndex const & index)
{
    if (empty(index.faiFilename))
        return false;  // Cannot write out if faiFilename member is empty.
    return save(index, toCString(index.faiFilename));
}

// ----------------------------------------------------------------------------
// Function getRecordInfo(Fastq);
// ----------------------------------------------------------------------------

template <typename TFwdIterator>
inline void getRecordInfo(FaiIndexEntry_ & entry, TFwdIterator & iter, Fasta)
{
    typedef EqualsChar<'>'> TFastaBegin;

    clear(entry);

    skipUntil(iter, TFastaBegin());     // forward to the next '>'
    entry.offset = position(iter);      // store file position
    skipOne(iter);                      // skip '>'

    readUntil(entry.name, iter, IsWhitespace());    // read Fasta id (up to first whitespace)
    skipLine(iter);
    entry.offset = position(iter);      // store offset

    FaiIndexEntry_ temp;
    FaiIndexEntry_ *entryPtr = &entry;
    FaiIndexEntry_ *cmpPtr = &entry;
    OrFunctor<IsNewline, CountFunctor<NotFunctor<IsWhitespace> > > countCharsPerLine;

    while (!atEnd(iter) && !TFastaBegin()(value(iter)))
    {
        // check for consistency
        if (cmpPtr->lineLength != entry.lineLength ||
            cmpPtr->overallLineLength != entry.overallLineLength)
        {
            SEQAN_THROW(ParseError("FastaIndex: Record has inconsistent line lengths or line endings"));
        }

        __uint64 start = position(iter);

        // skip to the end of line and count non-whitespace characters
        clear(countCharsPerLine.func2);
        skipUntil(iter, countCharsPerLine);
        entry.sequenceLength += value(countCharsPerLine.func2);

        // determine line length in bases
        entryPtr->lineLength = value(countCharsPerLine.func2);

        // skip linebreak
        skipLine(iter);

        // determine line length in bytes
        entryPtr->overallLineLength = (__uint64)position(iter) - start;
        cmpPtr = entryPtr;
        entryPtr = &temp;
    }
}

// ---------------------------------------------------------------------------
// Function build()
// ---------------------------------------------------------------------------

/*!
 * @fn FaiIndex#build
 * @brief Create a FaiIndex from FASTA file.
 *
 * @signature bool build(faiIndex, seqFileName[, faiFileName]);
 *
 * @param[out] faiIndex    The FaiIndex to build into.
 * @param[in]  seqFileName Path to the FASTA file to build an index for.  Type: <tt>char const *</tt>.
 * @param[in]  faiFileName Path to the FAI file to use as the index file.  Type: <tt>char const *</tt>.
 *                         Default: <tt>"${seqFileName}.fai"</tt>.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> otherwise.
 */

/**
.Function.FaiIndex#build
..class:Class.FaiIndex
..summary:Create @Class.FaiIndex@ from FASTA file.
..signature:build(faiIndex, fastaFilename[, faiFilename])
..description:The index can later be written out with @Function.FaiIndex#save@ and be loaded again using @Function.FaiIndex#read@.
..param.faiIndex:@Class.FaiIndex@ to write index to.
...type:Class.FaiIndex
..param.fastaFilename:Name of FASTA file to build an index for.
...type:Shortcut.CharString
..param.faiFilename:Name of FAI index file, stored in $faiIndex$. Optional.
...default:$fastaFilename + ".fai"$
...type:nolink:$char const *$
..returns:$int$, equal to 0 on success, != 0 otherwise.
..include:seqan/stream.h
 */

inline bool build(FaiIndex & index, char const * seqFilename, char const * faiFilename)
{
    index.fastaFilename = seqFilename;
    index.faiFilename = faiFilename;

    SeqFileIn seqFile(seqFilename);
    DirectionIterator<SeqFileIn, Input>::Type iter = directionIterator(seqFile, Input());
    
    if (!isEqual(format(seqFile), seqan::Fasta()))
        return false;  // Invalid format, not FASTA.

    // Clear everything.
    clear(index.seqNameStore);
    clear(index.seqNameStoreCache);
    clear(index.indexEntryStore);

    // Create FastaIndex
    FaiIndexEntry_ entry;
    while (!atEnd(iter))
    {
        getRecordInfo(entry, iter, Fasta());
        appendValue(index.seqNameStore, entry.name);
        appendValue(index.indexEntryStore, entry);
    }

    // Recreate name store cache.
    refresh(index.seqNameStoreCache);

    return true;
}

inline int build(FaiIndex & index, char const * seqFilename)
{
    CharString faiFilename(seqFilename);
    append(faiFilename, ".fai");
    return build(index, seqFilename, toCString(faiFilename));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_
