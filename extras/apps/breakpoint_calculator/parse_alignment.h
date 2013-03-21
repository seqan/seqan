// ==========================================================================
//                           Breakpoint Calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_
#define SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_

#include <seqan/stream.h>
#include <seqan/align.h>

#include "parse_alignment.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tags Xmfa, XmfaSwap, Maf
// ----------------------------------------------------------------------------

struct Xmfa_;
typedef Tag<Xmfa_> Xmfa;

struct XmfaSwap_;
typedef Tag<XmfaSwap_> XmfaSwap;

struct Maf_;
typedef Tag<Maf_> Maf;

// ----------------------------------------------------------------------------
// Class AlignmentBlockRow
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
struct AlignmentBlockRow
{
    TSize rowNum;             // row number
    CharString chromosomeId;  // chromosome identifier

    TPosition startPos;      // in source seq
    TPosition endPos;        // in source seq

    bool orientation;
    // startPos is always smaller or equal endPos also if orientation is false

    AlignmentBlockRow() {}

    AlignmentBlockRow(TSize id, CharString chr, TPosition start, TPosition end, bool ori)
    {
        rowNum = id;
        chromosomeId = chr;
        startPos = start;
        endPos = end;
        orientation = ori;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

template <typename TStream, typename TPassSpec, typename TTag>
int
skipHeader(RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           TTag const &)
{
    CharString buffer;
    int res = 0;

    // Skip comment lines
    while (value(recordReader) == '#')
    {
        res = skipLine(recordReader);
        if (res)
            return res;

        if (atEnd(recordReader))
            return 1;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                                 [Xmfa]
// ----------------------------------------------------------------------------

// Reads one collinear alignment block from a file in XMFA format

template <typename TSize, typename TStream, typename TPassSpec>
int
readRecord(std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           unsigned lineId,
           std::set<CharString> const & /*seqIds*/,
           bool swapPos,
           Xmfa const &)
{
    typedef AlignmentBlockRow<TSize, TSize> TRow;

    CharString buffer;
    int res = 0;

    unsigned startPos, endPos;
    CharString id, chromosomeId = "";
    bool orientation;

    while (!skipChar(recordReader, '>'))
    {
        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read seq id
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        id = buffer;

        res = skipChar(recordReader, ':');
        if (res)
            return res;

        // read start position
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        startPos = lexicalCast<unsigned>(buffer);

        res = skipChar(recordReader, '-');
        if (res)
            return res;

        // read end position
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        endPos = lexicalCast<unsigned>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read orientation
        clear(buffer);
        res = readNChars(buffer, recordReader, 1);
        if (res)
            return res;

        if (buffer  == "+")
            orientation = true;
        else if (buffer == "-")
            orientation = false;
        else
            return 1;

        if (swapPos && orientation == false)
        {
            unsigned help = startPos;
            startPos = endPos;
            endPos = help;
        }
        if (startPos != 0)
            --startPos;

        // skip rest of line
        res = skipLine(recordReader);
        if (res)
            return res;

        // skip gapped sequence
        while (value(recordReader) != '>' && value(recordReader) != '=')
        {
            res = skipLine(recordReader);
            if (res)
                return res;
        }

        if (endPos != startPos)
            appendValue(idToRowsMap[id], TRow(lineId, chromosomeId, startPos, endPos, orientation));
    }

    if (value(recordReader) != '=')
        return 1;

    skipUntilLineBeginsWithChar(recordReader, '>');

    return 0;
}

template <typename TSize, typename TStream, typename TPassSpec>
int
readRecord(std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           unsigned lineId,
           std::set<CharString> const & seqIds,
           Xmfa const &)
{
    return readRecord(idToRowsMap, recordReader, lineId, seqIds, false, Xmfa());
}

template <typename TSize, typename TStream, typename TPassSpec>
int
readRecord(std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           unsigned lineId,
           std::set<CharString> const & seqIds,
           XmfaSwap const &)
{
    return readRecord(idToRowsMap, recordReader, lineId, seqIds, true, Xmfa());
}

// ----------------------------------------------------------------------------
// Function readRecord()                                                  [Maf]
// ----------------------------------------------------------------------------

// Reads one collinear alignment block from a file in MAF format

template <typename TSize, typename TStream, typename TPassSpec>
int
readRecord(std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           unsigned lineId,
           std::set<CharString> const & seqIds,
           Maf const &)
{
    typedef AlignmentBlockRow<TSize, TSize> TRow;

    CharString buffer;
    int res = 0;

    unsigned startPos, len, seqLen;
    CharString seqId, chromosomeId = "";
    bool orientation;

    skipWhitespaces(recordReader);
    res = skipChar(recordReader, 'a');
    if (res)
        return res;

    res = skipLine(recordReader);
    if (res)
        return res;

    skipWhitespaces(recordReader);
    while (!skipChar(recordReader, 's'))
    {
        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read seq id
        clear(buffer);
        res = readUntilOneOf(buffer, recordReader, '.', ' ', '\t');
        if (res)
            return res;

        seqId = buffer;

        // skip line if seqId is not in set of specified seqIds
        if (seqIds.size() != 0 && seqIds.count(seqId) == 0)
        {
            skipLine(recordReader);
            continue;
        }

        res = skipChar(recordReader, '.');
        if (!res)
        {
            // read chromosome id
            clear(buffer);
            res = readUntilWhitespace(buffer, recordReader);
            if (res)
                return res;

            chromosomeId = buffer;
        }

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read start position
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        startPos = lexicalCast<unsigned>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read length
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        len = lexicalCast<unsigned>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // skip orientation
        clear(buffer);
        res = readNChars(buffer, recordReader, 1);
        if (res)
            return res;

        if (buffer  == '+')
            orientation = true;
        else if (buffer == '-')
            orientation = false;
        else
            return 1;

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // skip seq length
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        seqLen = lexicalCast<unsigned>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // skip gapped seq
        skipLine(recordReader);

        if (!orientation)
            startPos = seqLen - (startPos + len);

        if (len > 0)
            appendValue(idToRowsMap[seqId], TRow(lineId, chromosomeId, startPos, startPos + len, orientation));

        skipWhitespaces(recordReader);
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function parseAlignment()
// ----------------------------------------------------------------------------

template <typename TSize, typename TFile, typename TTag>
int
parseAlignment(String<std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > > & idToRowsMaps,
               TFile & file,
               std::set<CharString> const & seqIds,
               bool verbose,
               TTag const tag)
{
    typedef AlignmentBlockRow<TSize, TSize> TRow;
    typedef std::map<CharString, String<TRow> > TMap;

    RecordReader<TFile, SinglePass<> > recordReader(file);
    int res = 0;

    res = skipHeader(recordReader, tag);
    if (res)
        return res;

    unsigned lineId = 0;

    while (!atEnd(recordReader))
    {
        TMap idToRows;
        res = readRecord(idToRows, recordReader, lineId, seqIds, tag);
        if (res)
            return res;

        appendValue(idToRowsMaps, idToRows);

        skipWhitespaces(recordReader);
        while (!atEnd(recordReader) && value(recordReader) == '#')
            res = skipLine(recordReader);
        
        ++lineId;
    }
    if (!atEnd(recordReader))
        return 1;

    if (verbose)
    {
        std::cout << length(idToRowsMaps) << " alignment block";
        if (length(idToRowsMaps) != 1)
            std::cout << "s";
        std::cout << " loaded." << std::endl;
    }

    return 0;
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_
