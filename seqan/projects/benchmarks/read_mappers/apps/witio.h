/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
 ===========================================================================
  Copyright (C) 2007
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
 ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ===========================================================================
  I/O functions for WIT files.
 ===========================================================================*/

#ifndef WITIO_H_
#define WITIO_H_

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/misc/misc_parsing.h>

using namespace seqan;

// Represents a read match to a contig with an error.
//
// WeightedMatch(c, p, d) represents a match of a read at position p
// in contig c with distance d.
struct WeightedMatch
{
    size_t contigId;
    bool isForward;
    size_t pos;
    int distance;
    // Optional begin position, used in wit builder to smooth matches.
    // Is not written out or used in the less than operator.
    size_t beginPos;
  
    WeightedMatch() {}

    WeightedMatch(size_t _contigId, bool _isForward, size_t _pos, int _distance, size_t _beginPos)
            : contigId(_contigId), isForward(_isForward), pos(_pos), distance(_distance), beginPos(_beginPos)
    {}

    // Order lexicographically by (contigId, pos, -distance).
    bool operator<(WeightedMatch const & other) const
    {
        if (contigId < other.contigId) return true;
        if (contigId == other.contigId && isForward > other.isForward) return true;
        if (contigId == other.contigId && isForward == other.isForward &&
            pos < other.pos) return true;
        if (contigId == other.contigId && isForward == other.isForward &&
            pos == other.pos && distance > other.distance) return true;
        return false;
    }

    bool operator==(WeightedMatch const & other) const
    {
        if (contigId != other.contigId) return false;
        if (isForward != other.isForward) return false;
        if (pos != other.pos) return false;
        if (distance != other.distance) return false;
        if (beginPos != other.beginPos) return false;
        return true;
    }
};
    

// Stream output for WeightedMatch objects, for debugging.
template <typename TStream>
TStream &operator<<(TStream &out, const WeightedMatch &m)
{
    out << "(" << m.contigId << ", " << (m.isForward ? "F, " : "R, ") << m.pos << ", " << m.distance << ", " << m.beginPos << ")";
    return out;
}


struct WeightedMatchBeginPosNeqOrContigIdNeq : std::binary_function<WeightedMatch, WeightedMatch, bool> {
    bool operator()(WeightedMatch const & arg1, WeightedMatch const & arg2) {
//         std::cerr << "&arg1 == " << &arg1 << ", &arg2 == " << &arg2 << std::endl;
//         std::cerr << arg1.beginPos << " != " << arg2.beginPos << " == " << (arg1.beginPos != arg2.beginPos) << std::endl;
        return arg1.beginPos != arg2.beginPos || arg1.contigId != arg2.contigId || arg1.isForward != arg2.isForward;
    }
};
    

typedef String<WeightedMatch> TWeightedMatches;

/*
// Write WIT header line ("@HD\tVN:1.0") to stream.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.
template <typename TStream>
TStream & writeWitHeader(TStream &stream) {
  stream << "@HD\tVN:1.0" << std::endl;
  return stream;
}


// Write a WIT comment line ("# %s") to stream.
//
// TStream -- an output stream.
// TString -- the type of the comment string.
//
// stream -- stream to write to.
// str    -- string to write to stream.
template <typename TStream, typename TString>
TStream &writeWitComment(TStream &stream, const TString &str) {
  stream << "# " << str << std::endl;
  return stream;
}


// Write a WIT data line to stream.
//
// TStream -- an output stream.
// TString -- the type of the comment string.
// TScore  -- an integer type for the score.
// TPos    -- an integer type for the position.
//
// stream -- stream to write to.
// readName -- the name of the read.
// score  -- score of the match.
// referenceSequenceName -- name of the sequence
// first -- first entry value
// last  -- last entry value
template <typename TStream, typename TString, typename TScore, typename TPos>
TStream &writeWitRecord(TStream &stream, const TString &readName,
                        TScore score, const TString &referenceSequenceName,
                        TPos first, TPos last) {
  stream << readName << "\t" << score << "\t" << referenceSequenceName
         << "\t" << first << "\t" << last << std::endl;
  return stream;
}


// Read WIT header line ("@HD\tVN:1.0") from stream.
//
// TStream -- type of the stream.
// TChar   -- type of lookahead character.
//
// stream -- Stream to read from.
// c -- Lookahead character.
template <typename TStream, typename TChar>
    void readWitHeader(TStream &stream, TChar &c) {
    CharString tmp;
    // Read "@HD".
    c = _streamGet(stream);
    tmp = _parseReadWordUntilWhitespace(stream, c);
    if (tmp != CharString("@HD"))
        std::cerr << "WARNING: File did not begin with \"@HD\", was: \"" << tmp << "\"" << std::endl;
    // Skip "\t".
    _parseSkipWhitespace(stream, c);
    // Read "VN:1.0".
    tmp = _parseReadWordUntilWhitespace(stream, c);
    if (tmp != CharString("VN:1.0"))
        std::cerr << "WARNING: Version is not \"VN:1.0\"" << std::endl;
    // Skip to and after end of line.
    _parse_skipLine(stream, c);
}

// Represents an entry in a WIT file.
struct WitRecord {
    size_t readId;
    int score;
    size_t contigId;
    size_t first;
    size_t last;

    WitRecord() {
    }
    
    WitRecord(size_t readId_, int score_, size_t contigId_, size_t first_, size_t last_)
        : readId(readId_), score(score_), contigId(contigId_), first(first_), last(last_) {
    }

    // Lexicographic comparison.
    bool operator<(const WitRecord &other) const {
        if (readId < other.readId)
            return true;
        if (readId == other.readId and score < other.score)
            return true;
        if (readId == other.readId and score == other.score and
            contigId < other.contigId)
            return true;
        if (readId == other.readId and score == other.score and
            contigId == other.contigId and first < other.first)
            return true;
        if (readId == other.readId and score == other.score and
            contigId == other.contigId and first == other.first and
            last < other.last)
            return true;
        return false;
    }
};


template <typename TStream>
TStream &operator<<(TStream &stream, const WitRecord &witRecord) {
    stream << "WitRecord(" << witRecord.readId << ", " << witRecord.score
           << ", " << witRecord.contigId << ", " << witRecord.first
           << ", " << witRecord.last << ")";
    return stream;
}


// Read WIT record, skipping comments.
//
// TStream -- type of the stream.
// TString -- the type of the comment string.
// TScore  -- an integer type for the score.
// TPos    -- an integer type for the position.
// TChar   -- lookahead for the parser.
//
// stream -- Stream to read from.
// readName -- the name of the read.
// score  -- score of the match.
// referenceSequenceName -- name of the sequence
// first -- first entry value
// last  -- last entry value
// c -- lookahead for the parser.
//
// Returns true iff the record could be successfully read from the file.
template <typename TStream, typename TString, typename TScore, typename TPos, typename TChar>
bool readWitRecord(TStream &stream, TString &readName,
                   TScore &score, TString &referenceSequenceName,
                   TPos &first, TPos &last, TChar &c) {
    String<char> tmp;

    // Maybe skip comments.
    while (not _streamEOF(stream) and c == '#')
        _parse_skipLine(stream, c);
    if (_streamEOF(stream))
        return false;

    // Read line.
    _parseReadIdentifier(stream, readName, c);
    _parseSkipWhitespace(stream, c);
    score = _parseReadNumber(stream, c);
    _parseSkipWhitespace(stream, c);
    _parseReadIdentifier(stream, referenceSequenceName, c);
    _parseSkipWhitespace(stream, c);
    first = _parseReadNumber(stream, c);
    _parseSkipWhitespace(stream, c);
    last = _parseReadNumber(stream, c);
//     std::cout << "readName, score, first, last, " << readName << ", " << score << ", " << first << ", " << last << std::endl;
    
    // Skip to and after end of line.
    _parse_skipLine(stream, c);
    return true;
}
*/
#endif  // WITIO_H_
