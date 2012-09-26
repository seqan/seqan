/*
  Definition of the weighted interval type with supporting functions,
  i.e. for I/O.
 */

// TODO(holtgrew): Should maybe better be called LabeledInterval?

#ifndef WIT_BUILDER_INTERVALS_H_
#define WIT_BUILDER_INTERVALS_H_

#include <iostream>

#include <seqan/basic.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/sequence.h>

using namespace seqan;

// The columns in a WIT file.
#define WIT_COLUMN_NAMES "QNAME\tSCORE\tRNAME\tRDIR\tSPOS\tEPOS"


// A simple record class that stores the information of a weighted
// interval for WIT files.
struct WitRecord {
    // Name of the read in question.
    CharString readName;

    enum {
      FLAG_PAIRED = 0x01,
      FLAG_FIRST_MATE = 0x40,
      FLAG_SECOND_MATE = 0x80
    };

    // Flags, 0x01 - paired, 0x40 - first read in pair, 0x80 - second read in
    // pair.
    int flags;

    // Id of the read, possibly not set.
    size_t readId;

    // Distance associated with the interval.
    int distance;

    // Name of the contig the interval is defined on.
    CharString contigName;

    // Id of the contig, possibly not set.
    size_t contigId;

    // true iff the interval is on the forward strand.
    bool isForward;

    // First position on the interval.
    size_t firstPos;

    // Last position on the interval.
    size_t lastPos;

    // Default constructor.
    WitRecord() {}

    // Complete constructor for all properties.
    WitRecord(CharString const & _readName, int const & _flags,
              int const & _distance,
              CharString const & _contigName, bool const & _isForward,
              size_t const & _firstPos, size_t const & _lastPos)
            : readName(_readName), flags(_flags), readId(0), distance(_distance),
              contigName(_contigName), contigId(0), isForward(_isForward),
              firstPos(_firstPos), lastPos(_lastPos) {}

    // Lexicographic comparison.
    bool operator<(WitRecord const & other) const {
        if (readId < other.readId)
            return true;
        if (readId == other.readId && distance < other.distance)
            return true;
        if (readId == other.readId && distance == other.distance &&
            contigId < other.contigId)
            return true;
        if (readId == other.readId && distance == other.distance &&
            contigId == other.contigId && firstPos < other.firstPos)
            return true;
        if (readId == other.readId && distance == other.distance &&
            contigId == other.contigId && firstPos == other.firstPos &&
            lastPos < other.lastPos)
            return true;
        return false;
    }
};


struct WitRecord_Lt_ContigIdReadIdLastPos {
    bool operator()(WitRecord const & a, WitRecord const & b) {
        if (a.contigId < b.contigId)
            return true;
        if (a.contigId == b.contigId && a.readId < b.readId)
            return true;
        if (a.contigId == b.contigId && a.readId == b.readId && a.lastPos < b.lastPos)
            return true;
        return false;
    }
};


// Output-stream operator for WitRecord objects.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.
template <typename TStream>
TStream & operator<<(TStream & stream, WitRecord const & record) {
    stream << record.readName;
    if (record.flags & WitRecord::FLAG_PAIRED) {
      if (record.flags & WitRecord::FLAG_FIRST_MATE) {
        stream << "/0";
      } else if (record.flags & WitRecord::FLAG_SECOND_MATE) {
        stream << "/1";
      } else {
        stream << "/?";
      }
    }
    stream << "\t"
           << record.distance << "\t"
           << record.contigName << "\t"
           << (record.isForward ? "F" : "R") << "\t"
           << record.firstPos << "\t"
           << record.lastPos << "\t";
    return stream;
}


// Write WIT header line ("@WIT\tVN:1.0") to stream.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.
template <typename TStream>
TStream &writeWitHeader(TStream & stream) {
    stream << "@WIT\tVN:1.0" << std::endl;
    stream << "@MATES\tSEP:/\tTYPE:01" << std::endl;
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
TStream &writeWitComment(TStream &stream, TString const & str) {
    stream << "# " << str << std::endl;
    return stream;
}


// Write a WIT data line to stream.
//
// TStream -- an output stream.
//
// stream -- stream to write to.
// record -- the WitRecord to write out.
template <typename TStream, typename TString, typename TScore, typename TPos>
TStream &writeWitRecord(TStream & stream, WitRecord const & record) {
    return stream << record << std::endl;
}


// Read WIT header line ("@WIT\tVN:1.0") from stream.
//
// TStream -- type of the stream.
// TChar   -- type of lookahead character.
//
// stream -- Stream to read from.
// c -- Lookahead character.
template <typename TStream, typename TChar>
void readWitHeader(TStream &stream, TChar &c) {
    CharString tmp;
    // Read "@WIT".
    c = _streamGet(stream);
    tmp = _parseReadWordUntilWhitespace(stream, c);
    if (tmp != CharString("@WIT"))
        std::cerr << "WARNING: File did not begin with \"@WIT\", was: \"" << tmp << "\"" << std::endl;
    // Skip "\t".
    _parseSkipWhitespace(stream, c);
    // Read "VN:1.0".
    tmp = _parseReadWordUntilWhitespace(stream, c);
    if (tmp != CharString("VN:1.0"))
        std::cerr << "WARNING: Version is not \"VN:1.0\"" << std::endl;
    // Skip to and after end of line.
    _parseSkipLine(stream, c);
    // Maybe read/skip additional header lines.
    while (c == '@')
      _parseSkipLine(stream, c);
}


// Read WIT record, skipping comments.
//
// TStream -- type of the stream.
// TChar   -- lookahead for the parser.
//
// stream -- Stream to read from.
// record -- the Witrecord to read.
// c -- lookahead for the parser.
//
// Note that the WIT file contains 1-based entries but we subtract 1
// at this location.
//
// Returns true iff the record could be successfully read from the file.
template <typename TStream, typename TChar>
bool readWitRecord(TStream & stream, WitRecord & record, TChar & c) {
    String<char> tmp;

    // Maybe skip comments.
    while (not _streamEOF(stream) && c == '#')
        _parseSkipLine(stream, c);
    if (_streamEOF(stream))
        return false;

    // Read line.
    _parseReadIdentifier(stream, record.readName, c);
    _parseReadIdentifier(stream, record.readName, c);
    _parseSkipWhitespace(stream, c);
    record.distance = _parseReadNumber(stream, c);
    _parseSkipWhitespace(stream, c);
    _parseReadIdentifier(stream, record.contigName, c);
    _parseSkipWhitespace(stream, c);
    record.isForward = (_parseReadChar(stream, c) == 'F');
    _parseSkipWhitespace(stream, c);
    record.firstPos = _parseReadNumber(stream, c);
    _parseSkipWhitespace(stream, c);
    record.lastPos = _parseReadNumber(stream, c);
    
    // Skip to and after end of line.
    _parseSkipLine(stream, c);
    return true;
}

#endif  // WIT_BUILDER_INTERVALS_H_

