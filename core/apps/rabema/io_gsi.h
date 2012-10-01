// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// I/O functionality for GSI Records.
// ==========================================================================

#ifndef SEQAN_CORE_APPS_RABEMA_IO_GSI_H_
#define SEQAN_CORE_APPS_RABEMA_IO_GSI_H_

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// The columns in a WIT file.
#define GSI_COLUMN_NAMES "QNAME\tSCORE\tRNAME\tRDIR\tSPOS\tEPOS"

struct Gsi_;
typedef seqan::Tag<Gsi_> Gsi;

// ---------------------------------------------------------------------------
// Class GsiHeader
// ---------------------------------------------------------------------------

// Dummy only at the moment.
class GsiHeader
{
public:
    GsiHeader() {}
};

// ---------------------------------------------------------------------------
// Class GsiRecord
// ---------------------------------------------------------------------------

// A simple record class that stores the information of a weighted interval for GSI files.

struct GsiRecord
{
    // Name of the read in question.
    CharString readName;

    enum
    {
        FLAG_PAIRED = 0x01,
        FLAG_FIRST_MATE = 0x40,
        FLAG_SECOND_MATE = 0x80
    };

    // Flags, 0x01 - paired, 0x40 - first read in pair, 0x80 - second read in pair.
    int flags;

    // Id of the read, possibly not set.
    size_t readId;

    // Original distance of the interval before lowering.
    int originalDistance;

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
    GsiRecord() :
        flags(0), readId(0), originalDistance(0), distance(0), contigId(0), isForward(0), firstPos(0), lastPos(0)
    {}

    // Complete constructor for all properties.
    GsiRecord(CharString const & _readName,
              int const & _flags,
              int const & _distance,
              CharString const & _contigName,
              bool const & _isForward,
              size_t const & _firstPos,
              size_t const & _lastPos) :
        readName(_readName), flags(_flags), readId(0), originalDistance(_distance), distance(_distance),
        contigName(_contigName), contigId(0), isForward(_isForward), firstPos(_firstPos), lastPos(_lastPos)
    {}

    // Lexicographic comparison.
    bool operator<(GsiRecord const & other) const
    {
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

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

void clear(GsiRecord & record)
{
    clear(record.readName);
    record.flags = 0;
    record.readId = 0;
    record.originalDistance = 0;
    record.distance = 0;
    clear(record.contigName);
    record.contigId = 0;
    record.isForward = true;
    record.firstPos = 0;
    record.lastPos = 0;
}

// ---------------------------------------------------------------------------
// Function operator<<()                                           [GsiRecord]
// ---------------------------------------------------------------------------

// Output-stream operator for WitRecord objects.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.

template <typename TStream>
TStream & operator<<(TStream & stream, GsiRecord const & record)
{
    stream << record.readName;
    if (record.flags & GsiRecord::FLAG_PAIRED)
    {
        if (record.flags & GsiRecord::FLAG_FIRST_MATE)
            stream << "/0";
        else if (record.flags & GsiRecord::FLAG_SECOND_MATE)
            stream << "/1";
        else
            stream << "/?";
    }
    stream << "\t"
    << record.distance << "\t"
    << record.contigName << "\t"
    << (record.isForward ? "F" : "R") << "\t"
    << record.firstPos << "\t"
    << record.lastPos;
    return stream;
}

// ---------------------------------------------------------------------------
// Function readRecord()                                           [GsiHeader]
// ---------------------------------------------------------------------------

// Read GSI header line ("@GSI\tVN:1.1") from record reader.

template <typename TStream, typename TSpec>
int readRecord(GsiHeader & header, RecordReader<TStream, TSpec> & reader, Gsi const & /*tag*/)
{
    (void) header;

    CharString tmp;
    // Read "@GSI".
    if (readUntilTabOrLineBreak(tmp, reader) != 0)
        return 1;  // Could not read header.

    if (tmp != "@GSI")
        std::cerr << "WARNING: File did not begin with \"@GSI\", was: \"" << tmp << "\"" << std::endl;
    // Skip "\t".
    if (skipChar(reader, '\t') != 0)
        return 1;  // Next char was not TAB.

    // Read "VN:1.1".
    clear(tmp);
    if (readUntilTabOrLineBreak(tmp, reader) != 0)
        return 1;  // Could not read version.

    if (tmp != "VN:1.1")
        std::cerr << "WARNING: Version is not \"VN:1.1\", was: \"" << tmp << "\"" << std::endl;
    // Skip to and after end of line.
    skipLine(reader);
    // Maybe read/skip additional header lines.
    while (!atEnd(reader) && value(reader) == '@')
        if (skipLine(reader) != 0)
            return 1;

    // Skipy any trailing comment lines.
    while (!atEnd(reader) && value(reader) == '#')
        if (skipLine(reader) != 0)
            return 1;

    return 0;
}

// ---------------------------------------------------------------------------
// Function readRecord()                                           [GsiHeader]
// ---------------------------------------------------------------------------

// Read GSI record, skipping comments.
//
// TStream -- type of the stream.
// TChar   -- lookahead for the parser.
//
// stream -- Stream to read from.
// record -- the GsiRecord to read.
// c -- lookahead for the parser.
//
// Note that the GSI file contains 1-based entries but we subtract 1
// at this location.
//
// Returns true iff the record could be successfully read from the file.

template <typename TStream, typename TSpec>
int readRecord(GsiRecord & record, RecordReader<TStream, TSpec> & reader, Gsi const & /*tag*/)
{
    CharString buffer;

    // No more records in file.
    if (atEnd(reader))
        return 1;

    // Read read name.
    clear(record);
    if (readUntilTabOrLineBreak(record.readName, reader) != 0)
        return 1;

    if (value(reader) != '\t')
        return 1;

    skipChar(reader, '\t');

    // Interpret trailing characters for mate-pair identifier in read name.
    if (length(record.readName) >= 2u && record.readName[length(record.readName) - 2] == '/')
    {
        char c = back(record.readName);
        if (c == '0')
            record.flags = GsiRecord::FLAG_PAIRED | GsiRecord::FLAG_FIRST_MATE;
        else if (c == '1')
            record.flags = GsiRecord::FLAG_PAIRED | GsiRecord::FLAG_SECOND_MATE;
        else
            return 1;  // Could not interpret trailing mate indicator.

        resize(record.readName, length(record.readName) - 2);
    }

    // Read distance.
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;

    if (!lexicalCast2(record.distance, buffer))
        return 1;  // Could not convert distance.

    record.originalDistance = record.distance;
    if (value(reader) != '\t')
        return 1;

    skipChar(reader, '\t');

    // Read contig name.
    if (readUntilTabOrLineBreak(record.contigName, reader) != 0)
        return 1;

    if (value(reader) != '\t')
        return 1;

    skipChar(reader, '\t');

    // Read 'F'/'R'.
    clear(buffer);
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;

    if (buffer != "F" && buffer != "R")
        return 1;

    record.isForward = (buffer[0] == 'F');
    if (value(reader) != '\t')
        return 1;

    skipChar(reader, '\t');

    // Read first pos.
    clear(buffer);
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;

    if (!lexicalCast2(record.firstPos, buffer))
        return 1;

    if (value(reader) != '\t')
        return 1;

    skipChar(reader, '\t');

    // Read last pos.
    clear(buffer);
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;

    if (!lexicalCast2(record.lastPos, buffer))
        return 1;

    // We only allow '\t' here because output buggy.
    if (value(reader) != '\t' && value(reader) != '\r' && value(reader) != '\n')
        return 1;

    if (skipLine(reader) != 0)
        return 1;  // Skip line.

    // Skipy any trailing comment lines.
    while (!atEnd(reader) && value(reader) == '#')
        if (skipLine(reader) != 0)
            return 1;

    return 0;
}

// ---------------------------------------------------------------------------
// Function writeRecord()                                      [GsiHeader,Gsi]
// ---------------------------------------------------------------------------

// Write GSI header line ("@WIT\tVN:1.1") to stream.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.

template <typename TStream>
int writeRecord(TStream & stream, GsiHeader const & /*header*/, Gsi const & /*tag*/)
{
    if (streamPut(stream, "@GSI\tVN:1.1\n") != 0)
        return 1;
    if (streamPut(stream, "@MATES\tSEP:/\tTYPE:01\n") != 0)
        return 1;
    return 0;
}

// ---------------------------------------------------------------------------
// Function writeRecord()                                        [comment,Gsi]
// ---------------------------------------------------------------------------

// Write a GSI comment line ("# %s") to stream.
//
// TStream -- an output stream.
// TString -- the type of the comment string.
//
// stream -- stream to write to.
// str    -- string to write to stream.

template <typename TStream>
int writeRecord(TStream & stream, CharString const & str, Gsi const & /*tag*/)
{
    if (streamPut(stream, '#') != 0)
        return 1;
    if (streamPut(stream, str) != 0)
        return 1;
    if (streamPut(stream, '\n') != 0)
        return 1;
    return 0;
}

// ---------------------------------------------------------------------------
// Function writeRecord()
// ---------------------------------------------------------------------------

// Write a GSI data line to stream.
//
// TStream -- an output stream.
//
// stream -- stream to write to.
// record -- the GsiRecord to write out.

template <typename TStream>
int writeRecord(TStream & stream, GsiRecord const & record, Gsi const & /*tag*/)
{
    if (streamPut(stream, record.readName) != 0)
        return 1;
    if (record.flags & GsiRecord::FLAG_PAIRED)
    {
        if (record.flags & GsiRecord::FLAG_FIRST_MATE)
        {
            if (streamPut(stream, "/0") != 0)
                return 1;
        }
        else if (record.flags & GsiRecord::FLAG_SECOND_MATE)
        {
            if (streamPut(stream, "/1") != 0)
                return 1;
        }
        else
        {
            if (streamPut(stream, "/?") != 0)
                return 1;
        }
    }
    if (streamPut(stream, '\t') != 0)
        return 1;
    if (streamPut(stream, record.distance) != 0)
        return 1;
    if (streamPut(stream, '\t') != 0)
        return 1;
    if (streamPut(stream, record.contigName) != 0)
        return 1;
    if (streamPut(stream, '\t') != 0)
        return 1;
    if (streamPut(stream, (record.isForward ? 'F' : 'R')) != 0)
        return 1;
    if (streamPut(stream, '\t') != 0)
        return 1;
    if (streamPut(stream, record.firstPos) != 0)
        return 1;
    if (streamPut(stream, '\t') != 0)
        return 1;
    if (streamPut(stream, record.lastPos) != 0)
        return 1;
    if (streamPut(stream, '\n') != 0)
        return 1;

    return 0;
}

#endif  // #ifndef SEQAN_CORE_APPS_RABEMA_IO_GSI_H_
