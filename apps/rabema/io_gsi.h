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

#ifndef SEQAN_APPS_RABEMA_IO_GSI_H_
#define SEQAN_APPS_RABEMA_IO_GSI_H_

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
// Function readHeader()                                           [GsiHeader]
// ---------------------------------------------------------------------------

// Read GSI header line ("@GSI\tVN:1.1") from forward iterator.

template <typename TForwardIter>
void readHeader(GsiHeader & header, TForwardIter & iter, Gsi const & /*tag*/)
{
    (void) header;

    CharString tmp;
    // Read "@GSI".
    readUntil(tmp, iter, OrFunctor<IsTab, IsNewline>());

    if (tmp != "@GSI")
        std::cerr << "WARNING: File did not begin with \"@GSI\", was: \"" << tmp << "\"" << std::endl;
    // Skip "\t".
    skipOne(iter, IsTab());

    // Read "VN:1.1".
    clear(tmp);
    readUntil(tmp, iter, OrFunctor<IsTab, IsNewline>());

    if (tmp != "VN:1.1")
        std::cerr << "WARNING: Version is not \"VN:1.1\", was: \"" << tmp << "\"" << std::endl;
    // Skip to and after end of line.
    skipLine(iter);
    // Maybe read/skip additional header lines.
    while (!atEnd(iter) && *iter == '@')
        skipLine(iter);

    // Skipy any trailing comment lines.
    while (!atEnd(iter) && *iter == '#')
        skipLine(iter);
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

template <typename TForwardIter>
void readRecord(GsiRecord & record, TForwardIter & iter, Gsi const & /*tag*/)
{
    CharString buffer;

    // No more records in file.
    if (atEnd(iter))
        throw seqan::UnexpectedEnd();

    // Read read name.
    clear(record);
    readUntil(record.readName, iter, OrFunctor<IsTab, IsNewline>());

    skipOne(iter, IsTab());

    // Interpret trailing characters for mate-pair identifier in read name.
    if (length(record.readName) >= 2u && record.readName[length(record.readName) - 2] == '/')
    {
        char c = back(record.readName);
        if (c == '0')
            record.flags = GsiRecord::FLAG_PAIRED | GsiRecord::FLAG_FIRST_MATE;
        else if (c == '1')
            record.flags = GsiRecord::FLAG_PAIRED | GsiRecord::FLAG_SECOND_MATE;
        else
            throw seqan::ParseError("Could not interpret trailing mate indicator.");

        resize(record.readName, length(record.readName) - 2);
    }

    // Read distance and perform lexical casting.
    readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
    lexicalCastWithException(record.distance, buffer);

    record.originalDistance = record.distance;

    skipOne(iter, IsTab());

    // Read contig name.
    readUntil(record.contigName, iter, OrFunctor<IsTab, IsNewline>());

    skipOne(iter, IsTab());

    // Read 'F'/'R'.
    resize(buffer, 1, '?');
    readOne(buffer[0], iter, OrFunctor<EqualsChar<'F'>, EqualsChar<'R'> >());
    record.isForward = (buffer[0] == 'F');

    skipOne(iter, IsTab());

    // Read first pos and perform lexical casting.
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
    lexicalCastWithException(record.firstPos, buffer);

    skipOne(iter, IsTab());

    // Read last pos.
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());

    lexicalCastWithException(record.lastPos, buffer);

    skipLine(iter);

    // Skipy any trailing comment lines.
    while (!atEnd(iter) && *iter == '#')
        skipLine(iter);
}

// ---------------------------------------------------------------------------
// Function writeHeader()                                      [GsiHeader,Gsi]
// ---------------------------------------------------------------------------

// Write GSI header line ("@WIT\tVN:1.1") to stream.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.

template <typename TStream>
void writeHeader(TStream & stream, GsiHeader const & /*header*/, Gsi const & /*tag*/)
{
    stream << "@GSI\tVN:1.1\n"
           << "@MATES\tSEP:/\tTYPE:01\n";
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
void writeRecord(TStream & stream, CharString const & str, Gsi const & /*tag*/)
{
    stream << '#' << str << '\n';
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
void writeRecord(TStream & stream, GsiRecord const & record, Gsi const & /*tag*/)
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
    stream << '\t' << record.distance << '\t'
           << record.contigName << '\t'
           << (record.isForward ? 'F' : 'R') << '\t'
           << record.firstPos << '\t'
           << record.lastPos << '\n';
}

#endif  // #ifndef SEQAN_APPS_RABEMA_IO_GSI_H_
