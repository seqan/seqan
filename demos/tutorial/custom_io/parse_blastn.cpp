//![includes]
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;
//![includes]

//![tags]
struct BlastnTab_;
typedef Tag<BlastnTab_> BlastnTab;

struct BlastnTabComment_;
typedef Tag<BlastnTabComment_> BlastnTabComment;

struct BlastnTabAlignment_;
typedef Tag<BlastnTabAlignment_> BlastnTabAlignment;
//![tags]

//![record]
struct BlastnTabAlignmentRecord
{
    CharString queryName;
    CharString subjectName;
    double identity;
    unsigned alignmentLength;
    unsigned mismatches;
    unsigned gapOpens;
    unsigned queryBegin;
    unsigned queryEnd;
    unsigned subjectBegin;
    unsigned subjectEnd;
    double eValue;
    double bitScore;

    BlastnTabAlignmentRecord() :
        identity(0), alignmentLength(0), mismatches(0), gapOpens(0),
        queryBegin(0), queryEnd(0), subjectBegin(0), subjectEnd(0),
        eValue(0), bitScore(0)
    {}
};

template <typename TWriter>
inline void
writeRecord(TWriter & writer, BlastnTabAlignmentRecord & record)
{
    write(writer, "query name: ");
    write(writer, record.queryName);
    write(writer, "\nsubject name: ");
    write(writer, record.subjectName);
    write(writer, "\nidentity: ");
    appendNumber(writer, record.identity);
    write(writer, "\nalignment length: ");
    appendNumber(writer, record.alignmentLength);
    write(writer, "\nmismatches: ");
    appendNumber(writer, record.mismatches);
    write(writer, "\ngap opens: ");
    appendNumber(writer, record.gapOpens);
    write(writer, "\nquery begin: ");
    appendNumber(writer, record.queryBegin);
    write(writer, "\nquery end: ");
    appendNumber(writer, record.queryEnd);
    write(writer, "\nsubject begin: ");
    appendNumber(writer, record.subjectBegin);
    write(writer, "\nsubject end: ");
    appendNumber(writer, record.subjectEnd);
    write(writer, "\nevalue: ");
    appendNumber(writer, record.eValue);
    write(writer, "\nbit score: ");
    appendNumber(writer, record.bitScore);
    write(writer, "\n\n");
}
//![record]

//![next-is]
template <typename TReader>
inline bool
nextIs(TReader & reader, BlastnTabComment const & /*tag*/)
{
    return !atEnd(reader) && value(reader) == '#';
}

template <typename TReader>
inline bool
nextIs(TReader & reader, BlastnTabAlignment const & /*tag*/)
{
    return !atEnd(reader) && value(reader) != '#';
}
//![next-is]

//![read-record]
template <typename TCharSequence, typename TReader>
inline void
readRecord(TCharSequence & buffer, TReader const & reader, BlastnTabComment const & /*tag*/)
{
    SEQAN_ASSERT(nextIs(reader, BlastnTabComment()));
    clear(buffer);
    readLine(buffer, reader);
}

template <typename TReader>
inline void
readRecord(BlastnTabAlignmentRecord & record, CharString & buffer, TReader & reader, BlastnTabAlignment const & /*tag*/)
{
    SEQAN_ASSERT(nextIs(reader, BlastnTabAlignment()));

    // Read query name.
    clear(record.queryName);
    readUntil(record.queryName, reader, IsTab());
    skipOne(reader, IsTab());

    // Read subject name.
    clear(record.subjectName);
    readUntil(record.subjectName, reader, IsTab());
    skipOne(reader, IsTab());

    // Read identity.
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.identity, buffer);
    skipOne(reader, IsTab());

    // Read alignment length.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.alignmentLength, buffer);
    skipOne(reader, IsTab());

    // Read mismatches.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.mismatches, buffer);
    skipOne(reader, IsTab());

    // Read gap opens.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.gapOpens, buffer);
    skipOne(reader, IsTab());

    // Read query begin.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.queryBegin, buffer);
    skipOne(reader, IsTab());

    // Read query end.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.queryEnd, buffer);
    skipOne(reader, IsTab());

    // Read subject begin.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.subjectBegin, buffer);
    skipOne(reader, IsTab());

    // Read subject end.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.subjectEnd, buffer);
    skipOne(reader, IsTab());

    // Read evalue.
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.eValue, buffer);
    skipOne(reader, IsTab());

    // Read bit score, up to end of the line.
    clear(buffer);
    readLine(buffer, reader);
    lexicalCast(record.bitScore, buffer);
}
//![read-record]

//![skip-record]
template <typename TReader>
inline void
skipRecord(TReader & reader, BlastnTabComment const & /*tag*/)
{
    skipLine(reader);
}

template <typename TReader>
inline void
skipRecord(TReader & reader, BlastnTabAlignment const & /*tag*/)
{
    skipLine(reader);
}
//![skip-record]

//![batch-read]
template <typename TBlastnTabRecords, typename TReader>
inline void
readRecords(TBlastnTabRecords & records, TReader & reader, BlastnTab const & /*tag*/)
{
    BlastnTabAlignmentRecord record;
    CharString buffer;
    while (!atEnd(reader))
    {
        if (nextIs(reader, BlastnTabComment()))
        {
            skipRecord(reader, BlastnTabComment());
        }
        else if (nextIs(reader, BlastnTabAlignment()))
        {
            readRecord(record, buffer, reader, BlastnTabAlignment());
            appendValue(records, record);
        }
        else
        {
            SEQAN_THROW(ParseError("Unknown BLAST record type"));
        }
    }
}
//![batch-read]

//![main]
int main(int argc, char const * argv[])
{
    // Process command line arguments, open file.
    if (argc != 2)
    {
        std::cerr << "Incorrect argument count!" << std::endl;
        std::cerr << "USAGE: tutorial_parse_blastn INPUT.txt" << std::endl;
        return 1;
    }
    std::ifstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
    {
        std::cerr << "Could not open file " << argv[1] << std::endl;
        return 1;
    }

    // Read file.
    DirectionIterator<std::ifstream, Input>::Type reader = directionIterator(stream, Input());
    String<BlastnTabAlignmentRecord> records;
    readRecords(records, reader, BlastnTab());

    // Write read records.
    for (unsigned i = 0; i < length(records); ++i)
        writeRecord(std::cout, records[i]);

    return 0;
}
//![main]
