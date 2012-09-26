// FRAGMENT(includes)
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>

using namespace seqan;

// FRAGMENT(tags)
struct BlastnTab_;
typedef Tag<BlastnTab_> BlastnTab;

struct BlastnTabComment_;
typedef Tag<BlastnTabComment_> BlastnTabComment;

struct BlastnTabAlignment_;
typedef Tag<BlastnTabAlignment_> BlastnTabAlignment;

// FRAGMENT(record)
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

template <typename TStream>
int write(TStream & stream, BlastnTabAlignmentRecord & record)
{
    streamPut(stream, "query name: ");
    streamPut(stream, record.queryName);
    streamPut(stream, "\nsubject name: ");
    streamPut(stream, record.subjectName);
    streamPut(stream, "\nidentity: ");
    streamPut(stream, record.identity);
    streamPut(stream, "\nalignment length: ");
    streamPut(stream, record.alignmentLength);
    streamPut(stream, "\nmismatches: ");
    streamPut(stream, record.mismatches);
    streamPut(stream, "\ngap opens: ");
    streamPut(stream, record.gapOpens);
    streamPut(stream, "\nquery begin: ");
    streamPut(stream, record.queryBegin);
    streamPut(stream, "\nquery end: ");
    streamPut(stream, record.queryEnd);
    streamPut(stream, "\nsubject begin: ");
    streamPut(stream, record.subjectBegin);
    streamPut(stream, "\nsubject end: ");
    streamPut(stream, record.subjectEnd);
    streamPut(stream, "\nevalue: ");
    streamPut(stream, record.eValue);
    streamPut(stream, "\nbit score: ");
    streamPut(stream, record.bitScore);
    int res = streamPut(stream, "\n\n");
    return res;
}

void clear(BlastnTabAlignmentRecord & record)
{
    clear(record.queryName);
    clear(record.subjectName);
}

// FRAGMENT(next-is)
template <typename TStream, typename TSpec>
inline bool
nextIs(RecordReader<TStream, SinglePass<TSpec> > & reader, BlastnTabComment const & /*tag*/)
{
    return !atEnd(reader) && value(reader) == '#';
}

template <typename TStream, typename TSpec>
inline bool
nextIs(RecordReader<TStream, SinglePass<TSpec> > & reader, BlastnTabAlignment const & /*tag*/)
{
    return !atEnd(reader) && value(reader) != '#';
}

// FRAGMENT(read-record)
template <typename TCharSequence, typename TStream, typename TSpec>
inline int
readRecord(TCharSequence & buffer, RecordReader<TStream, SinglePass<TSpec> > const & reader, BlastnTabComment const & /*tag*/)
{
    SEQAN_ASSERT(nextIs(reader, BlastnTabComment()));
    clear(buffer);
    return readLine(buffer, reader);
}

template <typename TStream, typename TSpec>
inline bool
readRecord(BlastnTabAlignmentRecord & record, RecordReader<TStream, SinglePass<TSpec> > & reader, BlastnTabAlignment const & /*tag*/)
{
    SEQAN_ASSERT(nextIs(reader, BlastnTabAlignment()));
    int res = 0;
    CharString buffer;
    clear(record);

    // Read query name.
    res = readUntilChar(record.queryName, reader, '\t');
    if (res != 0)
        return res;
    goNext(reader);

    // Read subject name.
    res = readUntilChar(record.subjectName, reader, '\t');
    if (res != 0)
        return res;
    goNext(reader);
    
    // Read identity.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<double>(record.identity, buffer))
        return 1;  // Could not cast identity to double.
    goNext(reader);

    // Read alignment length.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.alignmentLength, buffer))
        return 1;  // Could not cast alignment length to unsigned.
    goNext(reader);

    // Read mismatches.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.mismatches, buffer))
        return 1;  // Could not cast mismatches to unsigned.
    goNext(reader);

    // Read gap opens.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.gapOpens, buffer))
        return 1;  // Could not cast gap opens to unsigned.
    goNext(reader);

    // Read query begin.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.queryBegin, buffer))
        return 1;  // Could not cast query begin to unsigned.
    goNext(reader);

    // Read query end.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.queryEnd, buffer))
        return 1;  // Could not cast query end to unsigned.
    goNext(reader);

    // Read subject begin.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.subjectBegin, buffer))
        return 1;  // Could not cast subject begin to unsigned.
    goNext(reader);

    // Read subject end.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.subjectEnd, buffer))
        return 1;  // Could not cast subject end to unsigned.
    goNext(reader);

    // Read evalue.
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<double>(record.eValue, buffer))
        return 1;  // Could not cast evalue to double.
    goNext(reader);

    // Read bit score, up to end of the line.
    clear(buffer);
    res = readLine(buffer, reader);
    if (res != 0)
        return res;
    if (!lexicalCast2<double>(record.bitScore, buffer))
        return 1;  // Could not cast bit score to double.
    
    return 0;
}

// FRAGMENT(skip-record)
template <typename TStream, typename TSpec>
inline int
skipRecord(RecordReader<TStream, SinglePass<TSpec> > & reader, BlastnTabComment const & /*tag*/)
{
    return skipLine(reader);
}

template <typename TStream, typename TPass>
inline int
skipRecord(RecordReader<TStream, TPass> & reader, BlastnTabAlignment const & /*tag*/)
{
    return skipLine(reader);
}

// FRAGMENT(batch-read)
template <typename TBlastnTabRecords, typename TStream, typename TSpec>
int read(TBlastnTabRecords & records, RecordReader<TStream, SinglePass<TSpec> > & reader, BlastnTab const & /*tag*/)
{
    BlastnTabAlignmentRecord record;
    while (!atEnd(reader))
    {
        if (nextIs(reader, BlastnTabComment()))
        {
            skipRecord(reader, BlastnTabComment());
            continue;
        }
        if (!nextIs(reader, BlastnTabAlignment()))
            return 1;
        if (readRecord(record, reader, BlastnTabAlignment()) != 0)
            return 1;
        appendValue(records, record);
    }
    return 0;
}
            

// FRAGMENT(main)
int main(int argc, char const ** argv)
{
    // Process command line arguments, open file.
    if (argc != 2)
    {
        std::cerr << "Incorrect argument count!" << std::endl;
        std::cerr << "USAGE: tutorial_parse_blastn INPUT.txt" << std::endl;
        return 1;
    }
    std::fstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
    {
        std::cerr << "Could not open file " << argv[1] << std::endl;
        return 1;
    }

    // Read file.
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    String<BlastnTabAlignmentRecord> records;
    int res = read(records, reader, BlastnTab());
    if (res != 0)
    {
        std::cerr << "Could not read BLASTN records." << std::endl;
        return 1;
    }

    // Write read records.
    for (unsigned i = 0; i < length(records); ++i)
        write(std::cout, records[i]);

    return 0;
}
