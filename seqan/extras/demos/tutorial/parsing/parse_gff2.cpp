// FRAGMENT(includes)
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

// FRAGMENT(tags-structs)
struct Gff2_;
typedef Tag<Gff2_> Gff2;

struct Gff2Record
{
    CharString seqName;
    CharString source;
    CharString feature;
    unsigned start;
    unsigned end;
    bool hasScore;
    double score;
    char strand;
    unsigned frame;
    CharString attributes;
    CharString comments;

    Gff2Record() : start(0), end(0), hasScore(false), score(0), strand('.'), frame(0)
    {}
};

void clear(Gff2Record & record)
{
    clear(record.seqName);
    clear(record.source);
    clear(record.feature);
    clear(record.attributes);
    clear(record.comments);
}

// FRAGMENT(read-record)
template <typename TStream, typename TPass>
int readRecord(Gff2Record & record, RecordReader<TStream, TPass> & reader, Gff2 const & /*tag*/)
{
    clear(record);
    CharString buffer;
    char c = '\0';
    int res = 0;

    // GFF2 records look like this:
    //
    // <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

    // <seqname>
    res = readUntilChar(record.seqName, reader, '\t');
    if (res != 0)
        return res;
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <source>
    res = readUntilChar(record.source, reader, '\t');
    if (res != 0)
        return res;
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <feature>
    res = readUntilChar(record.feature, reader, '\t');
    if (res != 0)
        return res;
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <start>
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.end, buffer))
        return 1;  // Could not cast!
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <end>
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (!lexicalCast2<unsigned>(record.end, buffer))
        return 1;  // Could not cast!
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <score>
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    record.hasScore = (buffer != '.');
    if (record.hasScore && !lexicalCast2<double>(record.score, buffer))
        return 1;  // Could not cast!
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <strand>
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (length(buffer) != 1u)
        return 1;  // More than one char or none.
    c = front(buffer);
    if (c != '.' && c != '+' && c != '-')
        return 1;  // Invalid strand.
    record.strand = c;
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <frame>
    clear(buffer);
    res = readUntilChar(buffer, reader, '\t');
    if (res != 0)
        return res;
    if (length(buffer) != 1u)
        return 1;  // More than one char or none.
    c = front(buffer);
    if (c != '.' && c != '0' && c != '1' && c != '2')
        return 1;  // Invalid frame.
    record.frame = c;
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;

    // <attributes>
    res = readUntilTabOrLineBreak(record.attributes, reader);
    if (res != 0 && res != EOF_BEFORE_SUCCESS)
        return res;
    if (atEnd(reader))
        return 0;
    if (value(reader) == '\t')
    {
        if (goNext(reader))
            return EOF_BEFORE_SUCCESS;
    }

    // <comment>
    res = readLine(record.seqName, reader);
    if (res != 0 && res != EOF_BEFORE_SUCCESS)
        return res;
    return 0;
}

// FRAGMENT(read-batch)
template <typename TGff2Records, typename TStream, typename TPass>
int read2(TGff2Records & records, RecordReader<TStream, TPass> & reader, Gff2 const & /*tag*/)
{
    Gff2Record record;
    while (!atEnd(reader))
    {
        clear(record);
        int res = readRecord(record, reader, Gff2());
        if (res != 0)
            return res;
        appendValue(records, record);
    }
    return 0;
}

// FRAGMENT(main)
int main(int argc, char const ** argv)
{
    // Handle command line arguments, open files.
    if (argc != 2)
        return 1;
    std::fstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
        return 1;

    // Read file.
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    String<Gff2Record> gffRecords;
    int res = read2(gffRecords, reader, Gff2());
    if (res != 0)
        return res;

    // Write out some of the data to stdout.
    for (unsigned i = 0; i < length(gffRecords); ++i)
        std::cout << gffRecords[i].seqName << "\t" << gffRecords[i].strand << "\t" << gffRecords[i].start << "\t"
                  << gffRecords[i].end << std::endl;
    
    return 0;
}
