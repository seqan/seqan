//![includes]
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;
//![includes]

//![tags-structs]
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

    Gff2Record() :
        start(0), end(0), hasScore(false), score(0), strand('.'), frame(0)
    {}
};
//![tags-structs]

//![read-record]
template <typename TReader>
inline void
readRecord(Gff2Record & record, CharString & buffer, TReader & reader, Gff2 const & /*tag*/)
{
    // GFF2 records look like this:
    //
    // <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

    // <seqname>
    clear(record.seqName);
    readUntil(record.seqName, reader, IsTab());
    skipOne(reader, IsTab());

    // <source>
    clear(record.source);
    readUntil(record.source, reader, IsTab());
    skipOne(reader, IsTab());

    // <feature>
    clear(record.feature);
    readUntil(record.feature, reader, IsTab());
    skipOne(reader, IsTab());

    // <start>
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.end, buffer);
    skipOne(reader, IsTab());

    // <end>
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    lexicalCast(record.end, buffer);
    skipOne(reader, IsTab());

    // <score>
    clear(buffer);
    readUntil(buffer, reader, IsTab());
    record.hasScore = (buffer != ".");
    if (record.hasScore)
        lexicalCast(record.score, buffer);
    skipOne(reader, IsTab());

    // <strand>
    readOne(record.strand, reader, OrFunctor<OrFunctor<EqualsChar<'-'>, EqualsChar<'+'> >, EqualsChar<'.'> >());
    skipOne(reader, IsTab());

    // <frame>
    readOne(record.frame, reader, OrFunctor<EqualsChar<'.'>, IsInRange<'0', '2'> >());
    skipOne(reader, IsTab());

    // <attributes>
    clear(record.attributes);
    clear(record.comments);
    readUntil(record.attributes, reader, OrFunctor<IsTab, IsNewline>());
    if (atEnd(reader) || IsNewline() (value(reader)))
    {
        skipLine(reader);
        return;
    }
    skipOne(reader, IsTab());

    // <comment>
    readLine(record.comments, reader);
}
//![read-record]

//![read-batch]
template <typename TGff2Records, typename TReader>
inline void
readRecords(TGff2Records & records, TReader & reader, Gff2 const & /*tag*/)
{
    Gff2Record record;
    CharString buffer;
    while (!atEnd(reader))
    {
        readRecord(record, buffer, reader, Gff2());
        appendValue(records, record);
    }
}
//![read-batch]

//![main]
int main(int argc, char const ** argv)
{
    // Handle command line arguments, open files.
    if (argc != 2)
        return 1;

    std::ifstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
        return 1;

    // Read file.
    DirectionIterator<std::ifstream, Input>::Type reader = directionIterator(stream, Input());
    String<Gff2Record> gffRecords;
    readRecords(gffRecords, reader, Gff2());

    // Write out some of the data to stdout.
    for (unsigned i = 0; i < length(gffRecords); ++i)
        std::cout << gffRecords[i].seqName << "\t" << gffRecords[i].strand << "\t" << gffRecords[i].start << "\t"
                  << gffRecords[i].end << std::endl;

    return 0;
}
//![main]
