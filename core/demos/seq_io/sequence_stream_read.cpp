#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    CharString path;

    if (argc > 1)
    {
        path = argv[1];
    }
    else
    {
        path = SEQAN_PATH_TO_ROOT();
        append(path, "/core/demos/seq_io/example.fa");
    }

//    typedef String<char, MMap<> >                   TString;
//    typedef Iterator<TString, Rooted>::Type         TIter;
//    TString string(toCString(path), OPEN_RDWR | OPEN_APPEND);
//    std::cout << length(string) << std::endl;
//    TIter it = begin(string, Rooted());

    typedef VirtualStream<char, Input>              TStream;
    typedef Iter<TStream, StreamIterator<Input> >   TIter;
    TStream stream(toCString(path), OPEN_RDONLY);
    TIter it(stream);

    CharString id;
    CharString seq;
//    Dna5String seq;
    CharString qual;

    double start, finish;

    start = sysTime();

    typename Size<CharString>::Type records = 0;
    typename Size<CharString>::Type bases = 0;

    while (!atEnd(it))
    {
        try
        {
            readRecord(id, seq, qual, it, Fastq());
        }
        catch (std::runtime_error & e)
        {
            std::cerr << "Record #" << records + 1 << ": " << e.what() << std::endl;
            continue;
        }

        records += 1;
        bases += length(seq);
    }

    finish = sysTime();

    std::cout << finish - start << " sec" << std::endl;
    std::cout << records << " records" << std::endl;
    std::cout << bases << " bases" << std::endl;

    return 0;
}
