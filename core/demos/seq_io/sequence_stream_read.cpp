#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

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
        append(path, "/core/tests/seq_io/test_dna.fq");
    }

    // TODO(esiragusa): define DefaultOpenMode<> on VirtualStream
    SequenceFile<Input> file(toCString(path), OPEN_RDONLY);

    CharString id;
    CharString seq;
//    Dna5String seq;
    CharString qual;

    double start, finish;

    start = sysTime();

    typename Size<CharString>::Type records = 0;
    typename Size<CharString>::Type bases = 0;

    while (!atEnd(file))
    {
        try
        {
            read(file, id, seq, qual);
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
