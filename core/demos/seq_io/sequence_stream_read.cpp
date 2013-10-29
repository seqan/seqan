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

    SequenceFile<Input> file(toCString(path));

    CharString id;
    CharString seq;
    CharString qual;

    typename Size<CharString>::Type records = 0;
    typename Size<CharString>::Type bases = 0;

    double start, finish;

    start = sysTime();

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

    close(file);

    std::cout << finish - start << " sec" << std::endl;
    std::cout << records << " records" << std::endl;
    std::cout << bases << " bases" << std::endl;

    return 0;
}
