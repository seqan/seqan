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
        append(path, "/tests/seq_io/test_dna.fq");
    }

    SeqFileIn file;
    if (!open(file, toCString(path)))
    {
        std::cerr << "Can't open the file." << std::endl;
        return 1;
    }

    CharString id;
    DnaString seq;
    CharString qual;

    Size<CharString>::Type records = 0;
    Size<CharString>::Type bases = 0;

    double start, finish;

    start = sysTime();

    while (!atEnd(file))
    {
        try
        {
            readRecord(id, seq, qual, file);
        }
        catch (UnexpectedEnd &)
        {
            break;
        }
        catch (ParseError & e)
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
