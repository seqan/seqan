#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

    CharString id;
    Dna5String seq;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    try
    {
        readRecord(id, seq, seqFileIn);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::cout << id << '\t' << seq << '\n';

    return 0;
}
