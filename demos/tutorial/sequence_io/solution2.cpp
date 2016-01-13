#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    CharString seqFileName = getAbsolutePath("demos/tutorial/sequence_io/example.fa");
    CharString id;
    Dna5String seq;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(seqFileName)))
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
