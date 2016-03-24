#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    CharString seqFileName = getAbsolutePath("demos/tutorial/sequence_io/example.fq");

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(seqFileName)))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    StringSet<CharString> quals;

    try
    {
        readRecords(ids, seqs, quals, seqFileIn);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(ids); ++i)
        std::cout << ids[i] << '\t' << seqs[i] << "\n+qual:\t" << quals[i] << '\n';

    return 0;
}
