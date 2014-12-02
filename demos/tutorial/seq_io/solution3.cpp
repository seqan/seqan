#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    try
    {
        readRecords(ids, seqs, seqFileIn);
    }
    catch (seqan::IOError const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(ids); ++i)
        std::cout << id[i] << '\t' << seq[i] << '\n';

    return 0;
}
