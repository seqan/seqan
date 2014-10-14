#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

    seqan::SeqFileOut seqFileOut;
    if (!open(seqFileOut, argv[1]))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    try
    {
        seqan::StringSet<seqan::CharString> ids;
        appendValue(ids, "seq1");
        appendValue(ids, "seq2");
        seqan::StringSet<seqan::Dna5String> seqs;
        appendValue(seqs, "CGAT");
        appendValue(seqs, "TTTT");

        writeRecords(seqFileOut, ids, seqs);
    }
    catch (std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
