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
        seqan::CharString id = "seq1";
        seqan::Dna5String seq = "CGAT";
        writeRecord(seqFileOut, id, seq);

        id = "seq1";
        seq = "TTTT";
        writeRecord(seqFileOut, id, seq);
    }
    catch (std::runtime_error &e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
