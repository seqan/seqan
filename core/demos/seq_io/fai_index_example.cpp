#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " SEQUENCE.fasta\n";
        return 1;
    }
    
    FaiIndex faiIndex;

    // Try to read the FAI index.
    bool readSuccess = (read(faiIndex, argv[1]) == 0);
    if (!readSuccess)
        std::cerr << "Could not read the FAI index.  Not fatal, we can just build it.\n";

    // Try to build the FAI index (in memory) if reading was unsuccessful.  If
    // building into memory succeeded, we try to write it out.
    if (!readSuccess)
    {
        if (build(faiIndex, argv[1]) != 0)
        {
            std::cerr << "FATAL: Could not build FAI index.\n";
            return 1;
        }

        if (write(faiIndex) != 0)
        {
            std::cerr << "FATAL: Could not write out FAI index after building.\n";
            return 1;
        }
    }

    // Now, read the first 1000 characters of chr1.
    unsigned idx = 0;
    if (!getIdByName(faiIndex, "chr1", idx))
    {
        std::cerr << "FATAL: chr1 not found in FAI index.\n";
        return 1;
    }
    CharString seq;
    if (readRegion(seq, faiIndex, idx, 0, 1000) != 0)
    {
        std::cerr << "FATAL: Problem reading FASTA file through FAI index.\n";
        return 1;
    }

    // Now print the first 1000 characters we just read.
    std::cerr << "chr1:1-1000 = " << seq << "\n";
    
    return 0;
}
