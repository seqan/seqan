#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/store.h>

using namespace seqan;

int main(int argc, char **argv) {
    // Check arguments.
    if (argc != 3) {
        std::cerr << "Usage:  sam_info GENOME.FASTA FILE.SAM" << std::endl;
        return 1;
    }

    // =======================================================================
    // Load files.
    // =======================================================================
    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;

    // Load contigs.
    std::cerr << "Reading FASTA contigs sequence file " << argv[1] << " ..." << std::endl;
    if (!loadContigs(fragments, argv[1])) {
        std::cerr << "Could not read contigs." << std::endl;
        return 1;
    }

    // Load SAM file.
    std::cerr << "Reading SAM file file " << argv[2] << " ..." << std::endl;
    {
        std::fstream fstrm(argv[2], std::ios_base::in | std::ios_base::binary);
        if (! fstrm.is_open()) {
            std::cerr << "Could not open SAM file." << std::endl;
            return 1;
        }
        read(fstrm, fragments, SAM());
    }
    std::cout << "length(fragments.alignedReadStore) == " << length(fragments.alignedReadStore) << std::endl;

    // =======================================================================
    // Perform evaluation.
    // =======================================================================
    std::cerr << "Evaluating..." << std::endl;
    performEvaluation(fragments);
    
    return 0;
}
