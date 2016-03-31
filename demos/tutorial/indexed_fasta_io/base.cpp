//![build_index1]
//![open_index1]
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    CharString pathToFile = getAbsolutePath("/demos/tutorial/indexed_fasta_io/example.fasta");
    FaiIndex faiIndex;

//![open_index1]
    if (!build(faiIndex, toCString(pathToFile)))
        std::cout << "ERROR: Could not build the index!\n";
//![build_index1]
    clear(faiIndex);
//![build_index2]
    CharString pathToFaiFile = pathToFile;
    append(pathToFaiFile, ".fai");
    if (!build(faiIndex, toCString(pathToFile), toCString(pathToFaiFile)))
        std::cout << "ERROR: Could not build the index!\n";
//![build_index2]

//![save_index]
    if (!save(faiIndex, toCString(pathToFaiFile)))
        std::cout << "ERROR: Could not save the index to file!\n";
//![save_index]

//![open_index1]
    if (!open(faiIndex, toCString(pathToFile)))
        std::cout << "ERROR: Could not load FAI index " << pathToFile << ".fai\n";
//![open_index1]

//![open_index2]
    if (!open(faiIndex, toCString(pathToFile), toCString(pathToFaiFile)))
        std::cout << "ERROR: Could not load FAI index " << pathToFaiFile << "\n";
//![open_index2]

//![idx]
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, "chr1"))
        std::cout << "ERROR: FAI index has no entry for chr1.\n";
//![idx]

//![example_functions]
    unsigned seqLength = sequenceLength(faiIndex, idx);

    // Load first 10 characters of chr1.
    CharString seqChr1Prefix;
    readRegion(seqChr1Prefix, faiIndex, idx, 0, 10);

    // Load all of chr1.
    CharString seqChr1;
    readSequence(seqChr1, faiIndex, idx);
//![example_functions]
    ignoreUnusedVariableWarning(seqLength);
//![example_functions]
    return 0;
}
//![example_functions]
