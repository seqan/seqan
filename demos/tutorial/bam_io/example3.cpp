#include <seqan/seq_io.h>  // for FaiIndex
#include <seqan/bam_io.h>  // for BamStream

using namespace seqan;

int main()
{
    // Open FASTA, FAI, and BAM file.  Build FAI if necessary.
    FaiIndex faiIndex;
    if (!open(faiIndex, "filename.fasta"))          // try to load
        if (!build(faiIndex, "filename.fasta"))     // try to build
            return 1;  // Error.
    BamFileIn bamFileIn("file.bam");
    BamHeader header;
    readRecord(header, bamFileIn);

    // Build mapping from bamSeqIds to fastaSeqIds;
    String<int> mapping;
    resize(mapping, length(nameStore(context(bamFileIn))), -1);
    for (unsigned i = 0; i < length(nameStore(context(bamFileIn))); ++i)
        if (!getIdByName(mapping[i], faiIndex, nameStore(context(bamFileIn))[i]))
        {
            std::cerr << "ERROR: Sequence "
                      << nameStore(context(bamFileIn))[i]
                      << "unknown in FASTA Index.\n";
            return 1; 
        }

    return 0;
}
