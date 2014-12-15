#include <seqan/seq_io.h>  // for FaiIndex
#include <seqan/bam_io.h>  // for BamStream

int main()
{
    // Open FASTA, FAI, and BAM file.  Build FAI if necessary.
    seqan::FaiIndex faiIndex;
    if (!open(faiIndex, "filename.fasta"))          // try to load
        if (!build(faiIndex, "filename.fasta"))     // try to build
            return 1;  // Error.
    seqan::BamFileIn bamFileIn("file.bam");
    seqan::BamHeader header;
    readRecord(header, bamFileIn);

    // Build mapping from bamSeqIds to fastaSeqIds;
    seqan::String<int> mapping;
    resize(mapping, length(contigNames(context(bamFileIn))), -1);
    for (unsigned i = 0; i < length(contigNames(context(bamFileIn))); ++i)
        if (!getIdByName(mapping[i], faiIndex, contigNames(context(bamFileIn))[i]))
        {
            std::cerr << "ERROR: Sequence "
                      << contigNames(context(bamFileIn))[i]
                      << "unknown in FASTA Index.\n";
            return 1; 
        }

    return 0;
}
