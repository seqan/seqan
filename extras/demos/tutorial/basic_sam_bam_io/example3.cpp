#include <seqan/seq_io.h>  // for FaiIndex
#include <seqan/bam_io.h>  // for BamStream

int main()
{
    // Open FASTA, FAI, and BAM file.  Build FAI if necessary.
    seqan::FaiIndex faiIndex;
    if (read(faiIndex, "filename.fasta") != 0)       // try to load
        if (build(faiIndex, "filename.fasta") != 0)  // try to build
            return 1;  // Error.
    seqan::BamStream bamIO("file.bam");

    // Build mapping from bamSeqIds to fastaSeqIds;
    seqan::String<unsigned> mapping;
    resize(mapping, length(bamIO.header.sequenceInfos), 0);
    for (unsigned i = 0; i < length(bamIO.header.sequenceInfos); ++i)
    {
        seqan::CharString seqName = bamIO.header.sequenceInfos[i].i1;
        if (!getIdByName(faiIndex, seqName, mapping[i]))
        {
            std::cerr << "ERROR: Sequene "
                      << bamIO.header.sequenceInfos[i].i1
                      << "unknown in FASTA Index.\n";
            return 1; 
        }
    }

    return 0;
}
