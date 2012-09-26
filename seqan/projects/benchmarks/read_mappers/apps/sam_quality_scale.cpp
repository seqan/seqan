#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/store.h>

using namespace seqan;

// Analyzes the aligned read from the fragment score to get the
// average quality value per alignment error.
template <typename TFragmentStore>
void computeQualityScale(const char * contigsFilename, const char * samFilename, TFragmentStore & store) {
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TGapAnchors TReadGapAnchors;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    typedef typename TContigStoreElement::TGapAnchors TContigGapAnchors;
    typedef Gaps<TReadSeq, AnchorGaps<TContigGapAnchors> > TReadGaps;
    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> > TContigGaps;
    typedef typename Iterator<TContigGaps, Standard>::Type TContigGapAnchorsIterator;
    typedef typename Iterator<TReadGaps, Standard>::Type TReadGapAnchorsIterator;

    size_t alignedReadId = 0;
    size_t totalWeightedDistance = 0;
    size_t totalUnweightedDistance = 0;
    TReadSeq readSeq;
    for (TAlignedReadsIter it = begin(store.alignedReadStore, Standard()); it != end(store.alignedReadStore, Standard()); ++it, ++alignedReadId) {
        // Get contig and read sequences.
        TContigSeq & contigSeq = store.contigStore[it->contigId].seq;
        readSeq = store.readSeqStore[it->readId];
        // Get gaps for contig and read.
        TContigGaps contigGaps(contigSeq, store.contigStore[it->contigId].gaps);
        TReadGaps readGaps(readSeq, store.alignedReadStore[alignedReadId].gaps);
        // Limit contig gaps to aligned read position.
        setBeginPosition(contigGaps, _min(it->beginPos, it->endPos));
        setEndPosition(contigGaps, _max(it->beginPos, it->endPos));

        // Can reverse complement readSeq in place since gaps store a holder.
        if (it->beginPos > it->endPos) {
//             std::cerr << "complementing..." << std::endl;
            reverseComplement(readSeq);
        }

//         std::cerr << "Read:      " << readSeq << std::endl;
//         std::cerr << "Qualities: ";
//         for (unsigned i = 0; i < length(readSeq); ++i)
//             std::cerr << getQualityValue(readSeq[i]) << ", ";
//         std::cerr << std::endl;

        size_t readIndex = 0;  // Index in read.
        TContigGapAnchorsIterator contigGapsIt = begin(contigGaps, Standard());
        for (TReadGapAnchorsIterator readGapsIt = begin(readGaps, Standard()); readGapsIt != end(readGaps, Standard()); ++contigGapsIt, ++readGapsIt) {
            if (isGap(readGapsIt)) {  // Insertion in read.
                if (readIndex == 0) {
                    // At beginning -- pay first base's quality.
                    totalWeightedDistance += getQualityValue(front(readSeq));
                } else if (readIndex == length(readSeq) - 1) {
                    // At end -- pay last base's quality.
                    totalWeightedDistance += getQualityValue(back(readSeq));
                } else {
                    // Otherwise, pay mean of neighbouring bases' qualities, but at least one.
                    int x = getQualityValue(readSeq[readIndex]);
                    int y = getQualityValue(readSeq[readIndex + 1]);
                    totalWeightedDistance += _max(1, static_cast<int>(ceil(0.5 * (x + y))));
                }
            } else if (isGap(contigGapsIt) || convert<Dna5>(*readGapsIt) != convert<Dna5>(*contigGapsIt)) {
                // Insertion in contig or mismatch -- pay quality of read base.
                totalWeightedDistance += getQualityValue(convert<Dna5Q>(*readGapsIt));
                readIndex += 1;  // Advance one in read.
            } else {
                readIndex += 1;  // Advance one in read.
                continue;
            }
            // Tally unweighted error.
            totalUnweightedDistance += 1;
        }
    }

    // Compute average weight per unit distance.
    double result = (totalUnweightedDistance == 0) ? 0 : 1.0 * totalWeightedDistance / totalUnweightedDistance;
//     std::cout << "# contigs file\tSAM file\ttotal unweighted distance\ttotal weighted distance\tweight factor" << std::endl;
    std::cout << "{ \"genome\": \"" << contigsFilename << "\", \"sam_filename\": \"" << samFilename << "\", \"total_unweighted\": " << totalUnweightedDistance << ", \"total_weighted\": " << totalWeightedDistance << ", \"factor\": " << result << " }" << std::endl;
//     std::cout << contigsFilename << "\t" << samFilename << "\t" << totalUnweightedDistance << "\t" << totalWeightedDistance << "\t" << result << std::endl;
}


int main(int argc, char **argv) {
    // =================================================================
    // Interpret command line arguments.
    // =================================================================
    if (argc != 3) {
        std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE:  sam_quality_scale CONTIGS.FASTA ALIGNMENT.Sam" << std::endl;
        return 1;
    }

    // =================================================================
    // Read contigs and Sam file.
    // =================================================================
    FragmentStore<> store;

    std::cerr << "Loading contigs from " << argv[1] << std::endl;
    if (!loadContigs(store, argv[1])) {
        std::cerr << "Could not read contigs." << std::endl;
            return 1;
    }
    
    {
        std::cerr << "Loading Sam file " << argv[2] << std::endl;
        std::fstream fstrm(argv[2], std::ios_base::in | std::ios_base::binary);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open Sam file." << std::endl;
            return 1;
        }
        read(fstrm, store, Sam());
    }

    // =================================================================
    // Compute the quality scale.
    // =================================================================
    std::cerr << "Computing quality scale..." << std::endl;
    computeQualityScale(argv[1], argv[2], store);

    return 0;
}

