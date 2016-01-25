//![includes]
#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/svg.h>

using namespace seqan;

int main()
{
    CharString fastaFileName = getAbsolutePath("demos/tutorial/fragment_store/example.fa");
    CharString samFileName = getAbsolutePath("demos/tutorial/fragment_store/example.sam");

    typedef FragmentStore<> TStore;

    TStore store;
    loadContigs(store, toCString(fastaFileName));
    BamFileIn file(toCString(samFileName));
    readRecords(store, file);

    SEQAN_ASSERT_GEQ(length(store.alignedReadStore), 5u);
//![includes]

//![typedefs]
    typedef Value<TStore::TContigStore>::Type                               TContig;
    typedef Value<TStore::TAlignedReadStore>::Type                          TAlignedRead;

    typedef Gaps<TContig::TContigSeq, AnchorGaps<TContig::TGapAnchors> >    TContigGaps;
    typedef Gaps<TStore::TReadSeq, AnchorGaps<TAlignedRead::TGapAnchors> >  TReadGaps;

    TStore::TReadSeq readSeq;
//![typedefs]

//![output]
    for (int i = 140; i < 160; i += 4)
    {
        TAlignedRead & ar = store.alignedReadStore[i];

        readSeq = store.readSeqStore[ar.readId];
        if (ar.endPos < ar.beginPos)
            reverseComplement(readSeq);

        TContigGaps contigGaps(
            store.contigStore[ar.contigId].seq,
            store.contigStore[ar.contigId].gaps);

        TReadGaps readGaps(
            readSeq,
            ar.gaps);

        setBeginPosition(contigGaps, std::min(ar.beginPos, ar.endPos));
        setEndPosition(contigGaps, std::max(ar.beginPos, ar.endPos));

        std::cout << "ALIGNMENT " << i << std::endl;
        std::cout << "\tcontig " << ar.contigId << ":\t" << contigGaps;
        std::cout << "     \t[" << beginPosition(contigGaps) << ".." << endPosition(contigGaps) << "[" << std::endl;
        std::cout << "\tread "   << ar.readId   << ":\t" << readGaps << std::endl;
        std::cout << std::endl;
    }
//![output]
//![appendix]
    return 0;
}
//![appendix]
