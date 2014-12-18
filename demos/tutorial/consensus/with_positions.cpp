//![includes]
#include <iostream>

#include <seqan/store.h>
#include <seqan/consensus.h>

using namespace seqan;

int main()
{
//![includes]
//![fill_store]
    FragmentStore<> store;
    // Resize contigStore and contigNameStore (required for printing the first layout).
    resize(store.contigStore, 1);
    appendValue(store.contigNameStore, "ref");

    // Actual read layout.
    //
    // AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTT
    //           AAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTT
    //                AGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTA
    //                               ACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATT
    //                                         AGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAG

    // Append reads (includes small errors).
    appendRead(store, "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTT");
    appendRead(store, "AAAGTAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTT");
    appendRead(store, "AGTTGTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAATTTTCAATACTGTA");
    appendRead(store, "ACATCTCTTAAAGAGCTTTGATGCTAATTTAGTCAAATT");
    appendRead(store, "AGAGCTTTGATGCTAATTTAGTCAAATTTTCAATACTGTACAATCTTCTCTAG");

    // The position used in the following are only approximate and would
    // not lead to the read layout above.
    appendAlignedRead(store, 0, 0, 0, (int)length(store.readSeqStore[0]));
    appendAlignedRead(store, 1, 0, 12, 12 + (int)length(store.readSeqStore[1]));
    appendAlignedRead(store, 2, 0, 14, 14 + (int)length(store.readSeqStore[2]));
    appendAlignedRead(store, 3, 0, 18, 18 + (int)length(store.readSeqStore[3]));
    appendAlignedRead(store, 4, 0, 25, 25 + (int)length(store.readSeqStore[4]));

    // Print the (wrong) alignment.
    std::cout << "Initial alignment\n\n";
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 80, 0, 30);
//![fill_store]

//![compute_consensus]
    ConsensusAlignmentOptions options;
    options.useContigID = true;
    consensusAlignment(store, options);
//![compute_consensus]

//![print_layout]
    std::cout << "Final alignment\n\n";
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 80, 0, 30);

    return 0;
}
//![print_layout]
