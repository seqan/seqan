#include <iostream>

#include <seqan/store.h>
#include <seqan/consensus.h>

using namespace seqan;

int main()
{
    FragmentStore<> store;

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

    ConsensusAlignmentOptions options;
    options.useContigID = false;
    consensusAlignment(store, options);

    std::cout << "Final alignment\n\n";
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ 80, 0, 30);

    return 0;
}
