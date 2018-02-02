//![include]
#include <iostream>
#include <seqan/stream.h>
#include <seqan/journaled_string_tree.h>
//![include]

//![match_printer]
template <typename TTraverser>
struct MatchPrinter
{
    TTraverser & trav;

    MatchPrinter(TTraverser & _trav) : trav(_trav)
    {}

    void
    operator()()
    {
        auto pos = position(trav);
        for (auto p : pos)
        {
            std::cout  << "Seq: " << p.i1 << " Pos: " << p.i2 << std::endl;
        }
    }
};
//![match_printer]

//![typedef]
using namespace seqan;

int main()
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef Pattern<DnaString, Horspool>   TPattern;
    typedef Traverser<TJst>::Type          TTraverser;
//![typedef]

//![init]
    DnaString seq = "AGATCGAGCGAGCTAGCGACTCAG";
    TJst jst(seq, 10);

    insert(jst, 1, 3, std::vector<unsigned>{1, 3, 5, 6, 7}, DeltaTypeDel());
    insert(jst, 8, "CGTA", std::vector<unsigned>{1, 2}, DeltaTypeIns());
    insert(jst, 10, 'C', std::vector<unsigned>{4, 9}, DeltaTypeSnp());
    insert(jst, 15, 2, std::vector<unsigned>{0, 4, 7}, DeltaTypeDel());
    insert(jst, 20, 'A', std::vector<unsigned>{0, 9}, DeltaTypeSnp());
    insert(jst, 20, Pair<unsigned, DnaString>(1, "CTC"), std::vector<unsigned>{1, 2, 3, 7}, DeltaTypeSV());
//![init]

//![prepare_search]
    DnaString ndl = "AGCGT";
    TTraverser trav(jst, length(ndl));

    TPattern pat(ndl);
    JstExtension<TPattern> ext(pat);
//![prepare_search]

//![search]
    MatchPrinter<TTraverser> delegate(trav);
    find(trav, ext, delegate);

    return 0;
}
//![search]
