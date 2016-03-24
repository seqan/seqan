#include <iostream>
#include <seqan/stream.h>
#include <seqan/journaled_string_tree.h>

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

using namespace seqan;

int main()
{
    typedef JournaledStringTree<DnaString>   TJst;
    typedef Pattern<DnaString, MyersUkkonen> TPattern;
    typedef Traverser<TJst>::Type            TTraverser;

    DnaString seq = "AGATCGAGCGAGCTAGCGACTCAG";
    TJst jst(seq, 10);

    insert(jst, 1, 3, std::vector<unsigned>{1, 3, 5, 6, 7}, DeltaTypeDel());
    insert(jst, 8, "CGTA", std::vector<unsigned>{1, 2}, DeltaTypeIns());
    insert(jst, 10, 'C', std::vector<unsigned>{4, 9}, DeltaTypeSnp());
    insert(jst, 15, 2, std::vector<unsigned>{0, 4, 7}, DeltaTypeDel());
    insert(jst, 20, 'A', std::vector<unsigned>{0, 9}, DeltaTypeSnp());
    insert(jst, 20, Pair<unsigned, DnaString>(1, "CTC"), std::vector<unsigned>{1, 2, 3, 7}, DeltaTypeSV());

    DnaString ndl = "CCTCCA";
    TTraverser trav(jst, length(ndl) + 2);

    TPattern pat(ndl, -2);
    JstExtension<TPattern> ext(pat);

    MatchPrinter<TTraverser> delegate(trav);
    find(trav, ext, delegate);

    return 0;
}
