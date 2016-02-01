//![includes]
#include <iostream>
#include <random>
#include <seqan/align.h>
#include <seqan/index.h>

using namespace seqan;
//![includes]

//![matrix_init]
template <typename TStringSet, typename TIndexSpec>
void qgramCounting(TStringSet & set, TIndexSpec)
{
    typedef Index<TStringSet, TIndexSpec> TIndex;
    typedef typename Fibre<TIndex, QGramCounts>::Type TCounts;
    typedef typename Fibre<TIndex, QGramCountsDir>::Type TCountsDir;
    typedef typename Value<TCountsDir>::Type TDirValue;
    typedef typename Iterator<TCounts, Standard>::Type TIterCounts;
    typedef typename Iterator<TCountsDir, Standard>::Type TIterCountsDir;

    TIndex index(set);
    indexRequire(index, QGramCounts());

    // initialize distance matrix
    int seqNum = countSequences(index);
    Matrix<int, 2> distMat;
    setLength(distMat, 0, seqNum);
    setLength(distMat, 1, seqNum);
    resize(distMat, 0);

    std::cout << std::endl << "Length of the CountsDir fibre: " << length(indexCountsDir(index)) << std::endl;
    TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
    TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
    TIterCounts itCountsBegin = begin(indexCounts(index), Standard());
//![matrix_init]

//![matrix_calculation]
    // for each bucket count common q-grams for each sequence pair
    TDirValue bucketBegin = *itCountsDir;
    for (++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir)
    {
        TDirValue bucketEnd = *itCountsDir;

        // q-gram must occur in at least 2 different sequences
        if (bucketBegin != bucketEnd)
        {
            TIterCounts itA = itCountsBegin + bucketBegin;
            TIterCounts itEnd = itCountsBegin + bucketEnd;
            for (; itA != itEnd; ++itA)
                for (TIterCounts itB = itA; itB != itEnd; ++itB)
                    distMat((*itA).i1, (*itB).i1) += _min((*itA).i2, (*itB).i2);
        }
        bucketBegin = bucketEnd;
    }

    std::cout << std::endl << "Common 5-mers for Seq_i, Seq_j" << std::endl;
    std::cout << distMat;
}
//![matrix_calculation]

//![initialization]
int main()
{
    //  for the sake of reproducibility
    std::mt19937 rng;

    // create StringSet of 3 random sequences
    StringSet<DnaString> stringSet;
    reserve(stringSet, 3);
    for (int seqNo = 0; seqNo < 3; ++seqNo)
    {
        DnaString tmp;
        int len = rng() % 100 + 10;
        for (int i = 0; i < len; ++i)
            appendValue(tmp, Dna(rng() % 4));
        appendValue(stringSet, tmp);
        std::cout << ">Seq" << seqNo << std::endl << tmp << std::endl;
    }

    qgramCounting(stringSet, IndexQGram<UngappedShape<5> >());
    qgramCounting(stringSet, IndexQGram<UngappedShape<5>, OpenAddressing>());
    return 0;
}
//![initialization]
