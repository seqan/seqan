#ifndef APPS_BS_TOOLS_BISAR_SCORE_H_
#define APPS_BS_TOOLS_BISAR_SCORE_H_


namespace seqan {

struct MyFragmentStoreConfig;

template<>
struct FragmentStoreConfig<MyFragmentStoreConfig> :
	public FragmentStoreConfig<>
{
	typedef double		TMappingQuality;    // ?
};

template <typename T = void>
struct TranslateTableDna5OrdValueToDna5OrdValueComplement_
{
    static int const VALUE[5];
};

template <typename T>
int const TranslateTableDna5OrdValueToDna5OrdValueComplement_<T>::VALUE[5] = {3, 2, 1, 0, 4};

template <typename TValue> struct FunctorDna5OrdValueComplement;

template <>
struct FunctorDna5OrdValueComplement<int> : public std::unary_function<int,int>
{
    inline int operator()(int x) const
    {
        return TranslateTableDna5OrdValueToDna5OrdValueComplement_<>::VALUE[x];
    }
};

////////////////////////////////////////////////////////////////////////////////////
/// Metafunction SequenceEntryForScore                       [Bs Score]
///////////////////////////////////////////////////////////////////////////////////


template <typename TSequence>
struct SequenceEntryForScore<Score<int, BsCaseCT >, TSequence>
{
    typedef typename Value<TSequence>::Type Type;   // For the beginning, base and quality are enough
};

template <typename TSequence>
struct SequenceEntryForScore<Score<int, BsCaseGA >, TSequence>
{
    typedef typename Value<TSequence>::Type Type;
};

template <typename TSequence, typename TPosition>
inline typename Value<TSequence>::Type
sequenceEntryForScore(Score<int, BsCaseCT > const & /*scoringScheme*/, TSequence const & seq, TPosition pos)
{
    return seq[pos];
}

template <typename TSequence, typename TPosition>
inline typename Value<TSequence>::Type
sequenceEntryForScore(Score<int, BsCaseGA > const & /*scoringScheme*/, TSequence const & seq, TPosition pos)
{
    return seq[pos];
}

//////////////////////////////////////////////////////////////////////
// Match/ Mismatch Bs scores
//////////////////////////////////////////////////////////////////////


template<typename TScore, typename TSubstMatrix, typename TSeqErrorFreqs>
inline void
computeBsScores(TScore &sc_data,
                TSubstMatrix const &s,
                TSeqErrorFreqs const &seqErrorFreqs,
                BsCaseCT const &,
                Left const &)
{
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned i = 0; i < 5; ++i)    // Ref base
        {
            resize(sc_data[q][i], 5, Exact());
            for (unsigned j = 0; j < 5; ++j)    // Read base
            {
                sc_data[q][i][j] = s[i*5 + j] * (1.0-e);   // No seq error

                for (unsigned o = 0; o < 5; ++o)
                {
                    if (o != j) sc_data[q][i][j] += s[i*5 + o] * e * seqErrorFreqs[o*5 + j]; // Seq error
                    //if (q == 62 && i == 1  && (j == 1 || j == 3) )
                        //std::cout << "Ref:C read:C/T: " << (Dna5)j << " o: " << (Dna5)o << "  subScore: " << (s[i*5 + o] * e * seqErrorFreqs[o*5 + j]) << " s: " << s[i*5 + o] << "  seqError:" << seqErrorFreqs[o*5 + j] << std::endl;
                }
            }
        }
    }
}

template<typename TScore, typename TSubstMatrix, typename TSeqErrorFreqs>
inline void
computeBsScores(TScore &sc_data,
                TSubstMatrix const &s,
                TSeqErrorFreqs const &seqErrorFreqs,
                BsCaseCT const &,
                Right const &)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned i = 0; i < 5; ++i)    // Ref base
        {
            resize(sc_data[q][i], 5, Exact());
            for (unsigned j = 0; j < 5; ++j)    // Read base
            {
                sc_data[q][i][j] = s[i*5 + j] * (1.0-e);   // No seq error

                for (unsigned o = 0; o < 5; ++o)
                {
                    if (o != j) sc_data[q][i][j] += s[i*5 + o] * e * seqErrorFreqs[fCompl(o)*5 + fCompl(j)]; // Seq error (regarding original read sequence)
                }
            }
        }
    }
}

template<typename TScore, typename TSubstMatrix, typename TSeqErrorFreqs>
inline void
computeBsScores(TScore &sc_data,
                TSubstMatrix const &s,
                TSeqErrorFreqs const &seqErrorFreqs,
                BsCaseGA const &,
                Left const &)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned i = 0; i < 5; ++i)    // Ref base
        {
            resize(sc_data[q][i], 5, Exact());
            for (unsigned j = 0; j < 5; ++j)    // Read base
            {
                sc_data[q][i][j] = s[i*5 + j] * (1.0-e);   // No seq error

                for (unsigned o = 0; o < 5; ++o)
                {
                    if (o != j) sc_data[q][i][j] += s[i*5 + o] * e * seqErrorFreqs[fCompl(o)*5 + fCompl(j)]; // Seq error (regarding original read sequence)
                }
            }
        }
    }
}


template<typename TScore, typename TSubstMatrix, typename TSeqErrorFreqs>
inline void
computeBsScores(TScore &sc_data,
                TSubstMatrix const &s,
                TSeqErrorFreqs const &seqErrorFreqs,
                BsCaseGA const &,
                Right const &)
{
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned i = 0; i < 5; ++i)    // Ref base
        {
            resize(sc_data[q][i], 5, Exact());
            for (unsigned j = 0; j < 5; ++j)    // Read base
            {
                sc_data[q][i][j] = s[i*5 + j] * (1.0-e);   // No seq error

                for (unsigned o = 0; o < 5; ++o)
                {
                    if (o != j) sc_data[q][i][j] += s[i*5 + o] * e * seqErrorFreqs[o*5 + j]; // Seq error
                    //if (q == 62 && i == 1  && (j == 1 || j == 3) )
                        //std::cout << "Ref:C read:C/T: " << (Dna5)j << " o: " << (Dna5)o << "  subScore: " << (s[i*5 + o] * e * seqErrorFreqs[o*5 + j]) << " s: " << s[i*5 + o] << "  seqError:" << seqErrorFreqs[o*5 + j] << std::endl;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// ReadGapScores
//////////////////////////////////////////////////////////////////////


template<typename TScores, typename TValue, typename TSubstMatrix, typename TDelErrorFreqs, typename TOptions>
inline void
computeReadGapScores(TScores &sc_data,
                     TValue const &gapScore,
                     TSubstMatrix const &s,
                     TDelErrorFreqs const &delErrorFreqs,
                     TOptions &options,
                     BsCaseCT const &,
                     Left const &)
{
    resize(sc_data, 5, Exact());
    TValue e = options.delErrorRate;
    for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be deleted
    {
        sc_data[b] = exp(gapScore*options.lambda) * (1.0-e);   // Genome deletion, but no seq error
        for (unsigned i = 0; i < 5; ++i)
        {
            sc_data[b] += s[b*5 + i] * e * delErrorFreqs[i] * options.scalingFactorDelErrors;  // No genome deletion, but deletion seq error -> Match/Mismatch
        }
    }
}

template<typename TScores, typename TValue, typename TSubstMatrix, typename TDelErrorFreqs, typename TOptions>
inline void
computeReadGapScores(TScores &sc_data,
                     TValue const &gapScore,
                     TSubstMatrix const &s,
                     TDelErrorFreqs const &delErrorFreqs,
                     TOptions &options,
                     BsCaseCT const &,
                     Right const &)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    resize(sc_data, 5, Exact());
    TValue e = options.delErrorRate;
    for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be deleted
    {
        sc_data[b] = exp(gapScore*options.lambda) * (1.0-e);   // Genome deletion, but no seq error
        for (unsigned i = 0; i < 5; ++i)
        {
            sc_data[b] += s[b*5 + i] * e * delErrorFreqs[fCompl(i)] * options.scalingFactorDelErrors;  // No genome deletion, but deletion seq error -> Match/Mismatch
        }
    }
}

template<typename TScores, typename TValue, typename TSubstMatrix, typename TDelErrorFreqs, typename TOptions>
inline void
computeReadGapScores(TScores &sc_data,
                     TValue const &gapScore,
                     TSubstMatrix const &s,
                     TDelErrorFreqs const &delErrorFreqs,
                     TOptions &options,
                     BsCaseGA const &,
                     Left const &)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    resize(sc_data, 5, Exact());
    TValue e = options.delErrorRate;
    for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be deleted
    {
        sc_data[b] = exp(gapScore*options.lambda) * (1.0-e);   // Genome deletion, but no seq error
        for (unsigned i = 0; i < 5; ++i)
        {
            sc_data[b] += s[b*5 + i] * e * delErrorFreqs[fCompl(i)] * options.scalingFactorDelErrors;  // No genome deletion, but deletion seq error -> Match/Mismatch
        }
    }
}

template<typename TScores, typename TValue, typename TSubstMatrix, typename TDelErrorFreqs, typename TOptions>
inline void
computeReadGapScores(TScores &sc_data,
                     TValue const &gapScore,
                     TSubstMatrix const &s,
                     TDelErrorFreqs const &delErrorFreqs,
                     TOptions &options,
                     BsCaseGA const &,
                     Right const &)
{
    resize(sc_data, 5, Exact());
    TValue e = options.delErrorRate;
    for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be deleted
    {
        sc_data[b] = exp(gapScore*options.lambda) * (1.0-e);   // Genome deletion, but no seq error
        for (unsigned i = 0; i < 5; ++i)
        {
            sc_data[b] += s[b*5 + i] * e * delErrorFreqs[i] * options.scalingFactorDelErrors;  // No genome deletion, but deletion seq error -> Match/Mismatch
        }
    }
}


//////////////////////////////////////////////////////////////////////
// RefGapScores
//////////////////////////////////////////////////////////////////////


template<typename TScores, typename TValue, typename TInsErrorFreqs, typename TSeqErrorFreqs, typename TOptions>
inline void
computeRefGapScores(TScores &sc_data,
                    TValue const &simpleScore,
                    TInsErrorFreqs * &insErrorFreqs,
                    TSeqErrorFreqs * &seqErrorFreqs,
                    TOptions &options,
                    BsCaseCT const &,
                    Left const &)
{
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be inserted
        {
            sc_data[q][b] = exp(simpleScore*options.lambda) * (1.0-e);   // Genome insertion, but no seq error
            sc_data[q][b] += 1.0 * e * insErrorFreqs[b];               // No genome insertion, but insertion seq error
            for (unsigned i = 0; i < 5; ++i)
            {
                if (i != b) sc_data[q][b] += exp(simpleScore*options.lambda) * e * seqErrorFreqs[i*5 + b];  // Genome insertion and seq error
            } // Assuming genome insertions are equally distributed
        }
    }
}

template<typename TScores, typename TValue, typename TInsErrorFreqs, typename TSeqErrorFreqs, typename TOptions>
inline void
computeRefGapScores(TScores &sc_data,
                    TValue const &simpleScore,
                    TInsErrorFreqs * &insErrorFreqs,
                    TSeqErrorFreqs * &seqErrorFreqs,
                    TOptions &options,
                    BsCaseCT const &,
                    Right const &)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be inserted
        {
            sc_data[q][b] = exp(simpleScore*options.lambda) * (1.0-e);   // Genome insertion, but no seq error
            sc_data[q][b] += 1.0 * e * insErrorFreqs[fCompl(b)];               // No genome insertion, but insertion seq error
            for (unsigned i = 0; i < 5; ++i)
            {
                if (i != b) sc_data[q][b] += exp(simpleScore*options.lambda) * e * seqErrorFreqs[fCompl(i)*5 + fCompl(b)];  // Genome insertion and seq error
            }
        }
    }
}

template<typename TScores, typename TValue, typename TInsErrorFreqs, typename TSeqErrorFreqs, typename TOptions>
inline void
computeRefGapScores(TScores &sc_data,
                    TValue const &simpleScore,
                    TInsErrorFreqs * &insErrorFreqs,
                    TSeqErrorFreqs * &seqErrorFreqs,
                    TOptions &options,
                    BsCaseGA const &,
                    Left const &)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be inserted
        {
            sc_data[q][b] = exp(simpleScore*options.lambda) * (1.0-e);   // Genome insertion, but no seq error
            sc_data[q][b] += 1.0 * e * insErrorFreqs[fCompl(b)];               // No genome insertion, but insertion seq error
            for (unsigned i = 0; i < 5; ++i)
            {
                if (i != b) sc_data[q][b] += exp(simpleScore*options.lambda) * e * seqErrorFreqs[fCompl(i)*5 + fCompl(b)];  // Genome insertion and seq error
            }
        }
    }
}

template<typename TScores, typename TValue, typename TInsErrorFreqs, typename TSeqErrorFreqs, typename TOptions>
inline void
computeRefGapScores(TScores &sc_data,
                    TValue const &simpleScore,
                    TInsErrorFreqs * &insErrorFreqs,
                    TSeqErrorFreqs * &seqErrorFreqs,
                    TOptions &options,
                    BsCaseGA const &,
                    Right const &)
{
    resize(sc_data, 63, Exact());
    for (unsigned q = 0; q < 63; ++q)
    {
        long double e = pow(10, -(long double)q/10.0);

        resize(sc_data[q], 5, Exact());
        for (unsigned b = 0; b < 5; ++b)    // Base which is assumed to be inserted
        {
            sc_data[q][b] = exp(simpleScore*options.lambda) * (1.0-e);   // Genome insertion, but no seq error
            sc_data[q][b] += 1.0 * e * insErrorFreqs[b];               // No genome insertion, but insertion seq error
            for (unsigned i = 0; i < 5; ++i)
            {
                if (i != b) sc_data[q][b] += exp(simpleScore*options.lambda) * e * seqErrorFreqs[i*5 + b];  // Genome insertion and seq error
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Bs score
//////////////////////////////////////////////////////////////////////


template<typename TBsCase, typename TModel, typename TSegment>
class Score<int, BsTagList<TBsCase, TModel, TSegment> >
{
public:
    String<String<String<double> > > data;  // precomputed score values for all possible qual values, bases etc.

    // The gap extension score.
    String<String<double> >     data_gap_extend_ref;    // Score dependent on read base and qual
    String<double>              data_gap_extend_read;   // Score dependent on ref base

    // The gap open score.
    String<String<double> >     data_gap_open_ref;      // Score dependent on read base and qual
    String<double>              data_gap_open_read;     // Score dependent on ref base

    double lambda;

    // Use default values apart of gap costs
    /*template<typename TOptions>
    Score(TValue simple_gap_extend, TValue simple_gap_open, TValue delErrorRate, TOptions &options)
    {
        BsSubstitutionMatrix<TValue, TBsCase, BsSimple> bsSubstitutionMatrix;
        TValue const * seqErrorFreqs = SeqErrorFreqs<TValue, BsNonSimple>::getData();
        TValue const * insErrorFreqs = InsErrorFreqs<TValue, BsNonSimple>::getData();
        TValue const * delErrorFreqs = DelErrorFreqs<TValue, BsNonSimple>::getData();
        computeBsScores((*this).data, bsSubstitutionMatrix.data_tab, seqErrorFreqs, TBsCase(), TSegment());
        computeReadGapScores((*this).data_gap_open_read, options.gapOpenScore, bsSubstitutionMatrix.data_tab, delErrorFreqs, options, TBsCase(), TSegment());
        computeReadGapScores((*this).data_gap_extend_read, options.gapExtendScore, bsSubstitutionMatrix.data_tab, delErrorFreqs, options, TBsCase(), TSegment());
        computeRefGapScores((*this).data_gap_open_ref, options.gapOpenScore, insErrorFreqs, seqErrorFreqs, options, TBsCase(), TSegment());
        computeRefGapScores((*this).data_gap_extend_ref, options.gapExtendScore, insErrorFreqs, seqErrorFreqs, options, TBsCase(), TSegment());
    }*/

    // User defined rates
    template<typename TOptions, typename TBsSubstitutionMatrix, typename TValue>
    Score(TOptions &options,
            TBsSubstitutionMatrix &bsSubstitutionMatrix,
            TValue const * &seqErrorFreqs,
            TValue const * &insErrorFreqs,
            TValue const * &delErrorFreqs) : lambda(0.0)
    {
        computeBsScores((*this).data, bsSubstitutionMatrix.data_tab, seqErrorFreqs, TBsCase(), TSegment());

        computeReadGapScores((*this).data_gap_open_read, options.gapOpenScore, bsSubstitutionMatrix.data_tab, delErrorFreqs, options, TBsCase(), TSegment());
        computeReadGapScores((*this).data_gap_extend_read, options.gapExtendScore, bsSubstitutionMatrix.data_tab, delErrorFreqs, options, TBsCase(), TSegment());

        computeRefGapScores((*this).data_gap_open_ref, options.gapOpenScore, insErrorFreqs, seqErrorFreqs, options, TBsCase(), TSegment());
        computeRefGapScores((*this).data_gap_extend_ref, options.gapExtendScore, insErrorFreqs, seqErrorFreqs, options, TBsCase(), TSegment());
    }
};

///////////////////////////////////////////////////////////////////////////////////
// Bs scores functions
///////////////////////////////////////////////////////////////////////////////////

template <typename TBsCase, typename TModel, typename TSegment, typename TSeqHValue, typename TSeqVValue>
inline int
scoreGapOpenVertical(
    Score<int, BsTagList<TBsCase, TModel, TSegment> > const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & seqVVal)
{
    unsigned int qual = getQualityValue(seqVVal);
    unsigned int j = (Dna5)seqVVal;  // Read base

    return (int)(std::log10(me.data_gap_open_ref[qual][j]) *10000 + 0.5);
}

template <typename TBsCase, typename TModel, typename TSegment, typename TSeqHValue, typename TSeqVValue>
inline int
scoreGapOpenHorizontal(
    Score<int, BsTagList<TBsCase, TModel, TSegment> > const & me,
    TSeqHValue const & seqHVal,
    TSeqVValue const & /*seqVVal*/)
{
    unsigned int i = (Dna5)seqHVal;  // Ref base
    return (int)(std::log10(me.data_gap_open_read[i]) *10000 + 0.5);
}

template <typename TBsCase, typename TModel, typename TSegment, typename TSeqHValue, typename TSeqVValue>
inline int
scoreGapExtendVertical(
    Score<int, BsTagList<TBsCase, TModel, TSegment> > const & me,
    TSeqHValue const & /*seqHVal*/,
    TSeqVValue const & seqVVal)
{
    unsigned int qual = getQualityValue(seqVVal);
    unsigned int j = (Dna5)seqVVal;  // Read base

    return (int)(std::log10(me.data_gap_extend_ref[qual][j]) *10000 + 0.5);
}

template <typename TBsCase, typename TModel, typename TSegment, typename TSeqHValue, typename TSeqVValue>
inline int
scoreGapExtendHorizontal(
    Score<int, BsTagList<TBsCase, TModel, TSegment> > const & me,
    TSeqHValue const & seqHVal,
    TSeqVValue const & /*seqVVal*/)
{
    unsigned int i = (Dna5)seqHVal;  // Ref base
    return (int)(std::log10(me.data_gap_extend_read[i]) *10000 + 0.5);
}

// Just for global alignment function to check, if open and extenc cost are the same
template <typename TBsCase, typename TModel, typename TSegment>
inline int
scoreGapExtend(Score<int, BsTagList<TBsCase, TModel, TSegment> > const & /*me*/) {
    //std::cout << "This should happen only once: scoreGapExtend()"<< std::endl;
    return (int)-100;
}
template <typename TBsCase, typename TModel, typename TSegment>
inline int
scoreGapOpen(Score<int, BsTagList<TBsCase, TModel, TSegment> > const & /*me*/) {
    //std::cout << "This should happen only once: scoreGapOpen()" << std::endl;
    return (int)-200;
}


template <typename TBsCase, typename TModel, typename TSegment, typename TVal1>
inline int
score(Score<int, BsTagList<TBsCase, TModel, TSegment> > const & sc, TVal1 val1, Dna5Q val2)
{
    unsigned int i = (Dna5) val1;
    unsigned int j = (Dna5) val2;
    int qual = getQualityValue(val2);

    //std::cout << "qual: " << qual << std::endl;
    return (int)(std::log10(sc.data[qual][i][j]) *10000 + 0.5);
}


}






#endif

