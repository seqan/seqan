#ifndef APPS_BS_TOOLS_BISAR_SCORE_DATA_H_
#define APPS_BS_TOOLS_BISAR_SCORE_DATA_H_


namespace seqan {

struct BsCaseCT_;
typedef Tag<BsCaseCT_>   BsCaseCT;

struct BsCaseGA_;
typedef Tag<BsCaseGA_>   BsCaseGA;

struct Left_;
typedef Tag<Left_>   Left;

struct Right_;
typedef Tag<Right_>  Right;



struct BsNonSimple_;        // Tag for assuming most simple model: simple substitution matrix, equal error distributions
typedef Tag<BsNonSimple_>   BsNonSimple;

struct BsSimple_;
typedef Tag<BsSimple_>      BsSimple;

template <typename TBsCase, typename TModel, typename TSegment, typename TCellDescriptor = InnerCell>
struct BsTagList
{
    typedef TBsCase             TypeBsCase;
    typedef TModel              TypeModel;
    typedef TSegment            TypeSeqment;
    typedef TCellDescriptor     TypeCellDescriptor;
};

// Substitution matrix
template <typename TValue = double, typename TSpec1 = BsCaseCT, typename TSpec2 = BsNonSimple>
class BsSubstitutionMatrix;


template <typename TValue>
class BsSubstitutionMatrix<TValue, BsCaseCT, BsSimple > {
public:
    // Static computation of the required array size.
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    // The data table.
    TValue data_tab[TAB_SIZE];

    BsSubstitutionMatrix(TValue &fMeth, TValue &fBsConversion, TValue &seqIdentity, TValue &fN)
    {
        TValue fA = 1.0/4.0 - fN/16.0;  // Uniform base distribution, taking reference Ns into account
        TValue fC = 1.0/4.0 - fN/16.0;
        TValue fG = 1.0/4.0 - fN/16.0;
        TValue fT = 1.0/4.0 - fN/16.0;

        TValue fMatch = (seqIdentity-fN)/4.0;
        TValue fSubst = (1.0-seqIdentity)/12.0;
        TValue fNtoX = fN/4.0;            // Ns in reference, X in read

        //std::cout << "Construct simple subst. matrix, CT" << std::endl;
        // Because of Bs treatment different rates in reads
        TValue fCRead = fC - fC *(1.0-fMeth) * fBsConversion;  // Compute expected frequency of Cs in Bs reads
        TValue fTRead = fT + fC *(1.0-fMeth) * fBsConversion;  // .. of Ts

        TValue pCC = fMatch*(1.0-fMeth)*(1.0-fBsConversion) + fMatch*fMeth;   // Prob. that C stays C: no substitution and (unmethylated and not converted or methylated)
        TValue pCT = fSubst + fMatch*(1.0-fMeth)*fBsConversion;         // Prob. that C turns into T: substitution or unmethylated and converted
        // add prob. that methylated and converted?
        static TValue const tab[25] = {
                      //A,                      C,                          G,                          T,                                  N (ref?)
                      fMatch/(fA*fA),           fSubst/(fA*fCRead),         fSubst/(fA*fG),             fSubst/(fA*fTRead),                         0,     // A
                      fSubst/(fC*fA),           pCC/(fC*fCRead),            fSubst/(fC*fG),             pCT/(fC*fTRead),                            0,     // C   Ref:C - read:T  ?
                      fSubst/(fG*fA),           fSubst/(fG*fCRead),         fMatch/(fG*fG),             fSubst/(fG*fTRead),                         0,     // G
                      fSubst/(fT*fA),           fSubst/(fT*fCRead),         fSubst/(fT*fG),             fMatch/(fT*fTRead),                         0,     // T
                      fNtoX/(fN*fA),            fNtoX/(fN*fCRead),          fNtoX/(fN*fG),              fNtoX/(fN*fTRead),                          0      // N in ref
        };
        arrayCopy(tab, tab + 25, (*this).data_tab);
    }
};


template <typename TValue>
class BsSubstitutionMatrix<TValue, BsCaseGA, BsSimple >
{
public:
    // Static computation of the required array size.
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    // The data table.
    TValue data_tab[TAB_SIZE];

    BsSubstitutionMatrix(TValue &fMeth, TValue &fBsConversion, TValue &seqIdentity, TValue &fN)
    {
        TValue fA = 1.0/4.0 - fN/16.0;
        TValue fC = 1.0/4.0 - fN/16.0;
        TValue fG = 1.0/4.0 - fN/16.0;
        TValue fT = 1.0/4.0 - fN/16.0;

        TValue fMatch = (seqIdentity-fN)/4.0;
        TValue fSubst = (1.0-seqIdentity)/12.0;
        TValue fNtoX = fN/4.0;            // Ns in reference, X in read

        // Because of Bs treatment different rates in reads
        TValue fGRead = fG - fG *(1.0-fMeth) * fBsConversion;  // Compute expected frequency of Cs in Bs reads
        TValue fARead = fA + fA *(1.0-fMeth) * fBsConversion;  // .. of Ts

        TValue pGG = fMatch*(1.0-fMeth)*(1.0-fBsConversion) + fMatch*fMeth;   // add prob. that methylated and converted?
        TValue pGA = fSubst + fMatch*(1.0-fMeth)*fBsConversion;

        static TValue const tab[25] = {
                      //A,                                      C,                          G,                                      T,                          N (ref?)
                      fMatch/(fA*fARead),                       fSubst/(fA*fC),             fSubst/(fA*fGRead),                     fSubst/(fA*fT),             0,     // A
                      fSubst/(fC*fARead),                       fMatch/(fC*fC),             fSubst/(fC*fGRead),                     fSubst/(fC*fT),             0,     // C   Ref:C - read:T  ?
                      pGA/(fG*fARead),                          fSubst/(fG*fC),             pGG/(fG*fGRead),                        fSubst/(fG*fT),             0,     // G
                      fSubst/(fT*fARead),                       fSubst/(fT*fC),             fSubst/(fT*fGRead),                     fMatch/(fT*fT),             0,     // T
                      fNtoX/(fN*fARead),                        fNtoX/(fN*fC),              fNtoX/(fN*fGRead),                      fNtoX/(fN*fT),              0      // N in ref
        };
        arrayCopy(tab, tab + 25, (*this).data_tab);
    }
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// error frequencies
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TValue = double, typename TSpec = BsNonSimple>
struct SeqErrorFreqs;

// Default sequencing error frequencies (substitutions only)
template <typename TValue>
struct SeqErrorFreqs<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {

                        // Column sum must be 1   (in the current model)
        static TValue const _data[TAB_SIZE] = {
                // To A, C, G, T, N (ref?)
                        0,            (0.14/0.38),   (0.05/0.16),    (0.05/0.21),    1.0/4.0,    // From A
                       (0.13/0.25),    0,            (0.02/0.16),    (0.04/0.21),    1.0/4.0,    // C
                       (0.04/0.25),   (0.08/0.38),    0,             (0.12/0.21),    1.0/4.0,
                       (0.08/0.25),   (0.15/0.38),   (0.09/0.16),     0,             1.0/4.0,
                       0,              0,             0,              0,             0      // N?
        };
        return _data;
    }
};

// Simple sequencing error frequencies (substitutions only)
// [ordValue(realBase) * 5 + ordValue(observed base)]
template <typename TValue>
struct SeqErrorFreqs<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {

        TValue fE = 1.0/3.0;                        // Column sum must be 1
        static TValue const _data[TAB_SIZE] = {
                // To A, C, G, T, N (ref?)
                      0,    fE,   fE,   fE, 1.0/4.0,    // From A
                      fE,   0,    fE,   fE, 1.0/4.0,    // C
                      fE,   fE,   0,    fE, 1.0/4.0,
                      fE,   fE,   fE,   0,  1.0/4.0,
                      0,    0,    0,    0,  0
        };
        return _data;
    }
};

template <typename TValue = double, typename TSpec = BsNonSimple>
struct SeqErrorFreqsTo;

template <typename TValue>
struct SeqErrorFreqsTo<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
        TValue f = 1.0/5.0;
                                                // A,    C,     G,     T,    N
        static TValue const _data[VALUE_SIZE] = {f, f, f, f, f};      // siehe Dohm et.al. 2011

        return _data;
    }
};


template <typename TValue>
struct SeqErrorFreqsTo<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
                                                // A,    C,     G,     T,    N
        static TValue const _data[VALUE_SIZE] = {0.25, 0.38, 0.16, 0.21, 0.2};      // siehe Dohm et.al. 2011 // TODO N?

        return _data;
    }
};

template <typename TValue = double, typename TSpec = BsNonSimple>
struct BaseErrorFreqsFrom;

template <typename TValue>
struct BaseErrorFreqsFrom<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE
    };

    static inline TValue const * getData() {
                                                // A,  C,    G,    T     N
        static TValue const _data[VALUE_SIZE] =  {0.25, 0.19, 0.24, 0.33, 0.0};      // siehe Dohm et.al. 2008
        return _data;                                                                // N  stays N
    }
};

template <typename TValue>
struct BaseErrorFreqsFrom<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE
    };

    static inline TValue const * getData() {
        TValue f = 1.0/5.0;
                                                // A,  C,    G,    T     N
        static TValue const _data[VALUE_SIZE] = {f, f, f, f, f};      // siehe Dohm et.al. 2008
        return _data;
    }
};




template <typename TValue = double, typename TSpec = BsNonSimple>
struct InsErrorFreqs;

template <typename TValue>
struct InsErrorFreqs<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
                                                // A,    C,     G,     T,    N
        static TValue const _data[VALUE_SIZE] = {0.43, 0.065, 0.065, 0.43, 0.01};      // siehe Dohm et.al. 2011 // TODO N?

        return _data;
    }
};

template <typename TValue>
struct InsErrorFreqs<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
        TValue f = 1.0/5.0;
        static TValue const _data[VALUE_SIZE] = {f, f, f, f, f};
        return _data;
    }
};


template <typename TValue = double, typename TSpec = BsNonSimple>
struct DelErrorFreqs;

template <typename TValue>
struct DelErrorFreqs<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
                                                // A,    C,     G,     T,    N
        static TValue const _data[VALUE_SIZE] = {0.42, 0.075, 0.075, 0.42, 0.00};      // siehe Dohm et.al. 2011

        return _data;
    }
};

template <typename TValue>
struct DelErrorFreqs<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE
    };

    static inline TValue const * getData() {
        TValue f = 1.0/5.0;
        static TValue const _data[VALUE_SIZE] = {f, f, f, f, f};
        return _data;
    }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// No simulated data for these non-simple cases yet...
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class BsSubstitutionMatrix<TValue, BsCaseCT, BsNonSimple > {
public:
    // Static computation of the required array size.
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    // The data table.
    TValue data_tab[TAB_SIZE];

    // Frequencies:
    TValue fMatch;
    TValue fSubst;
    TValue fNtoX;
    TValue fTransi;     // Transition rate: A <-> G, C <-> T
    TValue fTransv;     // Transversion rate: everything else
    TValue fMeth;
    TValue fBsConversion;
    // Base frequencies in genome
    TValue fA;
    TValue fC;
    TValue fG;
    TValue fT;
    TValue fN;

    BsSubstitutionMatrix(TValue &fMeth, TValue &fBsConversion):
        fMatch(0.89/4.0),
        fSubst(0.09/12.0),
        fNtoX(0.02/4.0),
        fTransi(0.5/4.0),
        fTransv(0.5/8.0),
        fA(0.9/4.0),
        fC(0.9/4.0),
        fG(0.9/4.0),
        fT(0.9/4.0),
        fN(0.4/4.0)
    {
        std::cout << "Construct non-simple subst. matrix, CT" << std::endl;
        // Because of Bs treatment different rates in reads
        TValue fCRead = fC - fC *(1.0-fMeth) * fBsConversion;  // Compute expected frequency of Cs in Bs reads
        TValue fTRead = fT + fC *(1.0-fMeth) * fBsConversion;  // .. of Ts

        TValue pCC = fMatch*(1.0-fMeth)*(1.0-fBsConversion) + fMatch*fMeth;   // Prob. that C stays C: no substitution and (unmethylated and not converted or methylated)
        TValue pCT = fSubst*fTransi + fMatch*(1.0-fMeth)*fBsConversion;         // Prob. that C turns into T: substitution or unmethylated and converted
        // add prob. that methylated and converted?
        static TValue const tab[25] = {
                      //A,                      C,                                  G,                          T,                                                  N (ref?)
                      fMatch/(fA*fA),           fSubst*fTransv/(fA*fCRead),         fSubst*fTransi/(fA*fG),     fSubst*fTransv/(fA*fTRead),                         0,     // A
                      fSubst*fTransv/(fC*fA),   pCC/(fC*fCRead),                    fSubst*fTransv/(fC*fG),     pCT/(fC*fTRead),                                    0,     // C   Ref:C - read:T  ?
                      fSubst*fTransi/(fG*fA),   fSubst*fTransv/(fG*fCRead),         fMatch/(fG*fG),             fSubst*fTransv/(fG*fTRead),                         0,     // G
                      fSubst*fTransv/(fT*fA),   fSubst*fTransi/(fT*fCRead),         fSubst*fTransv/(fT*fG),     fMatch/(fT*fTRead),                                 0,     // T
                      fNtoX/(fN*fA),            fNtoX/(fN*fCRead),                  fNtoX/(fN*fG),              fNtoX/(fN*fTRead),                                  0
        };
        arrayCopy(tab, tab + 25, (*this).data_tab);
    }
};

template <typename TValue>
class BsSubstitutionMatrix<TValue, BsCaseGA, BsNonSimple >
{
public:
    // Static computation of the required array size.
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    // The data table.
    TValue data_tab[TAB_SIZE];

    TValue fMatch;
    TValue fSubst;
    TValue fNtoX;
    TValue fTransi;
    TValue fTransv;
    TValue fA;
    TValue fC;
    TValue fG;
    TValue fT;
    TValue fN;

    BsSubstitutionMatrix(TValue &fMeth, TValue &fBsConversion):
        fMatch(0.89/4.0),
        fSubst(0.09/12.0),
        fNtoX(0.02/4.0),
        fTransi(0.5/4.0),
        fTransv(0.5/8.0),
        fA(0.9/4.0),
        fC(0.9/4.0),
        fG(0.9/4.0),
        fT(0.9/4.0),
        fN(0.4/4.0)
    {
        // Because of Bs treatment different rates in reads
        TValue fGRead = fG - fG *(1.0-fMeth) * fBsConversion;  // Compute expected frequency of Cs in Bs reads
        TValue fARead = fA + fA *(1.0-fMeth) * fBsConversion;  // .. of Ts

        TValue pGG = fMatch*(1.0-fMeth)*(1.0-fBsConversion) + fMatch*fMeth;   // add prob. that methylated and converted?
        TValue pGA = fSubst*fTransi + fMatch*(1.0-fMeth)*fBsConversion;

        static TValue const tab[25] = {
                      //A,                                              C,                          G,                                              T,                          N (ref?)
                      fMatch/(fA*fARead),                               fSubst*fTransv/(fA*fC),     fSubst*fTransi/(fA*fGRead),                     fSubst*fTransv/(fA*fT),     0,     // A
                      fSubst*fTransv/(fC*fARead),                       fMatch/(fC*fC),             fSubst*fTransv/(fC*fGRead),                     fSubst*fTransi/(fC*fT),     0,     // C   Ref:C - read:T  ?
                      pGA/(fG*fARead),                                  fSubst*fTransv/(fG*fC),     pGG/(fG*fGRead),                                fSubst*fTransv/(fG*fT),     0,     // G
                      fSubst*fTransv/(fT*fARead),                       fSubst*fTransi/(fT*fC),     fSubst*fTransv/(fT*fGRead),                     fMatch/(fT*fT),             0,     // T
                      fNtoX/(fN*fARead),                                fNtoX/(fN*fC),              fNtoX/(fN*fGRead),                              fNtoX/(fN*fT),              0
        };
        arrayCopy(tab, tab + 25, (*this).data_tab);
    }
};


}

#endif
