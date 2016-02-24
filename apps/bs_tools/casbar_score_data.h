#ifndef APPS_BS_TOOLS_CASBAR_SCORE_DATA_H__
#define APPS_BS_TOOLS_CASBAR_SCORE_DATA_H__

namespace seqan {

struct Bs_;
typedef Tag<Bs_>   Bs;


// normalized for columns
template <typename TValue = double, typename TSpec = BsNonSimple>
struct SeqErrorFreqsN;

// Default sequencing error frequencies (substitutions only)
template <typename TValue>
struct SeqErrorFreqsN<TValue, BsNonSimple> {    // TODO -> BsNonSimple
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {

                        // Column sum must be 1   (in the current model)
        static TValue const _data[TAB_SIZE] = {
                // To A, C, G, T, N (ref?)
                        0,            (0.14/0.38),   (0.05/0.16),    (0.05/0.21),    0,    // From A
                       (0.13/0.25),    0,            (0.02/0.16),    (0.04/0.21),    0,    // C
                       (0.04/0.25),   (0.08/0.38),    0,             (0.12/0.21),    0,
                       (0.08/0.25),   (0.15/0.38),   (0.09/0.16),     0,             0,
                       0,              0,             0,              0,             0      // N?
        };
        return _data;
    }
};

// Simple sequencing error frequencies (substitutions only)
// [ordValue(realBase) * 5 + ordValue(observed base)]
template <typename TValue>
struct SeqErrorFreqsN<TValue, BsSimple> {
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
struct ReadBaseFreqs;
// Corresponding to original top and bottom
// strands
template <typename TValue>
struct ReadBaseFreqs<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {

        TValue f = 0.9/4.0;
        static TValue const _data[TAB_SIZE] = {f,    f,    f,    f,   0.4/4.0};    // TODO bs case
        return _data;
    }
};

template <typename TValue>
struct ReadBaseFreqs<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {

        TValue f = 0.9/4.0;
        static TValue const _data[TAB_SIZE] = {f,    f,    f,    f,   0.4/4.0};    // TODO bs case
        return _data;
    }
};


template <typename TValue = double, typename TSpec = BsNonSimple>
struct RefBaseFreqs;

template <typename TValue>
struct RefBaseFreqs<TValue, BsNonSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {

        TValue f = 0.9/4.0;
        static TValue const _data[TAB_SIZE] = {f,    f,    f,    f,   0.4/4.0};
        return _data;
    }
};

template <typename TValue>
struct RefBaseFreqs<TValue, BsSimple> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline TValue const * getData() {

        TValue f = 0.9/4.0;
        static TValue const _data[TAB_SIZE] = {f,    f,    f,    f,   0.4/4.0};
        return _data;
    }
};



}

#endif
