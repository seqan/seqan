/*==========================================================================



==========================================================================*/

#ifndef APPS_BS_TOOLS_CASBAR_ALPHABETS_H_
#define APPS_BS_TOOLS_CASBAR_ALPHABETS_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

//#include "snp_meth_store.h"
//#include "meths.h"


namespace seqan
{


struct DnaM_ {};
typedef SimpleType<unsigned char, DnaM_> DnaM;

template <>
struct ValueSize<DnaM>
{
    typedef uint8_t Type;
    static const Type VALUE = 6;
};

template <>
struct BitsPerValue<DnaM>
{
    typedef uint8_t Type;
    static const Type VALUE = 3;
};




// ============================================================================
// Translate Tables
// ============================================================================

template <typename T = void>
struct TranslateTableDnaMToChar_
{
    static char const VALUE[6];
};

template <typename T>
char const TranslateTableDnaMToChar_<T>::VALUE[6] = {'A', 'C', 'G', 'T', 'D', 'H'};

template <typename T = void>
struct TranslateTableCharToDnaM_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableCharToDnaM_<T>::VALUE[256] =
{
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

    0,   0,   0,   1,   4,   0,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0, //4
//   ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,

    0,   0,   0,   0,   3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//  P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,

    0,   0,   0,   1,   4,   0,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

    0,   0,   0,   0,   3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,

    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

template <typename T = void>
struct TranslateTableDnaMToDna_
{
    static char const VALUE[6];
};

template <typename T>
char const TranslateTableDnaMToDna_<T>::VALUE[6] =
{
    0,   1,   2,   3,   1,   2      //
};

template <typename T = void>
struct TranslateTableByteToDnaM_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableByteToDnaM_<T>::VALUE[256] =
{
    0,   1,   2,   3,   4,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

// ============================================================================
// Assignment / Conversion Functions
// ============================================================================

inline void assign(char & c_target,
                   DnaM const & source)
{
    c_target = TranslateTableDnaMToChar_<>::VALUE[source.value];
}

template <>
struct CompareType<DnaM, uint8_t>
{
    typedef DnaM Type;
};

inline void assign(DnaM & target, uint8_t c_source)
{
    target.value = TranslateTableByteToDnaM_<>::VALUE[c_source];
}

template <>
struct CompareType<DnaM, char>
{
    typedef DnaM Type;
};

inline void assign(DnaM & target, char c_source)
{
    target.value = TranslateTableCharToDnaM_<>::VALUE[(unsigned char) c_source];
}


template <>
struct CompareType<DnaM, Dna>
{
    typedef Dna Type;
};

inline void assign(DnaM & target, Dna const & c_source)
{
    target.value = c_source.value;
}

inline void assign(Dna & target, DnaM const & c_source)
{
    target.value = TranslateTableDnaMToDna_<>::VALUE[(unsigned char)c_source];
;
}

////////////////////////////////////////////////////////////////////////////
// DnaMR (for profile: A, C, G, T, (top strand)
//                     a, c, g, t, (bottom strand)
//                     R (ref., stores ord value), V (count mapped based)
//                     X (gaps forward)
//                    ( Y (gaps reverse), implicit )
////////////////////////////////////////////////////////////////////////////

// TODO get rid of 'R' and 'V'

struct DnaMR_ {};
typedef SimpleType<unsigned char, DnaMR_> DnaMR;

template <>
struct ValueSize<DnaMR>
{
    typedef uint8_t Type;
    static const Type VALUE = 11;
};

template <>
struct BitsPerValue<DnaMR>
{
    typedef uint8_t Type;
    static const Type VALUE = 4;
};


// ============================================================================
// Translate Tables
// ============================================================================

template <typename T = void>
struct TranslateTableDnaMRToChar_
{
    static char const VALUE[11];
};

template <typename T>
char const TranslateTableDnaMRToChar_<T>::VALUE[11] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'R', 'V', 'X'};



// ProfileChar<DnaMR, double>
inline bool
empty(ProfileChar<DnaMR, double> const & source)
{
    for (unsigned i = 0; i < 8; ++i)
        if (source.count[i] > 0.00001)
            return false;
    if ((int)source.count[8] >= 0 && (int)source.count[8] <= 3) return false;
    return true;
}

}
#endif
