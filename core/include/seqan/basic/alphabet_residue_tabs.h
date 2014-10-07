// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Conversion tables for residue SimpleType specializations.
//
// Dna and Dna5 share their tables for conversion to char.  So do Rna and
// Rna5.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_RESIDUE_TABS_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_RESIDUE_TABS_H_


namespace seqan {

// --------------------------------------------------------------------------
// Dna and Dna5
// --------------------------------------------------------------------------

template <typename T = void>
struct TranslateTableDna5ToAscii_
{
    static char const VALUE[5];
};

template <typename T>
char const TranslateTableDna5ToAscii_<T>::VALUE[5] = {'A', 'C', 'G', 'T', 'N'};


template <typename T = void>
struct TranslateTableDna5ToIupac_
{
    static char const VALUE[5];
};

template <typename T>
char const TranslateTableDna5ToIupac_<T>::VALUE[5] = {0x01, 0x02, 0x04, 0x08, 0x0f};

template <typename T = void>
struct TranslateTableAsciiToDna_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableAsciiToDna_<T>::VALUE[256] = 
{
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

    0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //4
//   ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

    0,   0,   0,   0,   3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//  P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

    0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //6
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
struct TranslateTableAsciiToDna5_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableAsciiToDna5_<T>::VALUE[256] = 
{
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //0
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //1
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //2
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //3

    4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //4
//   ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

    4,   4,   4,   4,   3,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //5
//  P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

    4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

    4,   4,   4,   4,   3,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //8
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //9
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //10
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //11
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //12
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //13
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //14
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4  //15
};

template <typename T = void>
struct TranslateTableByteToDna_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableByteToDna_<T>::VALUE[256] = 
{
    0,   1,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
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

template <typename T = void>
struct TranslateTableByteToDna5_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableByteToDna5_<T>::VALUE[256] =
{
    0,   1,   2,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //0
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //1
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //2
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //3
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //4
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //5
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //6
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //7
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //8
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //9
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //10
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //11
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //12
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //13
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //14
    4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4  //15
};

// --------------------------------------------------------------------------
// Rna and Rna5
// --------------------------------------------------------------------------

template <typename T = void>
struct TranslateTableRna5ToAscii_
{
    static char const VALUE[5];
};

template <typename T>
char const TranslateTableRna5ToAscii_<T>::VALUE[5] = {'A', 'C', 'G', 'U', 'N'}; 

// other tables identical to Dna(5)

// --------------------------------------------------------------------------
// Iupac
// --------------------------------------------------------------------------

template <typename T = void>
struct TranslateTableIupacToAscii_
{
    static char const VALUE[16];
};

template <typename T>
char const TranslateTableIupacToAscii_<T>::VALUE[16] = 
{        //TGCA
    '=', //0000=0 = or U
    'A', //0001=1
    'C', //0010=2
    'M', //0011=3 AC
    'G', //0100=4
    'R', //0101=5 AG (purine)
    'S', //0110=6 CG
    'V', //0111=7 non-T
    'T', //1000=8
    'W', //1001=9 TA
    'Y', //1010=A TC (pyrimidine)
    'H', //1011=B not-G
    'K', //1100=C TG
    'D', //1101=D not-C
    'B', //1110=E non-A
    'N'  //1111=F any
};

template <typename T = void>
struct TranslateTableIupacToDna_
{
    static char const VALUE[16];
};

template <typename T>
char const TranslateTableIupacToDna_<T>::VALUE[16] = 
{      //TGCA
    0, //0000=0 = or U
    0, //0001=1
    1, //0010=2
    0, //0011=3 AC
    2, //0100=4
    0, //0101=5 AG (purine)
    1, //0110=6 CG
    0, //0111=7 non-T
    3, //1000=8
    0, //1001=9 TA
    1, //1010=A TC (pyrimidine)
    0, //1011=B not-G
    2, //1100=C TG
    0, //1101=D not-C
    1, //1110=E non-A
    0  //1111=F any
};

template <typename T = void>
struct TranslateTableIupacToDna5_
{
    static char const VALUE[16];
};

template <typename T>
char const TranslateTableIupacToDna5_<T>::VALUE[16] = 
{      //TGCA
    3, //0000=0 = or U
    0, //0001=1
    1, //0010=2
    4, //0011=3 AC
    2, //0100=4
    4, //0101=5 AG (purine)
    4, //0110=6 CG
    4, //0111=7 non-T
    3, //1000=8
    4, //1001=9 TA
    4, //1010=A TC (pyrimidine)
    4, //1011=B not-G
    4, //1100=C TG
    4, //1101=D not-C
    4, //1110=E non-A
    4  //1111=F any
};

template <typename T = void>
struct TranslateTableAsciiToIupac_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableAsciiToIupac_<T>::VALUE[256] = 
{
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    //                                                                =
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,   0,  15,  15,
    //    A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,
    15,   1,  14,   2,  13,  15,  15,   4,  11,  15,  15,  12,  15,   3,  15,  15,
    //    Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z
    15,  15,   5,   6,   8,   0,   7,   9,  15,  10,  15,  15,  15,  15,  15,  15,
    //    a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,
    15,   1,  14,   2,  13,  15,  15,   4,  11,  15,  15,  12,  15,   3,  15,  15,
    //    q,   r,   s,   t,   u,   v,   w,   x,   y,   z
    15,  15,   5,   6,   8,   0,   7,   9,  15,  10,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15
};

template <typename T = void>
struct TranslateTableByteToIupac_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableByteToIupac_<T>::VALUE[256] = 
{
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15, //0
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //1
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //2
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //3
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //4
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //5
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //6
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //7
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //8
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //9
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //10
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //11
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //12
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //13
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //14
    15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15  //15
};

// --------------------------------------------------------------------------
// Amino Acid
// --------------------------------------------------------------------------

template <typename T = void>
struct TranslateTableAAToAscii_
{
    static char const VALUE[24];
};
template <typename T>
char const TranslateTableAAToAscii_<T>::VALUE[24] = 
{
    'A', // 0 Ala Alanine                 
    'R', // 1 Arg Arginine                
    'N', // 2 Asn Asparagine              
    'D', // 3 Asp Aspartic Acid           
    'C', // 4 Cys Cystine                 
    'Q', // 5 Gln Glutamine               
    'E', // 6 Glu Glutamic Acid           
    'G', // 7 Gly Glycine                 
    'H', // 8 His Histidine               
    'I', // 9 Ile Isoleucine              
    'L', //10 Leu Leucine                 
    'K', //11 Lys Lysine                  
    'M', //12 Met Methionine              
    'F', //13 Phe Phenylalanine           
    'P', //14 Pro Proline                 
    'S', //15 Ser Serine                  
    'T', //16 Thr Threonine               
    'W', //17 Trp Tryptophan              
    'Y', //18 Tyr Tyrosine                
    'V', //19 Val Valine                  
    'B', //20 Aspartic Acid, Asparagine   
    'Z', //21 Glutamic Acid, Glutamine    
    'X', //22 Unknown                     
    '*'  //23 Terminator                  
};

template <typename T = void>
struct TranslateTableAsciiToAA_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableAsciiToAA_<T>::VALUE[256] = 
{
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //0
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //1
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  23,  22,  22,  22,  22,  22, //2
//                                                     *    
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //3
    22,   0,  20,   4,   3,   6,  13,   7,   8,   9,  22,  11,  10,  12,   2,  22, //4
//    ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,

    14,   5,   1,  15,  16,  22,  19,  17,  22,  18,  21,  22,  22,  22,  22,  22, //5
//   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    ,

    22,   0,  20,   4,   3,   6,  13,   7,   8,   9,  22,  11,  10,  12,   2,  22, //6
//    ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

    14,   5,   1,  15,  16,  22,  19,  17,  22,  18,  21,  22,  22,  22,  22,  22, //7
//   p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,    ,

    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //8
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //9
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //10
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //11
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //12
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //13
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //14
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22  //15
};

template <typename T = void>
struct TranslateTableByteToAA_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableByteToAA_<T>::VALUE[256] = 
{
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, //0
    16,  17,  18,  19,  20,  21,  22,  23,  22,  22,  22,  22,  22,  22,  22,  22, //1
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //2
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //3
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //4
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //5
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //6
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //7
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //8
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //9
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //10
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //11
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //12
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //13
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //14
    22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22  //15
};

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_BASIC_ALPHABET_RESIDUE_TABS_H_
