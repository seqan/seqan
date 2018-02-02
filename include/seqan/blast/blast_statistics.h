// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// bitScore and e-Value calculation
// ==========================================================================
// contains empiric values from blast_stat.c as well as small portions of
// code from blast_stat.c written by Tom Madden and released into the public
// domain by the NCBI.
// ==========================================================================

#ifndef __BLAST_STATISTICS_H__
#define __BLAST_STATISTICS_H__

#ifndef NCBI_INT2_MAX
#define NCBI_INT2_MAX    32767
#endif

namespace seqan
{

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class KarlinAltschulValues
// ----------------------------------------------------------------------------

//NOTE(h-2): removed to clean up interface
/*
 * @class KarlinAltschulValues
 * @headerfile <seqan/blast.h>
 * @brief An object that holds Statistical parameters
 * @signature struct KarlinAltschulValues<TScore, TSpec> { ... };
 *
 * This is an adapter around the SeqAn scoring scheme, for use with BLAST and
 * other e-value statistics. If you only use high-level Blast IO, you do not
 * have to worry about this, except that you have to know that not every scoring
 * scheme supported by SeqAn and not every combination of gap scores is supported.
 * [@link BlastIOContext#setBlastScoringScheme @endlink/@link BlastIOContext#setScoringScheme @endlink will fail if
 * a non-supported scoringScheme is passed]
 *
 * @section Details
 *
 * This data structure holds Karlin-Altschul statistical scoring
 * parameters (&lambda;, &Kappa;, H, &alpha;, &beta;, &theta;) for the given
 * scoring scheme. These values are empirically derived. The ones in
 * SeqAn are imported from the NCBI blast source code (blast-2.2.26).
 * For an explanation of the role of &lambda; and &Kappa; please see these
 * slides: <a href="https://www.cs.umd.edu/class/fall2011/cmsc858s/Local_Alignment_Statistics.pdf">
 * https://www.cs.umd.edu/class/fall2011/cmsc858s/Local_Alignment_Statistics.pdf</a>.
 *
 * <b>Warning</b>, not every scoring scheme supported by SeqAn and not every
 * combination of gap scores is supported.
 */

/// GENERIC
template <typename TScore, typename TSpec = void>
struct KarlinAltschulValues
{
    typedef uint8_t TSize;

    /* statics */
    static constexpr TSize nParams = 0;
    static constexpr TSize nParamSets = 0;
    static double const VALUE[1][8]; // size 8 to silence warnings
};

template <typename TScore, typename TSpec>
double const KarlinAltschulValues<TScore, TSpec>::VALUE[1][8] =
{
    {0,0,0,0,0,0,0,0}
};

/// BLOSUM30
// not implemented in BLAST

/// BLOSUM45
template <typename TSpec>
struct KarlinAltschulValues<Blosum45, TSpec>
{
    typedef uint8_t TSize;

    /* statics */
    static constexpr TSize nParams = 8;
    static constexpr TSize nParamSets = 16;
    static double const VALUE[nParamSets][nParams];
};

template <typename TSpec>
double const KarlinAltschulValues<Blosum45, TSpec>::VALUE
  [KarlinAltschulValues<Blosum45, TSpec>::nParamSets]
  [KarlinAltschulValues<Blosum45, TSpec>::nParams] =
{
    {(double) NCBI_INT2_MAX, (double) NCBI_INT2_MAX,
            (double) NCBI_INT2_MAX, 0.2291, 0.0924, 0.2514, 0.9113, -5.7},
    {13, 3, (double) NCBI_INT2_MAX, 0.207,  0.049,  0.14,   1.5, -22},
    {12, 3, (double) NCBI_INT2_MAX, 0.199,  0.039,  0.11,   1.8, -34},
    {11, 3, (double) NCBI_INT2_MAX, 0.190,  0.031,  0.095,  2.0, -38},
    {10, 3, (double) NCBI_INT2_MAX, 0.179,  0.023,  0.075,  2.4, -51},
    {16, 2, (double) NCBI_INT2_MAX, 0.210,  0.051,  0.14,   1.5, -24},
    {15, 2, (double) NCBI_INT2_MAX, 0.203,  0.041,  0.12,   1.7, -31},
    {14, 2, (double) NCBI_INT2_MAX, 0.195,  0.032,  0.10,   1.9, -36},
    {13, 2, (double) NCBI_INT2_MAX, 0.185,  0.024,  0.084,  2.2, -45},
    {12, 2, (double) NCBI_INT2_MAX, 0.171,  0.016,  0.061,  2.8, -65},
    {19, 1, (double) NCBI_INT2_MAX, 0.205,  0.040,  0.11,   1.9, -43},
    {18, 1, (double) NCBI_INT2_MAX, 0.198,  0.032,  0.10,   2.0, -43},
    {17, 1, (double) NCBI_INT2_MAX, 0.189,  0.024,  0.079,  2.4, -57},
    {16, 1, (double) NCBI_INT2_MAX, 0.176,  0.016,  0.063,  2.8, -67},
};


/// BLOSUM50
// not implemented in seqan

/// BLOSUM62
template <typename TSpec>
struct KarlinAltschulValues<Blosum62, TSpec>
{
    typedef uint8_t TSize;

    /* statics */
    static constexpr TSize nParams = 8;
    static constexpr TSize nParamSets = 12;
    static double const VALUE[nParamSets][nParams];
};

template <typename TSpec>
double const KarlinAltschulValues<Blosum62, TSpec>::VALUE
  [KarlinAltschulValues<Blosum62, TSpec>::nParamSets]
  [KarlinAltschulValues<Blosum62, TSpec>::nParams] =
{
        {(double) NCBI_INT2_MAX, (double) NCBI_INT2_MAX,
                (double) NCBI_INT2_MAX, 0.3176, 0.134, 0.4012, 0.7916, -3.2},
        {11, 2, (double) NCBI_INT2_MAX, 0.297,  0.082, 0.27,   1.1,    -10} ,
        {10, 2, (double) NCBI_INT2_MAX, 0.291,  0.075, 0.23,   1.3,    -15} ,
        { 9, 2, (double) NCBI_INT2_MAX, 0.279,  0.058, 0.19,   1.5,    -19} ,
        { 8, 2, (double) NCBI_INT2_MAX, 0.264,  0.045, 0.15,   1.8,    -26} ,
        { 7, 2, (double) NCBI_INT2_MAX, 0.239,  0.027, 0.10,   2.5,    -46} ,
        { 6, 2, (double) NCBI_INT2_MAX, 0.201,  0.012, 0.061,  3.3,    -58} ,
        {13, 1, (double) NCBI_INT2_MAX, 0.292,  0.071, 0.23,   1.2,    -11} ,
        {12, 1, (double) NCBI_INT2_MAX, 0.283,  0.059, 0.19,   1.5,    -19} ,
        {11, 1, (double) NCBI_INT2_MAX, 0.267,  0.041, 0.14,   1.9,    -30} ,
        {10, 1, (double) NCBI_INT2_MAX, 0.243,  0.024, 0.10,   2.5,    -44} ,
        { 9, 1, (double) NCBI_INT2_MAX, 0.206,  0.010, 0.052,  4.0,    -87} ,
};
/// BLOSUM80
template <typename TSpec>
struct KarlinAltschulValues<Blosum80, TSpec>
{
    typedef uint8_t TSize;

    /* statics */
    static constexpr TSize nParams = 8;
    static constexpr TSize nParamSets = 10;
    static double const VALUE[nParamSets][nParams];
};

template <typename TSpec>
double const KarlinAltschulValues<Blosum80, TSpec>::VALUE
 [KarlinAltschulValues<Blosum80, TSpec>::nParamSets]
 [KarlinAltschulValues<Blosum80, TSpec>::nParams] =
{
    {(double) NCBI_INT2_MAX, (double) NCBI_INT2_MAX,
            (double) NCBI_INT2_MAX, 0.3430,  0.177, 0.6568, 0.5222, -1.6},
    {25, 2, (double) NCBI_INT2_MAX, 0.342,   0.17,  0.66,   0.52,   -1.6},
    {13, 2, (double) NCBI_INT2_MAX, 0.336,   0.15,  0.57,   0.59,   -3},
    { 9, 2, (double) NCBI_INT2_MAX, 0.319,   0.11,  0.42,   0.76,   -6},
    { 8, 2, (double) NCBI_INT2_MAX, 0.308,   0.090, 0.35,   0.89,   -9},
    { 7, 2, (double) NCBI_INT2_MAX, 0.293,   0.070, 0.27,   1.1,   -14},
    { 6, 2, (double) NCBI_INT2_MAX, 0.268,   0.045, 0.19,   1.4,   -19},
    {11, 1, (double) NCBI_INT2_MAX, 0.314,   0.095, 0.35,   0.90,   -9},
    {10, 1, (double) NCBI_INT2_MAX, 0.299,   0.071, 0.27,   1.1,   -14},
    { 9, 1, (double) NCBI_INT2_MAX, 0.279,   0.048, 0.20,   1.4,   -19},
};

/// BLOSUM90
// not implemented in seqan

/// PAM30
// not implemented in seqan

/// PAM40
// not implemented in blast

/// PAM70
// not implemented in seqan

/// PAM120
// not implemented in blast

/// PAM200
// not implemented in blast

/// PAM250
template <typename TSpec>
struct KarlinAltschulValues<Pam250, TSpec>
{
    typedef uint8_t TSize;

    /* statics */
    static constexpr TSize nParams = 8;
    static constexpr TSize nParamSets = 16;
    static double const VALUE[nParamSets][nParams];
};

template <typename TSpec>
double const KarlinAltschulValues<Pam250, TSpec>::VALUE
  [KarlinAltschulValues<Pam250, TSpec>::nParamSets]
  [KarlinAltschulValues<Pam250, TSpec>::nParams] =
{
    {(double) NCBI_INT2_MAX, (double) NCBI_INT2_MAX,
            (double) NCBI_INT2_MAX, 0.2252, 0.0868, 0.2223, 0.98,  -5.0},
    {15, 3, (double) NCBI_INT2_MAX, 0.205,  0.049,  0.13,   1.6,  -23},
    {14, 3, (double) NCBI_INT2_MAX, 0.200,  0.043,  0.12,   1.7,  -26},
    {13, 3, (double) NCBI_INT2_MAX, 0.194,  0.036,  0.10,   1.9,  -31},
    {12, 3, (double) NCBI_INT2_MAX, 0.186,  0.029,  0.085,  2.2,  -41},
    {11, 3, (double) NCBI_INT2_MAX, 0.174,  0.020,  0.070,  2.5,  -48},
    {17, 2, (double) NCBI_INT2_MAX, 0.204,  0.047,  0.12,   1.7,  -28},
    {16, 2, (double) NCBI_INT2_MAX, 0.198,  0.038,  0.11,   1.8,  -29},
    {15, 2, (double) NCBI_INT2_MAX, 0.191,  0.031,  0.087,  2.2,  -44},
    {14, 2, (double) NCBI_INT2_MAX, 0.182,  0.024,  0.073,  2.5,  -53},
    {13, 2, (double) NCBI_INT2_MAX, 0.171,  0.017,  0.059,  2.9,  -64},
    {21, 1, (double) NCBI_INT2_MAX, 0.205,  0.045,  0.11,   1.8,  -34},
    {20, 1, (double) NCBI_INT2_MAX, 0.199,  0.037,  0.10,   1.9,  -35},
    {19, 1, (double) NCBI_INT2_MAX, 0.192,  0.029,  0.083,  2.3,  -52},
    {18, 1, (double) NCBI_INT2_MAX, 0.183,  0.021,  0.070,  2.6,  -60},
    {17, 1, (double) NCBI_INT2_MAX, 0.171,  0.014,  0.052,  3.3,  -86},
};

/// VTML200
// not implemented in blast

/// NUCLEOTIDES
template <typename TSpec>
struct KarlinAltschulValues<Score<int, Simple>, TSpec>
{
    typedef uint8_t TSize;

    /* statics */
    static constexpr TSize nParams = 11;
    static constexpr TSize nParamSets = 60;
    /*  0 = matchScore
        1 = mismatch cost
        2 = Gap opening cost,
        3 = Gap extension cost,
        4 = Lambda,
        5 = K,
        6 = H,
        7 = Alpha,
        8 = Beta,
        9 = Theta
        10 = only applicable to even score, uneven need to be rounded down */
    static double const VALUE[nParamSets][nParams];
};

template <typename TSpec>
double const KarlinAltschulValues<Score<int, Simple>, TSpec>::VALUE
 [KarlinAltschulValues<Score<int, Simple>, TSpec>::nParamSets]
 [KarlinAltschulValues<Score<int, Simple>, TSpec>::nParams] =
{
    { 1, 5, 0, 0, 1.39, 0.747, 1.38, 1.00,  0, 100, 0 } ,
    { 1, 5, 3, 3, 1.39, 0.747, 1.38, 1.00,  0, 100, 0 } ,

    { 1, 4, 0, 0, 1.383, 0.738, 1.36, 1.02,  0, 100, 0 } ,
    { 1, 4, 1, 2,  1.36,  0.67,  1.2,  1.1,  0,  98, 0 } ,
    { 1, 4, 0, 2,  1.26,  0.43, 0.90,  1.4, -1,  91, 0 } ,
    { 1, 4, 2, 1,  1.35,  0.61,  1.1,  1.2, -1,  98, 0 } ,
    { 1, 4, 1, 1,  1.22,  0.35, 0.72,  1.7, -3,  88, 0 } ,

    { 2, 7, 0, 0,  0.69, 0.73, 1.34, 0.515,  0, 100, 1 } ,
    { 2, 7, 2, 4,  0.68, 0.67,  1.2,  0.55,  0,  99, 1 } ,
    { 2, 7, 0, 4,  0.63, 0.43, 0.90,   0.7, -1,  91, 1 } ,
    { 2, 7, 4, 2, 0.675, 0.62,  1.1,   0.6, -1,  98, 1 } ,
    { 2, 7, 2, 2,  0.61, 0.35, 0.72,   1.7, -3,  88, 1 } ,

    { 1, 3, 0, 0, 1.374, 0.711, 1.31, 1.05,  0, 100, 0 } ,
    { 1, 3, 2, 2,  1.37,  0.70,  1.2,  1.1,  0,  99, 0 } ,
    { 1, 3, 1, 2,  1.35,  0.64,  1.1,  1.2, -1,  98, 0 } ,
    { 1, 3, 0, 2,  1.25,  0.42, 0.83,  1.5, -2,  91, 0 } ,
    { 1, 3, 2, 1,  1.34,  0.60,  1.1,  1.2, -1,  97, 0 } ,
    { 1, 3, 1, 1,  1.21,  0.34, 0.71,  1.7, -2,  88, 0 } ,

    { 2, 5, 0, 0, 0.675, 0.65,  1.1,  0.6, -1, 99, 1 },
    { 2, 5, 2, 4,  0.67, 0.59,  1.1,  0.6, -1, 98, 1 },
    { 2, 5, 0, 4,  0.62, 0.39, 0.78,  0.8, -2, 91, 1 },
    { 2, 5, 4, 2,  0.67, 0.61,  1.0, 0.65, -2, 98, 1 },
    { 2, 5, 2, 2,  0.56, 0.32, 0.59, 0.95, -4, 82, 1 },

    { 1, 2, 0, 0, 1.28, 0.46, 0.85, 1.5, -2, 96, 0 },
    { 1, 2, 2, 2, 1.33, 0.62,  1.1, 1.2,  0, 99, 0 },
    { 1, 2, 1, 2, 1.30, 0.52, 0.93, 1.4, -2, 97, 0 },
    { 1, 2, 0, 2, 1.19, 0.34, 0.66, 1.8, -3, 89, 0 },
    { 1, 2, 3, 1, 1.32, 0.57,  1.0, 1.3, -1, 99, 0 },
    { 1, 2, 2, 1, 1.29, 0.49, 0.92, 1.4, -1, 96, 0 },
    { 1, 2, 1, 1, 1.14, 0.26, 0.52, 2.2, -5, 85, 0 },

    { 2, 3, 0, 0,  0.55, 0.21, 0.46,  1.2, -5, 87, 1 },
    { 2, 3, 4, 4,  0.63, 0.42, 0.84, 0.75, -2, 99, 1 },
    { 2, 3, 2, 4, 0.615, 0.37, 0.72, 0.85, -3, 97, 1 },
    { 2, 3, 0, 4,  0.55, 0.21, 0.46,  1.2, -5, 87, 1 },
    { 2, 3, 3, 3, 0.615, 0.37, 0.68,  0.9, -3, 97, 1 },
    { 2, 3, 6, 2,  0.63, 0.42, 0.84, 0.75, -2, 99, 1 },
    { 2, 3, 5, 2, 0.625, 0.41, 0.78,  0.8, -2, 99, 1 },
    { 2, 3, 4, 2,  0.61, 0.35, 0.68,  0.9, -3, 96, 1 },
    { 2, 3, 2, 2, 0.515, 0.14, 0.33, 1.55, -9, 81, 1 },

    { 3, 4, 6, 3, 0.389, 0.25, 0.56, 0.7, -5, 95, 0},
    { 3, 4, 5, 3, 0.375, 0.21, 0.47, 0.8, -6, 92, 0},
    { 3, 4, 4, 3, 0.351, 0.14, 0.35, 1.0, -9, 86, 0},
    { 3, 4, 6, 2, 0.362, 0.16, 0.45, 0.8, -4, 88, 0},
    { 3, 4, 5, 2, 0.330, 0.092, 0.28, 1.2, -13, 81, 0},
    { 3, 4, 4, 2, 0.281, 0.046, 0.16, 1.8, -23, 69, 0},

    { 4, 5, 0, 0, 0.22, 0.061, 0.22, 1.0, -15, 74, 0 },
    { 4, 5, 6, 5, 0.28,  0.21, 0.47, 0.6 , -7, 93, 0 },
    { 4, 5, 5, 5, 0.27,  0.17, 0.39, 0.7,  -9, 90, 0 },
    { 4, 5, 4, 5, 0.25,  0.10, 0.31, 0.8, -10, 83, 0 },
    { 4, 5, 3, 5, 0.23, 0.065, 0.25, 0.9, -11, 76, 0 },

    { 1, 1, 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99, 0 },
    { 1, 1, 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97, 0 },
    { 1, 1, 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92, 0 },
    { 1, 1, 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72, 0 },
    { 1, 1, 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98, 0 },
    { 1, 1, 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96, 0 },
    { 1, 1, 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90, 0 },

    { 3, 2,  5, 5, 0.208, 0.030, 0.072, 2.9, -47, 77, 0},
    { 5, 4, 10, 6, 0.163, 0.068, 0.16,  1.0, -19, 85, 0 },
    { 5, 4,  8, 6, 0.146, 0.039, 0.11,  1.3, -29, 76, 0 }

};

// ----------------------------------------------------------------------------
// Class BlastScoringScheme
// ----------------------------------------------------------------------------

/*!
 * @class BlastScoringScheme
 * @brief Wrapper around a SeqAn scoring scheme
 * @signature struct BlastScoringScheme<TScoringScheme> { ... };
 * @headerfile <seqan/blast.h>
 *
 * This class wraps around a SeqAn scoring scheme and add functions for the retrievel of pre-computed
 * Karlin-Altschul-Values (required for e-value statistics) and it modifies the gap-scoring slightly.
 * If you use high-level IO (@link BlastReport @endlink, @link BlastTabular @endlink) you need not worry about the
 * specifics, however you should always make sure that the scoring scheme @link BlastScoringScheme#isValid @endlink
 * after you modify it (or check the return values of the set* functions).
 *
 * To retrieve a score-object compatible with regular SeqAn-functions, call
 * @link BlastScoringScheme#seqanScheme @endlink.
 *
 * @section Gap cost computation
 *
 * Gap-costs are computed differently in SeqAn and in BLAST.
 * Blast (and many other tools) compute scores of a stretch of gaps as
 * <tt>s = gO + n * gE</tt>
 * where gO is the gapOpen score, gE is the gap extend score and n ist the
 * total number of gap characters.
 *
 * SeqAn, however, normally computes them as as <tt>s = gO + (n-1) * gE</tt>.
 *
 * Effectively this results only in a different gap open score. To make dealing with this easier, there are
 * convenience functions for setting and getting the gap open score in blast convention:
 * @link BlastScoringScheme#scoreGapOpenBlast @endlink and @link BlastScoringScheme#setScoreGapOpenBlast @endlink.
 *
 * @section Karlin-Altschul values
 *
 * E-Value statistics require certain pre-computed constants. The BlastScoringScheme picks these automatically
 * whenever you modify it. Unfortunately there aren't parameters for all combinations of scoring schemes and gap
 * costs, so you have to check @link BlastScoringScheme#isValid @endlink after modifying it or verify the return
 * value of the set* functions which is bool here (instead of void).
 *
 * More details on the scoring parameters is in these slides:
 * <a href="https://www.cs.umd.edu/class/fall2011/cmsc858s/Local_Alignment_Statistics.pdf">
 * https://www.cs.umd.edu/class/fall2011/cmsc858s/Local_Alignment_Statistics.pdf</a>.
 *
 * The constants used in SeqAn were imported from the NCBI blast source code (blast-2.2.26)
 *
 * @tparam TScoringScheme A SeqAn Score type
 */

/* NOTE(h-2): Maybe this should specialize Score<> instead of being a wrapper.
 * This would however involve adding a new template parameter to score and
 * changing lots of code. Also it would become confusing as the shortcuts
 * like Blosum62 would no longer work...
 * TODO for SeqAn3: include this files implications when redesigning scoring schemes
 * and reduce complexity
 */
template <typename TScore>
struct BlastScoringScheme
{
    typedef typename Value<TScore>::Type TValue;
    typedef typename KarlinAltschulValues<TScore, void>::TSize TNumValues;

    /* related to ScoringScheme conversion */
    TValue _gapOpenBlast = 0;
    TScore _internalScheme;

    /* parameter selection */
    TNumValues parameterIndex = std::numeric_limits<TNumValues>::max();

    /* hacks for the dynamic matrix */
    double const * _m;
    TNumValues _nParams;
};

// ============================================================================
// Functions
// ============================================================================

/*!
 * @fn BlastScoringScheme#seqanScheme
 * @headerfile <seqan/blast.h>
 * @brief Retrieve a @link Score @endlink object for use with e.g. @link localAlignment @endlink.
 * @signature TScoringScheme const & seqanScheme(blastScoringScheme);
 * @return TScoringScheme_const_& the internal scoring object for use with other SeqAn functions.
 *
 * @section Remarks
 *
 * The reference returned is always const, because you should not modify the internal score object
 * directly. Use @link BlastScoringScheme @endlink's interface instead.
 */

template <typename TScore>
inline TScore const &
seqanScheme(BlastScoringScheme<TScore> const & bscheme)
{
    return bscheme._internalScheme;
}

// ----------------------------------------------------------------------------
// Function overloads for getScore*
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#scoreMatch
 * @headerfile <seqan/blast.h>
 * @brief The Match score (only for SimpleScore specialization).
 * @signature int scoreMatch(blastScoringScheme);
 */

inline int
scoreMatch(BlastScoringScheme<Score<int, Simple>> & scheme)
{
    return scoreMatch(scheme._internalScheme);
}

/*!
 * @fn BlastScoringScheme#scoreMismatch
 * @headerfile <seqan/blast.h>
 * @brief The Mismatch score (only for SimpleScore specialization).
 * @signature int scoreMismatch(blastScoringScheme);
 */

inline int
scoreMismatch(BlastScoringScheme<Score<int, Simple>> & scheme)
{
    return scoreMismatch(scheme._internalScheme);
}

/*!
 * @fn BlastScoringScheme#scoreGapOpen
 * @headerfile <seqan/blast.h>
 * @brief The Gap Open score in SeqAn convention.
 * @signature TValue scoreGapOpen(blastScoringScheme);
 */

template <typename TScore>
inline typename Value<TScore>::Type
scoreGapOpen(BlastScoringScheme<TScore> const & scheme)
{
    return scoreGapOpen(scheme._internalScheme);
}

/*!
 * @fn BlastScoringScheme#scoreGapOpenBlast
 * @headerfile <seqan/blast.h>
 * @brief The Gap Open score in Blast convention.
 * @signature TValue scoreGapOpenBlast(blastScoringScheme);
 */

template <typename TScore>
inline typename Value<TScore>::Type
scoreGapOpenBlast(BlastScoringScheme<TScore> const & scheme)
{
    return scheme._gapOpenBlast;
}

/*!
 * @fn BlastScoringScheme#scoreGapExtend
 * @headerfile <seqan/blast.h>
 * @brief The Gap Extend score.
 * @signature TValue scoreGapExtend(blastScoringScheme);
 */

template <typename TScore>
inline typename Value<TScore>::Type
scoreGapExtend(BlastScoringScheme<TScore> const & scheme)
{
    return scoreGapExtend(scheme._internalScheme);
}

// ----------------------------------------------------------------------------
// Function isValid
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#isValid
 * @headerfile <seqan/blast.h>
 * @brief check whether valid KarlinAltschulValues are selected
 * @signature bool isValid(blastScoringScheme);
 */

template <typename TScore>
inline bool
isValid(BlastScoringScheme<TScore> const & scheme)
{
    typedef KarlinAltschulValues<TScore> TKAValues;
    return scheme.parameterIndex != std::numeric_limits<typename TKAValues::TSize>::max();
}

// ----------------------------------------------------------------------------
// Function _selectSet
// ----------------------------------------------------------------------------

template <typename TScore, typename TKAValues>
inline bool
_selectSet(BlastScoringScheme<TScore> & scheme)
{
    for (typename TKAValues::TSize i = 0; i < TKAValues::nParamSets; ++i)
    {
        if ((TKAValues::VALUE[i][0] == -scoreGapOpenBlast(scheme)) &&
            (TKAValues::VALUE[i][1] == -scoreGapExtend(scheme)))
        {
            scheme.parameterIndex = i;
            return true;
        }
    }
    // no suitable adapter
    scheme.parameterIndex = std::numeric_limits<typename TKAValues::TSize>::max();
    return false;
}

template <typename TValue, typename TSpec>
inline bool
_selectSet(BlastScoringScheme<Score<TValue, ScoreMatrix<AminoAcid, TSpec>>> & scheme)
{
    using TScore = Score<TValue, ScoreMatrix<AminoAcid, TSpec>>;
    using TKAValues = KarlinAltschulValues<TScore>;
    return _selectSet<TScore, TKAValues>(scheme);
}

template <typename TValue>
inline bool
_selectSet(BlastScoringScheme<Score<TValue, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>> & scheme)
{
    using TScoreOrig = Score<TValue, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>;
    bool ret = false;
    impl::score::matrixTagDispatch(impl::score::MatrixTags(),
                                   getScoreMatrixId(seqanScheme(scheme)),
                                   [&] (auto const & tag)
    {
        using TScoreMod = Score<TValue, ScoreMatrix<AminoAcid, std::decay_t<decltype(tag)>>>;
        using TKAValues = KarlinAltschulValues<TScoreMod>;
        ret =  _selectSet<TScoreOrig, TKAValues>(scheme);
        // save some KAV data in scheme for retrievel without tag dispatching
        scheme._m = &TKAValues::VALUE[0][0];
        scheme._nParams = TKAValues::nParams;
    });
    return ret;
}

inline bool
_selectSet(BlastScoringScheme<Score<int, Simple>> & scheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;

    for (typename TKAValues::TSize i = 0; i < TKAValues::nParamSets; ++i)
    {
        if ((TKAValues::VALUE[i][0] ==  scoreMatch(scheme))  &&
            (TKAValues::VALUE[i][1] == -scoreMismatch(scheme)) &&
            (TKAValues::VALUE[i][2] == -scoreGapOpenBlast(scheme)) &&
            (TKAValues::VALUE[i][3] == -scoreGapExtend(scheme)))
        {
            scheme.parameterIndex = i;
            return true;
        }
    }
    // no suitable adapter
    scheme.parameterIndex = std::numeric_limits<typename TKAValues::TSize>::max();
    return false;
}

// ----------------------------------------------------------------------------
// Function overloads for setScore*
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#setScoreMatch
 * @brief Set match score.
 * @signature bool setScoreMatch(blastScoringScheme, value);
 * @param[in,out] blastScoringScheme The @link BlastScoringScheme @endlink to modify.
 * @param[in]     value              The new value.
 * @return isValid(blastScoringScheme) whether Karlin-Altschul values for the new scheme exist.
 */

inline bool
setScoreMatch(BlastScoringScheme<Score<int, Simple>> & scheme, int const val)
{
    setScoreMatch(scheme._internalScheme, val);
    return _selectSet(scheme);
}

/*!
 * @fn BlastScoringScheme#setScoreMismatch
 * @brief Set mismatch score.
 * @signature bool setScoreMismatch(blastScoringScheme, value);
 * @param[in,out] blastScoringScheme The @link BlastScoringScheme @endlink to modify.
 * @param[in]     value              The new value.
 * @return isValid(blastScoringScheme) whether Karlin-Altschul values for the new scheme exist.
 */

inline bool
setScoreMismatch(BlastScoringScheme<Score<int, Simple>> & scheme, int const val)
{
    setScoreMismatch(scheme._internalScheme, val);
    return _selectSet(scheme);
}

/*!
 * @fn BlastScoringScheme#setScoreGapOpenBlast
 * @headerfile <seqan/blast.h>
 * @brief Set gap open score in BLAST convention.
 * @signature bool setScoreGapOpenBlast(blastScoringScheme, value);
 * @param[in,out] blastScoringScheme The @link BlastScoringScheme @endlink to modify.
 * @param[in]     value              The new value.
 * @return isValid(blastScoringScheme) whether Karlin-Altschul values for the new scheme exist.
 */

template <typename TScore>
inline bool
setScoreGapOpenBlast(BlastScoringScheme<TScore> & scheme, typename Value<TScore>::Type const val)
{
    scheme._gapOpenBlast = val;
    setScoreGapOpen(scheme._internalScheme, scoreGapOpenBlast(scheme) + scoreGapExtend(scheme));

    return _selectSet(scheme);
}

/*!
 * @fn BlastScoringScheme#setScoreGapOpen
 * @brief Set gap open score in SeqAn convention.
 * @signature bool setScoreGapOpen(blastScoringScheme, value);
 * @param[in,out] blastScoringScheme The @link BlastScoringScheme @endlink to modify.
 * @param[in]     value              The new value.
 * @return isValid(blastScoringScheme) whether Karlin-Altschul values for the new scheme exist.
 */

template <typename TScore>
inline bool
setScoreGapOpen(BlastScoringScheme<TScore> & scheme, typename Value<TScore>::Type const val)
{
    return setScoreGapOpenBlast(scheme, val - scoreGapExtend(scheme));
}

/*!
 * @fn BlastScoringScheme#setScoreGapExtend
 * @brief Set gap extend score.
 * @signature bool setScoreGapExtend(blastScoringScheme, value);
 * @param[in,out] blastScoringScheme The @link BlastScoringScheme @endlink to modify.
 * @param[in]     value              The new value.
 * @return isValid(blastScoringScheme) whether Karlin-Altschul values for the new scheme exist.
 */

template <typename TScore>
inline bool
setScoreGapExtend(BlastScoringScheme<TScore> & scheme, typename Value<TScore>::Type const val)
{
    setScoreGapExtend(scheme._internalScheme, val);
    setScoreGapOpen(scheme._internalScheme, scoreGapOpenBlast(scheme) + scoreGapExtend(scheme));

    return _selectSet(scheme);
}

// ----------------------------------------------------------------------------
// Function getLambda
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#getLambda
 * @headerfile <seqan/blast.h>
 * @brief Get the &lambda; value of the Karlin-Altschul parameters.
 * @signature double getLambda(blastScoringScheme);
 */

inline double
getLambda(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>> const & scoringScheme)
{
    SEQAN_ASSERT(isValid(scoringScheme));
    return *(scoringScheme._m + scoringScheme.parameterIndex * scoringScheme._nParams + 3);
}

template <typename TMatrixSpec>
inline double
getLambda(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][3];
}

inline double
getLambda(BlastScoringScheme<Score<int, Simple>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][4];
}

// ----------------------------------------------------------------------------
// Function getKappa
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#getKappa
 * @headerfile <seqan/blast.h>
 * @brief Get the &Kappa; value of the Karlin-Altschul parameters.
 * @signature double getKappa(blastScoringScheme);
 */

inline double
getKappa(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>> const & scoringScheme)
{
    SEQAN_ASSERT(isValid(scoringScheme));
    return *(scoringScheme._m + scoringScheme.parameterIndex * scoringScheme._nParams + 4);
}

template <typename TMatrixSpec>
inline double
getKappa(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int,
        ScoreMatrix< AminoAcid, TMatrixSpec>>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][4];
}

inline double
getKappa(BlastScoringScheme<Score<int, Simple>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][5];
}

// ----------------------------------------------------------------------------
// Function getH
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#getH
 * @headerfile <seqan/blast.h>
 * @brief Get the H value of the Karlin-Altschul parameters.
 * @signature double getH(blastScoringScheme);
 */

inline double
getH(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>> const & scoringScheme)
{
    SEQAN_ASSERT(isValid(scoringScheme));
    return *(scoringScheme._m + scoringScheme.parameterIndex * scoringScheme._nParams + 5);
}

template <typename TMatrixSpec>
inline double
getH(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int,
        ScoreMatrix< AminoAcid, TMatrixSpec>>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][5];
}

inline double
getH(BlastScoringScheme<Score<int, Simple>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][6];
}

// ----------------------------------------------------------------------------
// Function getAlpha
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#getAlpha
 * @headerfile <seqan/blast.h>
 * @brief Get the &alpha; value of the Karlin-Altschul parameters.
 * @signature double getAlpha(blastScoringScheme);
 */

inline double
getAlpha(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>> const & scoringScheme)
{
    SEQAN_ASSERT(isValid(scoringScheme));
    return *(scoringScheme._m + scoringScheme.parameterIndex * scoringScheme._nParams + 6);
}


template <typename TMatrixSpec>
inline double
getAlpha(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][6];
}

inline double
getAlpha(BlastScoringScheme<Score<int, Simple>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][7];
}

// ----------------------------------------------------------------------------
// Function getBeta
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#getBeta
 * @headerfile <seqan/blast.h>
 * @brief Get the &beta; value of the Karlin-Altschul parameters.
 * @signature double getBeta(blastScoringScheme);
 */

inline double
getBeta(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable>>> const & scoringScheme)
{
    SEQAN_ASSERT(isValid(scoringScheme));
    return *(scoringScheme._m + scoringScheme.parameterIndex * scoringScheme._nParams + 7);
}

template <typename TMatrixSpec>
inline double
getBeta(BlastScoringScheme<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][7];
}

inline double
getBeta(BlastScoringScheme<Score<int, Simple>> const & scoringScheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    SEQAN_ASSERT(isValid(scoringScheme));
    return TKAValues::VALUE[scoringScheme.parameterIndex][8];
}

// ----------------------------------------------------------------------------
// Function _lengthAdjustment
// ----------------------------------------------------------------------------

// Computes the length adjustment for E-value computation
// Based on the public domain NCBI BLAST code by Tom Madden.
template <typename TSize, typename TSize2, typename TScore>
inline TSize
_lengthAdjustment(TSize     const & dbLength,
                  TSize2    const & queryLength,
                  BlastScoringScheme<TScore> const & scheme)
{
    const double lambda = getLambda(scheme);
    const double K      = getKappa(scheme);
    const double logK   = std::log(K);
    const double alpha  = getAlpha(scheme);
    const double alphaByLambda = alpha/lambda;
    const double beta   = getBeta(scheme);
    const TSize maxIterations = 20;

    double n = (double)dbLength;
    double m = (double)queryLength;
    double totalLen;

    double val = 0, val_min = 0, val_max;
    bool converged = false;

     /* Choose val_max to be the largest nonnegative value that satisfies
      *    K * (m - val) * (n - N * val) > max(m,n)
      * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */

    { // scope of mb, and c, the coefficients in the quadratic formula (the variable mb is -b, a=1 ommited)
        double mb = m + n;
        double c  = n * m - std::max(m, n) / K;

        if(c < 0) {
            return 0;
        } else {
            val_max = 2 * c / (mb + std::sqrt(mb * mb - 4 * c));
        }
    } // end scope of mb and c

    for(TSize i = 1; i <= maxIterations; i++) {
        totalLen = (m - val) * (n - val);
        double val_new  = alphaByLambda * (logK + std::log(totalLen)) + beta;  // proposed next value of val
        if(val_new >= val) { // val is no bigger than the true fixed point
            val_min = val;
            if(val_new - val_min <= 1.0) {
                converged = true;
                break;
            }
            if(val_min == val_max) { // There are no more points to check
                break;
            }
        } else { // val is greater than the true fixed point
            val_max = val;
        }
        if(val_min <= val_new && val_new <= val_max) { // ell_new is in range. Accept it.
            val = val_new;
        } else { // else val_new is not in range. Reject it.
            val = (i == 1) ? val_max : (val_min + val_max) / 2;
        }
    }

    if(converged) { // the iteration converged
        // If val_fixed is the (unknown) true fixed point, then we wish to set lengthAdjustment to floor(val_fixed).
        // We assume that floor(val_min) = floor(val_fixed)
        return (TSize) val_min;

        // But verify that ceil(val_min) != floor(val_fixed)
        val = std::ceil(val_min);
        if( val <= val_max ) {
          totalLen = (m - val) * (n - val);
          if(alphaByLambda * (logK + std::log(totalLen)) + beta >= val) {
            // ceil(val_min) == floor(val_fixed)
            return (TSize) val;
          }
        }
    } else { // the iteration did not converge
        // Use the best value seen so far.
        return (TSize) val_min;
    }
}

// ----------------------------------------------------------------------------
// Function computeAlignmentStats
// ----------------------------------------------------------------------------

/*!
 * @fn BlastMatch#computeAlignmentStats
 * @headerfile <seqan/blast.h>
 * @brief Compute the @link BlastMatch::alignStats @endlink member of a @link BlastMatch @endlink.
 * @signature void computeAlignmentStats(blastMatch, context);
 *
 * @param[in,out]   blastMatch  A @link BlastMatch @endlink that has a valid align member.
 * @param[in]       context     A @link BlastIOContext @endlink with parameters and buffers.
 *
 * The alignRow-members (@link BlastMatch::alignRow0 @endlink, @link BlastMatch::alignRow1 @endlink) are used as 
 * in-parameter to compute the @link BlastMatch::alignStats @endlink member. This includes the raw
 * score, amount of gaps, mismatches et cetera. This is a prerequisite for printing a match to file or computing it's
 * e-value.
 */

template <typename TBlastMatch,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
computeAlignmentStats(TBlastMatch & match,
                      BlastIOContext<TScore, p, h> const & context)
{
    computeAlignmentStats(match.alignStats, match.alignRow0, match.alignRow1, seqanScheme(context.scoringScheme));
}

// ----------------------------------------------------------------------------
// Function computeBitScore
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringScheme#computeBitScore
 * @brief Compute the bit-score for a given raw score.
 * @signature double computeBitScore(rawScore, blastScoringScheme);
 *
 * @param[in]   rawScore            The raw score of the alignment.
 * @param[in]   blastScoringScheme The @link BlastScoringScheme @endlink.
 *
 * @return double The bit-score computed.
 * @headerfile <seqan/blast.h>
 */

template <typename TScore>
inline double
computeBitScore(double const rawScore, BlastScoringScheme<TScore> const & scheme)
{
    return (getLambda(scheme) * rawScore - std::log(getKappa(scheme))) / std::log(2);
}

/*!
 * @fn BlastMatch#computeBitScore
 * @brief Compute the bit-score for a @link BlastMatch @endlink.
 * @signature double computeBitScore(blastMatch, context);
 *
 * @param[in,out]   blastMatch  A @link BlastMatch @endlink that has a valid align member.
 * @param[in]       context     A @link BlastIOContext @endlink with parameters and buffers.
 *
 * @return double blastMatch.@link BlastMatch::bitScore @endlink after computation
 * @headerfile <seqan/blast.h>
 *
 * @section Remarks
 *
 * This  will compute the bit-score for a @link BlastMatch @endlink. At least the following
 * members need to be set before calling this function:
 * <li> blastMatch.@link BlastMatch::alignStats @endlink.@link AlignmentStats::alignmentScore @endlink (if you have valid
 * alignRow-members (@link BlastMatch::alignRow0 @endlink, @link BlastMatch::alignRow1 @endlink), you can call
 * @link BlastMatch#computeAlignmentStats @endlink to compute the stats member). </li>
 */

template <typename TBlastMatch,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline double
computeBitScore(TBlastMatch & match,
                BlastIOContext<TScore, p, h> const & context)
{
    match.bitScore = computeBitScore(match.alignStats.alignmentScore, context.scoringScheme);
    return match.bitScore;
}

/*!
 * @fn BlastScoringScheme#computeEValue
 * @brief Compute the E-Value for a given raw score and sequence lengths.
 * @signature double computeEValue(rawScore, queryLength, dbLength, blastScoringScheme);
 *
 * @param[in]   rawScore            The raw score of the alignment.
 * @param[in]   queryLength         The length of the query sequence.
 * @param[in]   dbLength            The length of the database in case of many subject sequences, otherwise the
 * length of the subject sequence.
 * @param[in]   blastScoringScheme The @link BlastScoringScheme @endlink.
 *
 * @return double The e-value computed.
 * @headerfile <seqan/blast.h>
 */

template <typename TScore>
inline double
_computeEValue(double const rawScore,
               double const adjustedQueryLength,
               double const adjustedDbLength,
               BlastScoringScheme<TScore> const & scheme)
{
    return getKappa(scheme) * adjustedQueryLength * adjustedDbLength * std::exp(-getLambda(scheme) * rawScore);
}

template <typename T, typename TScore>
inline void
_conditionalDec(T &, BlastScoringScheme<TScore> const &)
{}

template <typename T>
inline void
_conditionalDec(T & val, BlastScoringScheme<Score<int, Simple>> const & scheme)
{
    typedef KarlinAltschulValues<Score<int, Simple>> TKAValues;
    if (TKAValues::VALUE[scheme.parameterIndex][10])
        --val;
}

template <typename TScore>
inline double
computeEValue(uint64_t rawScore,
              uint64_t const queryLength,
              uint64_t const dbLength,
              BlastScoringScheme<TScore> const & scheme)
{
    // for some parameters the score has to be "rounded down" to being even
    _conditionalDec(rawScore, scheme);

    uint64_t adj = _lengthAdjustment(dbLength, queryLength, scheme);
    return _computeEValue(rawScore, queryLength - adj, dbLength - adj, scheme);
}

/*!
 * @fn BlastMatch#computeEValue
 * @brief Compute the E-Value for a @link BlastMatch @endlink.
 * @signature double computeEValue(blastMatch, qLength, context);
 * [[deprecated]] double computeEValue(blastMatch, context);
 *
 * @param[in,out]   blastMatch  A @link BlastMatch @endlink that has a valid align member.
 * @param[in]       qLength     The length of the query sequence (pass @link BlastRecord::qLength @endlink).
 * @param[in,out]   context     A @link BlastIOContext @endlink with parameters and buffers.
 *
 * @return double blastMatch.@link BlastMatch::eValue @endlink after computation
 * @headerfile <seqan/blast.h>
 *
 * @section Remarks
 *
 * This  will compute the e-value for a @link BlastMatch @endlink. At least the following
 * members need to be set before calling this function:
 * <li> blastMatch.@link BlastMatch::alignStats @endlink.@link AlignmentStats::alignmentScore @endlink (if you have valid
 * alignRow-members (@link BlastMatch::alignRow0 @endlink, @link BlastMatch::alignRow1 @endlink), you can call
 * @link BlastMatch#computeAlignmentStats @endlink to compute the stats member). </li>
 * <li> blastMatch.@link BlastMatch::qLength @endlink (only in the deprecated interface)</li>
 * <li> context.@link BlastIOContext::dbTotalLength @endlink </li>
 *
 * Note, that in contrast to the general interface (@link BlastScoringScheme#computeEValue @endlink), this interface
 * caches some computation steps in the context and is therefore faster.
 */

template <typename TBlastMatch,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline double
computeEValue(TBlastMatch & match,
              uint64_t ql,
              BlastIOContext<TScore, p, h> & context)
{
    SEQAN_ASSERT_GT(ql, 0ull);

    // convert to 64bit and divide for translated sequences
    ql = ql / (qIsTranslated(context.blastProgram) ? 3 : 1);
    // length adjustment not yet computed
    if (context._cachedLengthAdjustments.find(ql) == context._cachedLengthAdjustments.end())
        context._cachedLengthAdjustments[ql] = _lengthAdjustment(context.dbTotalLength, ql, context.scoringScheme);

    uint64_t adj = context._cachedLengthAdjustments[ql];

    match.eValue = _computeEValue(match.alignStats.alignmentScore,
                                  ql - adj,
                                  context.dbTotalLength - adj,
                                  context.scoringScheme);
    return match.eValue;
}

template <typename TBlastMatch,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
[[deprecated("Use the interface with an explicit query length parameter instead (use the record's member).")]]
inline double
computeEValue(TBlastMatch & match,
              BlastIOContext<TScore, p, h> & context)
{
    return computeEValue(match, match.qLength, context);
}

}

#endif // ndef __BLAST_STATISTICS_H__
