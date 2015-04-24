// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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

    /////////////////////////////////////////////////////////////////////////
    //                        BLAST NUMBER VOODOO                          //
    /////////////////////////////////////////////////////////////////////////

/** Karlin-Altschul parameter values taken from blast_stat.c */

#ifndef __BLAST_STATISTICS_H__
#define __BLAST_STATISTICS_H__

#ifndef NCBI_INT2_MAX
#define NCBI_INT2_MAX    32767
#endif

namespace SEQAN_NAMESPACE_MAIN
{

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class BlastScoringAdapter
// ----------------------------------------------------------------------------

/*!
 * @class BlastScoringAdapter
 * @headerfile <seqan/blast.h>
 * @brief An object that holds Statistical parameters in addition to the
 * scoring Scheme
 * @signature struct BlastScoringAdapter<TScore, TSpec> { ... };
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
 * combination of gap scores is supported. Run
 * @link BlastScoringAdapter#isValid @endlink after construction.
 *
 * @tparam TScoringScheme A SeqAn Score type
 * @tparam TSpec This only prevents instantiation of tables that aren't used
 * (defaults to void)
 *
 * @fn BlastScoringAdapter::BLASTScoringAdapter(TScore const & scheme)
 * @signature BlastScoringAdapter(TScore const & scheme)
 * @param scheme A scoring scheme.
 * @brief The Constructor. Upon construction the statistical values are
 * automatically selected for the given scheme. Run
 * @link BlastScoringAdapter#isValid @endlink to see, if this was successful.
 */

/// GENERIC
template <typename TScore, typename TSpec>
struct BlastScoringAdapter;

// ============================================================================
// Forwards
// ============================================================================

template <typename TScoringScheme>
bool
assignScoreScheme(BlastScoringAdapter<TScoringScheme> &,
                  TScoringScheme const &);

// ============================================================================
// Classes
// ============================================================================

/// BLOSUM30
// not implemented in BLAST

/// BLOSUM45
template <typename TSpec>
struct BlastScoringAdapter<Blosum45, TSpec>
{
    typedef uint8_t size_t;

    /* statics */
    static constexpr size_t nParams = 8;
    static constexpr size_t nParamSets = 16;
    static double const karlinAltschulValues[nParamSets][nParams];

    /* members */
    Blosum45 scheme;
    size_t index = std::numeric_limits<size_t>::max();
    // member signifying the paramSet belonging to scheme

    BlastScoringAdapter() {};
    BlastScoringAdapter(Blosum45 const & _scheme)
    {
        assignScoreScheme(*this, _scheme);
    }
};

template <typename TSpec>
double const BlastScoringAdapter<Blosum45, TSpec>::karlinAltschulValues
  [BlastScoringAdapter<Blosum45, TSpec>::nParamSets]
  [BlastScoringAdapter<Blosum45, TSpec>::nParams] =
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
struct BlastScoringAdapter<Blosum62, TSpec>
{
    typedef uint8_t size_t;

    /* statics */
    static constexpr size_t nParams = 8;
    static constexpr size_t nParamSets = 12;
    static double const karlinAltschulValues[nParamSets][nParams];

    /* members */
    Blosum62 scheme;
    size_t index = std::numeric_limits<size_t>::max();
    // member signifying the paramSet belonging to scheme

    BlastScoringAdapter() {};
    BlastScoringAdapter(Blosum62 const & _scheme)
    {
        assignScoreScheme<Blosum62>(*this, _scheme);
    }
};

template <typename TSpec>
double const BlastScoringAdapter<Blosum62, TSpec>::karlinAltschulValues
  [BlastScoringAdapter<Blosum62, TSpec>::nParamSets]
  [BlastScoringAdapter<Blosum62, TSpec>::nParams] =
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
struct BlastScoringAdapter<Blosum80, TSpec>
{
    typedef uint8_t size_t;

    /* statics */
    static constexpr size_t nParams = 8;
    static constexpr size_t nParamSets = 10;
    static double const karlinAltschulValues[nParamSets][nParams];


    /* members */
    Blosum80 scheme;
    size_t index = std::numeric_limits<size_t>::max();
    // member signifying the paramSet belonging to scheme

    BlastScoringAdapter() {};
    BlastScoringAdapter(Blosum80 const & _scheme)
    {
        assignScoreScheme(*this, _scheme);
    }
};

template <typename TSpec>
double const BlastScoringAdapter<Blosum80, TSpec>::karlinAltschulValues
 [BlastScoringAdapter<Blosum80, TSpec>::nParamSets]
 [BlastScoringAdapter<Blosum80, TSpec>::nParams] =
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
struct BlastScoringAdapter<Pam250, TSpec>
{
    typedef uint8_t size_t;

    /* statics */
    static constexpr size_t nParams = 8;
    static constexpr size_t nParamSets = 16;
    static double const karlinAltschulValues[nParamSets][nParams];


    /* members */
    Pam250 scheme;
    size_t index = std::numeric_limits<size_t>::max();
    // member signifying the paramSet belonging to scheme

    BlastScoringAdapter() {};
    BlastScoringAdapter(Pam250 const & _scheme)
    {
        assignScoreScheme(*this, _scheme);
    }
};

template <typename TSpec>
double const BlastScoringAdapter<Pam250, TSpec>::karlinAltschulValues
  [BlastScoringAdapter<Pam250, TSpec>::nParamSets]
  [BlastScoringAdapter<Pam250, TSpec>::nParams] =
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
struct BlastScoringAdapter<Score<int, Simple>, TSpec>
{
    typedef uint8_t size_t;

    /* statics */
    static constexpr size_t nParams = 11;
    static constexpr size_t nParamSets = 60;
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
    static double const karlinAltschulValues[nParamSets][nParams];

    /* members */
    Score<int, Simple> scheme;
    size_t index = std::numeric_limits<size_t>::max();
    // member signifying the paramSet belonging to scheme
    BlastScoringAdapter() {};
    BlastScoringAdapter(Score<int, Simple> const & _scheme)
    {
        assignScoreScheme(*this, _scheme);
    }
};

template <typename TSpec>
double const BlastScoringAdapter<Score<int, Simple>, TSpec>::karlinAltschulValues
 [BlastScoringAdapter<Score<int, Simple>, TSpec>::nParamSets]
 [BlastScoringAdapter<Score<int, Simple>, TSpec>::nParams] =
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
// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isValid
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#isValid
 * @headerfile <seqan/blast.h>
 * @brief check whether object holds correct Karlin-Altschul values for its scheme
 * @signature bool isValid(blastScoringAdapter);
 */

template <typename TScoringScheme>
inline bool
isValid(BlastScoringAdapter<TScoringScheme> const & adapter)
{
    typedef BlastScoringAdapter<TScoringScheme> TAdapter;
    return adapter.index != std::numeric_limits<typename TAdapter::size_t>::max();
}

// ----------------------------------------------------------------------------
// Function _selectSet
// ----------------------------------------------------------------------------

template <typename TMatrixSpec>
inline bool
_selectSet(BlastScoringAdapter<
                Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> & adapter)
{
    int scO = scoreGapOpen(adapter.scheme) - scoreGapExtend(adapter.scheme);
    typedef BlastScoringAdapter<
                Score<int, ScoreMatrix<AminoAcid, TMatrixSpec> > > TAdapter;

    for (typename TAdapter::size_t i = 0;
         i < TAdapter::nParamSets; //karlinAltschulValues.size();
         ++i)
    {
        if ((TAdapter::karlinAltschulValues[i][0] == -scO)  &&
            (TAdapter::karlinAltschulValues[i][1] == -scoreGapExtend(adapter.scheme)))
        {
            adapter.index = i;
            return true;
        }
    }
    // no suitable adapter
    adapter.index = std::numeric_limits<typename TAdapter::size_t>::max();
    return false;
}

inline bool
_selectSet(BlastScoringAdapter<Score<int, Simple>> & adapter)
{
    int scO = scoreGapOpen(adapter.scheme) - scoreGapExtend(adapter.scheme);
    typedef BlastScoringAdapter<Score<int, Simple>> TAdapter;

    for (typename TAdapter::size_t i = 0;
         i < TAdapter::nParamSets; //karlinAltschulValues.size();
         ++i)
    {
        if ((TAdapter::karlinAltschulValues[i][0] == scoreMatch(adapter.scheme))  &&
            (TAdapter::karlinAltschulValues[i][1] == -scoreMismatch(adapter.scheme)) &&
            (TAdapter::karlinAltschulValues[i][2] == -scO) &&
            (TAdapter::karlinAltschulValues[i][3] == -scoreGapExtend(adapter.scheme)))
        {
            adapter.index = i;
            return true;
        }
    }
    // no suitable adapter
    adapter.index = std::numeric_limits<typename TAdapter::size_t>::max();
    return false;
}

// ----------------------------------------------------------------------------
// Function assignScoreScheme
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#assignScoreScheme
 * @headerfile <seqan/blast.h>
 * @brief assign a new SeqAn @link Score @endlink to the adapter
 * @signature bool assignScoreScheme(blastScoringAdapter, scoringScheme);
 *
 * @param blastScoringAdapter the adapter to modify
 * @param scoringScheme the new scoring scheme
 *
 * @return whether the adapter @link BlastScoringAdapter#isValid @endlink after assignment
 *
 */

template <typename TScoringScheme>
inline bool
assignScoreScheme(BlastScoringAdapter<TScoringScheme> & adapter,
                  TScoringScheme const & _scheme)
{
    adapter.scheme = _scheme;
    return _selectSet(adapter);
}

// ----------------------------------------------------------------------------
// Function getScoreScheme
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#getScoreScheme
 * @headerfile <seqan/blast.h>
 * @brief access the scoring scheme inside the adapter
 * @signature TScoringScheme getScoreScheme(blastScoringAdapter);
 *
 * @param blastScoringAdapter the adapter
 *
 * @return const Reference to the @link Score @endlink object inside the adapter
 */

template <typename TScoringScheme>
TScoringScheme const &
getScoreScheme(BlastScoringAdapter<TScoringScheme> const & adapter)
{
    return adapter.scheme;
}

// ----------------------------------------------------------------------------
// Function getLambda
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#getLambda
 * @headerfile <seqan/blast.h>
 * @brief get the &lambda; value of the Karlin-Altschul parameters.
 * @signature double getLambda(blastScoringAdapter);
 */

template <typename TMatrixSpec>
inline double
getLambda(BlastScoringAdapter<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec>>> const & adapter)
{
    typedef BlastScoringAdapter<Score<int,
        ScoreMatrix< AminoAcid, TMatrixSpec>>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][3];
}

inline double
getLambda(BlastScoringAdapter<Score<int, Simple>> const & adapter)
{
    typedef BlastScoringAdapter<Score<int, Simple>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][4];
}

// ----------------------------------------------------------------------------
// Function getKappa
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#getKappa
 * @headerfile <seqan/blast.h>
 * @brief get the &Kappa; value of the Karlin-Altschul parameters.
 * @signature double getKappa(blastScoringAdapter);
 */

template <typename TMatrixSpec>
inline double
getKappa(BlastScoringAdapter<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec> > > const & adapter)
{
    typedef BlastScoringAdapter<Score<int,
        ScoreMatrix< AminoAcid, TMatrixSpec>>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][4];
}

inline double
getKappa(BlastScoringAdapter<Score<int, Simple>> const & adapter)
{
    typedef BlastScoringAdapter<Score<int, Simple>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][5];
}

// ----------------------------------------------------------------------------
// Function getH
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#getH
 * @headerfile <seqan/blast.h>
 * @brief get the H value of the Karlin-Altschul parameters.
 * @signature double getH(blastScoringAdapter);
 */

template <typename TMatrixSpec>
inline double
getH(BlastScoringAdapter<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec> > > const & adapter)
{
    typedef BlastScoringAdapter<Score<int,
        ScoreMatrix< AminoAcid, TMatrixSpec>>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][5];
}

inline double
getH(BlastScoringAdapter<Score<int, Simple>> const & adapter)
{
    typedef BlastScoringAdapter<Score<int, Simple>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][6];
}

// ----------------------------------------------------------------------------
// Function getAlpha
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#getAlpha
 * @headerfile <seqan/blast.h>
 * @brief get the &alpha; value of the Karlin-Altschul parameters.
 * @signature double getAlpha(blastScoringAdapter);
 */

template <typename TMatrixSpec>
inline double
getAlpha(BlastScoringAdapter<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec> > > const & adapter)
{
    typedef BlastScoringAdapter<Score<int,
        ScoreMatrix< AminoAcid, TMatrixSpec>>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][6];
}

inline double
getAlpha(BlastScoringAdapter<Score<int, Simple>> const & adapter)
{
    typedef BlastScoringAdapter<Score<int, Simple>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][7];
}

// ----------------------------------------------------------------------------
// Function getBeta
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#getBeta
 * @headerfile <seqan/blast.h>
 * @brief get the &beta; value of the Karlin-Altschul parameters.
 * @signature double getBeta(blastScoringAdapter);
 */

template <typename TMatrixSpec>
inline double
getBeta(BlastScoringAdapter<Score<int, ScoreMatrix<AminoAcid, TMatrixSpec> > > const & adapter)
{
    typedef BlastScoringAdapter<Score<int,
        ScoreMatrix< AminoAcid, TMatrixSpec>>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][7];
}

inline double
getBeta(BlastScoringAdapter<Score<int, Simple>> const & adapter)
{
    typedef BlastScoringAdapter<Score<int, Simple>> TAdapter;
    SEQAN_ASSERT(isValid(adapter));
    return TAdapter::karlinAltschulValues[adapter.index][8];
}

// ----------------------------------------------------------------------------
// Function for converting scoringSchemes
// ----------------------------------------------------------------------------

// see the documentation for setBlastScoringScheme() in blast_io_context.h

template <typename TValue, typename TSpec>
inline void
seqanScoringScheme2blastScoringScheme(Score<TValue, TSpec> & scheme)
{
    setScoreGapOpen(scheme, scoreGapOpen(scheme) - scoreGapExtend(scheme));
}

template <typename TValue, typename TSpec>
inline void
blastScoringScheme2seqanScoringScheme(Score<TValue, TSpec> & scheme)
{
    setScoreGapOpen(scheme, scoreGapOpen(scheme) + scoreGapExtend(scheme));
}

// ----------------------------------------------------------------------------
// Function _lengthAdjustment
// ----------------------------------------------------------------------------

// Computes the length adjustment for E-value computation
// Based on the public domain NCBI BLAST code by Tom Madden.
template <typename TSize, typename TSize2, typename TScore>
inline
TSize
_lengthAdjustment(TSize     const & dbLength,
                  TSize2    const & queryLength,
                  BlastScoringAdapter<TScore> const & adapter)
{
    const double lambda = getLambda(adapter);
    const double K      = getKappa(adapter);
    const double logK   = std::log(K);
    const double alpha  = getAlpha(adapter);
    const double alphaByLambda = alpha/lambda;
    const double beta   = getBeta(adapter);
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
        double c  = n * m - _max(m, n) / K;

        if(c < 0) {
            return 0;
        } else {
            val_max = 2 * c / (mb + sqrt(mb * mb - 4 * c));
        }
    } // end scope of mb and c

    for(TSize i = 1; i <= maxIterations; i++) {
        totalLen = (m - val) * (n - val);
        double val_new  = alphaByLambda * (logK + log(totalLen)) + beta;  // proposed next value of val
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
          if(alphaByLambda * (logK + log(totalLen)) + beta >= val) {
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
 * The @link BlastMatch::align @endlink member is used as in-parameter to compute the
 * @link BlastMatch::alignStats @endlink member. This includes the raw
 * score, amount of gaps, mismatches et cetera. This is a prerequisite for printing a match to file or computing it's
 * e-value.
 */

template <typename TBlastMatch,
          typename TScore,
          typename TConString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
computeAlignmentStats(TBlastMatch & match,
                      BlastIOContext<TScore, TConString, p, h> const & context)
{
    computeAlignmentStats(match.alignStats, match.align, context.scoringAdapter.scheme);
}

// ----------------------------------------------------------------------------
// Function computeBitScore
// ----------------------------------------------------------------------------

/*!
 * @fn BlastScoringAdapter#computeBitScore
 * @brief Compute the bit-score for a given raw score.
 * @signature double computeBitScore(rawScore, blastScoringAdapter);
 *
 * @param[in]   rawScore            The raw score of the alignment.
 * @param[in]   blastScoringAdapter The @link BlastScoringAdapter @endlink.
 *
 * @return double The bit-score computed.
 * @headerfile <seqan/blast.h>
 */

template <typename TScore>
inline double
computeBitScore(double const rawScore, BlastScoringAdapter<TScore> const & adapter)
{
    return (getLambda(adapter) * rawScore - std::log(getKappa(adapter))) / std::log(2);
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
 * <li> blastMatch.@link BlastMatch::alignStats @endlink.@link AlignmentStats::alignmentScore @endlink (if you have a valig
 * @link BlastMatch::align @endlink member set, you can call
 * @link BlastMatch#computeAlignmentStats @endlink on to compute the stats member). </li>
 */

template <typename TBlastMatch,
          typename TScore,
          typename TConString,
          BlastProgram p,
          BlastTabularSpec h>
inline double
computeBitScore(TBlastMatch & match,
                BlastIOContext<TScore, TConString, p, h> const & context)
{
    match.bitScore = computeBitScore(match.alignStats.alignmentScore, context.scoringAdapter);
    return match.bitScore;
}

/*!
 * @fn BlastScoringAdapter#computeEValue
 * @brief Compute the E-Value for a given raw score and sequence lengths.
 * @signature double computeEValue(rawScore, queryLength, dbLength, blastScoringAdapter);
 *
 * @param[in]   rawScore            The raw score of the alignment.
 * @param[in]   queryLength         The length of the query sequence.
 * @param[in]   dbLength            The length of the database in case of many subject sequences, otherwise the
 * length of the subject sequence.
 * @param[in]   blastScoringAdapter The @link BlastScoringAdapter @endlink.
 *
 * @return double The e-value computed.
 * @headerfile <seqan/blast.h>
 */

template <typename TScore>
inline double
_computeEValue(double const rawScore,
               double const adjustedQueryLength,
               double const adjustedDbLength,
               BlastScoringAdapter<TScore> const & adapter)
{
    return getKappa(adapter) * adjustedQueryLength * adjustedDbLength * std::exp(-getLambda(adapter) * rawScore);
}

template <typename TScore>
inline double
computeEValue(uint64_t rawScore,
              uint64_t const queryLength,
              uint64_t const dbLength,
              BlastScoringAdapter<TScore> const & adapter)
{
    typedef BlastScoringAdapter<TScore> TAdapter;
    // for some parameters the score has to be "rounded down" to being even
    if ((TAdapter::nParams >= 11) && (TAdapter::karlinAltschulValues[adapter.index][10]))
        --rawScore;

    uint64_t adj = _lengthAdjustment(dbLength, queryLength, adapter);
    return _computeEValue(rawScore, queryLength - adj, dbLength - adj, adapter);
}

/*!
 * @fn BlastMatch#computeEValue
 * @brief Compute the E-Value for a @link BlastMatch @endlink.
 * @signature double computeEValue(blastMatch, context);
 *
 * @param[in,out]   blastMatch  A @link BlastMatch @endlink that has a valid align member.
 * @param[in,out]   context     A @link BlastIOContext @endlink with parameters and buffers.
 *
 * @return double blastMatch.@link BlastMatch::eValue @endlink after computation
 * @headerfile <seqan/blast.h>
 *
 * @section Remarks
 *
 * This  will compute the e-value for a @link BlastMatch @endlink. At least the following
 * members need to be set before calling this function:
 * <li> blastMatch.@link BlastMatch::alignStats @endlink.@link AlignmentStats::alignmentScore @endlink (if you have a valig
 * @link BlastMatch::align @endlink member set, you can call
 * @link BlastMatch#computeAlignmentStats @endlink to compute the stats member). </li>
 * <li> blastMatch.@link BlastMatch::qLength @endlink </li>
 * <li> context.@link BlastIOContext::dbTotalLength @endlink </li>
 *
 * Note, that in contrast to the general interface (@link BlastScoringAdapter#computeEValue @endlink), this interface
 * caches some computation steps in the context and is therefore faster.
 */

template <typename TBlastMatch,
          typename TScore,
          typename TConString,
          BlastProgram p,
          BlastTabularSpec h>
inline double
computeEValue(TBlastMatch & match,
              BlastIOContext<TScore, TConString, p, h> & context)
{
    uint64_t ql = match.qLength; // convert to 64bit, once
    // length adjustment not yet computed
    if (context.cachedLengthAdjustments.find(ql) == context.cachedLengthAdjustments.end())
        context.cachedLengthAdjustments[ql] = _lengthAdjustment(context.dbTotalLength, ql, context.scoringAdapter);

    uint64_t adj = context.cachedLengthAdjustments[ql];
    match.eValue = _computeEValue(match.alignStats.alignmentScore,
                                  ql - adj,
                                  context.dbTotalLength - adj,
                                  context.scoringAdapter);
    return match.eValue;
}

}

#endif // ndef __BLAST_STATISTICS_H__
