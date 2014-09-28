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
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================
// Definition of all AFScore and the specialisations, D2, D2star, D2z
// and N2.
// ==========================================================================

// TODO(holtgrew): Make struct a class here.

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_BASE_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_BASE_H_

namespace seqan {

/*!
 * @class AFScore
 * @headerfile <seqan/alignment_free.h>
 * @brief Used to specify parameters and methods for alignment-free sequence comparison.
 * 
 * @signature template <typename TSpec>
 *            struct AFScore;
 * 
 * @tparam TSpec Tag for specialization.
 * 
 * @see alignmentFreeComparison
 */

/**
.Class.AFScore:
..cat:Alignment Free
..summary:Used to specify parameters and methods for alignment-free sequence comparison
..signature:AFScore<TScoreType>
..param.TScoreType:Method to use: N2, D2, D2Star, D2z
..include:seqan/alignment_free.h
..see:Function.alignmentFreeComparison
*/

template <typename TSpec>
struct AFScore;

/*!
 * @class D2AFScore D2 AFScore
 * @extends AFScore
 * @headerfile <seqan/alignment_free.h>
 * 
 * @brief D2 computes the inner product of the kmer count vectors.
 * 
 * @signature template <>
 *            struct AFScore<D2>;
 * 
 * To be used for alignment free comparison.
 *
 * @section References
 *
 * Lippert RA, et al. Distributional regimes for the number of k-word matches between two random sequences. Proc. Natl
 * Acad. Sci. USA 2002.
 *
 * @see alignmentFreeComparison
 *
 * @var unsigned D2AFScore::kmerSize;
 * @brief Size of the kmers.
 * 
 * @var bool D2AFScore::verbose;
 * @brief <tt>true</tt> to enable verbose debug output.
 * 
 * @fn D2AFScore::AFScore
 * 
 * @brief Constructor
 * 
 * @signature AFScore::AFScore(kmerSize, verbose);
 * 
 * @param[in] verbose  This option will report progress to standard output (<tt>bool</tt).
 * @param[in] kmerSize Size of kmer (<tt>unsigned</tt>).
 * 
 * @see alignmentFreeComparison
 */

/**
.Spec.D2 AFScore
..general:Class.AFScore
..cat:Alignment Free
..summary:D2 computes the inner product of the kmer count vectors
..signature:AFScore<D2>
..remarks:D2 can be used for alignment-free sequence comparison
..cite:Lippert RA, et al. Distributional regimes for the number of k-word matches between two random sequences. Proc. Natl Acad. Sci. USA 2002.
..include:seqan/alignment_free.h

.Memvar.D2 AFScore#kmerSize:
..class:Spec.D2 AFScore
..summary:Size of the kmers

.Memfunc.D2 AFScore#AFScore
..class:Spec.D2 AFScore
..summary:Constructor
..signature:AFScore(kmerSize, verbose)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..param.verbose:This option will report progress to standard output.
...type:nolink:bool
...default:false
..remarks:
...text:No remarks
..see:Function.alignmentFreeComparison
*/

struct D2_;     // Inner product of k-mer counts, d2 score
typedef Tag<D2_> const D2;

template <>
struct AFScore<D2>
{
    unsigned kmerSize;
    bool verbose;
    AFScore(unsigned k, bool verbose_ = false) : kmerSize(k), verbose(verbose_)
    {}
};

/*!
 * @class D2StarAFScore D2Star AFSScore
 * @extends AFScore
 * @headerfile <seqan/alignment_free.h>
 * @brief D2Star computes the inner product of the standardised kmer count vectors.
 * 
 * @signature template <>
 *            struct AFScore<D2Star>;
 * 
 * D2Star can be used for alignment-free sequence comparison, this version calculates the background model on the
 * concatenation of both sequences
 *
 * @section References
 *
 * Reinert, G.; Chew, D.; Sun, F.,  Waterman, M. S. Alignment-Free Sequence Comparison (I): Statistics and Power. J
 * Comput Biol, 2009.
 *
 * @see alignmentFreeComparison
 * 
 * @fn D2StarAFScore::AFScore
 * @brief Constructor
 * @signature AFScore::AFScore(kmerSize, bgModelOrder, verbose);
 * 
 * @param[in] kmerSize      Size of kmer, <tt>unsigned</tt>.
 * @param[in] bgModelOrder Order of the background Markov model, <tt>unsigned</tt>.
 * @param[in] verbose      This option will report progress to standard output, <tt>bool</tt>.
 * 
 * @var unsigned D2StarAFScore::kmerSize
 * @brief Size of the kmers.
 * 
 * @var CharString D2StarAFScore::outputFile
 * @brief When specified, all kmerWeights will be written to this file, for every sequence, and for every sequence
 *        comparison.
 * 
 * @var unsigned D2StarAFScore::bgModelOrder
 * @brief Order of the background model.
 */

/**
.Spec.D2Star AFScore
..general:Class.AFScore
..cat:Alignment Free
..summary:D2Star computes the inner product of the standardised kmer count vectors
..signature:AFScore<D2Star>
..remarks:D2Star can be used for alignment-free sequence comparison, this version calculates the background model on the concatenation of both sequences
..cite:Reinert, G.; Chew, D.; Sun, F. & Waterman, M. S. Alignment-Free Sequence Comparison (I): Statistics and Power. J Comput Biol, 2009.
..include:seqan/alignment_free.h

.Memvar.D2Star AFScore#kmerSize:
..class:Spec.D2Star AFScore
..summary:Size of the kmers

.Memvar.D2Star AFScore#bgModelOrder:
..class:Spec.D2Star AFScore
..summary:Order of the background model

.Memvar.D2Star AFScore#outputFile:
..class:Spec.D2Star AFScore
..summary: When specified, all kmerWeights will be written to this file, for every sequence, and for every sequence comparison

.Memfunc.D2Star AFScore#AFScore:
..class:Spec.D2Star AFScore
..summary:Constructor
..signature:AFScore(kmerSize, bgModelOrder, verbose)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..param.bgModelOrder: Order of the background Markov model
...type:nolink:unsigned
..param.verbose:This option will report progress to standard output.
...type:nolink:bool
...default:false
..remarks:
...text:No remarks
..see:Function.alignmentFreeComparison
*/

struct D2Star_;        // Reinert and Waterman, D2 with centralised and standardised counts
typedef Tag<D2Star_> const D2Star;

template <>
struct AFScore<D2Star>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    bool verbose;


    AFScore<D2Star>(unsigned k, unsigned m, bool verbose_ = false) :
        kmerSize(k), bgModelOrder(m), verbose(verbose_)
    {}
};

/*!
 * @class N2AFScore N2 AFScore
 * @extends AFScore
 * @headerfile <seqan/alignment_free.h>
 * @brief N2 computes the inner product of the standardised neighbourhood kmer count vectors.
 * 
 * @signature template <>
 *            struct AFScore<N2>;
 * 
 * N2 can be used for alignment-free sequence comparison.
 *
 * @section References
 *
 * Jonathan Goeke, Marcel H. Schulz, Julia Lasserre, and Martin Vingron.
 Estimation of Pairwise Sequence Similarity of Mammalian Enhancers with Word Neighbourhood Counts. Bioinformatics
 (2012). 
 * 
 * @fn N2AFScore::AFScore
 * 
 * @brief Constructor
 * 
 * @signature AFScore::AFScore(kmerSize, bgModelOrder, outputFile, verbose);
 * @signature AFScore::AFScore(kmerSize, bgModelOrder, revCom, outputFile, verbose);
 * @signature AFScore::AFScore(kmerSize, bgModelOrder, revCom, mismatches, mismatchWeight, outputFile, verbose);
 * 
 * @param[in] kmerSize       Size of kmer, <tt>unsigned</tt>.
 * @param[in] bgModelOrder   Order of the background Markov model, <tt>unsigned</tt>.
 * @param[in] outputFile     When specified, all normalised and standardised kmer neighbourhood counts will be written
 *                           to this file for every sequence, @link CharString @endlink.
 * @param[in] revCom         Scoring of reverse complements words [''/'max'/'min'/'mean'/'both_strands'/],
 *                           @link CharString @endlink.
 * @param[in] verbose        This option will report progress to standard output, <tt>bool</tt>, defaults to
 *                           <tt>false</tt>.
 * @param[in] mismatches     Includes words with one mismatch into the word neighbourhood, <tt>unsigned</tt>, 0 or 1.
 * @param[in] mismatchWeight Weight of word counts with one mismatch, double.
 * 
 * @see alignmentFreeComparison
 * 
 * @var unsigned N2AFScore::kmerSize;
 * @brief Size of the kmers.
 * 
 * @var double N2AFScore::mismatchWeight;
 * @brief Weight for approximate word matches
 * 
 * @var CharString N2AFScore::revCom;
 * @brief Scoring of reverse complements words, @link CharString @endlink [''/'max'/'min'/'mean'/'both_strands'/].
 * 
 * @var CharString N2AFScore::outputFile;
 * @brief When specified, all kmerWeights for every sequence will be written to this file.
 * 
 * @var unsigned N2AFScore::bgModelOrder;
 * @brief Order of the background model
 * 
 * @var unsigned N2AFScore::mismatches;
 * @brief Approximate word matches [0(exact)/1(one mismatch)]
 *
 * @var bool N2AFScore::verbose;
 * @brief <tt>true</tt> to enable verbose debug output.
 */

/**
.Spec.N2 AFScore
..general:Class.AFScore
..cat:Alignment Free
..summary:N2 computes the inner product of the standardised neighbourhood kmer count vectors
..signature:AFScore<N2>
..remarks:N2 can be used for alignment-free sequence comparison.
..cite:Jonathan Goeke et al. (to appear).
..include:seqan/alignment_free.h

.Memvar.N2 AFScore#kmerSize:
..class:Spec.N2 AFScore
..summary:Size of the kmers

.Memvar.N2 AFScore#bgModelOrder:
..class:Spec.N2 AFScore
..summary:Order of the background model

.Memvar.N2 AFScore#revCom:
..class:Spec.N2 AFScore
..summary:Scoring of reverse complements words [''/'max'/'min'/'mean'/'both_strands'/]

.Memvar.N2 AFScore#mismatches:
..class:Spec.N2 AFScore
..summary:Approximate word matches [0(exact)/1(one mismatch)]

.Memvar.N2 AFScore#mismatchWeight:
..class:Spec.N2 AFScore
..summary:Weight for approximate word matches

.Memvar.N2 AFScore#outputFile:
..class:Spec.N2 AFScore
..summary: When specified, all kmerWeights for every sequence will be written to this file

.Memfunc.N2 AFScore#AFScore:
..class:Spec.N2 AFScore
..summary:Constructor
..signature:AFScore(kmerSize, bgModelOrder, outputFile, verbose)
..signature:AFScore(kmerSize, bgModelOrder, revCom, outputFile, verbose)
..signature:AFScore(kmerSize, bgModelOrder, revCom, mismatch, mismatchWeight, outputFile, verbose)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..param.bgModelOrder: Order of the background Markov model
...type:nolink:unsigned
..param.revCom: Scoring of reverse complements words [''/'max'/'min'/'mean'/'both_strands'/]
...type:nolink:String<char>
..param.mismatch: Includes words with one mismatch into the word neighbourhood. [0/1]
...type:nolink:unsigned
..param.mismatchWeight: Weight of word counts with one mismatch
...type:nolink:double
..param.outputFile:When specified, all normalised and standardised kmer neighbourhood counts will be written to this file for every sequence
...type:nolink:String<char>
...default:""
..param.verbose:This option will report progress to standard output.
...type:nolink:bool
...default:false
..remarks:
...text:No remarks
..see:Function.alignmentFreeComparison
*/

struct N2_;     // Reinert and Waterman, D2 with centralised and standardised counts
typedef Tag<N2_> const N2;

template <>
struct AFScore<N2>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    String<char> revCom;    // Count reverse complement words?
                            // revCom="";"mean","max","both_strands"
    unsigned mismatches;    // Currently 0 or 1
    double mismatchWeight;  // Weight of words in the mismatch neighbourhood
    bool verbose;
    bool norm;              // Normalize score? Needed to provide a proper similarity measure
    String<char> outputFile;  // Output of all kmer weights for every sequence into this file

    // Constructor for the simple case with only exact word counts (N2*)
    AFScore(unsigned k, unsigned m, String<char> kmerWeightsFile = "", bool verbose_ = false)
    {
        kmerSize = k;
        bgModelOrder = m;
        outputFile = kmerWeightsFile;
        verbose = verbose_;
        revCom = "";
        mismatches = 0;
        mismatchWeight = 1.0;
        norm = true;

    };

    // Constructor for the case with exact word counts and reverse complement (N2rc)
    AFScore(unsigned k, unsigned m, String<char> revCom_, String<char> kmerWeightsFile = "", bool verbose_ = false)
    {
        kmerSize = k;
        bgModelOrder = m;
        revCom = revCom_;
        outputFile = kmerWeightsFile;
        verbose = verbose_;
        mismatches = 0;
        mismatchWeight = 1.0;
        norm = true;
    };

    // Constructor for the case with mismatch-neighbourhood word counts and reverse complement (N2mmrc)
    AFScore(unsigned k,
            unsigned m,
            String<char> revCom_,
            unsigned mm, double mmw,
            String<char> kmerWeightsFile = "",
            bool verbose_ = false)
    {
        kmerSize = k;
        bgModelOrder = m;
        revCom = revCom_;
        mismatches = mm;
        mismatchWeight = mmw;
        outputFile = kmerWeightsFile;
        verbose = verbose_;
        norm = true;
    };
};

/*!
 * @class D2zAFScore D2zAFScore
 * @extends AFScore
 * @headerfile <seqan/alignment_free.h>
 * @brief D2z computes a z-score of the inner product of kmer count vectors
 * 
 * @signature template <>
 *            struct AFScore<D2z>;
 * 
 * D2z can be used for alignment-free sequence comparison. The algorithm differs from the original implementation by
 * the way masked sequences are handled
 * 
 * @section References
 * 
 * Kantorovitz, M. R.; Robinson, G. E., Sinha, S. A statistical method for alignment-free comparison of regulatory
 * sequences. Bioinformatics, 2007.
 * 
 * @fn D2zAFScore::AFScore
 * 
 * @brief Constructor
 * 
 * @signature AFScore::AFScore(kmerSize, bgModelOrder[, verbose]);
 * 
 * @param[in] kmerSize     Size of kmer, <tt>unsigned</tt>.
 * @param[in] bgModelOrder Order of the background Markov model, <tt>unsigned</tt>.
 * @param[in] verbose      This option will report progress to standard output; <tt>bool</tt>, defaults to
 *                         <tt>false</tt>.
 * 
 * @var unsigned D2zAFScore::bgModelOrder;
 * @brief Order of the background model
 * 
 * @var unsigned D2zAFScore::kmerSize;
 * @brief Size of the kmers
 *
 * @var bool D2zAFScore::verbose;
 * @brief <tt>true</tt> to enable verbose debug output.
] */ 

/**
.Spec.D2z AFScore
..cat:Alignment Free
..summary:D2z computes a z-score of the inner product of kmer count vectors
..signature:AFScore<D2z>
..general:Class.AFScore
..remarks:D2z can be used for alignment-free sequence comparison. The algorithm differs from the original implementation by the way masked sequences are handled
..cite:Kantorovitz, M. R.; Robinson, G. E. & Sinha, S. A statistical method for alignment-free comparison of regulatory sequences. Bioinformatics, 2007.
..include:seqan/alignment_free.h

.Memvar.D2z AFScore#kmerSize:
..class:Spec.D2z AFScore
..summary:Size of the kmers

.Memvar.D2z AFScore#bgModelOrder:
..class:Spec.D2z AFScore
..summary:Order of the background model

.Memfunc.D2z AFScore#AFScore:
..class:Spec.D2z AFScore
..summary:Constructor
..signature:AFScore(kmerSize, bgModelOrder, verbose)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..param.bgModelOrder: Order of the background Markov model
...type:nolink:unsigned
..param.verbose:This option will report progress to standard output.
...type:nolink:bool
...default:false
..remarks:
...text:nolink:No remarks
*/

struct D2z_;  // Inner product of k-mer counts, d2 score with z-score
typedef Tag<D2z_> const D2z;

template <>
struct AFScore<D2z>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    bool verbose;
    AFScore<D2z>(unsigned k, unsigned m, bool verbose_ = false) :
        kmerSize(k), bgModelOrder(m), verbose(verbose_)
    {}
};

}  // namespace seqan

#endif  // SEQAN_EXTRAS_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_BASE_H_
