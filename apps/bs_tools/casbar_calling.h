/*==========================================================================

Methylation Calling: Version 2

Calculate snp and methylation probability in one step
        Pr(Dj | G = A:T)
    Pr(Dj | G = C^(1-beta)C_M^(beta):T)
    ...

Add additional cases: e.g. prob. converted and sequencing error
Add underconversion rates...

==========================================================================*/

#ifndef __APPS_BS_TOOLS_CASBAR_CALLING_H__
#define __APPS_BS_TOOLS_CASBAR_CALLING_H__

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

using namespace std;
using namespace seqan;

// TODO:
// precalculate probs and get rid of thousand different cases to check
// Compute prob. to observe given base i under the assumption that the underlying haplotype h
// For reads mapped to the forward strand
template<typename TProb, typename TErrorProb, typename TMethOptions>
inline void
getSingleBaseProbHaploF(TProb &singleProb, Dna i, DnaM h, TErrorProb &e, bool & origin, TMethOptions &methOptions)
{
    if (methOptions.uniformSeqErrorsCalling)
    {
        if ( h == 'C')                              // Haplotype C
        {
            if (i == 'C')                                               // correct and not converted + converted and error
                singleProb = (1.0-e)*(1.0-methOptions.convRate) ; //+ methOptions.convRate*(e/3.0);
            else if (i == 'T')                                          // error + correct and converted
                singleProb = e/3.0 + (1.0-e)*(methOptions.convRate);
            else                                                        // error
                singleProb = e/3.0;
        }
        else if ( h == 'D')                         // Haplotype Cm
        {
            if (i == 'C')                                               // correct and not converted + converted and error
                singleProb = (1.0-e)*(1.0-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
            else if (i == 'T')                                          // error + correct and converted
                singleProb = e/3.0 + (1.0-e)*methOptions.methConvRate;
            else                                                        // error
                singleProb = e/3.0;
        }
        else if ( h == 'G')                         // Haplotype G
        {
            if (i == 'G')
                singleProb = 1.0 - e;
            else
                singleProb = e/3.0;
        }
            else if ( h == 'H')                         // Haplotype Gm
        {
            if (i == 'G')
                singleProb = 1.0 - e;
            else
                singleProb = e/3.0;
        }
        else                                        // Haplotype T, A: no bs
        {
            if (i == h)
                singleProb = 1.0 - e;
            else
                singleProb = e/3.0;
        }
    }
    else
    {
        // temporary, just to try
        //
        TProb const *seqErrorFreqs = SeqErrorFreqsN<TProb, BsNonSimple>::getData();
        //TProb const *seqErrorFreqsTo = SeqErrorFreqsTo<TProb, BsNonSimple>::getData();

        if (origin)
        {
            if ( h == 'C')                              // Haplotype C
            {
                if (i == 'C')
                    singleProb = (1.0-e)*(1-methOptions.convRate) ;
                else if (i == 'T')
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] + (1.0-e)*(methOptions.convRate);
                else
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)];
            }
            else if ( h == 'D')                         // Haplotype Cm
            {
                if (i == 'C')
                    singleProb = (1.0-e)*(1.0-methOptions.methConvRate);
                else if (i == 'T')
                    singleProb = e*seqErrorFreqs[1*5 + ordValue(i)]  + (1.0-e)*methOptions.methConvRate;
                else
                    singleProb = e*seqErrorFreqs[1*5 + ordValue(i)] ;
            }
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] ;
            }
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[2*5 + ordValue(i)] ;
            }
            else                                        // Haplotype T, A: no bs
            {
                if (i == h)
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] ;
            }
        }
        else    // use sequencing error probs from RC
        {
            FunctorDna5OrdValueComplement<int> fCompl;
             if ( h == 'C')                              // Haplotype C
            {
                if (i == 'C')
                    singleProb = (1.0-e)*(1-methOptions.convRate) ;
                else if (i == 'T')
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] + (1.0-e)*(methOptions.convRate);
                else
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))];
            }
            else if ( h == 'D')                         // Haplotype Cm
            {
                if (i == 'C')
                    singleProb = (1.0-e)*(1.0-methOptions.methConvRate);
                else if (i == 'T')
                    singleProb = e*seqErrorFreqs[fCompl(1)*5 + fCompl(ordValue(i))]  + (1.0-e)*methOptions.methConvRate;
                else
                    singleProb = e*seqErrorFreqs[fCompl(1)*5 + fCompl(ordValue(i))] ;
            }
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] ;
            }
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[fCompl(2)*5 + fCompl(ordValue(i))] ;
            }
            else                                        // Haplotype T, A: no bs
            {
                if (i == h)
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] ;
            }
        }
    }
}
// For reads mapped to the reverse strand
template<typename TProb, typename TErrorProb, typename TMethOptions>
inline void
getSingleBaseProbHaploR(TProb &singleProb, Dna i, DnaM h, TErrorProb &e, bool & origin, TMethOptions &methOptions)
{
    if (methOptions.uniformSeqErrorsCalling)
    {
        // i is read base corresponding to forward strand
        if ( h == 'C')                              // Haplotype C
        {
            if (i == 'C')
                singleProb = 1.0-e;
            else
                singleProb = e/3.0;
        }
        else if ( h == 'D')                         // Haplotype Cm
        {
            if (i == 'C')
                singleProb = 1.0-e;
            else
                singleProb = e/3.0;
        }
        else if ( h == 'G')                         // Haplotype G
        {
            if (i == 'G')
                singleProb = (1.0-e)*(1-methOptions.convRate); //+ methOptions.convRate*(e/3.0);    // TODO: why not included???
            else if (i == 'A')
                singleProb = e/3.0 + (1.0-e)*methOptions.convRate;
            else
                singleProb = e/3.0;
        }
        else if ( h == 'H')                         // Haplotype Gm
        {
            if (i == 'G')
                singleProb = (1.0-e)*(1-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
            else if (i == 'A')
                singleProb = e/3.0 + (1.0-e)*methOptions.methConvRate;
            else
                singleProb = e/3.0;
         }
        else                                        // Haplotype T, A: no bs
        {
            if (i == h)
                singleProb = 1.0 - e;
            else
                singleProb = e/3.0;
        }
    }
    else
    {
        // temporary, just to try
        TProb const *seqErrorFreqs = SeqErrorFreqsN<TProb, BsNonSimple>::getData();
        if (origin)
        {
            FunctorDna5OrdValueComplement<int> fCompl;
           // i is read base corresponding to forward strand
            if ( h == 'C')                              // Haplotype C
            {
                if (i == 'C')
                    singleProb = 1.0-e;
                else
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))];
            }
            else if ( h == 'D')                         // Haplotype Cm
            {
                if (i == 'C')
                    singleProb = 1.0-e;
                else
                    singleProb = e*seqErrorFreqs[fCompl(1)*5 + fCompl(ordValue(i))];
            }
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')
                    singleProb = (1.0-e)*(1-methOptions.convRate); //+ methOptions.convRate*(e/3.0);    // TODO: why not included???
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))] + (1.0-e)*methOptions.convRate;
                else
                    singleProb = e/3.0;
            }
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')
                    singleProb = (1.0-e)*(1-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[fCompl(2)*5 + fCompl(ordValue(i))] + (1.0-e)*methOptions.methConvRate;
                else
                    singleProb = e*seqErrorFreqs[fCompl(2)*5 + fCompl(ordValue(i))];
             }
            else                                        // Haplotype T, A: no bs
            {
                if (i == h)
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[fCompl(ordValue(h))*5 + fCompl(ordValue(i))];
            }
        }
        else
        {
            // This read was projected to original BS strand, thus this is the original base
            if ( h == 'C')                              // Haplotype C
            {
                if (i == 'C')
                    singleProb = 1.0-e;
                else
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)];
            }
            else if ( h == 'D')                         // Haplotype Cm
            {
                if (i == 'C')
                    singleProb = 1.0-e;
                else
                    singleProb = e*seqErrorFreqs[1*5 + ordValue(i)];
            }
            else if ( h == 'G')                         // Haplotype G
            {
                if (i == 'G')
                    singleProb = (1.0-e)*(1-methOptions.convRate); //+ methOptions.convRate*(e/3.0);    // TODO: why not included???
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)] + (1.0-e)*methOptions.convRate;
                else
                    singleProb = e/3.0;
            }
            else if ( h == 'H')                         // Haplotype Gm
            {
                if (i == 'G')
                    singleProb = (1.0-e)*(1-methOptions.methConvRate); // + methOptions.methConvRate*(e/3.0);
                else if (i == 'A')
                    singleProb = e*seqErrorFreqs[2*5 + ordValue(i)] + (1.0-e)*methOptions.methConvRate;
                else
                    singleProb = e*seqErrorFreqs[2*5 + ordValue(i)];
             }
            else                                        // Haplotype T, A: no bs
            {
                if (i == h)
                    singleProb = 1.0 - e;
                else
                    singleProb = e*seqErrorFreqs[ordValue(h)*5 + ordValue(i)];
            }
        }
    }
}


// For NaiveMult multiplication of lHoods
template<typename TConstants, typename TValue>
inline void
addFactorToConstants(TConstants &constants, TValue &pC, TValue &pCm, TValue &pOther, unsigned &r, NaiveMult const &)
{
    SEQAN_ASSERT_NEQ(pC, 0.0);
    SEQAN_ASSERT_NEQ(pCm, 0.0);

    // Compute constant values:
    if (pOther > 0)
    {
        constants[0][r] = 0.5*(pC + pOther);
        constants[1][r] = 0.5*(-pC + pCm);
    }
    else
    {
        constants[0][r] = pC;
        constants[1][r] = -pC + pCm;
    }
}


// For Log function
template<typename TConstants, typename TValue>
inline void
addFactorToConstants(TConstants &constants, TValue &pC, TValue &pCm, TValue &pOther, unsigned &r, LogFunction const &)
{
    addFactorToConstants(constants, pC, pCm, pOther, r, NaiveMult());
}

template<typename TConstantSet, typename TStrand>
inline void
adjustConstantsSize(TConstantSet &constantSet, TStrand strand, NaiveMult const &)
{
    if (strand == 'F')
    {
        for (unsigned i = 0; i <= 1; ++i)
        {
            eraseBack(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')][i]);
            eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')][i]);
            eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')][i]);
            eraseBack(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')][i]);
        }
    }
    else
    {
        for (unsigned i = 0; i <= 1; ++i)
        {
            eraseBack(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')][i]);
            eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')][i]);
            eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')][i]);
            eraseBack(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')][i]);
        }
    }
}


template<typename TConstantSet, typename TStrand>
inline void
adjustConstantsSize(TConstantSet &constantSet, TStrand strand, LogFunction const &)
{
    adjustConstantsSize(constantSet, strand, NaiveMult());
}

template<typename TCoeffs>
inline TCoeffs polynomialMultipl(TCoeffs &coeffs1, TCoeffs &coeffs2)
{
    int n = length(coeffs1)-1;
    int m = length(coeffs2)-1;
    TCoeffs result;
    resize(result, n+m+1, 0.0, Exact());
    for (int k = 0; k <= m+n; ++k)
    {
        for (int i = 0; (i <= n && i <=k); ++i)
        {
            if (k-i <= m)                                   // clean up
                result[k] += coeffs1[i]*coeffs2[k-i];
        }
    }
    return result;
}


// Get coeffs for meth cases and (partial) likelihoods for other cases
template<typename TConstantSet, typename TLHoods, typename TQStrings, typename TMapqs, typename TOriginString, typename TCounts, typename TMethOptions, typename TRefContext, typename TMethod>
inline void
constructConstantsAndLHoods(TConstantSet &constantSet,
                            TLHoods &lHoods,
                            TQStrings &qualF, TQStrings &qualR,
                            TMapqs &mapqsF, TMapqs &mapqsR,
                            TOriginString & originStringF, TOriginString & originStringR,
                            TCounts &countF, TCounts &countR,
                            TMethOptions &methOptions,
                            TRefContext &refContext,
                            TMethod const &)
{
    unsigned minCountCT = 1;
    unsigned countF_CT = countF[ordValue((Dna)'C')] + countF[ordValue((Dna)'T')];
    unsigned countR_CT = countR[ordValue((Dna)'G')] + countR[ordValue((Dna)'A')];
    unsigned rF = 0;
    unsigned rR = 0;

    // ATTENTIONE: values for C are stored in constantSet[CG] and betas[CG]; values for G are stored in constantSet[GC] and betas[GC]

    // Get constants for reads on corresponding strand
    // For reads on other strand: caclulate already likelihood, since this is independent on beta
    // Multiply later
    // for each observed base type
    for (unsigned i = 0; i < 4; ++i)
    {
        // for all reads mapped on forward strand
        for (unsigned j = 0; j < length(qualF[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualF[i][j])-33);
            /*
            if (candidatePos + startCoord == 908640 || candidatePos + startCoord == 985089 )
            {
                std::cout << qual << std::endl;
            }
            */

            // If quality is below threshold, ignore read for all further calculations
            if (qual < 1 || (methOptions.useMapq && mapqsF[i][j] < 1))
            {
                adjustConstantsSize(constantSet, 'F', TMethod());
                continue;
            }
            long double e = pow(10.0L, (long double)(-qual/10.0));
            // for each possible candidate genotype
            // likelihood to observe single base under assumption of given genotype
            String<long double> singleProbs;
            resize(singleProbs, 6);
            for (unsigned h = 0; h < 6; ++h)
            {
                getSingleBaseProbHaploF(singleProbs[h], (Dna)i, (DnaM)h, e, originStringF[i][j], methOptions);
                if (methOptions.useMapq)    // somehow not practical ...? why? base quals? read specific error count?
                {
                    //singleProbs[h] *= mapqsF[i][j];
                    //if (countF_CT == 7) std::cerr << "Use mapq : " << mapqsF[i][j]  << std::endl;
                    long double e_m = pow(10.0L, (long double)(- mapqsF[i][j]/10.0));
                    singleProbs[h] = (1.0 -  e_m) * singleProbs[h];
                    double sim = 0.98;  // expected similarity
                    //  TODO how to deal with prob. to observe this base in wrong mapping ?
                    if ((Dna)i == 'C' && (Dna)refContext.refAllele == 'C') singleProbs[h] += e_m * (sim*0.5*0.25);
                    else if ((int)i == refContext.refAllele) singleProbs[h] += e_m * (sim*0.25);
                    else if ((Dna)i == 'T' && (Dna)refContext.refAllele == 'C') singleProbs[h] += e_m * (sim*0.5*0.25);
                    else if ((Dna)i == 'C') singleProbs[h] += e_m * ((1.0-sim)*0.5*0.25);
                    else if ((Dna)i == 'T') singleProbs[h] += e_m * ((1.0-sim)*1.5*0.25);
                    else singleProbs[h] += e_m * ((1.0-sim)*1.0*0.25);
                }
            }

            // Calculate likelihood for diploid genotypes
            for (unsigned h1 = 0; h1 < 4; ++h1)
            {
                for (unsigned h2 = h1; h2 < 4; ++h2)
                {
                    //std::cout << " test F" << (Dna)h1 << (Dna)h2 << " rF: " << rF << std::endl;
                    // Build up polynom coeffs
                    if ((Dna)h1 == 'C' && (Dna)h2 == 'C')   // CC
                    {
                        if (countF_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'C')];
                            long double pCm = singleProbs[ordValue((DnaM)'D')];
                            long double pOther = 0.0;
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rF, TMethod());
                        }
                    }
                    else if ( (Dna)h1 == 'C' && (Dna)h2 == 'G') // CG TODO same as CX
                    {
                        if ((countF_CT >= minCountCT) && (countR_CT >= minCountCT))
                        {
                            long double pC = singleProbs[ordValue((Dna)'C')];  // Beta independent of meth Level of Cs on other strand, hence to work with G is enough (?)
                            long double pCm = singleProbs[ordValue((DnaM)'D')];
                            long double pOther = singleProbs[ordValue((Dna)'G')];
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rF, TMethod());
                        }
                    }
                    else if ((Dna)h1 == 'C' || (Dna)h2 == 'C')  // CX
                    {
                        if (countF_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'C')];
                            long double pCm = singleProbs[ordValue((DnaM)'D')];
                            long double pOther = singleProbs[((Dna)h1=='C')? h2:h1];
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rF, TMethod());
                        }
                    }
                    else            // XX, no beta to maximize
                    {
                        long double p = 0.5*singleProbs[h1] + 0.5*singleProbs[h2];
                        lHoods[(h1<<2)|h2] *= p;
                    }
                }
            }
            ++rF;
        }
        // for all reads mapped on reverse strand
        for (unsigned j = 0; j < length(qualR[i]); ++j)
        {
            long double qual =  static_cast<long double>(ordValue(qualR[i][j])-33);
            // If quality is below threshold, ignore read for all further calculations
            if (qual < 1 || (methOptions.useMapq && mapqsR[i][j] < 1))
            {
                adjustConstantsSize(constantSet, 'R', TMethod());
                continue;
            }

            long double e = pow(10.0L, (long double)(-qual/10.0));
            // likelihood to observe single base under assumption of given genotype
            String<long double> singleProbs;
            resize(singleProbs, 6);
            // for each possible candidate genotype
            for (unsigned h = 0; h < 6; ++h)
            {
                getSingleBaseProbHaploR(singleProbs[h], (Dna)i, (DnaM)h, e, originStringR[i][j], methOptions);
                if (methOptions.useMapq)
                {
                    long double e_m = pow(10.0L, (long double)(- mapqsR[i][j]/10.0));
                    singleProbs[h] = (1.0 -  e_m) * singleProbs[h];
                    double sim = 0.98;
                    if ((Dna)i == 'G' && (Dna)refContext.refAllele == 'C') singleProbs[h] += e_m * (sim*0.5*0.25);
                    else if ((int)i == refContext.refAllele) singleProbs[h] += e_m * (sim*0.25);
                    else if ((Dna)i == 'A' && (Dna)refContext.refAllele == 'C') singleProbs[h] += e_m * (sim*0.5*0.25);
                    else if ((Dna)i == 'G') singleProbs[h] += e_m * ((1.0-sim)*0.5*0.25);
                    else if ((Dna)i == 'A') singleProbs[h] += e_m * ((1.0-sim)*1.5*0.25);
                    else singleProbs[h] += e_m * ((1.0-sim)*1.0*0.25);
              }
            }

            // Calculate likelihood for diploid genotypes
            for (int h1 = 0; h1 < 4; ++h1)
            {
                for (int h2 = h1; h2 < 4; ++h2)
                {
                    // Build up polynom coeffs
                    if ( (Dna)h1 == 'C' && (Dna)h2 == 'G') // CG
                    {
                        if ( (countF_CT >= minCountCT) && (countR_CT >= minCountCT))
                        {
                            // G
                            long double pC = singleProbs[ordValue((Dna)'G')];
                            long double pCm = singleProbs[ordValue((DnaM)'H')];
                            long double pOther = singleProbs[ordValue((DnaM)'C')];
                            addFactorToConstants(constantSet[(h2<<2)|h1], pC, pCm, pOther, rR, TMethod());
                        }
                    }
                    else if ((Dna)h1 == 'G' && (Dna)h2 == 'G')  // GG
                    {
                        if (countR_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'G')];
                            long double pCm = singleProbs[ordValue((DnaM)'H')];
                            long double pOther = 0.0;
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rR, TMethod());
                        }
                    }
                    else if ((Dna)h1 == 'G' || (Dna)h2 == 'G')   // GX
                    {
                        if (countR_CT >= minCountCT)
                        {
                            long double pC = singleProbs[ordValue((Dna)'G')];
                            long double pCm = singleProbs[ordValue((DnaM)'H')];
                            long double pOther = singleProbs[((Dna)h1=='G')? h2:h1];
                            addFactorToConstants(constantSet[(h1<<2)|h2], pC, pCm, pOther, rR, TMethod());
                        }
                    }
                    else            // XX, no beta to maximize
                    {
                        long double p = 0.5*singleProbs[h1] + 0.5*singleProbs[h2];
                        lHoods[(h1<<2)|h2] *= p;
                    }
                }
            }
            ++rR;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////


// Check if borders are better (for newton method)
template<typename TLHood, typename TBeta, typename TFunctor, typename TEvalMethod>
inline TBeta
verifyBeta(TLHood &lHood, TBeta &beta, TFunctor &functor, TEvalMethod const &)
{
    bool eBViolated = false;
    TLHood fBeta = functor(eBViolated, beta);
    if (eBViolated) fBeta = 0.0;

    TBeta b_0 = 0.0;
    TLHood f_0 = functor(eBViolated, b_0);
    if (eBViolated) f_0 = 0.0;

    TBeta b_1 = 1.0;
    TLHood f_1 = functor(eBViolated, b_1);
    if (eBViolated) f_1 = 0.0;

    if (fBeta >= f_0 && fBeta >= f_1)
    {
        lHood = fBeta;
        return beta;
    }
    else if (f_1 >= f_0)
    {
        lHood = f_1;
        return 1.0;
    }
    else
    {
        lHood = f_0;
        return 0.0;
    }
}


// NSpace
template<typename TBeta, typename TLHood, typename TConstants, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, TLHood &lHood, TConstants &constants, TBeta /*&guess*/, TMethOptions &/*methOptions*/, Sampling const &, NaiveMult const &)
{
#ifdef CALL_PROFILE
    double timeStamp = sysTime();
#endif
    FctNaive_0N<long double> fctNaiveMult(constants);

    TLHood maxlHood = 0.0;
    TBeta betaMax = 0.0;
    for (TBeta currBeta = 0.0; currBeta <= 1.0; currBeta+=0.01)       // TODO: choose betas dependend guess and not equally distributed
    {
        TLHood currlHood =  fctNaiveMult(currBeta);
        //std::cout << "curr Beta: " << currBeta << "  CurrlHood: " << currlHood << std::endl;
        if (currlHood > maxlHood)
        {
            betaMax = currBeta;
            maxlHood = currlHood;
        }
    }
    beta = betaMax;
    lHood = maxlHood;
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif
}


// NSpace, Log function
template<typename TBeta, typename TLHood, typename TConstants, typename TMethOptions>
inline void
getMaximizingBeta(TBeta &beta, TLHood &lHood, TConstants &constants, TBeta &guess, TMethOptions &methOptions, Newton const &, LogFunction const &)
{
#ifdef CALL_PROFILE
    double timeStamp = sysTime();
#endif
    typedef typename boost::math::tuple<TBeta, TBeta> TTuple;

    boost::uintmax_t maxIter = 10;
    //std::cout << " Guess: " << guess << std::endl;
    int digits = std::numeric_limits<long double>::digits;
    FctLog_12N<long double> myF12(constants);
    beta = boost::math::tools::newton_raphson_iterate(myF12, guess, (TBeta)0.0, (TBeta)1.0, digits, maxIter);
    // Check if f2 < 0 and get lHood
    FctLog_02N<long double> myF02(constants);
    TTuple tupleR = myF02(beta);
    long double f_2R = boost::math::get<1>(tupleR);
    if (f_2R >= 0)
    {
        getMaximizingBeta(beta, lHood, constants, guess, methOptions, Sampling(), NaiveMult());
        ++methOptions.countPlanB;
    }
    else
    {
        TLHood logLHood = boost::math::get<0>(tupleR);
        lHood = pow(10, logLHood);
        ++methOptions.countNoPlanB;
    }

    if (methOptions.helpPrint)   // certain pos
    {
        std::cout << "Log function: " << std::endl;
        for (TBeta b =0.0; b <= 1.0; b+=0.05)
        {
            TTuple t02 = myF02(b);
            TTuple t12 = myF12(b);
            std::cout << "beta: " << b << "\t f: " << boost::math::get<0>(t02) << "\t f': " << boost::math::get<0>(t12) << "\t f'': " <<  boost::math::get<1>(t12)  << std::endl;
        }
        std::cout << "beta: " << beta << " lHood: " << lHood << std::endl;
    }
#ifdef CALL_PROFILE
    Times::instance().time_optimization += (sysTime() - timeStamp);
#endif

}


template<typename TBetas, typename TLHoods, typename TConstantSet, typename TCounts, typename TMethOptions, typename TMethod, typename TEvalMethod>
inline void
getBetasAndLHoods(TBetas &betas, TLHoods &lHoods, TConstantSet &constantSet, TCounts &countF, TCounts &countR, TMethOptions &methOptions, TMethod const &, TEvalMethod const &, unsigned &pos)
{
    typedef typename Value<TLHoods>::Type   TLHood;

    unsigned minCountCT = 1;    // We do not have to limit this here, since this is taken into account in prob. calculations later anyway
    resize(betas, 4*4, 0.666, Exact());
    long double guess;
    unsigned countF_CT = countF[ordValue((Dna)'C')] + countF[ordValue((Dna)'T')];
    unsigned countR_CT = countR[ordValue((Dna)'G')] + countR[ordValue((Dna)'A')];

    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            // Compute beta values and assign corresponding probs to lHoods
            if ( ((Dna)h1 == 'C' && (Dna)h2 == 'G')) // CG
            {
                if ((countF_CT >= minCountCT) && (countR_CT >= minCountCT))                // TODO C/T threshold to take meth level into account ?? problem: could bias score, result
                {
                    // C: lHood from froward strand reads
                    guess = (long double)countF[ordValue((Dna)'C')]/(long double)(countF_CT);
                    TLHood lHood;
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    lHoods[(h1<<2)|h2] = lHood;
                    // G: lHood from reverse strand reads
                    guess = (long double)countR[ordValue((Dna)'G')]/(long double)(countR_CT);
                    getMaximizingBeta(betas[h2<<2|h1], lHood, constantSet[(h2<<2)|h1], guess, methOptions, TMethod(), TEvalMethod());
                    lHoods[(h1<<2)|h2] *= lHood;
                }
                else        // if not enough Cs and Ts: do not iterate and set prob to 0
                {
                    betas[h1<<2|h2] = 666.0;
                    betas[h2<<2|h1] = 666.0;
                    lHoods[(h1<<2)|h2] = 0.0;
               }
            }
            else if ( ((Dna)h1 == 'C' && (Dna)h2 == 'T'))
            {
                if (countF_CT >= minCountCT)
                {
					int diff = (int)countF[ordValue((Dna)'T')] - (int)countF_CT / 2;
                    if (diff > 0) guess = (long double)countF[ordValue((Dna)'C')]/(long double)(countF[ordValue((Dna)'C')] + diff);
                    else guess = 1.0;   // Set to 0.9 to avoid algorithm to get stuck?
                    TLHood lHood;
                    if (pos == 104)
                    {
                        methOptions.helpPrint = true;
                        std::cout << (Dna)h1 << (Dna)h2 << std::endl;
                    }
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    methOptions.helpPrint = false;
                    lHoods[(h1<<2)|h2] *= lHood;
               }
                else        // if not enough Cs and Ts: do not iterate and set prob to 0
                {
                    betas[h1<<2|h2] = 666.0;
                    lHoods[(h1<<2)|h2] = 0.0;
                }
            }
            else if ( ((Dna)h1 == 'A' && (Dna)h2 == 'G'))
            {
                if (countR_CT >= minCountCT)
                {
                    int diff = (int)countR[ordValue((Dna)'A')] - (int)countR_CT / 2;
                    if (diff > 0) guess = (long double)countR[ordValue((Dna)'G')]/(long double)(countR[ordValue((Dna)'G')] + diff);
                    else guess = 1.0;   // Set to 0.9 to avoid algorithm to get stuck?
                    TLHood lHood;
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    lHoods[(h1<<2)|h2] *= lHood;
               }
                else
                {
                    betas[h1<<2|h2] = 666.0;
                    lHoods[(h1<<2)|h2] = 0.0;
                }
            }
            else if ((Dna)h1 == 'C' || (Dna)h2 == 'C')
            {
                if (countF_CT >= minCountCT)
                {
                    guess = (long double)countF[ordValue((Dna)'C')]/(long double)(countF_CT);
                    TLHood lHood;
                    if (pos == 104 && (Dna)h1 == 'C' && (Dna)h2 == 'C')
                    {
                        methOptions.helpPrint = true;
                        std::cout << (Dna)h1 << (Dna)h2 << std::endl;
                    }
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    methOptions.helpPrint = false;
                    lHoods[(h1<<2)|h2] *= lHood;
               }
                else
                {
                    betas[h1<<2|h2] = 666.0;
                    lHoods[(h1<<2)|h2] = 0.0;
                }
            }
            else if ((Dna)h1 == 'G' || (Dna)h2 == 'G')
            {
                if (countR_CT >= minCountCT)
                {
                    guess = (long double)countR[ordValue((Dna)'G')]/(long double)(countR_CT);
                    TLHood lHood;
                    getMaximizingBeta(betas[h1<<2|h2], lHood, constantSet[(h1<<2)|h2], guess, methOptions, TMethod(), TEvalMethod());
                    lHoods[(h1<<2)|h2] *= lHood;
              }
                else
                {
                    betas[h1<<2|h2] = 666.0;
                    lHoods[(h1<<2)|h2] = 0.0;
               }
            }
        }
    }
}


template<typename TConstantSet, typename TLHoods, typename TCounts>
inline void
setUpConstants(TConstantSet &constantSet, TLHoods &lHoods, TCounts &countF, TCounts &countR, NaiveMult const &)
{
    clear(lHoods);
    resize(lHoods, 4*4, 1.0, Exact());

    clear(constantSet);
    resize(constantSet, 4*4, Exact());
    unsigned covF = countF[0] + countF[1] + countF[2] + countF[3];
    unsigned covR = countR[0] + countR[1] + countR[2] + countR[3];

    // ATTENTIONE: values for C are stored in constantSet[CG] and betas[CG]; values for G are stored in constantSet[GC] and betas[GC]

    // resize for as and bs
    for (unsigned i = 0; i <= 1; ++i)
    {
        // Forward strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')], 2, Exact());
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')], 2, Exact());
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')], 2, Exact());
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')], 2, Exact());
        // Reverse strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')], 2, Exact());
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')], 2, Exact());
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')], 2, Exact());
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')], 2, Exact());
    }
    // resize for number of reads
    for (unsigned j = 0; j <= 1; ++j)
    {
        // Forward strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'C')][j], covF, 0.0, Exact());
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'C')][j], covF, 0.0, Exact());
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'G')][j], covF, 0.0, Exact());
        resize(constantSet[ordValue((Dna)'C')<<2|ordValue((Dna)'T')][j], covF, 0.0, Exact());
        // Reverse strand methylations
        resize(constantSet[ordValue((Dna)'A')<<2|ordValue((Dna)'G')][j], covR, 0.0, Exact());
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'C')][j], covR, 0.0, Exact());
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'G')][j], covR, 0.0, Exact());
        resize(constantSet[ordValue((Dna)'G')<<2|ordValue((Dna)'T')][j], covR, 0.0, Exact());
    }
}

template<typename TConstantSet, typename TLHoods, typename TCounts>
inline void
setUpConstants(TConstantSet &constantSet, TLHoods &lHoods, TCounts &countF, TCounts &countR, LogFunction const &)
{
    setUpConstants(constantSet, lHoods, countF, countR, NaiveMult());
}


template<typename TPostProbs, typename TLHoods, typename TRefContext, typename TOptions, typename TMethOptions>
inline void
computePostProbs(TPostProbs &postProbs, TLHoods &lHoods, TRefContext &refContext, TOptions &options, TMethOptions &methOptions)
{
    // Calculate Pr(D)
    long double obsBasesProb = 0.0;
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            obsBasesProb += methOptions.genPriors[refContext.refAllele<<4| ordValue((Dna)h1)<<2| ordValue((Dna)h2)] * lHoods[(h1<<2)| h2];
         }
    }
    // Calculate posterior probs. for possible genotypes
    // Pr(G|D) = pi(G) * Pr(D|G) / Pr(D9)
    // For bs genotypes take prior methylation/non-methylation probability into account, dependent on context
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            if(options._debugLevel > 1)
                std::cout << (Dna)refContext.refAllele << (Dna)h1 << (Dna)h2 << "  " << std::setprecision (25) << "genPrior: " << std::setprecision (25) << (long double)methOptions.genPriors[ refContext.refAllele<<4| ordValue((Dna)((DnaM)h1))<<2 | ordValue((Dna)((DnaM)h2))] << "  lHood  " << (long double)lHoods[(h1<<2)| h2] << std::endl;

            // future: maybe take bsPriors into account, take context into account !!!
            postProbs[(h1<<2)| h2] = methOptions.genPriors[ refContext.refAllele<<4| h1<<2 | h2] * lHoods[(h1<<2)| h2] / obsBasesProb;

            if(options._debugLevel > 1)
                 std::cout << std::setprecision (25) << "candidateProb: " << (Dna)refContext.refAllele << (Dna)h1 << (Dna)h2 << "  " << postProbs[(h1<<2)| h2] << "  context" << refContext.contextF << refContext.contextR << std::endl;
        }
    }
}


template<typename TProbs, typename TBetas, typename TMethOptions, typename TOptions, typename TQStrings, typename TMapqs, typename TOriginString, typename TCounts, typename TRefContext>
inline void
getCandidateProbs(TProbs &postProbs, TBetas &betas,
                  TMethOptions &methOptions, TOptions &options,
                  TQStrings &qualF, TQStrings &qualR,
                  TMapqs &mapqsF, TMapqs &mapqsR,
                  TOriginString & originStringF,
                  TOriginString & originStringR,
                  TCounts &countF, TCounts &countR,
                  TRefContext &refContext)
{
    if (methOptions.ignoreBs) // Snp calling without bs conversions
        methOptions.convRate = 0.0;

    String<long double> lHoods;  // likelihoods to observe observed data under assumption of given genotypes
    String<String<String<long double> > > constantSet;

    // NaiveMult and Sampling Method
    if (methOptions.betaSampling)           // sampling, NaiveMult, NSpace
    {
        setUpConstants(constantSet, lHoods, countF, countR, NaiveMult());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, refContext, NaiveMult());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Sampling(), NaiveMult(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
    else  // Newton, logFunction, NSpace
    {
        //std::cout << " Run with LogFunction ...................................................!" << std::endl;
        setUpConstants(constantSet, lHoods, countF, countR, LogFunction());
        constructConstantsAndLHoods(constantSet, lHoods, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, methOptions, refContext, LogFunction());
        getBetasAndLHoods(betas, lHoods, constantSet, countF, countR, methOptions, Newton(), LogFunction(), refContext.pos);
        computePostProbs(postProbs, lHoods, refContext, options, methOptions);
    }
}

///////////////////////////////////////////////////////////////////////
// SNP and meth calling in one step
template<typename TCounts, typename TQualities, typename TMapqs, typename TOriginString, typename TRefContext, typename TMethOptions, typename TOptions, typename TMethylVariant>
inline bool
doBsCalling(TCounts & countF,
          TCounts & countR,
          TQualities & qualF,
          TQualities & qualR,
          TMapqs & mapqsF,
          TMapqs & mapqsR,
          TOriginString & originStringF,
          TOriginString & originStringR,
          TRefContext & refContext,
          TMethOptions &methOptions,
          TOptions & options,
          TMethylVariant &meth
          )
{
    int genotypeRef = (refContext.refAllele<<2) | refContext.refAllele;

    String<long double> candidateProbs;
    String<long double> betas;
    resize(candidateProbs, 4*4); // for simplicity; not all are used

    getCandidateProbs(candidateProbs, betas, methOptions, options, qualF, qualR, mapqsF, mapqsR, originStringF, originStringR, countF, countR, refContext);
    //Choose genotype which maximizes the posterior prob.

    int genotype1 = genotypeRef;
    int allele1 = 666;
    int allele2 = 666;
    int genotype2 = genotypeRef;
    long double maxProb1 = 0.0;
    long double maxProb2 = 0.0;
    for (int h1 = 0; h1 < 4; ++h1)
    {
        for (int h2 = h1; h2 < 4; ++h2)
        {
            /*if (refContext.pos == 30032)
            {
                std::cout << std::setprecision (50) << "current genotype: " << (Dna)h1 << (Dna)h2 <<  "probs: " << candidateProbs[(h1<<2)| h2] << std::endl;
            }*/
            if (candidateProbs[(h1<<2)| h2] >= maxProb1)
            {
                maxProb2 = maxProb1;
                genotype2 = (allele1<<2)|allele2;
                maxProb1 = candidateProbs[(h1<<2)| h2];
                allele1 = h1;
                allele2 = h2;
            }
            else if (candidateProbs[(h1<<2)| h2] >= maxProb2)
            {
                maxProb2 = candidateProbs[(h1<<2)| h2];
                genotype2 = (h1<<2)|h2;
            }
        }
    }
    genotype1 = (allele1<<2)|allele2;

    //unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
                           //+ countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
    meth.genotype = genotype1;
    meth.score = log( (long double)candidateProbs[genotype1]/ (long double)candidateProbs[genotype2]);     // genotype calling score (only underlying genotypes with best beta, no bs types)?

    if (meth.score <= methOptions.minScoreToCallSnp)
    {
        meth.genotypeCalled = false;
        ++methOptions.countScoreTooLow;
    }
    else
        meth.genotypeCalled = true;

    /*if (refContext.pos == 693)
    {
        std::cout << " Pos 693" << std::endl;
        std::cout << "Score: " << meth.score << " msc: " << methOptions.minScoreToCallSnp << std::endl;
        std::cout << " Called: " << (meth.genotypeCalled ? "true":"false") << std::endl;
    }*/

    meth.genotypeProb = candidateProbs[genotype1];
    meth.methLevel1 = betas[genotype1];
    if ((Dna)allele1 == 'C' && (Dna)allele2 == 'G')
        meth.methLevel2 = betas[(ordValue((Dna)'G')<<2)|ordValue((Dna)'C')];

    // If bs case:
    // Genotype was called = score was good enough to call
    if (meth.genotypeCalled && ((Dna)allele1 == 'C' || (Dna)allele1 == 'G' || (Dna)allele2 == 'C' || (Dna)allele2 == 'G') )     // TODO: think how to hand call threshold
    {
        meth.bsCalled = true;
    }
    else
        meth.bsCalled = false;


    /*if (refContext.pos == 30032 )
    {
        std::cout << "Pos: " << refContext.pos << std::endl;
        std::cout << std::setprecision (25) << " prob genotype1.." << (long double)candidateProbs[genotype1] << " prob genotype2.." << (long double)candidateProbs[genotype2]  <<  std::endl;
        std::cout << " genotype1  allele1: "<< (Dna)(genotype1>>2) << "  allele2: " << (Dna)(genotype1 % 4) << "  beta: " << betas[genotype1] << std::endl;
        std::cout << " genotype2  allele1: "<< (Dna)(genotype2>>2) << "  allele2: " << (Dna)(genotype2 % 4) << "  beta: " << betas[genotype2] << std::endl;
        std::cout << std::setprecision (50) << " prob genotype CG.." << (long double)candidateProbs[ordValue((Dna)'C')<<2|ordValue((Dna)'G')]  << "  beta: " << betas[ordValue((Dna)'C')<<2|ordValue((Dna)'G')] <<  std::endl;
        //std::cout << std::setprecision (25) << " prob genotype3.." << (long double)candidateProbs[0>>2|2]  <<  std::endl;
    }*/

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Output
////////////////////////////////////////////////////////////////////////////////////////////
// VCF + bed
template<typename TVcfStream, typename TBedStream, typename TMethylVariant, typename TQualities, typename TRefContext, typename TMethOptions, typename TOptions>
inline bool
writeMeth(TVcfStream &vcfStream,
       TBedStream &bedStream,
       TMethylVariant &meth,
       TQualities &qualityStringF,
       TQualities &qualityStringR,
       TRefContext &refContext,
       unsigned realCoverage,
       TMethOptions &methOptions,
       TOptions &/*options*/)
{
#ifdef CALL_PROFILE
    double timeStamp = sysTime();
#endif

    // VCF
     if(meth.genotypeCalled && ((meth.genotype>>2) != refContext.refAllele || (meth.genotype%4) != refContext.refAllele) )
     {
        VcfRecord record;
        record.rID = 0; // regarding temp vcfStream // refContext.genomeID;
        record.beginPos = refContext.pos;
        record.id = "snp";
        resize(record.ref, 1);
        record.ref[0] = (Dna)refContext.refAllele;
        CharString alt;
        if ((meth.genotype>>2) != refContext.refAllele && (meth.genotype%4) != refContext.refAllele)
        {
            appendValue(alt, (Dna)(meth.genotype>>2));
            appendValue(alt, ',');
            appendValue(alt, (Dna)(meth.genotype%4));
        }
        else if ((meth.genotype>>2) != refContext.refAllele)
        {
            appendValue(alt, (Dna)(meth.genotype>>2));
        }
        else
        {
            appendValue(alt, (Dna)(meth.genotype%4));
        }

        record.alt = alt;
        record.qual = '.';
        record.filter = "PASS";

        record.info = "FB:";
        std::stringstream ss;
        ss << length(qualityStringF[0]);   // Top: A
        CharString str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ",");
        ss <<  length(qualityStringF[1]);   // Top: C
        str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ",");
        ss <<  length(qualityStringF[2]);   // Top: G
        str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ",");
        ss <<  length(qualityStringF[3]);   // Top: T
        str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ";");
        append(record.info, "RB=");
        ss <<  length(qualityStringR[0]);   // Bot: A
        str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ",");
        ss <<  length(qualityStringR[1]);   // Bot: C
        str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ",");
        ss <<  length(qualityStringR[2]);   // Bot: G
        str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ",");
        ss <<  length(qualityStringR[3]);   // Bot: T
        str = ss.str();
        ss.str("");
        ss.clear();
        append(record.info, str);
        append(record.info, ";");
        // TODO vv qualities

        record.format =  "GT:DP:GS:MP";
        CharString genotypeInfos;
        if ((meth.genotype>>2) != refContext.refAllele && (meth.genotype%4) != refContext.refAllele)
        {
            append(genotypeInfos, "1|2:");
        }
        else
        {
            append(genotypeInfos, "0|1:");
        }
        ss <<  realCoverage;
        str = ss.str();
        ss.str("");
        ss.clear();
        append(genotypeInfos, str);
        append(genotypeInfos, ":");
        ss <<  meth.score;
        str = ss.str();
        ss.str("");
        ss.clear();
        append(genotypeInfos, str);
        append(genotypeInfos, ":");
        ss <<  meth.genotypeProb;
        str = ss.str();
        ss.str("");
        ss.clear();
        append(genotypeInfos, str);
        appendValue(record.genotypeInfos, genotypeInfos);

        writeRecord(vcfStream, record);
    }

    // BED
    if (meth.bsCalled)
    {
        BedRecord<Bed6> bedRecord;
        bedRecord.ref = contigNames(context(vcfStream))[0];
        bedRecord.beginPos = refContext.pos;
        bedRecord.endPos = refContext.pos + 1;
        bedRecord.name = ".";
        std::stringstream ss;
        ss <<  meth.score;
        CharString str = ss.str();
        ss.str("");
        ss.clear();
        bedRecord.score = str;

        if ((Dna)(meth.genotype>>2)  == 'C' && (Dna)(meth.genotype%4) == 'G')
        {
            bedRecord.strand = '+';
            ss <<  meth.methLevel1;
            str = ss.str();
            ss.str("");
            ss.clear();
            bedRecord.data = str;
            writeRecord(bedStream, bedRecord);
            bedRecord.strand = '-';
            ss << meth.methLevel2;
            str = ss.str();
            ss.str("");
            ss.clear();
            bedRecord.data = str;
            writeRecord(bedStream, bedRecord);

            // Some stats...
            if (refContext.contextF == 0)
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel1;
            }
            else if (refContext.contextF == 1)
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel1;
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel1;
            }
            if (refContext.contextR == 0)
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel2;
            }
            else if (refContext.contextR == 1)
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel2;
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel2;
            }
        }
        else if ((Dna)(meth.genotype>>2)  == 'C')
        {
            bedRecord.strand = '+';
            ss <<  meth.methLevel1;
            str = ss.str();
            ss.str("");
            ss.clear();
            bedRecord.data = str;
            writeRecord(bedStream, bedRecord);

            if (refContext.contextF == 0)
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel1;
            }
            else if (refContext.contextF == 1)
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel1;
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel1;
            }
        }
        else if ((Dna)(meth.genotype>>2)  == 'G')
        {
            bedRecord.strand = '-';
            ss <<  meth.methLevel1;
            str = ss.str();
            ss.str("");
            ss.clear();
            bedRecord.data = str;
            writeRecord(bedStream, bedRecord);

            if (refContext.contextR == 0)
            {
                ++methOptions.countCG;
                methOptions.statsCGMethylated += meth.methLevel1;
            }
            else if (refContext.contextR == 1)
            {
                ++methOptions.countCHG;
                methOptions.statsCHGMethylated += meth.methLevel1;
            }
            else
            {
                ++methOptions.countCHH;
                methOptions.statsCHHMethylated += meth.methLevel1;
            }
        }
    }

#ifdef CALL_PROFILE
    Times::instance().time_IO += (sysTime() - timeStamp);
#endif
    return true;
}



#endif  // #ifndef SANDBOX_KRAKAU_APPS_SNP_METH_STORE_BS_ONE_CALLING_H_
