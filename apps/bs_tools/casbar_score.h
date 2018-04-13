#ifndef APPS_BS_TOOLS_CASBAR_SCORE_H__
#define APPS_BS_TOOLS_CASBAR_SCORE_H__

namespace seqan {

struct BsTop_;
typedef Tag<BsTop_> const BsTop;

struct BsBottom_;
typedef Tag<BsBottom_> const BsBottom;

// TODO Rename to top and bottom
struct BsProfileScoreCT_;
typedef Tag<BsProfileScoreCT_> const BsProfileScoreCT;

struct BsProfileScoreCTRight_;
typedef Tag<BsProfileScoreCTRight_> const BsProfileScoreCTRight;    // Top strand, right mate

struct BsProfileScoreGA_;
typedef Tag<BsProfileScoreGA_> const BsProfileScoreGA;

struct BsProfileScoreGARight_;
typedef Tag<BsProfileScoreGARight_> const BsProfileScoreGARight;    // Bottom strand, right mate

struct BsProfileScoreRef_;
typedef Tag<BsProfileScoreRef_> const BsProfileScoreRef;

//////////////////////////////////////////////////////////////////////////////
// Convenience
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

//////////////////////////////////////////////////////////////////////////////


template <typename TValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor1>
class Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor1> >
{
public:
    String<String<TValue> > targetFreqs;    // New target frequencies for each position in profile
    String<TValue>          gapFreqs;    // New gap frequencies for each position in profile

    TValue const *readBaseFreqs;
    TValue const *refBaseFreqs;

    TValue const *seqErrorFreqs;
    TValue const *delErrorFreqs;
    TValue const *insErrorFreqs;

    TValue delErrorRate;    // sequencing del error rate (over whole read)
    TValue delRate;         // deletion rate (from reference, e.g. corresponding mason)
    TValue insErrorRate;     // percentage of seq. errors being insertions (assuming seq error)
    TValue endGapScore;
    TValue scoreLimit;

    TValue scalingFactorDelErrors;

    template <typename TOptions>
	Score(TOptions &options,
	        TValue const * &_seqErrorFreqs,
            TValue const * &_insErrorFreqs,
            TValue const * &_delErrorFreqs,
            TValue _scalingFactorDelErrors)
	{
        (*this).delRate = options.delRate;
        (*this).delErrorRate = options.delErrorRate;
        (*this).insErrorRate = options.insErrorRate;
        (*this).scalingFactorDelErrors = _scalingFactorDelErrors;
        (*this).endGapScore = options.endGapScore;
        (*this).scoreLimit = options.scoreLimit;

	    (*this).readBaseFreqs = ReadBaseFreqs<TValue, TModel>::getData();
        (*this).refBaseFreqs = RefBaseFreqs<TValue, TModel>::getData();

        (*this).seqErrorFreqs = _seqErrorFreqs;
        (*this).delErrorFreqs = _insErrorFreqs;
        (*this).insErrorFreqs = _delErrorFreqs;

    }

    template <typename TCellDescriptor2>
	Score(Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor2> > const &_score)
	{
	    (*this).targetFreqs = _score.targetFreqs;
	    (*this).gapFreqs = _score.gapFreqs;

        (*this).delRate = _score.delRate;
        (*this).delErrorRate = _score.delErrorRate;
        (*this).insErrorRate = _score.insErrorRate;
        (*this).scalingFactorDelErrors = _score.scalingFactorDelErrors;
        (*this).endGapScore = _score.endGapScore;
        (*this).scoreLimit = _score.scoreLimit;

	    (*this).readBaseFreqs = _score.readBaseFreqs;
        (*this).refBaseFreqs = _score.refBaseFreqs;

        (*this).seqErrorFreqs = _score.seqErrorFreqs;
        (*this).delErrorFreqs = _score.insErrorFreqs;
        (*this).insErrorFreqs = _score.delErrorFreqs;
    }
};


// For aligning read against top strand
template <typename TValue, typename TString, typename TBsProfileScore, typename TModel, typename TCellDescriptor, typename TMethOptions>
inline void
assignTargetFreqs(Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> > &me,
			      TString const & profile, TMethOptions &/*options*/, BsTop const &)
{
	typedef typename Size<TString>::Type TSize;
	resize(me.targetFreqs, length(profile));
	resize(me.gapFreqs, length(profile));

	for(TSize i = 0; i < length(profile); ++i) {
	    resize(me.targetFreqs[i], 5);
		
	    double sum = 0; // For normalization, since the sum of all freqs should be 1
        for (TSize j = 0; j < 4; ++j)
            sum += profile[i].count[j];


        if (sum != 0 || (int)profile[i].count[8] != -1)
        {
            sum += 1;    // For reference assuming max qual and max mapq, can be base or gap
            sum += profile[i].count[10];   // Take top gap count into account
            me.targetFreqs[i][0] = profile[i].count[0]/sum; // A
            me.targetFreqs[i][1] = profile[i].count[1]/sum; // C
            me.targetFreqs[i][2] = profile[i].count[2]/sum; // G
            me.targetFreqs[i][3] = profile[i].count[3]/sum; // T
        }
        else    // Avoid target frequencies beeing 0 (on other strand are reads mapped, otherwise this column would have been removed)
        {
            sum += 1;    // For reference assuming max qual and max mapq, can be base or gap
            sum += profile[i].count[10];   // Take top gap count into account
            me.targetFreqs[i][0] = 1.0/4.0;
            me.targetFreqs[i][1] = 1.0/4.0;
            me.targetFreqs[i][2] = 1.0/4.0;
            me.targetFreqs[i][3] = 1.0/4.0;
        }

        // Note: Seems impossible to take information from bases mapped to bottom strand into account
        // since we can't assume that ref base is correct at this position
        // and we don't know if these are not caused by bs conversions on bottom strand
        // (otherwise we could take e.g. count of C(G)s, T(A)s into account to get an idea, if Ts are caused by s conversion, problem with A(T)s etc. then)
        // Sinse realign mainly orders gaps new, and gaps are taken from both strands into account, this should be such a big problem (?)
        // Use only top or bottom reads for new target frequencies

        // Add reference values to corresponding freqs:
        TValue estMethLevel = 0;

        if (profile[i].count[1]+profile[i].count[3] > 0) estMethLevel = 0.5 * profile[i].count[1]/(profile[i].count[1]+profile[i].count[3]);   // bs conversion rate, to avoid C punished too much, if only Ts observed

        TValue refGap = 0;
        switch ((int)profile[i].count[8])
        {
            case 0: // Ref A
                me.targetFreqs[i][0] += 1.0/sum;
                break;
            case 1: // Ref C
                me.targetFreqs[i][1] += estMethLevel/sum;           // -> C (dep. on estimated meth. level)
                me.targetFreqs[i][3] += (1.0 - estMethLevel)/sum;   // -> T
                break;
            case 2: // Ref G
                me.targetFreqs[i][2] += 1.0/sum;
                break;
            case 3: // Ref T
                me.targetFreqs[i][3] += 1.0/sum;
                break;
            case 4: // Ref N
                me.targetFreqs[i][0] += (1.0/sum)/4.0;                  // A
                me.targetFreqs[i][1] += (estMethLevel/sum)/4.0;         // C
                me.targetFreqs[i][2] += (1.0/sum)/4.0;                  // G
                me.targetFreqs[i][3] += ((1.0-estMethLevel)/sum)/4.0;   // T
                SEQAN_FALLTHROUGH
            case -1: // Ref gap
                refGap = 1;     // Helper to take ref gap into account
                break;
        }
        me.gapFreqs[i] = (profile[i].count[10] + profile[i].count[11] + refGap)/(sum +       // sum: bases and gaps on F + ref
                                                                                 profile[i].count[4] + profile[i].count[5] + profile[i].count[6] + profile[i].count[7] +     // bases on R
                                                                                 profile[i].count[11]     // gaps on R
                                                                                ); // gaps rate with information from both strands
    }
}


// For aligning read against bottom strand
// We look at reverse complements of reads aligned against reverse strand
// -> deal with GA case
template <typename TValue, typename TString, typename TBsProfileScore, typename TModel, typename TCellDescriptor, typename TMethOptions>
inline void
assignTargetFreqs(Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> > &me,
			      TString const & profile, TMethOptions &/*options*/, BsBottom const &)
{
	typedef typename Size<TString>::Type TSize;
	resize(me.targetFreqs, length(profile));
	resize(me.gapFreqs, length(profile));

	for(TSize i = 0; i < length(profile); ++i) {
	    resize(me.targetFreqs[i], 5);
		
	    double sum = 0; // For normalization, since the sum of all freqs should be 1
        for (TSize j = 4; j < 8; ++j)
            sum += profile[i].count[j];


        if (sum != 0 || (int)profile[i].count[8] != -1)
        {
            sum += 1;    // For reference assuming max qual and max mapq, can be base or gap
            sum += profile[i].count[11];   // Take bottom gap count into account
            me.targetFreqs[i][0] = profile[i].count[4]/sum; // A
            me.targetFreqs[i][1] = profile[i].count[5]/sum; // C
            me.targetFreqs[i][2] = profile[i].count[6]/sum; // G
            me.targetFreqs[i][3] = profile[i].count[7]/sum; // T
       }
        else    // Avoid target frequencies beeing 0
        {
            sum += 1;    // For reference assuming max qual and max mapq, can be base or gap
            sum += profile[i].count[11];   // Take bottom gap count into account
            me.targetFreqs[i][0] = 1.0/4.0;
            me.targetFreqs[i][1] = 1.0/4.0;
            me.targetFreqs[i][2] = 1.0/4.0;
            me.targetFreqs[i][3] = 1.0/4.0;
        }

        // Add reference values to corresponding freqs:
        TValue estMethLevel = 0;
        if (profile[i].count[6]+profile[i].count[4] > 0) estMethLevel = 0.5 * profile[i].count[6]/(profile[i].count[6]+profile[i].count[4]);   // bs conversion rate, to avoid C punished too much, if only Ts observed
        TValue refGap = 0;
        switch ((int)profile[i].count[8])
        {
            case 0: // Ref A
                me.targetFreqs[i][0] += 1.0/sum;
                break;
            case 1: // Ref C
                me.targetFreqs[i][1] += 1.0/sum;
                break;
            case 2: // Ref G
                me.targetFreqs[i][2] += estMethLevel/sum;       // -> G
                me.targetFreqs[i][0] += (1.0-estMethLevel)/sum; // -> A
                break;
            case 3: // Ref T
                me.targetFreqs[i][3] += 1.0/sum;
                break;
            case 4: // Ref N
                me.targetFreqs[i][0] += ((1.0-estMethLevel)/sum)/4.0;   // A
                me.targetFreqs[i][1] += (1.0/sum)/4.0;                  // C
                me.targetFreqs[i][2] += (estMethLevel/sum)/4.0;         // G
                me.targetFreqs[i][3] += (1.0/sum)/4.0;                  // T
                SEQAN_FALLTHROUGH
            case -1: // Ref gap
                refGap = 1;
                break;
        }
        me.gapFreqs[i] = (profile[i].count[10] + profile[i].count[11] + refGap)/(sum +       // sum: bases and gaps on R + ref
                                                                                 profile[i].count[0] + profile[i].count[1] + profile[i].count[2] + profile[i].count[3] +     // bases on F
                                                                                 profile[i].count[10]     // gaps on F
                                                                                ); // gaps rate with information from both strands
   }
}


// For aligning ref
template <typename TValue, typename TString, typename TModel, typename TCellDescriptor, typename TMethOptions>
inline void
assignTargetFreqs(Score<TValue, BsTagList<BsProfileScoreRef, TModel, TCellDescriptor> > &me,
			      TString const & profile, TMethOptions &/*Options*/)
{
	typedef typename Size<TString>::Type TSize;
	resize(me.targetFreqs, length(profile));
	resize(me.gapFreqs, length(profile));

	for(TSize i = 0; i < length(profile); ++i) {
	    resize(me.targetFreqs[i], 2*5); // TODO size not necessary anymore
		
	    double sumF = 0;
	    double sumR = 0;
	    sumF += profile[i].count[0];    // TODO think: sum over all?
        sumF += profile[i].count[2];
        sumR += profile[i].count[5];
        sumR += profile[i].count[7];
        // Do not sum up to 1 anymore ! -> /2
        if (sumF != 0 || sumR != 0)
        {
            sumF += profile[i].count[10];     // Gap count top
            sumR += profile[i].count[11];     // Gap count bottom
            me.targetFreqs[i][0] = (profile[i].count[0]/sumF)/2.0; // A
            me.targetFreqs[i][2] = (profile[i].count[2]/sumF)/2.0; // G
            me.targetFreqs[i][5] = (profile[i].count[5]/sumR)/2.0; // C
            me.targetFreqs[i][7] = (profile[i].count[7]/sumR)/2.0; // T
        }
        else    // Avoid target frequencies beeing 0 (some read are mapped, otherwise this column would have been removed, but at C/T on F or G/A on R)
        {
            sumF += profile[i].count[10];     // Gap count top
            sumR += profile[i].count[11];     // Gap count bottom
            me.targetFreqs[i][0] = 1.0/4.0;
            me.targetFreqs[i][2] = 1.0/4.0;
            me.targetFreqs[i][5] = 1.0/4.0;
            me.targetFreqs[i][7] = 1.0/4.0;
        }
        // A bit dirty
        me.gapFreqs[i] = (profile[i].count[10] + profile[i].count[11])/(sumF + sumR); // gaps rate with information from both strands
    }
}


// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore
// --------------------------------------------------------------------------

template <typename TValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor, typename TSequence>
struct SequenceEntryForScore<Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> >, TSequence>    // To avoid conflict with other bs scores, use specific one here
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template <typename TValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor, typename TSequence>
struct SequenceEntryForScore<Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> > const, TSequence> :
       SequenceEntryForScore<Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> >, TSequence>
{};

template <typename TScoreValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> > const &, TSequence const & seq, TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// se-reads: GA
// seqErrorFreqs <- complement
// delErrorFreqs <- complement
// insErrorFreqs <- complement
// readBaseFreqs <- complement

// pe-reads: right mates
// top, mapped against reverse strand
// readbaseFreq: /
// seqErrorFreq: fCompl
// bottom, mapped against forward strand
// readbaseFreq: fCompl
// seqErrorFreq: /


// Modify to use different scoring function at first and last row
// (end gaps score different)
// Computes the score and tracks it if enabled.
template <typename TDPScout, typename TTraceMatrixNavigator,
          typename TDPCell,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoreValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor2,
          typename TColumnDescriptor,
          typename TDPProfile>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TDPCell & c,
             TDPCell & d,
             TDPCell const & h,
             TDPCell & v,
             TSequenceHValue const & seqHVal,
             TSequenceVValue const & seqVVal,
             Score<TScoreValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor2> > const & scoringScheme,
             TColumnDescriptor const &,
             FirstCell const &,   // One of FirstCell, InnerCell or LastCell.
             TDPProfile const &)
{
    typedef FirstCell TCellDescriptor;
    typedef DPMetaColumn_<TDPProfile, TColumnDescriptor> TMetaColumn;

    Score<TScoreValue, BsTagList<TBsProfileScore, TModel, FirstCell> > scoringSchemeDummy(scoringScheme);
    assignValue(traceMatrixNavigator,
                _computeScore(c, d, h, v, seqHVal, seqVVal,
                              scoringSchemeDummy, typename RecursionDirection_<TMetaColumn, TCellDescriptor>::Type(),
                              TDPProfile()));
//	std::cout << "("<< activeCell._score << "," << previousDiagonal._score << "," << previousHorizontal._score << "," << previousVertical._score << ") ";
    if (TrackingEnabled_<TMetaColumn, TCellDescriptor>::VALUE)
    {
        bool isLastColumn = IsSameType<typename TColumnDescriptor::TColumnProperty, DPFinalColumn>::VALUE;
        bool isLastRow = And<IsSameType<TCellDescriptor, LastCell>,
                             Or<IsSameType<typename TColumnDescriptor::TLocation, PartialColumnBottom>,
                                IsSameType<typename TColumnDescriptor::TLocation, FullColumn> > >::VALUE;
        _setVerticalScoreOfCell(c, _verticalScoreOfCell(v));
        _scoutBestScore(scout, c, traceMatrixNavigator, isLastColumn, isLastRow);
    }
}

// Modify to use different scoring function at first and last row
// (end gaps score different)
// Computes the score and tracks it if enabled.
template <typename TDPScout, typename TTraceMatrixNavigator,
          typename TDPCell,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoreValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor2,
          typename TColumnDescriptor,
          typename TDPProfile>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TDPCell & c,
             TDPCell & d,
             TDPCell const & h,
             TDPCell & v,
             TSequenceHValue const & seqHVal,
             TSequenceVValue const & seqVVal,
             Score<TScoreValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor2> > const & scoringScheme,
             TColumnDescriptor const &,
             LastCell const &,   // One of FirstCell, InnerCell or LastCell.
             TDPProfile const &)
{
    typedef LastCell TCellDescriptor;
    typedef DPMetaColumn_<TDPProfile, TColumnDescriptor> TMetaColumn;

    Score<TScoreValue, BsTagList<TBsProfileScore, TModel, LastCell> > scoringSchemeDummy(scoringScheme);
    assignValue(traceMatrixNavigator,
                _computeScore(c, d, h, v, seqHVal, seqVVal,
                              scoringSchemeDummy, typename RecursionDirection_<TMetaColumn, TCellDescriptor>::Type(),
                              TDPProfile()));
//	std::cout << "("<< activeCell._score << "," << previousDiagonal._score << "," << previousHorizontal._score << "," << previousVertical._score << ") ";
    if (TrackingEnabled_<TMetaColumn, TCellDescriptor>::VALUE)
    {
        bool isLastColumn = IsSameType<typename TColumnDescriptor::TColumnProperty, DPFinalColumn>::VALUE;
        bool isLastRow = And<IsSameType<TCellDescriptor, LastCell>,
                             Or<IsSameType<typename TColumnDescriptor::TLocation, PartialColumnBottom>,
                                IsSameType<typename TColumnDescriptor::TLocation, FullColumn> > >::VALUE;
        _setVerticalScoreOfCell(c, _verticalScoreOfCell(v));
        _scoutBestScore(scout, c, traceMatrixNavigator, isLastColumn, isLastRow);
    }
}


////////////////////////////////////////////////////////////////////////////
// scoreGapExtendHorizontal
////////////////////////////////////////////////////////////////////////////
// New gap in read
// After entry2

// For end gaps
template <typename TValue, typename TBsProfileScore, typename TModel, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, BsTagList<TBsProfileScore, TModel, FirstCell> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
	return std::log10(me.endGapScore);
}
// For end gaps
template <typename TValue, typename TBsProfileScore, typename TModel, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, BsTagList<TBsProfileScore, TModel, LastCell> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
	return std::log10(me.endGapScore);
}

// Top strand, original (forward strand)
template <typename TValue, typename TModel, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, BsTagList<BsProfileScoreCT, TModel, InnerCell> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    TValue sc = (me.gapFreqs[position(entry1)] / me.delRate) * (1.0-me.delErrorRate);                                   // P(X-)/f(-) * (1-delErrorRate)
    for (unsigned i = 0; i < 4; ++i)
    {
        sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[i]) * me.delErrorRate * me.delErrorFreqs[i]; //  * me.scalingFactorDelErrors;      // P(Xb)/f(b) * derErrorRate * seqError(b -> -)
    }
	return  ((std::log10(sc/2.0) > -10)? std::log10(sc/2.0):-10);
}

// Top strand, rev. compl. of original (projected from reverse to forward strand)
template <typename TValue, typename TModel, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, BsTagList<BsProfileScoreCTRight, TModel, InnerCell> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    TValue sc = (me.gapFreqs[position(entry1)] / me.delRate) * (1.0-me.delErrorRate);
    for (unsigned i = 0; i < 4; ++i)
    {
        sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[i]) * me.delErrorRate * me.delErrorFreqs[fCompl(i)]; // * me.scalingFactorDelErrors;     No scaling, because we sum up!
    }
	return  ((std::log10(sc/2.0) > -10)? std::log10(sc/2.0):-10);
}

// Bottom strand, original (projected from reverse to forward strand)
template <typename TValue, typename TModel, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, BsTagList<BsProfileScoreGA, TModel, InnerCell> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    TValue sc = (me.gapFreqs[position(entry1)] / me.delRate) * (1.0-me.delErrorRate);
    for (unsigned i = 0; i < 4; ++i)
    {
        sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[fCompl(i)]) * me.delErrorRate * me.delErrorFreqs[fCompl(i)]; //  * me.scalingFactorDelErrors;
    }
	return  ((std::log10(sc/2.0) > -10)? std::log10(sc/2.0):-10);

}

// Bottom strand, rev. compl. of original (forward strand)
template <typename TValue, typename TModel, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, BsTagList<BsProfileScoreGARight, TModel, InnerCell> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    FunctorDna5OrdValueComplement<int> fCompl;
    TValue sc = (me.gapFreqs[position(entry1)] / me.delRate) * (1.0-me.delErrorRate);
    for (unsigned i = 0; i < 4; ++i)
    {
        sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[fCompl(i)]) * me.delErrorRate * me.delErrorFreqs[i]; //  * me.scalingFactorDelErrors;
    }
	return  ((std::log10(sc/2.0) > -10)? std::log10(sc/2.0):-10);
}


template <typename TValue, typename TModel, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, BsTagList<BsProfileScoreRef, TModel, InnerCell> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
    // TODO free ends?
    TValue sc = me.gapFreqs[position(entry1)] / me.delRate;                                   // P(X-)/f(-)  // insRate?
	return  ((std::log10(sc/2.0) > -10)? std::log10(sc/2.0):-10);
}


////////////////////////////////////////////////////////////////////////////
// scoreGapOpendHorizontal (same as extending)
////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
	Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
    return scoreGapExtendHorizontal(me, entry1, entry2);
}


////////////////////////////////////////////////////////////////////////////
// scoreGapExtendVertical
////////////////////////////////////////////////////////////////////////////
// New gap in current profile
// After entry1

// Top strand, original (forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, BsTagList<BsProfileScoreCT, TModel, TCellDescriptor> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,  // Curr. column profile
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)      // Curr. read base
{

	TValue e =  pow(10, -(long double)value(entry2).count[1]/10.0);
	TValue pseudoMinVal = 0; /* 0.00000005/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                            value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                value(entry1).count[10] + value(entry1).count[11] +
                                1     // for ref
	                            );  */  // Use count of bases from position before, the more reads, the lower the prob to insert gap into profile
    TValue sc = (pseudoMinVal/me.readBaseFreqs[(unsigned)value(entry2).count[0]]) * (1.0-e);          // prob. to insert gap into profile (pseudo, assume we observed one with low quality)
    sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[(unsigned)value(entry2).count[0]];                               // 1 * e * seqError(- -> a)
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}

// Top strand, rev. compl. of original (projected from reverse to forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, BsTagList<BsProfileScoreCTRight, TModel, TCellDescriptor> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,  // Curr. column profile
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)      // Curr. read base
{
    FunctorDna5OrdValueComplement<int> fCompl;
	TValue e =  pow(10, -(long double)value(entry2).count[1]/10.0);
	TValue pseudoMinVal = 0; /* 0.00000005/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                            value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                value(entry1).count[10] + value(entry1).count[11] +
                                1     // for ref
	                            ); */   // Use count of bases from position before, the more reads, the lower the prob to insert gap into profile
    TValue sc = (pseudoMinVal/me.readBaseFreqs[(unsigned)value(entry2).count[0]]) * (1.0-e);          // prob. to insert gap into profile (pseudo, assume we observed one with low quality)
    sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[fCompl((int)value(entry2).count[0])];                               // 1 * e * seqError(- -> a)
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}

// Bottom strand, original (projected from reverse to forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, BsTagList<BsProfileScoreGA, TModel, TCellDescriptor> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,  // Curr. column profile
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)      // Curr. read base
{
    FunctorDna5OrdValueComplement<int> fCompl;
	TValue e =  pow(10, -(long double)value(entry2).count[1]/10.0);
	TValue pseudoMinVal = 0; /*0.00000005/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                            value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                value(entry1).count[10] + value(entry1).count[11] +
                                1     // for ref
	                            ); */ // Use count of bases from position before, the more reads, the lower the prob to insert gap into profile
	TValue sc = (pseudoMinVal/me.readBaseFreqs[fCompl((int)value(entry2).count[0])]) * (1.0-e);          // prob. to insert gap into profile (pseudo, assume we observed one with low quality); read base freq. regarding original strand
    sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[fCompl((int)value(entry2).count[0])];          // 1 * e * seqError(- -> a)
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}

// Bottom strand, rev. compl. of original (forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, BsTagList<BsProfileScoreGARight, TModel, TCellDescriptor> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,  // Curr. column profile
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)      // Curr. read base
{
    FunctorDna5OrdValueComplement<int> fCompl;
	TValue e =  pow(10, -(long double)value(entry2).count[1]/10.0);
	TValue pseudoMinVal = 0; /*0.00000005/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                            value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                value(entry1).count[10] + value(entry1).count[11] +
                                1     // for ref
	                            );*/  // Use count of bases from position before, the more reads, the lower the prob to insert gap into profile
	TValue sc = (pseudoMinVal/me.readBaseFreqs[fCompl((int)value(entry2).count[0])]) * (1.0-e);          // prob. to insert gap into profile (pseudo, assume we observed one with low quality); read base freq. regarding original strand
    sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[(unsigned)value(entry2).count[0]];          // 1 * e * seqError(- -> a)
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}


template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, BsTagList<BsProfileScoreRef, TModel, TCellDescriptor> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,      // Curr. column profile
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)      // Curr. read base
{
 	TValue pseudoMinVal = 0; /*0.00000005/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                            value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                value(entry1).count[10] + value(entry1).count[11] +
                                1     // for ref
	                            ); */
	TValue sc = pseudoMinVal/me.refBaseFreqs[(unsigned)value(entry2).count[0]];                       // prob. to insert gap into profile (pseudo, assume we observed one with low quality)
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}


////////////////////////////////////////////////////////////////////////////
// scoreGapOpenVertical
////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TBsProfileScore, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
	Score<TValue, BsTagList<TBsProfileScore, TModel, TCellDescriptor> > const & me,
    ConsensusScoreSequenceEntry<TSeq1> const & entry1,
    ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
   return scoreGapExtendVertical(me, entry1, entry2);
}

////////////////////////////////////////////////////////////////////////////
// Score
////////////////////////////////////////////////////////////////////////////

// Top strand, original (forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, BsTagList<BsProfileScoreCT, TModel, TCellDescriptor> > const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,    // Curr. column profile
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)    // Curr. read base
{
    TValue e =  pow(10, -value(entry2).count[1]/10.0);
    TValue sc = (me.targetFreqs[position(entry1)][(unsigned)value(entry2).count[0]] / me.readBaseFreqs[(unsigned)value(entry2).count[0]]) * (1.0-e);
    for (unsigned i = 0; i < 4; ++i)
    {
        if (i != (unsigned)value(entry2).count[0])
            sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[i]) * e * me.seqErrorFreqs[i*5 + (unsigned)value(entry2).count[0]];
    }
    // Check if only gaps (only base was removed by this read)
    // Similar pseudo count as for insertion of gap into profile/ insertion of read base
    if (me.targetFreqs[position(entry1)][0] + me.targetFreqs[position(entry1)][1] +
        me.targetFreqs[position(entry1)][2] + me.targetFreqs[position(entry1)][3] < 0.001)     {
        TValue pseudoMinVal = 0; /* 0.00001/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                                value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                    value(entry1).count[10] + value(entry1).count[11] +
                                    1) ; */                                     // Use count of gaps, the more reads, the lower the prob to insert
        sc = (pseudoMinVal/me.readBaseFreqs[(unsigned)value(entry2).count[0]]) * (1.0-e);           // prob. to insert gap into profile (pseudo, assume we observed one with low quality)
        sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[(unsigned)value(entry2).count[0]];                         // 1 * e * seqError(- -> a)
    }
    // TODO:   for pe reads: if right mate: we need to use seqErrorFreqs from complements, since read was projected on rev.compl. strand
    // other freqs are the same
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}

// Top strand, rev. compl. of original (projected from reverse to forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, BsTagList<BsProfileScoreCTRight, TModel, TCellDescriptor> > const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,    // Curr. column profile
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)    // Curr. read base
{
    FunctorDna5OrdValueComplement<int> fCompl;
    TValue e =  pow(10, -value(entry2).count[1]/10.0);
    TValue sc = (me.targetFreqs[position(entry1)][(unsigned)value(entry2).count[0]] / me.readBaseFreqs[(unsigned)value(entry2).count[0]]) * (1.0-e);
    for (unsigned i = 0; i < 4; ++i)
    {
        if (i != (unsigned)value(entry2).count[0])
            sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[i]) * e * me.seqErrorFreqs[fCompl(i)*5 + fCompl((int)value(entry2).count[0])];
    }
    // Check if only gaps (only base was removed by this read)
    // Similar pseudo count as for insertion of gap into profile/ insertion of read base
    if (me.targetFreqs[position(entry1)][0] + me.targetFreqs[position(entry1)][1] +
        me.targetFreqs[position(entry1)][2] + me.targetFreqs[position(entry1)][3] < 0.001)     {
        TValue pseudoMinVal = 0; /* 0.00001/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                                value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                    value(entry1).count[10] + value(entry1).count[11] +
                                    1) ;    */                                  // Use count of gaps, the more reads, the lower the prob to insert
        sc = (pseudoMinVal/me.readBaseFreqs[(unsigned)value(entry2).count[0]]) * (1.0-e);           // prob. to insert gap into profile (pseudo, assume we observed one with low quality)
        sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[fCompl((int)value(entry2).count[0])];                         // 1 * e * seqError(- -> a)
    }
    // TODO:   for pe reads: if right mate: we need to use seqErrorFreqs from complements, since read was projected on rev.compl. strand
    // other freqs are the same
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}


// Bottom strand, original (projected from reverse to forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, BsTagList<BsProfileScoreGA, TModel, TCellDescriptor> > const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,    // Curr. column profile
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)    // Curr. read base
{
    FunctorDna5OrdValueComplement<int> fCompl;
    TValue e =  pow(10, -(long double)value(entry2).count[1]/10.0);
    TValue sc = (me.targetFreqs[position(entry1)][(unsigned)value(entry2).count[0]] / me.readBaseFreqs[fCompl((int)value(entry2).count[0])]) * (1.0-e);  //
    for (unsigned i = 0; i < 4; ++i)
    {
        if (i != (unsigned)value(entry2).count[0])
            sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[fCompl(i)]) * e * me.seqErrorFreqs[fCompl(i)*5 + fCompl((int)value(entry2).count[0])];
            // We are dealing with reverse complements here: to get real rates: look at compl. of curr. base
    }
    // for pe reads: if right mate: we need to use seqErrorFreqs from complements

    // Check if only gaps (only base was removed by this read)
    // Similar pseudo count as for insertion of gap into profile/ insertion of read base
    if (me.targetFreqs[position(entry1)][0] + me.targetFreqs[position(entry1)][1] +
        me.targetFreqs[position(entry1)][2] + me.targetFreqs[position(entry1)][3] < 0.001)     {
        TValue pseudoMinVal = 0; /*0.00001/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                                  value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                      value(entry1).count[10] + value(entry1).count[11] +
                                      1) ;  */                                    // Use count of gaps, the more reads, the lower the prob to insert
        sc = (pseudoMinVal/me.readBaseFreqs[fCompl((int)value(entry2).count[0])]) * (1.0-e);           // prob. to insert gap into profile (pseudo, assume we observed one with low quality)
        sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[fCompl((int)value(entry2).count[0])];                         // 1 * e * seqError(- -> a)
    }
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}

// Bottom strand, rev. compl. of original (forward strand)
template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, BsTagList<BsProfileScoreGARight, TModel, TCellDescriptor> > const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,    // Curr. column profile
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)    // Curr. read base
{
    FunctorDna5OrdValueComplement<int> fCompl;
    TValue e =  pow(10, -(long double)value(entry2).count[1]/10.0);
    TValue sc = (me.targetFreqs[position(entry1)][(unsigned)value(entry2).count[0]] / me.readBaseFreqs[fCompl((int)value(entry2).count[0])]) * (1.0-e);  //
    for (unsigned i = 0; i < 4; ++i)
    {
        if (i != (unsigned)value(entry2).count[0])
            sc += (me.targetFreqs[position(entry1)][i] / me.readBaseFreqs[fCompl(i)]) * e * me.seqErrorFreqs[i*5 + (unsigned)value(entry2).count[0]];
            // We are dealing with reverse complements here: to get real rates: look at compl. of curr. base
    }
    // for pe reads: if right mate: we need to use seqErrorFreqs from complements

    // Check if only gaps (only base was removed by this read)
    // Similar pseudo count as for insertion of gap into profile/ insertion of read base
    if (me.targetFreqs[position(entry1)][0] + me.targetFreqs[position(entry1)][1] +
        me.targetFreqs[position(entry1)][2] + me.targetFreqs[position(entry1)][3] < 0.001)     {
        TValue pseudoMinVal = 0; /*0.00001/(value(entry1).count[0] + value(entry1).count[1] + value(entry1).count[2] + value(entry1).count[3] +
	                                  value(entry1).count[4] + value(entry1).count[5] + value(entry1).count[6] + value(entry1).count[7] +
                                      value(entry1).count[10] + value(entry1).count[11] +
                                      1) ;    */                                  // Use count of gaps, the more reads, the lower the prob to insert
        sc = (pseudoMinVal/me.readBaseFreqs[fCompl((int)value(entry2).count[0])]) * (1.0-e);           // prob. to insert gap into profile (pseudo, assume we observed one with low quality)
        sc += 1.0 * e * me.insErrorRate * me.insErrorFreqs[(unsigned)value(entry2).count[0]];                         // 1 * e * seqError(- -> a)
    }
	return  ((std::log10(sc) > -10)? std::log10(sc):-10);

}

template <typename TValue, typename TModel, typename TCellDescriptor, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, BsTagList<BsProfileScoreRef, TModel, TCellDescriptor> > const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,    // Curr. column profile
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)    // Curr. read base
{
    TValue sc = 0;
	if (value(entry2).count[0] == ordValue(Dna5('A')))
        sc = me.targetFreqs[position(entry1)][0];       // For the beginning: use profile value from non bs (regarding this ref. base) strand
	else if (value(entry2).count[0] == ordValue(Dna5('C')))
        sc = me.targetFreqs[position(entry1)][5];
	else if (value(entry2).count[0] == ordValue(Dna5('G')))
        sc = me.targetFreqs[position(entry1)][2];
	else if (value(entry2).count[0] == ordValue(Dna5('T')))
        sc = me.targetFreqs[position(entry1)][7];

	else if (value(entry2).count[0] == ordValue(Dna5('N')))
        sc = 1.0/4.0;   // ?

    sc = sc/me.refBaseFreqs[(unsigned)value(entry2).count[0]];

	return  ((std::log10(sc) > -10)? std::log10(sc):-10);
}


}

#endif
