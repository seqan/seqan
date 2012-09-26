// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_BANDED_ALIGN_H
#define SEQAN_HEADER_BANDED_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
//needleman wunsch alignment
template <typename TScoreValue, unsigned DIMENSION, typename TString, typename TValue, typename TSpecSeed>
TScoreValue
_bandedNeedlemanWunsch(Matrix<TScoreValue, DIMENSION> & matrix_,
						 Seed<TValue, TSpecSeed> const &seed,
						 TValue k,
						 TString const & str1_,
						 TString const & str2_,
						 Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT
	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Standard>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables
	TSize height = leftDiagonal(seed) - rightDiagonal(seed)+1;
	TSize str1_length = height+2*k;
	TSize str2_length = length(str2_);

	TStringIterator x_begin = begin(str1_);
	TStringIterator x_end = end(str1_) -1;
	TStringIterator y_begin = begin(str2_);
	--y_begin;
	TStringIterator y_end = end(str2_)-1;
	
	TValue up_height = leftDiagonal(seed) - startDiagonal(seed) + k; //equals length of empty lower triangle
	//TValue up_width = str1_length - up_height;
	TValue down_height = endDiagonal(seed) - rightDiagonal(seed) + k; //equals length of empty lower triangle
	TValue down_width = str1_length - down_height;
	TSize length_right_diag = str2_length - down_height;
	TSize length_left_diag = str2_length - up_height;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap = scoreGapExtend(score_);

	TScoreValue horizontalValue = 0;
	TScoreValue border_ = score_gap;
	TScoreValue verticalValue = border_;

	setDimension(matrix_, 2);
	setLength(matrix_, 0, str1_length + 2);
	setLength(matrix_, 1, str2_length + 1);
	resize(matrix_);
	
	TMatrixIterator col_ = begin(matrix_);
	TMatrixIterator finger1;
	TMatrixIterator finger2;
	TMatrixIterator finger3;

	TStringIterator x = x_begin;
	TStringIterator y = y_begin;

	//-------------------------------------------------------------------------
	// init

	TValue inf = -100000;
	finger2 = col_;

	setPosition(finger2, length(matrix_, 0)-1);
	for (unsigned int i = 1; i != length_right_diag; ++i){
		*finger2 = inf;
		goNext(finger2, 1);
	}

	border_ = down_height * score_gap;

	*finger2 = inf;
	for (int i = -1; i != down_height; ++i){
		goPrevious(finger2, 0);
		goNext(finger2,1);
		*finger2 = border_;
		border_ -= score_gap;
	}
		
	border_ = score_gap;
	finger3 = finger2;

	for (int i = 0; i != down_width; ++i){
		goPrevious(finger2, 0);
		*finger2 = border_;
		border_ += score_gap;
	}

	*finger2 = inf;
	TScoreValue tmp;
	//-------------------------------------------------------------------------
	//first section
	
	y=y_end;
	TSize run_length = down_width;
	TSize measure = 0;

	for (int i = -1; i != down_height; ++i){	
		verticalValue = *finger3;
		finger1 = finger3;
		goPrevious(finger1,0);
		goPrevious(finger3,1);
		finger2 = finger3;
		goNext(finger3,0);
		horizontalValue = *finger3;
		x = x_end;
		for (unsigned int j = 0; j != run_length; ++j){
			if (*x==*y){
				*finger2 = verticalValue + score_match;
			} else {
				tmp = *finger1;
				TScoreValue s1 = verticalValue + score_mismatch;
				TScoreValue s2 = score_gap + ((horizontalValue > tmp) ? horizontalValue : tmp);
				*finger2 = (s1 > s2) ? s1 : s2;
			}
			horizontalValue = *finger2;
			goPrevious(finger2);
			verticalValue = *finger1;
			goPrevious(finger1);
			--x;
		}
		*finger2 = inf;
		++measure;
		if (measure < length_left_diag)
			++run_length;
		--y;
	}
	--run_length;

	goPrevious(finger3);

	while (y!= y_begin)
	{	
		verticalValue = *finger3;
		finger1 = finger3;
		goPrevious(finger1,0);
		goPrevious(finger3,1);
		finger2 = finger3;
		x = --x_end;
		horizontalValue = inf;
		for (unsigned int j = 0; j != run_length; ++j){
			if (*x==*y){
				*finger2 = verticalValue + score_match;
			} else {
				tmp = *finger1;
				TScoreValue s1 = verticalValue + score_mismatch;
				TScoreValue s2 = score_gap + ((horizontalValue > tmp) ? horizontalValue : tmp);
				*finger2 = (s1 > s2) ? s1 : s2;
			}
			horizontalValue = *finger2;
			goPrevious(finger2);
			verticalValue = *finger1;
			goPrevious(finger1);
			--x;
		}
		*finger2 = inf;
		++measure;
		if (measure >= length_left_diag)
			--run_length;
		--y;
	}
	return *(++finger2);
}

//////////////////////////////////////////////////////////////////////////////
//traceback through needleman wunsch matrix


//Position berechnen!
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION>
void
_bandedNeedlemanWunschTrace(Align<TTargetSource, TTargetSpec> & target_,
						Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_,
						Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT

	typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
	TTargetIterator target_0 = iter(row(target_, 0), 0);
	TTargetIterator target_1 = iter(row(target_, 1), 0);

	typedef typename Iterator<TTargetSourceSegment, Standard>::Type TStringIterator;
	TStringIterator it_0 = begin(str_0);
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = begin(str_1);
	TStringIterator it_1_end = end(str_1);

	TScoreValue score_diff = scoreMismatch(score_) - scoreGapExtend(score_);

	//-------------------------------------------------------------------------
	//follow the trace until the border is reached
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		bool gv;
		bool gh;
		if (*it_0 == *it_1)
		{
			gv = gh = true;
		}
		else
		{
			TMatrixIterator it_ = source_;

			goNext(it_, 0);
			TScoreValue h = *it_;

			it_ = source_;
			goNext(it_, 1);
			TScoreValue d = *it_;
		
			goPrevious(it_, 0);
			TScoreValue v = *it_;
			
			gv = (v >= h) | (d + score_diff >= h);
			gh = (h >  v) | (d + score_diff >= v);
		}
		if (gv){
			++it_1;
			goNext(source_, 1);
			goPrevious(source_, 0);
		}else{
			insertGap(target_1);
		}

		if (gh) {
			++it_0;
			goNext(source_, 0);
		}else{
			insertGap(target_0);
		}
		++target_0;
		++target_1;
	}
}



///////////////////////////////////////////////////////////////////////////////////////
//Gotoh
//Banded alignment with affine gap costs
template <typename TScoreValue, unsigned DIMENSION, typename TString, typename TValue, typename TSpecSeed>
TScoreValue
_bandedGotoh(Matrix<TScoreValue, DIMENSION> & diag_matrix_,
	   Matrix<TScoreValue, DIMENSION> & vert_matrix_,
	   Matrix<TScoreValue, DIMENSION> & hori_matrix_,
	   Seed<TValue, TSpecSeed> const &seed,
	   TValue k,
	   TString const & str1_,
	   TString const & str2_,
	   Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT


	typedef Matrix<TScoreValue, DIMENSION> TMatrix;

	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Iterator<TMatrix, Rooted>::Type TMatrixIterator;

	typedef typename Iterator<TString const, Rooted>::Type TStringIterator;
	typedef typename Value<TString const>::Type TAlphabet;

	//-------------------------------------------------------------------------
	//define some variables

	TSize height = leftDiagonal(seed) - rightDiagonal(seed)+1;
	TSize str1_length = height+2*k;
	TSize str2_length = length(str2_);

	TValue up_height = leftDiagonal(seed) - startDiagonal(seed) + k; //equals length of empty lower triangle
	//TValue up_width = str1_length - up_height;
	TValue down_height = endDiagonal(seed) - rightDiagonal(seed) + k; //equals length of empty lower triangle
	TValue down_width = str1_length - down_height;
	TValue length_right_diag = str2_length - down_height;
	TValue length_left_diag = str2_length - up_height;


	// TStringIterator x_begin = begin(str1_) - 1;
	TStringIterator x_end = end(str1_) - 1;
	TStringIterator y_begin = begin(str2_) - 1;
	TStringIterator y_end = end(str2_) - 1;

	TStringIterator x = x_end;
	TStringIterator y;

	TScoreValue score_match = scoreMatch(score_);
	TScoreValue score_mismatch = scoreMismatch(score_);
	TScoreValue score_gap_open = scoreGapOpen(score_);
	TScoreValue score_gap_extend = scoreGapExtend(score_);

	//TScoreValue v;

	setDimension(diag_matrix_, 2);
	setLength(diag_matrix_, 0, str1_length + 2);
	setLength(diag_matrix_, 1, str2_length + 1);
	resize(diag_matrix_);
	setDimension(vert_matrix_, 2);
	setLength(vert_matrix_, 0, str1_length + 2);
	setLength(vert_matrix_, 1, str2_length + 1);
	resize(vert_matrix_);
	setDimension(hori_matrix_, 2);
	setLength(hori_matrix_, 0, str1_length + 2);
	setLength(hori_matrix_, 1, str2_length + 1);
	resize(hori_matrix_);

	TSize pos = length(diag_matrix_, 0)-1;
	// TMatrixIterator diag_col_ = begin(diag_matrix_);
	TMatrixIterator diag_finger1(diag_matrix_,pos);
	TMatrixIterator diag_finger2;
	TMatrixIterator diag_finger3;

	// TMatrixIterator vert_col_ = begin(vert_matrix_);
	TMatrixIterator vert_finger1(vert_matrix_,pos);
	TMatrixIterator vert_finger2;
	TMatrixIterator vert_finger3;
	// TMatrixIterator hori_col_ = begin(hori_matrix_);
	TMatrixIterator hori_finger1(hori_matrix_,pos);
	TMatrixIterator hori_finger2;
	TMatrixIterator hori_finger3;

	
	//-------------------------------------------------------------------------
	// init
	
	TScoreValue inf = -1000000;
	for (int i = 1; i != length_right_diag; ++i){
		*diag_finger1 = inf;
		*hori_finger1 = inf;
		*vert_finger1 = inf;
		goNext(diag_finger1, 1);
		goNext(hori_finger1, 1);
		goNext(vert_finger1, 1);
	}

	*diag_finger1 = inf;
	*hori_finger1 = inf;
	*vert_finger1 = inf;

	TScoreValue border_ = score_gap_open + (down_height-1) * score_gap_extend;

	for (int i = -1; i != down_height; ++i){
		goPrevious(diag_finger1, 0);
		goNext(diag_finger1,1);
		goPrevious(vert_finger1, 0);
		goNext(vert_finger1,1);
		goPrevious(hori_finger1, 0);
		goNext(hori_finger1,1);
		*diag_finger1 = border_;
		*vert_finger1 = border_;
		*hori_finger1 = inf;
		border_ -= score_gap_extend;
	}

	*diag_finger1 = 0;
	*vert_finger1 = inf;
	*hori_finger1 = inf;
	
	diag_finger3 = diag_finger1;
	hori_finger3 = hori_finger1;
	vert_finger3 = vert_finger1;

	for (int i = 0; i < down_width; ++i){
		goPrevious(diag_finger1, 0);
		goPrevious(hori_finger1, 0);
		goPrevious(vert_finger1, 0);
		border_ += score_gap_extend;
		*diag_finger1 = border_;
		*hori_finger1 = border_;
		*vert_finger1 = inf;
	}

	*diag_finger1 = inf;
	*hori_finger1 = inf;
	
	//-------------------------------------------------------------------------
	//first section
	TScoreValue hori_value, vert_value, diag_value, tmp_diag, tmp_vert, diag_front, diag_match, diag_back;
	y=y_end;
	TSize run_length = down_width;
	TSize measure = 0;

	for (int i = -1; i != down_height; ++i){	
		goPrevious(hori_finger3,1);
		hori_finger1 = hori_finger3;
		goNext(hori_finger3);
		hori_value = *hori_finger3;

		vert_finger2 = vert_finger3;
		goPrevious(vert_finger2);
		goPrevious(vert_finger3,1);
		vert_finger1 = vert_finger3;
		goNext(vert_finger3);

		diag_finger2 = diag_finger3;
		diag_match=*diag_finger2;
		goPrevious(diag_finger2);
		goPrevious(diag_finger3,1);
		diag_finger1 = diag_finger3;
		goNext(diag_finger3);
		diag_back=*diag_finger3;

		x = x_end;
		for (unsigned int j = 0; j != run_length; ++j){
			diag_front= *diag_finger2;
			hori_value = (hori_value+score_gap_extend > diag_back+score_gap_open) ? hori_value+score_gap_extend : diag_back+score_gap_open;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+score_gap_extend > diag_front+score_gap_open) ? tmp_vert+score_gap_extend : diag_front+score_gap_open;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			
			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? score_match : score_mismatch);
			if (diag_value < tmp_diag)
				diag_value = tmp_diag;
			diag_back = diag_value;
			*diag_finger1 = diag_value;
			goPrevious(diag_finger1);
			diag_match = *diag_finger2;
			goPrevious(diag_finger2);
			--x;
		}
		*diag_finger1 = inf;
		*vert_finger1 = inf;
		++measure;
		if (static_cast<TValue>(measure) < length_left_diag)
			++run_length;
		--y;
	}
	--run_length;

	while(y!= y_begin)//for (int i = 0; i != main; ++i){
	{
		goPrevious(diag_finger3);
		goPrevious(hori_finger3);
		goPrevious(vert_finger3);

		goPrevious(hori_finger3,1);
		hori_finger1 = hori_finger3;
		goNext(hori_finger3);
		hori_value = inf;

		vert_finger2 = vert_finger3;
		goPrevious(vert_finger2);
		goPrevious(vert_finger3,1);
		vert_finger1 = vert_finger3;
		goNext(vert_finger3);

		diag_finger2 = diag_finger3;
		diag_match=*diag_finger2;
		goPrevious(diag_finger2);
		goPrevious(diag_finger3,1);
		diag_finger1 = diag_finger3;
		goNext(diag_finger3);
		diag_back= inf;
		x = --x_end;
		for (unsigned int j = 0; j != run_length; ++j){
			diag_front= *diag_finger2;
			hori_value = (hori_value+score_gap_extend > diag_back+score_gap_open) ? hori_value+score_gap_extend : diag_back+score_gap_open;
			*hori_finger1 = hori_value;
			goPrevious(hori_finger1);

			tmp_vert = *vert_finger2;
			vert_value = (tmp_vert+score_gap_extend > diag_front+score_gap_open) ? tmp_vert+score_gap_extend : diag_front+score_gap_open;
			*vert_finger1 = vert_value;
			goPrevious(vert_finger1);
			goPrevious(vert_finger2);

			tmp_diag = (vert_value > hori_value) ? vert_value  : hori_value;
			diag_value = diag_match + ((*x == *y) ? score_match : score_mismatch);
			if (diag_value < tmp_diag)
				diag_value = tmp_diag;
			diag_back = diag_value;
			*diag_finger1 = diag_value;
			goPrevious(diag_finger1);
			diag_match = *diag_finger2;
			goPrevious(diag_finger2);
			--x;
		}
		*vert_finger1 = inf;
		*diag_finger1 = inf;
		++measure;
		if (static_cast<TValue>(measure) >= length_left_diag)
			--run_length;
		--y;
	}
	--run_length;
	++diag_finger1;

	return *diag_finger1;
}


//////////////////////////////////////////////////////////////////////////////
//gotoh trace
template <typename TTargetSource, typename TTargetSpec, typename TScoreValue, unsigned DIMENSION, typename TValue>
void
_bandedGotohTrace(Align<TTargetSource, TTargetSpec> & target_,
					Matrix<TScoreValue, DIMENSION> & diag_matrix_,
					Matrix<TScoreValue, DIMENSION> & vert_matrix_,
					Matrix<TScoreValue, DIMENSION> & hori_matrix_,
					TValue position_)
{
SEQAN_CHECKPOINT
	typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

	typedef Align<TTargetSource, TTargetSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Standard>::Type TTargetIterator;
    typedef typename Infix<TTargetSource>::Type TTargetSourceSegment;
	typedef typename Iterator<TTargetSource, Standard>::Type TStringIterator;

	TTargetSourceSegment str_0 = sourceSegment(row(target_, 0));
	TTargetSourceSegment str_1 = sourceSegment(row(target_, 1));
	typename Size<TTargetSourceSegment>::Type dim_0_len = length(hori_matrix_,0);

	TPosition pos = position_;

	TTargetIterator target_0 = iter(row(target_, 0), 0, Standard());
	TTargetIterator target_1 = iter(row(target_, 1), 0, Standard());

	TStringIterator it_0 = iter(str_0, 0, Standard());
	TStringIterator it_0_end = end(str_0);

	TStringIterator it_1 = iter(str_1, 0, Standard());
	TStringIterator it_1_end = end(str_1);
	
	//-------------------------------------------------------------------------
	//follow the trace until the border is reached
	while ((it_0 != it_0_end) && (it_1 != it_1_end))
	{
		if (getValue(diag_matrix_,pos) > getValue(hori_matrix_,pos))
		{
			if (getValue(diag_matrix_,pos) > getValue(vert_matrix_,pos))
			{
				++it_0;
				++it_1;
				pos += dim_0_len;
			} else{
				++it_1;
				insertGap(target_0);
				pos += dim_0_len-1;
			}
		}
		else
		{
			if (getValue(hori_matrix_,pos) > getValue(vert_matrix_,pos)) 
			{					
				++it_0;
				insertGap(target_1);
				++pos;
			} 
			else
			{
				++it_1;
				insertGap(target_0);
				pos += dim_0_len-1;
			}
		}
		++target_0;
		++target_1;
	}
}

//////////////////////////////////////////////////////////////////////////////
//needleman wunsch alignment


template <typename TSource, typename TSpec, typename TScoreValue, typename TValue, typename TSpecSeed>
TScoreValue
bandedAlignment(Align<TSource, TSpec> & align_,
				Seed<TValue, TSpecSeed> const & seed,
				TValue k,
				Score<TScoreValue, Simple> const & score_,
				NeedlemanWunsch)
{
SEQAN_CHECKPOINT

	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

	TScoreValue ret;
	Matrix<TScoreValue> matr;
	ret = _bandedNeedlemanWunsch(matr, seed, k, sourceSegment(row(align_, 0)), sourceSegment(row(align_, 1)), score_);
	Iter< Matrix<TScoreValue>, PositionIterator > iter_ = begin(matr);
	setPosition(iter_, 1+k + leftDiagonal(seed) - startDiagonal(seed));
	_bandedNeedlemanWunschTrace(align_, iter_, score_);
	return ret;
}

template <typename TSource, typename TSpec, typename TScoreValue, typename TValue, typename TSpecSeed>
TScoreValue
bandedAlignment(Align<TSource, TSpec> & align_,
				Seed<TValue, TSpecSeed> const & seed,
				TValue k,
				Score<TScoreValue, Simple> const & score_,
				Gotoh)
{
SEQAN_CHECKPOINT
	clearGaps(row(align_,0));
	clearGaps(row(align_,1));

	TScoreValue ret;

	Matrix<TScoreValue> d_matr;
	Matrix<TScoreValue> v_matr;
	Matrix<TScoreValue> h_matr;
	ret = _bandedGotoh(d_matr, v_matr, h_matr, seed, k, sourceSegment(row(align_, 0)), sourceSegment(row(align_, 1)), score_);
	_bandedGotohTrace(align_, d_matr, v_matr, h_matr, 1+k + leftDiagonal(seed) - startDiagonal(seed));	
	return ret;
}


/**
.Function.bandedAlignment:
..class:Class.Seed
..summary:Calculates a banded alignment around a Seed.
..cat:Seed Handling
..signature:bandedAlignmnet(align, seed, k, score)
..param.align:An alignment dataStructure containing two sequences (sections). The sequence (parts) must be equal to the part which is covered by the seed.
...type:Class.Align
..param.seed: The seed for wich the banded alignmnent shall be constructed.
...type:Class.Seed
..param.k: A value describing the additional extension of the diagonal band.
..param.score: The score matrix used.
...type:Spec.Simple Score
..returns:The score of the optimal banded alignment given in align.
..remarks:Use the function @Function.globalAlignment@ with the tag $BandedNeedlemanWunsch$ or $BandedGotoh$ for more general banded alignment without a seed.
..see:Function.globalAlignment
..see:Tag.Pairwise Global Alignment Algorithms
..include:seqan/seeds.h
*/
template <typename TSource, typename TSpec, typename TScoreValue, typename TValue, typename TSpecSeed>
TScoreValue
bandedAlignment(Align<TSource, TSpec> & align_,
				Seed<TValue, TSpecSeed> const & seed,
				TValue k,
				Score<TScoreValue, Simple> const & score_)
{
SEQAN_CHECKPOINT

	if(scoreGapOpen(score_)==scoreGapExtend(score_))
	{//linear gap costs
		return bandedAlignment(align_, seed, k, score_, NeedlemanWunsch());
	}
	else
	{//affine gap costs
		return bandedAlignment(align_, seed, k, score_, Gotoh());
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

