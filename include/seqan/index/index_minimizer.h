#ifndef SEQAN_HEADER_INDEX_MINIMIZER_H
#define SEQAN_HEADER_INDEX_MINIMIZER_H

namespace seqan 
{
    //overloading _qgramFillsuffixArray 
	template <
		typename TSA,
		typename TText, 
		typename TValue, 
        unsigned TSPAN,
        unsigned TWEIGHT,
		typename TDir, 
		typename TBucketMap, 
		typename TWithConstraints, 
		typename TStepSize >
	inline void
	_qgramFillSuffixArray(
		TSA &sa, 
		TText const &text, 
		Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > shape, 
		TDir &dir, 
		TBucketMap &bucketMap, 
		TStepSize stepSize,
		TWithConstraints const)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type						TSize;

		if (empty(shape) || length(text) < length(shape)) return;

		TSize num_qgrams = length(text) - length(shape) + 1;
		TIterator itText = begin(text, Standard());
        String<TSize> delta;
        resize(delta, num_qgrams);

		if (TWithConstraints::VALUE) {
			TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;			// first hash
			if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = 0;						// if bucket is enabled
		} else
			sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = 0;			// first hash

		if (stepSize == 1)
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				if (TWithConstraints::VALUE) {
					TSize bktNo = getBucket(bucketMap, hashNext(shape, itText)) + 1;	// next hash
					if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;					// if bucket is enabled
				} else
					sa[dir[getBucket(bucketMap, hashNext(shape, itText)) + 1]++] = i;	// next hash
			}
		else
			for(TSize i = stepSize; i < num_qgrams; i += stepSize)
			{
				itText += stepSize;
				if (TWithConstraints::VALUE) {
					TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;		// next hash (we mustn't use hashNext here)
					if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;					// if bucket is enabled
				} else
					sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = i;		// next hash
			}

   
        compareHash<Shape<TValue, MinimizerShape<TSPAN,TWEIGHT> >, TIterator, HashNorm> compHash(begin(text));
        for (unsigned i = 0; i < length(dir) - 1; i++) 
            std::sort(begin(sa) + dir[i], begin(sa) + dir[i + 1], compHash);
	}
/*
	// multiple sequences
	template <
		typename TSA, 
		typename TString, 
		typename TSpec, 
		typename TShape, 
		typename TDir,
		typename TBucketMap,
		typename TStepSize,
		typename TWithConstraints >
	inline void
	_qgramFillSuffixArray(
		TSA &sa, 
		StringSet<TString, TSpec> const &stringSet,
		TShape shape, 
		TDir &dir,
		TBucketMap &bucketMap,
		TStepSize stepSize,
		TWithConstraints const)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type							TSize;

		if (empty(shape)) return;

		if (stepSize == 1)
			for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;
				TSize num_qgrams = length(sequence) - length(shape) + 1;

				typename Value<TSA>::Type localPos;
				assignValueI1(localPos, seqNo);
				assignValueI2(localPos, 0);

				TIterator itText = begin(sequence, Standard());
				if (TWithConstraints::VALUE) {
					TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;					// first hash
					if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;						// if bucket is enabled
				} else
					sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;			// first hash

				for(TSize i = 1; i < num_qgrams; ++i)
				{
					++itText;
					assignValueI2(localPos, i);
					if (TWithConstraints::VALUE) {
						TSize bktNo = getBucket(bucketMap, hashNext(shape, itText)) + 1;			// next hash
						if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;					// if bucket is enabled
					} else
						sa[dir[getBucket(bucketMap, hashNext(shape, itText)) + 1]++] = localPos;	// next hash
				}
			}
		else
			for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
			{
				TString const &sequence = value(stringSet, seqNo);
				if (length(sequence) < length(shape)) continue;
				TSize num_qgrams = length(sequence) - length(shape) + 1;

				typename Value<TSA>::Type localPos;
				assignValueI1(localPos, seqNo);
				assignValueI2(localPos, 0);

				TIterator itText = begin(sequence, Standard());
				for(TSize i = 0; i < num_qgrams; i += stepSize)
				{
					assignValueI2(localPos, i);
					if (TWithConstraints::VALUE) {
						TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;				// hash
						if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;					// if bucket is enabled
					} else
						sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;		// hash
					itText += stepSize;
				}
			}
	}
*/

//overloading getOccurrences
    template < 
        typename TText, 
        unsigned TSPAN, 
        unsigned TWEIGHT, 
        typename TSpec, 
        unsigned TSPAN2, 
        unsigned TWEIGHT2, 
        typename TValue >
    inline typename Infix< typename Fibre< Index< TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, TSpec> >, FibreSA>::Type const >::Type
    getOccurrences(
        Index< TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, TSpec> > const &index, 
        Shape< TValue, MinimizerShape<TSPAN2, TWEIGHT2> > const &shape)
    {
        typedef typename Size<typename Fibre< Index< TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, TSpec> >, FibreDir>::Type>::Type TDirSize;
        typedef typename Iterator<TText const, Standard>::Type  TIter; 
        TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
        compareHash<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> >, TIter, HashNormLess> compHashL(begin(indexText(index)));
        compareHash<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> >, TIter, HashNormGreat> compHashG(begin(indexText(index)));
        unsigned _lower_bound = std::lower_bound(begin(indexSA(index)) + indexDir(index)[bucket], begin(indexSA(index)) + indexDir(index)[bucket + 1], shape.nhValue, compHashL) - begin(indexSA(index));
        unsigned _upper_bound = std::upper_bound(begin(indexSA(index)) + indexDir(index)[bucket], begin(indexSA(index)) + indexDir(index)[bucket + 1], shape.nhValue, compHashG) - begin(indexSA(index));
        return infix(indexSA(index), begin(indexSA(index)) + _lower_bound, begin(indexSA(index)) + _upper_bound);
    } 
}

#endif //#ifndef SEQAN HEADER_...:
