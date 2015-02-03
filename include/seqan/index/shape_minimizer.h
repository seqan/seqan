#ifndef SEQAN_HEADER_SHAPE_MINIMIZER_H
#define SEQAN_HEADER_SHAPE_MINIMIZER_H

namespace seqan
{

    struct HashNorm_;
    struct HashMini_;
    struct HashNormLess_;
    struct HashNormGreat_; 

    typedef Tag<HashNorm_> const  HashNorm;
    typedef Tag<HashMini_> const    HashMini;     
    typedef Tag<HashNormLess_> const HashNormLess;
    typedef Tag<HashNormGreat_> const HashNormGreat;

    template <unsigned TSPAN, unsigned TWEIGHT>
    struct MinimizerShape{};

	template <typename TValue, unsigned TSPAN, unsigned TWEIGHT>
	class Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> >
	{
	public:
        typedef typename Value<Shape>::Type THashValue;
        
        unsigned span; 
        unsigned weight; 
        THashValue hValue;     //current minimizer
        THashValue nhValue;     //normal hash value;
        THashValue leftChar;   //leftmost character

        Shape():
            span(TSPAN),
            weight(TWEIGHT),
            leftChar(-1){} 
	};
    
    //compare hash values used as a function object in sorting function like std::sort, std::lower_bound etc.

    template <typename TShape, typename TIter, typename TSpec = Default>
    struct compareHash{
        TShape shape;
        TIter it;        
        bool operator()(unsigned i, unsigned j){
            return (hash(shape, it + i) < hash(shape, it + j));
        }
        compareHash(TShape _Shape, TIter _It ):
            it(_It){} 
    };
   
    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter>   
    struct compareHash< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> >, TIter, HashNorm>{
        Shape<TValue, UngappedShape<TSPAN> > shape; 
        TIter it;
        bool operator()(unsigned i, unsigned j){
            return (hash(shape, it + i) < hash(shape, it + j));
        }
        compareHash(TIter _It):
            it(_It){}
    };

    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter>   
    struct compareHash< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> >, TIter, HashNormLess>{
        typedef typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type THashValue;
        Shape<TValue, UngappedShape<TSPAN> > shape; 
        TIter it;
        bool operator()(unsigned i, THashValue _hashValue){
        return (hash(shape, it + i) < _hashValue);
        }
        compareHash(TIter _It):
            it(_It){}
    };

    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter>   
    struct compareHash< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> >, TIter, HashNormGreat>{
        typedef typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type THashValue;
        Shape<TValue, UngappedShape<TSPAN> > shape; 
        TIter it;
        bool operator()(THashValue _hashValue, unsigned i){
        return (_hashValue < hash(shape, it + i)) ;
        }
        compareHash(TIter _It):
            it(_It){}
    };

    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT>
    struct LENGTH< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >
    {
		enum { VALUE = TSPAN };
	};

    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT>
	struct WEIGHT< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >
	{
		enum { VALUE = TWEIGHT };
	};

    //return  lexicographically smaller of S as the minimizer
    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter> 
    inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type 
    getMinimizer(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &me, const TIter it)
    {
        typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type THValue;
      
        SEQAN_ASSERT_GT((unsigned)me.span, 0u);
        SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 

        Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
        TIter leftIt = it;
        THValue miniTmp = hash(tmpShape, leftIt);
        for (unsigned i = 1; i < me.span - me.weight + 1; i++)
        {
            if(miniTmp > hashNext(tmpShape, leftIt + i))
                miniTmp = tmpShape.hValue;        
        }
        me.hValue = miniTmp;
        return me.hValue;

    }
    
    //return  lexicographically smaller of S and the reverse complement of S as the minimizer 
    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter> 
    inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type 
    getMinimizer_2(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &me, const TIter it)
    {
        typedef typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type THValue;
 
        SEQAN_ASSERT_GT((unsigned)me.span, 0u);
        SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 
        TIter leftIt = it;
        Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
        hashInit(tmpShape, leftIt);
        THValue miniTmp = hash(tmpShape, leftIt);

        for (int i = 1 ; i < tmpShape.span; i++)
        {
            if(miniTmp > hashNext(tmpShape, leftIt + i))
                miniTmp = tmpShape.hValue;        
        }

        String<TValue> rc;
        resize(rc, leftIt -it + me.span - 1 );
        TIter rcIt = begin(rc);
        arrayCopy(it, leftIt + length(me), rcIt);
        reverseComplement(rc);
        leftIt = begin(rc);

        THValue miniTmpRC = hash(tmpShape, leftIt);
        for (int i = 1 ; i < tmpShape.span; i++)
        {
            if(miniTmpRC > hashNext(tmpShape, leftIt + i))
                miniTmpRC = me.hValue;
        }
        me.hValue = (miniTmp < miniTmpRC ? miniTmp : miniTmpRC);
        return me.hValue;
    }

    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT>
    inline SEQAN_HOST_DEVICE
    typename Size< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type
    weight(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > const &me)
    {
    SEQAN_CHECKPOINT
        return me.weight;
    }

	template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter>
    inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type
	hash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &me, TIter it)
	{
        Shape<TValue, UngappedShape<TSPAN> > _shape;
        me.nhValue = hash(_shape, it);
        return getMinimizer(me, it);
	}

    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TIter>
    inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > >::Type
    hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &me, TIter const &it)
    {
        Shape<TValue, UngappedShape<TSPAN> > _shape;
        me.nhValue = hash(_shape, it); 
        return getMinimizer(me, it);
    }


}	// namespace seqan

#endif
