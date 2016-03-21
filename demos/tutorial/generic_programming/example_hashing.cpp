#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

//![hashAll]
template <typename TShape, typename TString> 
void hashAll(TShape & shape, TString & str)
{ 
    typedef typename Iterator<TString>::Type TIterator;
    TIterator it = begin(str);
    TIterator it_end = end(str) - span(shape);
    while (it != it_end) 
    {
        unsigned int hash_value = hash(shape, it);
        /* do some things with the hash value */ ++it;
//![hashAll]
        ignoreUnusedVariableWarning(hash_value);
//![hashAll]
    }
}
//![hashAll]

//![classShape]
struct SimpleShape_;
typedef Tag<SimpleShape_> SimpleShape;

template <typename TValue, typename TSpec = SimpleShape> 
class Shape;
//![classShape]

//![classSimpleShape]
template <typename TValue> 
class Shape< TValue, SimpleShape >
{
    public:
        unsigned int span;
};
//![classSimpleShape]

//![classUngappedShape]
template <unsigned int q = 0> 
struct UngappedShape;

template <typename TValue, unsigned int q> 
class Shape< TValue, UngappedShape<q> >
{
  public:
     static unsigned int const span = q;
};
//![classUngappedShape]

//![hashNext]
template <typename TValue, unsigned int q, typename TIterator>
inline unsigned int
hashNext(Shape< TValue, UngappedShape<q> > const & shape, TIterator it, unsigned int prev)
{
    unsigned int val = prev * ValueSize<TValue>::VALUE - *it * shape.fac
                            + *(it + shape.span);
    return val;
// shape.fac stores |Σ|^q
}
//![hashNext]

//![specializedHashAll]
template <typename TValue, unsigned int q, typename TString>
void specializedHashAll(Shape< TValue, UngappedShape<q> > & shape, TString & str)
{
    typedef typename Iterator<TString>::Type TIterator;
    TIterator it = begin(str); 
    TIterator it_end = end(str) - span(shape);
    unsigned int hash_value = hash(shape, it);
    /* do some things with the hash value */

//![specializedHashAll]
    ignoreUnusedVariableWarning(hash_value);
//![specializedHashAll]
    while (++it != it_end)
    {
        unsigned int hash_value = hashNext(shape, it, hash_value);
        /* do some things with the hash value */
    }
}
//![specializedHashAll]

int main()
{
    return 0;
}
