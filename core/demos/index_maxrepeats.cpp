///A tutorial about finding maximal repeats.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
///We begin with a @Class.String@ to store our sequence.
	String<char> myString = "How many wood would a woodchuck chuck.";

///Then we create an @Class.Index@ of this @Class.StringSet@.
	typedef Index< String<char> > TMyIndex;
	TMyIndex myIndex(myString);

///To find maximal repeats, we use SeqAn's @Spec.MaxRepeats Iterator@
///and set the minimum repeat length to 3.
	typedef Iterator< TMyIndex, MaxRepeats >::Type TMaxRepeatIterator;
	TMaxRepeatIterator myRepeatIterator(myIndex, 3);

	while (!atEnd(myRepeatIterator)) 
	{
///A repeat can be represented by its length and positions it occurs at.
///$myRepeatIterator$ iterates over all repeat strings.
///Please note that in contrast to supermaximal repeats, given a maximal repeat string,
///not all pairs of its occurrences are maximal repeats.
///So we need an iterator to iterate over all maximal pairs of this repeat string.
///The @Spec.MaxRepeats Iterator@ can be seen as a container and be iterated for itself.
		Iterator<TMaxRepeatIterator>::Type myRepeatPair(myRepeatIterator);
		while (!atEnd(myRepeatPair)) {
			std::cout << *myRepeatPair << ", ";
			++myRepeatPair;
		}

///@Function.repLength@ returns the length of the repeat string.
		std::cout << repLength(myRepeatIterator) << "   ";

///The repeat string itself can be determined with @Function.representative@
		std::cout << "\t\"" << representative(myRepeatIterator) << '\"' << std::endl;

		++myRepeatIterator;
	}

	return 0;
}
