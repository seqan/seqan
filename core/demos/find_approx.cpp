///A tutorial about the use of approximate find algorithms.
#include <iostream>
#include <seqan/find.h>

using namespace seqan;

///Example 1: This program finds all occurrences of $CCT$ in $AACTTAACCTAA$ with $\leq 1$ error using the @Shortcut.MyersUkkonen@ approximate search algorithm.
int main() 
{
	String<char> haystk("AACTTAACCTAA");
	String<char> ndl("CCT");

	Finder<String<char> > fnd(haystk);
	Pattern<String<char>, MyersUkkonen> pat(ndl);
///The function @Function.setScoreLimit@ sets the limit score an occurrence must reach.
///Since the used scoring scheme is a distance measure (edit distance), all scores are negative.
///A score limit of $\geq -1$ therefore means an edit distance $\leq 1$.
///Note that @Function.position@ returns the position of the last found occurrence.
	setScoreLimit(pat, -1);
	while (find(fnd, pat)) {
		std::cout << position(fnd) << ": " << getScore(pat) << "\n";
	}

///Example 2: Finding all start and endpositions
	String<char> t = "babybanana";
	String<char> p = "babana";
	Finder<String<char> > finder(t);
	Pattern<String<char>, Myers<FindInfix> > pattern(p);
///Instead of using @Function.setScoreLimit@, we pass the score limit $-2$ as a third argument to find
	while (find(finder, pattern, -2)) {
		std::cout << "end: " << endPosition(finder) << std::endl;
///In order to find the begin position, we have to call @Function.findBegin@.
///Note that the third argument of @Function.findBegin@ is optional.
///The default is the score limit that was used during the last call of @Function.find@ (i.e. -2 in this example).
		while (findBegin(finder, pattern, getScore(pattern))) {
			std::cout << "begin: " << beginPosition(finder) << std::endl;
			std::cout << infix(finder) << " matches with score ";
			std::cout << getBeginScore(pattern) << std::endl;
		}
	}
	return 0;
}

