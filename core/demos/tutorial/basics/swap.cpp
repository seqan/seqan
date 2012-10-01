// FRAGMENT(swap-headers)
#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>


using namespace seqan;
using namespace std;

// FRAGMENT(swap-declaration)
template <typename T> void swap(T& container, int i, int j, int k)
{

// FRAGMENT(swap-metafunction)
	// define helper variable
	T help;
	resize(help,k);
	
	for (int x=0; x<k; ++x) 
		value(help,x) = container[i+x];

// FRAGMENT(swap-work)	
	for (int x=0; x<k; ++x) 
		value(container,i+x) = value(container,j+x);
	for (int x=0; x<k; ++x) 
		value(container,j+x) = help[x];

	return;
}


// FRAGMENT(swap-main)
int main()
{
	typedef String<Dna> TDnaString;
	TDnaString dna = "AAAATTTT";
	
	typedef String<int> TIntString;
	typedef Iterator<String<int>, Rooted >::Type TIntIterator;
	
	TIntString numbers;
    appendValue(numbers,1);   appendValue(numbers,1);   appendValue(numbers,1);
	appendValue(numbers,1);   appendValue(numbers,1);   appendValue(numbers,1);
	appendValue(numbers,3);	  appendValue(numbers,3);	appendValue(numbers,3);
	appendValue(numbers,3);   appendValue(numbers,3);	appendValue(numbers,3);
	
// FRAGMENT(swap-apply)
	swap(dna,1,4,2);
	cout << dna << endl;
	
    swap(numbers,1,7,2);
	for (TIntIterator it=begin(numbers); !atEnd(it); goNext(it)) {
		std::cout << *it;
	}
	cout << endl;
	
	return 0;
}

