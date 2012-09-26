#define SEQAN_PROFILE

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "STNode.h"
#include "Tree.h"
#include <seqan.h>
#include "JournalString.h"
#include "JournalIterator.h"
#include "Utility.h"

using namespace std;
using namespace seqan;

class Base{
public:

	Base( int f ) : m( f ) { };
	virtual ~Base() {
		cout << "Base Destroy! " << m << endl;
	};

	int m;
};

class Derived : public Base {
public:
	Derived( int f ) : Base( f ) { };
	virtual ~Derived() {
		cout << "Derived Destroy!" << endl;
	}
};



int main(){

	typedef String< char > TSequence;
	typedef String< char , Journal<Alloc<>, Basic > > TJournal;
	typedef String< TJournal > TJournalString;

	TSequence test1 = "Test1";
	TSequence test2 = "Test2";

//	String< char > * ptr_wombat = &wombat;
	TJournal journal1( test1 );
	TJournal jounral2( test2 );

	TJournalString testString;
	appendValue(testString, jounral2);
	appendValue(testString, journal1);


	Iterator< TJournal >::Type it = begin(journal1);
	Iterator< TJournal >::Type irt_end = end(journal1);

	while (it != irt_end) {
		cout << *it;
		++it;
	}
	cout << endl;

	insert(journal1, 4, "Test");
	erase(journal1, 1,4);

	Iterator< TJournal >::Type it_2 = begin(journal1);
	Iterator< TJournal >::Type irt_end_2 = end(journal1);

	cout << "After that nothing happens" << endl;

	while (it_2 != irt_end_2) {
		cout << *it_2;
		++it_2;
	}
	cout << endl;

//	size_t t = length(testString);
//	cout << "Length of testString: " << t << endl;

//	size_t t = length(grombat);
//	cout << t << endl;
//	String < String <char, Journal< Alloc<>, Basic > > > string_grombat;
//	resize(string_grombat, 2);
//	string_grombat[0] = grombat;
//	string_grombat[1] = grombat_cpy;
//	cout << ( string_grombat[1]).underlying() << endl;

//	underlying ( grombat);
//	SeqTree< STNode< JournalNode> > tree = (journal_tree (grombat));
//	TreeIterator< SeqTree< STNode< JournalNode> > > root_it = root( tree);
//
//	TreeIterator< SeqTree< STNode< JournalNode> > > it_tree = begin( tree );
//	TreeIterator< SeqTree< STNode< JournalNode> > > end_it = end( tree );

//	insertion_string (grombat);
//	Iterator< String< char> >::Type it = begin(wombat);
//	Iterator< String< char> >::Type it_end = end(wombat);
//
//	while(it != it_end) {
//		cout << *it;
//		++it;
//	}
//	cout << endl;
//	String< char> test ;
//	String<char> ins1 = "logopdie";
//	unsigned int pos = 3;
//	insert( grombat, pos, ins1);
//	JournalIterator< char, Alloc<>, Basic > itBeg_grom = begin( grombat);
//	JournalIterator< char, Alloc<>, Basic > itEnd_grom = end ( grombat );
//
//	while (itBeg_grom != itEnd_grom) {
//		cout << *itBeg_grom;
//		++itBeg_grom;
//	}
//	cout << underlying(grombat) << endl;

//	Iterator< String < char, Journal <Alloc<>, Basic > > > itEnd_grom = end(grombat);

//
//	String< char > wombat2 = "Wombat2";
//	String< char, Journal< Alloc<>, Basic > > bla( wombat2 );
//
//	String<char> test = bla.underlying();
//
//	cout << test << endl;
//
//	Test t(10);
//
//	cout << length(t) << endl;


//	std::cout << grombat << std::endl;
//	erase( grombat, 1, 4 );
//	std::cout << "erasing 1..3:" << std::endl;
//	std::cout << grombat << std::endl;
//	std::cout << "flattening . . ." << std::endl;
//	flatten( grombat );
//	std::cout << "journal state:" << std::endl;
//	std::cout << grombat << std::endl;
//	std::cout << "underlying string:" << std::endl;
//	std::cout << wombat << std::endl;
//
//	seqan::util::generate_random_string( 36, wombat, "abcdexx" );
//	String< char, Journal< Alloc<>, Basic > > schlombat( wombat );
//
//	std::cout << schlombat << std::endl;
//	std::cout << "replacing 'x' with '_':" << std::endl;
//	String< char > ins = "_";
//
//	for( unsigned int q = 0; q < length( schlombat ); ++ q ){
//		if( value( schlombat, q ) == 'x' ){
//			erase( schlombat, q, q+1 );
//			insert( schlombat, q, ins );
//		}
//	}
//
//	std::cout << schlombat << std::endl;
//	std::cout << "flattening. . ." << std::endl;
//	flatten( schlombat );
//	std::cout << "journal state:" << std::endl;
//	std::cout << schlombat << std::endl;
//	std::cout << "underlying string:" << std::endl;
//	std::cout << wombat << std::endl;
//
	return 23;
//	// BEGIN: Tests that are not run anymore at the moment
//	cout << "Length\tDeletions\tJournal\tAlloc" << endl;
//	SEQAN_PROTIMESTART(_alloc);
//	SEQAN_PROTIMESTART(_journal);
//
//	String< char > str_alloc;
//	String< unsigned int > positions;
//	for( unsigned int i = 100000; i <= 1600000; i *= 2 ){
//
//		for( unsigned int k = 1; k <= 20; ++k ){
//			for( unsigned int l = 0; l <= 10; ++l ){
//				seqan::util::generate_random_string( i, str_alloc, "acgt" );
//				String< char, Journal< Alloc<>, Basic > > str_journal( str_alloc );
//
//				for( unsigned int k = 0; k < 1000; ++k ){
//					appendValue( positions, seqan::util::randomNumber( seqan::length( str_alloc ) - 1001 ), Generous() );
//				}
//
//				cout << i << "\t" << 20*l << "\t";
//
//				for( unsigned int j = 0; j < 20*l; ++j ){
//					erase( str_journal, positions[j], positions[j] + 1 );
//				}
//
//				char tmp;
//				SEQAN_PROTIMEUPDATE(_journal);
//				for( unsigned int k = 0; k < 100; ++k ){
//					for( unsigned int j = 0; j < length( positions ); ++j ){
//						tmp = value( str_journal, positions[j] );
//					}
//				}
//				cout << SEQAN_PROTIMEUPDATE(_journal) << "\t";
//
//				SEQAN_PROTIMEUPDATE(_alloc);
//				for( unsigned int k = 0; k < 100; ++k ){
//					for( unsigned int j = 0; j < length( positions ); ++j ){
//						tmp = value( str_alloc, positions[j] );
//					}
//				}
//				cout << SEQAN_PROTIMEUPDATE(_alloc) << "\n";
//
//				clear( positions );
//				clear( str_alloc );
//			}
//		}
//	}
//
//	cout << endl;
//	cout << " Fin " << endl;
}

