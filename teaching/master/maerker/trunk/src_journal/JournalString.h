/*
 * JournalString.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#include <seqan.h>
#include "Tree.h"
#include "STNode.h"

#ifndef JOURNALSTRING_H_
#define JOURNALSTRING_H_

namespace seqan {


	template< typename TStringSpec, typename TJournalSpec >
	struct Journal;

	template <typename TValue, typename TStringSpec, typename TJournalSpec>
	struct Value< String< TValue, Journal< TStringSpec, TJournalSpec > > > {
		typedef String<TValue> Type;
	};

	struct Operation{
		Operation() : length( 0 ) { };
		Operation( size_t len ) : length( len ) { };
		Operation( Operation const & other ) : length( other.length ) { };
		virtual ~Operation() { };
		virtual bool isDel(){
			return false;
		}

		virtual bool isIns(){
			return false;
		}

		size_t length;
	};

	struct Insertion : public Operation{
		Insertion( size_t len ) : Operation( len ) {};
		virtual ~Insertion() { };
		virtual bool isIns(){
			return true;
		}
	};

	struct Deletion : public Operation{
		Deletion( size_t len ) : Operation( len ) {};
		virtual ~Deletion() { };
		virtual bool isDel(){
			return true;
		}

	};

	inline bool isDeletion( Operation* op ){
		return op->isDel();
	}

	inline bool isInsertion( Operation* op ){
		return op->isIns();
	}
//	TODO: SCHMARRN
//	inline bool isIns( Operation * op ){
//		return false;
//	}
//
//	inline bool isIns( Insertion * op ){
//		return true;
//	}

	struct JournalNode{

		JournalNode() : physical_position( 0 ), virtual_position( 0 ), blocksize( 0 ), op() { };
		JournalNode( size_t pos, Operation* op ) : virtual_position( pos ), op( op ) { };

		size_t physical_position;
		size_t virtual_position;
		size_t blocksize;
		Operation* op;

	};

	/**
	 * Prints the content of the JournalNode
	 */
	template<>
	inline std::string stringify( JournalNode const & cargo ){
		std::stringstream sstr;
		sstr << "<" << ( isInsertion( cargo.op ) ? "Insertion " : ( isDeletion( cargo.op ) ? "Deletion " : "Original " ) ) << "size=" << (*cargo.op).length << " vpos=" << cargo.virtual_position << " ppos=" << cargo.physical_position << " blocksize=" << cargo.blocksize << " />";
		return sstr.str();
	}

	/**
	 * Defines the String-template specialization for journal strings.
	 *
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	class String< TValue, Journal< TStringSpec, TJournalSpec > > {

	public:
		String() {

		};
		String( String< TValue, TStringSpec > & other ) : m_holder( other ) {
			JournalNode root( 0, new Operation( seqan::length( other ) ) );
			root.physical_position = 0;
			root.blocksize = seqan::length( other );
			seqan::insert( m_journal_tree, root );
			len = seqan::length( other );
			JournalNode dummy( len, new Insertion( 1 ) );
			dummy.physical_position = 0;
			dummy.blocksize = 1;
			seqan::insert( m_journal_tree, dummy );
			seqan::appendValue( m_insertion_string, TValue(), Generous() );
		};
		String( String< TValue, TStringSpec > const & other ) : m_holder( other ) {
			JournalNode root( 0, new Operation( seqan::length( other ) ) );
			root.physical_position = 0;
			root.blocksize = seqan::length( other );
			seqan::insert( m_journal_tree, root );
			len = seqan::length( other );
			JournalNode dummy( len, new Insertion( 1 ) );
			dummy.physical_position = 0;
			dummy.blocksize = 1;
			seqan::insert( m_journal_tree, dummy );
			seqan::appendValue( m_insertion_string, TValue(), Generous() );
		};

		String< TValue, TStringSpec> & underlying() const {
			return value( m_holder );
		}

		seqan::SeqTree< seqan::STNode< JournalNode > > & journal_tree() {
			return m_journal_tree;
		}

		seqan::String< TValue > & insertion_string() {
			return m_insertion_string;
		}

//		size_t length() const {
//			return len;
//		}

		void setLength( size_t l ){
			len = l;
		}

		template< typename TPosition >
		inline TValue operator[]( TPosition const & pos ){
			return value( *this, pos );
		}

		seqan::SeqTree< seqan::STNode< JournalNode > > m_journal_tree;
		typename DefaultStringSpec< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type m_insertion_string;
		Holder< String< TValue, TStringSpec > > m_holder;

		size_t len;
	};

	/**
	 * Returns the journal tree of the given journal
	 */
//	template< typename TValue, typename TStringSpec, typename TJournalSpec >
//	seqan::SeqTree< seqan::STNode< JournalNode > > &
//	journal_tree( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal ){
//		return journal.m_journal_tree;
//	}

	/**
	 * Returns the journal tree of the given journal
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	seqan::SeqTree< seqan::STNode< JournalNode > > const &
	journal_tree( String< TValue, Journal< TStringSpec, TJournalSpec > > const & journal ){
		return journal.m_journal_tree;
	}

	/**
	 * Returns the journal tree of the given journal
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	seqan::SeqTree< seqan::STNode< JournalNode > > &
	journal_tree( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal ){
		return journal.m_journal_tree;
	}

	/**
	 * Returns the inserted string of the given journal
	 */
//	template< typename TValue, typename TStringSpec, typename TJournalSpec >
//	typename DefaultStringSpec< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type &
//	insertion_string( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal ){
//		return journal.m_insertion_string;
//	}

	/**
	 * Returns the inserted string of the given journal.
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	typename DefaultStringSpec< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type &
	insertion_string( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal ) {
		return journal.m_insertion_string;
	}

	/**
	 * Returns the underlying sequence of the given journal.
	 */
//	template< typename TValue, typename TStringSpec, typename TJournalSpec >
//	String< TValue, TStringSpec> &
//	underlying( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal ){
//		return value( journal.m_holder );
//	}

	/**
	 * Returns the underlying sequence of the given journal.
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	String< TValue, TStringSpec> &
	underlying( String< TValue, Journal< TStringSpec, TJournalSpec > > const & journal ){
		return value( journal.m_holder );
	}

	/**
	 * Overloads the operator <<, which returns the output stream containing the
	 * value of the Iterator over a journal
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	std::ostream & operator << ( std::ostream& os, String< TValue, Journal< TStringSpec, TJournalSpec > > const & journal ) {
		typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type it_journal = begin( journal );
		typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type it_journal_end = end( journal );
		while( it_journal != it_journal_end ){
			os << *it_journal;
			++it_journal;
		}
		return os;
	}

	/**
	 * Sets the given length to the given journal
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	inline void _setLength( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal, size_t len ){
		journal.setLength( len );
	}

	/**
	 * Defines the default string type of the insertion string
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	struct 	DefaultStringSpec< String< TValue, Journal< TStringSpec, TJournalSpec > > >{
		typedef String< TValue > Type;
	};

	/**
	 *
	 */
	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TPosition >
	inline TValue value( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal, TPosition pos ){
		JournalNode entry = cargo( (*seqan::find( journal.journal_tree(), pos )) );
		if( !isInsertion( entry.op ) ){
			return journal.underlying()[ entry.physical_position + pos - entry.virtual_position ];
		}else{
			return journal.insertion_string()[ entry.physical_position + pos - entry.virtual_position ];
		}
	}

	template< typename TPosition >
	bool isLess( JournalNode const & jn, TPosition const & pos ){
		return jn.virtual_position + jn.blocksize <= pos;
	}

	template< typename TPosition >
	bool isGreater( JournalNode const & jn, TPosition const & pos ){
		return jn.virtual_position > pos;
	}

	template< typename TPosition >
	bool isLess( TPosition const & pos, JournalNode const & jn ){
		return pos < jn.virtual_position;
	}

	template< typename TPosition >
	bool isGreater( TPosition const & pos, JournalNode const & jn ){
		return  pos >= jn.virtual_position + jn.blocksize;
	}

	bool isLess( JournalNode const & lhs, JournalNode const & rhs ){
		return lhs.virtual_position + lhs.blocksize <= rhs.virtual_position;
	}

	bool isGreater( JournalNode const & lhs, JournalNode const & rhs ){
		return lhs.virtual_position >= rhs.virtual_position + rhs.blocksize;
	}

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	size_t length( String< TValue, Journal< TStringSpec, TJournalSpec > > const & journal ){
		return journal.len;
	}

	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TPosition, typename TString >
	void insert( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal, TPosition pos, TString str ){

		// New node will represent the insertion in the journal tree
		JournalNode new_node;
		new_node.virtual_position = pos;
		new_node.op = new Insertion( seqan::length( str ) );
		new_node.physical_position = seqan::length( journal.insertion_string() ); //block of data will start at current end of insertion string
		new_node.blocksize = seqan::length( str );

		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > node_iter = seqan::find( journal.journal_tree(), pos );
		JournalNode old_node = cargo(*node_iter); //store old node that contains insertion position

		JournalNode left_split( old_node.virtual_position, old_node.op ); //part of the old node that lies to the right of the insertion
		left_split.physical_position = old_node.physical_position;
		left_split.blocksize = pos - old_node.virtual_position;

		JournalNode right_split( pos + seqan::length( str ), old_node.op ); // part of the old node that lies to the left of the insertion
		right_split.physical_position = old_node.physical_position + pos - old_node.virtual_position;
		right_split.blocksize = old_node.virtual_position + old_node.blocksize - pos;

		append( journal.insertion_string(), str, Generous() ); //append inserted values to insertion string

		seqan::erase( journal.journal_tree(), node_iter ); // delete old node from journal tree

		/* Shift everything right of the insertion accordingly: */
		node_iter = seqan::find( journal.journal_tree(), old_node.virtual_position + old_node.blocksize ); //find successor of old node in journal tree
		do{
			cargo(*node_iter).virtual_position += seqan::length( str ); //shift nodes by insertion length
		}while( seqan::goNext( node_iter ) );
		node_iter = seqan::find( journal.journal_tree(), old_node.virtual_position + old_node.blocksize );
		assert( seqan::insert( journal.journal_tree(), new_node ) && "new_node insertion failed!");
		if( right_split.blocksize != 0 )
			assert( seqan::insert( journal.journal_tree(), right_split ) && "right_split insertion failed!");
		if( left_split.blocksize != 0 )
			assert( seqan::insert( journal.journal_tree(), left_split ) && "left_split isnertion failed!");

		_setLength( journal, seqan::length( journal ) + seqan::length( str ) ); //update length
	}

	//TODO must be implemented
//	void eraseFromTo( SeqTree< STNode< JournalNode > > & tree, size_t pos_begin, size_t pos_end ){
//
//	}

	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TPosition >
	void erase( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal, TPosition pos, TPosition pos_end ){


		// Initial Assertions to check for correct parameters
		assert( pos < pos_end && "Deletion of negative Length!" );
		assert( pos_end <= length(journal) && "Deletion exceeds string limits!" );

		// The new Node will represent the Deletion in the Tree
		JournalNode new_node;
		new_node.virtual_position = pos;
		new_node.op = new Deletion( pos_end - pos ); //Deletion of length pos_end - pos

		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > node_iter = seqan::find( journal.journal_tree(), pos ); //find the node that contains the starting position of the deletion
		JournalNode old_node = cargo( (*node_iter) ); //store for further use

		//left_split represents the part of the original node that preceedes the deletion
		JournalNode left_split( old_node.virtual_position, new Operation( *old_node.op ) );
		left_split.physical_position = old_node.physical_position;
		left_split.blocksize = pos - old_node.virtual_position;
		(*left_split.op).length = left_split.blocksize;	//adjusting the operations length to the new blocksize

		JournalNode right_split; //Node to hold the part of the node containing the deletion-end-position that is left after deleting from 'pos' to 'pos_end'

		// The String positions will hold the virtual positions of all nodes n, that need to be removed.
		// For those nodes T[pos..pos_end] \ T[n.virtual_position..n.virtual_position+n.blocksize] is not empty,
		// which means that their block of data contains values that will be deleted
		seqan::String< size_t > positions;
		while( cargo( (*node_iter) ).virtual_position < pos_end ){ //find all nodes that need to be removed
			appendValue( positions, cargo( (*node_iter) ).virtual_position );
			right_split = cargo( (*node_iter) ); //the last node that still has a virtual_position smaller than pos_end will contains the right end of the deletion
			if( !seqan::goNext( node_iter )){
				std::cout << std::endl << seqan::info( *node_iter ) << std::endl;
				assert( false && "Attempting to delete trailing dummy-node . . . this should not have happened!" );
			}
		};

		node_iter = seqan::find( journal.journal_tree(), right_split.virtual_position ); //get iterator to node containing pos_end

		//Adapt right_split to correctly represent the succeeding part of the node ocntaining pos_end
		right_split.physical_position += pos_end - right_split.virtual_position;
		right_split.blocksize -= pos_end - right_split.virtual_position;
		right_split.virtual_position = pos;
		(*right_split.op).length = right_split.blocksize;

		/* Shift everything right of the deletion accordingly: */
		seqan::goNext( node_iter );
		do{
			cargo( (*node_iter) ).virtual_position -= pos_end - pos; //virtual positions of nodes succeeding the deletion are diminshed by deletio nlength
		}while( seqan::goNext( node_iter ) );

		for( unsigned int k = 0; k < seqan::length( positions ); ++k ){
			erase( journal.journal_tree(), positions[k] ); //erase all nodes that are spanned by the deletion as determined above
		}

		seqan::insert( journal.journal_tree(), new_node ); //insert the new_node
		if( right_split.blocksize != 0 ){ //is the right split not of length zero? (i.e. the deletion did not end exactly at the border between two blocks)
			if( !seqan::insert( journal.journal_tree(), right_split )){
				std::cout << "Right-Split:" << seqan::stringify( right_split ) << std::endl;
				assert( false && "Insertion failed!");
			};
		}
		if( left_split.blocksize != 0 ){ //is the left split not of length zero (i.e. the deletion did not begin exactly at the border between two blocks)
			if( !seqan::insert( journal.journal_tree(), left_split )) {
				std::cout << "Left-Split:" << seqan::stringify( left_split ) << std::endl;
				assert( false && "Insertion failed!");
			};
		}

		_setLength( journal, length( journal ) - pos_end + pos );

	}

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	void flatten( String< TValue, Journal< TStringSpec, TJournalSpec > > & journal ){
		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > it_tree = end( journal.journal_tree() ); // get tree-iterator to last node
		resize( journal.underlying(), length( journal.underlying() ) + length( journal.insertion_string() ), Generous() ); // resize underlying string to make sufficient space for flattening

		typename Iterator< String< TValue, TStringSpec > >::Type it_underlying_write = end( journal.underlying() ); // get write-iterator to underlying string
		typename Iterator< String< TValue, TStringSpec > >::Type it_underlying_read = begin( journal.underlying() ); // get read-iterator to underlying string
		typename Iterator< typename DefaultStringSpec< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type >::Type it_insertion_read = begin( journal.insertion_string() ); // get read-iterator to insertion string

		--it_underlying_write; //we want to start at the last position that can hold values

		while( goPrev( it_tree ) ){ // while there are still nodes to process

			it_underlying_write -= cargo( *it_tree ).blocksize; //position write iterator

			if( isInsertion( cargo( *it_tree ).op ) ){
				it_insertion_read = begin( journal.insertion_string() ) + cargo( *it_tree ).physical_position; //position read-iterator
				for( unsigned int i = 0; i < cargo( *it_tree ).blocksize; ++i ){
					*it_underlying_write = *it_insertion_read; //copy values
					++it_underlying_write;
					++it_insertion_read;
				}
			}else if( !isDeletion( cargo( *it_tree ).op ) ){
				it_underlying_read = begin( journal.underlying() ) + cargo( *it_tree ).physical_position; //position read-iterator
				for( unsigned int i = 0; i < cargo( *it_tree ).blocksize; ++i ){
					*it_underlying_write = *it_underlying_read; //copy values
					++it_underlying_write;
					++it_underlying_read;
				}
			}
			it_underlying_write -= cargo( *it_tree ).blocksize; //reposition write iterator to point to the position left of the block we just wrote
		};

		SeqTree< STNode< JournalNode > > journal_tree; //new empty journal string
		JournalNode root( 0, new Operation( length( journal.underlying() ) ) ); //make rot node represent whole underlying string (including garbage at the beginning of the underlying string)
		int pos_end = length( journal.underlying() ) - length( journal ) - 1; //position left of the first value of the flattened string
		root.physical_position = 0;
		root.blocksize = length( journal.underlying() );
		insert( journal_tree, root ); //insert root
		resize( journal.insertion_string(), 0, Generous() ); //empty insertion string
		JournalNode dummy( length( journal.underlying() ), new Insertion( 1 ) ); //trailing dummy node
		dummy.physical_position = 0;
		dummy.blocksize = 1;
		insert( journal_tree, dummy ); //insert dummy
		appendValue( journal.insertion_string(), TValue(), Generous() ); //append zero value for dummy value to point to
		journal.journal_tree() = journal_tree; //assign new journal tree to journal
		_setLength( journal, length( journal.underlying() ) ); //assign correct length
		erase( journal, 0, pos_end ); //erase the first 0..pos_end values that are left over from the flattening process but do not contain values in the flattened string
	}

//	template< typename TValue, typename TStringSpec, typename TJournalSpec >
//	void
//	synchronize(
//		String< TValue, Journal< TStringSpec, TJournalSpec > > & journal,
//		Index< String< TValue, Journal< TStringSpec, TJournalSpec > >, Index_ESA > & index,
//		String< SAValue< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type > & inverse_suffix_array
//	){
//		assert( indexSupplied( index, ESA_SA ) && "Suffix-Array not present in index, use indexCreate instead of synchronize!" );
//		assert( indexSupplied( index, ESA_LCP ) && "lcp-Table not present in index, use indexCreate instead of synchronize!" );
//
//		String< SAValue< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type > suffix_array = indexSA( index );
//		String< size_t > lcp_table = indexLCP( index );
//
//		assert( length( lcp_table ) == length( suffix_array ) && length( suffix_array ) == length( inverse_suffix_array ) && "length mismatch between index tables!" );
//
//		String< SAValue< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type > suffixes;
//		String< size_t > deletions; //will hold the indexes (within the SA) of all suffixes that have to be removed from the suffix array
//
//		TreeIterator< SeqTree< STNode< JournalNode > > > it_tree = begin( journal.journal_tree() ); //get first node
//
//		SAValue< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type block_start = 0;
//
//		do{ //loop over all nodes in the tree
//			JournalNode jnode = cargo( *it_tree );
//			if( isInsertion( jnode.op ) ){ // node represents Insertion
//				for( unsigned int i = 0; i < jnode.blocksize; ++i ){
//					appendValue( suffixes, jnode.virtual_position + i, Generous() ); //indexes of suffixes that have to be inserted into the sa
//				}
//			}else if( isDeletion( jnode.op ) ){ // node represents Deletion
//				for( unsigned int i = 0; i < jnode.blocksize; ++i ){
//					appendValue( deletions, jnode.virtual_position + i, Generous() ); //indexes of suffixes that have to be deleted from the sa
//				}
//			}
//
//			if( jnode.virtual_position == 0 )
//				continue; //first node cannot influence any suffixes
//
//			unsigned int offset = 1; //use this to go backwards up to the beginning of the string
//			SAValue< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type suffix_index;
//
//			do{
//				suffix_index = jnode.virtual_position - offset; //current suffix
//				if( _max( // is the suffix influenced? <=> is the operation within the sensitive interval?
//							lcp_table[ inverse_suffix_array[ suffix_index ] ],
//							( inverse_suffix_array[ suffix_index ] > 0 ) ?
//							lcp_table[ inverse_suffix_array[ suffix_index ] - 1 ] : 0
//						) >= offset
//				){
//					appendValue( deletions, suffix_index, Generous() ); //mark suffix for deletion
//					appendValue( suffixes, suffix_index, Generous() ); //mark suffix for re-insertion
//				}
//			}while( suffix_index > block_start );
//
//			block_start = jnode.vitual_position;
//		}while( goNext( it_tree ) );
//
//		// start by erasing from the suffix-array, lcp-table and inverse suffix-array all suffixes in the 'deletions'-list
//
//	}
};


#endif /* JOURNALSTRING_H_ */
