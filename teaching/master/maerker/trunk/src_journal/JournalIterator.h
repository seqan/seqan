/*
 * JournalIterator.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#include "JournalString.h"

#ifndef JOURNALITERATOR_H_
#define JOURNALITERATOR_H_

namespace seqan {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//JournalIterator class

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	struct JournalIterator{

		typedef typename seqan::Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type Type;

		JournalIterator( seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & str )
			: m_string( str ),
			m_iter_underlying( seqan::begin( underlying( str ) ) ),
			m_iter_insertion( seqan::begin( insertion_string( str ) ) ),
			m_iter_tree( seqan::begin( journal_tree( str ) ) )
		{
			while( isDeletion( cargo(*m_iter_tree).op ) && seqan::goNext( m_iter_tree ) ) { };
			m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
			m_iter_underlying += cargo(*m_iter_tree).physical_position;
			m_iter_insertion += cargo(*m_iter_tree).physical_position;
		}

 		JournalIterator(
 				seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & str,
 				seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > it_tree
 			) : m_string( str ),
				m_iter_underlying( seqan::begin( str.underlying() ) + cargo(*it_tree).physical_position ),
				m_iter_insertion( seqan::begin( str.insertion_string() ) + cargo(*it_tree).physical_position ),
				m_iter_tree( it_tree )
		{
			m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
		}

		inline void operator= ( Type & other ){
			m_string = other.string();
			m_iter_underlying = other.iter_underlying();
			m_iter_insertion = other.iter_insertion();
			m_iter_tree = other.iter_tree();
			m_remaining_blocksize = other.remaining_blocksize();
		}

		operator void* ()
		{
			if( !isInsertion( cargo(*m_iter_tree).op ) ){
				return m_iter_underlying;
			}else{
				return m_iter_insertion;
			}
		}

		inline bool operator== ( Type & other ){
			return ( m_iter_underlying == other.iter_underlying() && m_iter_insertion == other.iter_insertion() && m_remaining_blocksize == other.remaining_blocksize() );
			//m_string == other.string() && && m_iter_tree == other.iter_tree() && m_remaining_blocksize == other.remaining_blocksize() );
		}

		inline bool operator!= ( Type & other ){
			return !( operator==(other) );
		}

		inline TValue operator*(){
			if( !isInsertion( cargo(*m_iter_tree).op ) ){
				return *m_iter_underlying;
			}else{
				return *m_iter_insertion;
			}
		}

		inline Type & operator++(){
			if( m_remaining_blocksize > 0 ){
				++m_iter_insertion;
				++m_iter_underlying;
				--m_remaining_blocksize;
			}else{
				while( seqan::goNext( m_iter_tree ) && isDeletion( cargo(*m_iter_tree).op ) ) { }; //TODO: specify behavior for end of string reached -> is reset to tree root regardless of root-node-type atm.
				m_iter_insertion = seqan::begin( m_string.insertion_string() ) + cargo(*m_iter_tree).physical_position;
				m_iter_underlying = seqan::begin( m_string.underlying() ) + cargo(*m_iter_tree).physical_position;
				m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
			}
			return *this;
		}

		inline Type operator++( int ){
			Type temp_it = *this;
			++(*this);
			return temp_it;
		}

		inline Type & operator--(){
			minus_one();
			return *this;
		}

		inline Type operator--( int ){
			Type temp_it = *this;
			minus_one();
			return temp_it;
		}

		template <typename TSize>
		inline Type & operator+=(TSize other) {
			for (; other > 0; --other)
				++(*this);
			return *this;
		}

		template <typename TSize>
		inline Type const operator+( TSize const & other ) {
			Type result = *this;
			result += other;
			return result;
		}

		template <typename TSize>
		inline Type & operator-=(TSize other)
		{
			for (; other > 0; --other) {
				--(*this);
			}
		}
		template <typename TSize>
		inline Type const operator-(TSize const & other) {
			Type result = *this;
			result -= other;
			return result;
		}

		int remaining_blocksize() {
			return m_remaining_blocksize;
		}

		typename seqan::Iterator< String< TValue, TStringSpec > >::Type iter_underlying(){
			return m_iter_underlying;
		}

		typename seqan::Iterator< String< TValue > >::Type iter_insertion(){
			return m_iter_insertion;
		}

		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > iter_tree(){
			return m_iter_tree;
		}

		seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & string(){
			return m_string;
		}

	private:

		inline void minus_one(){
			if( m_remaining_blocksize < (int)cargo(*m_iter_tree).blocksize - 1 ){
				--m_iter_insertion;
				--m_iter_underlying;
				++m_remaining_blocksize;
			}else{
				while( seqan::goPrev( m_iter_tree ) && isDeletion( cargo(*m_iter_tree).op ) ) { }; //TODO: specify behavior for end of string reached -> is reset to tree root regardless of root-node-type atm.
				m_iter_insertion = seqan::begin( m_string.insertion_string() ) + cargo(*m_iter_tree).physical_position + cargo(*m_iter_tree).blocksize - 1;
				m_iter_underlying = seqan::begin( m_string.underlying() ) + cargo(*m_iter_tree).physical_position + cargo(*m_iter_tree).blocksize - 1;
				m_remaining_blocksize = 0;
			}
		}

		seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & m_string;
		typename seqan::Iterator< String< TValue, TStringSpec > >::Type m_iter_underlying;
		typename seqan::Iterator< String< TValue > >::Type m_iter_insertion;
		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > m_iter_tree;
		int m_remaining_blocksize;
	};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ConstJournalIterator class

#if 0
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	struct ConstJournalIterator{

		typedef typename seqan::Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > const >::Type Type;

		ConstJournalIterator( seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > const & str )
			: m_string( str ),
			m_iter_underlying( seqan::begin( underlying( str ) ) ),
			m_iter_insertion( seqan::begin( insertion_string( str ) ) ),
			m_iter_tree( seqan::begin( journal_tree( str ) ) )
		{
			while( isDeletion( cargo(*m_iter_tree).op ) && seqan::goNext( m_iter_tree ) ) { };
			m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
			m_iter_underlying += cargo(*m_iter_tree).physical_position;
			m_iter_insertion += cargo(*m_iter_tree).physical_position;
		}

		ConstJournalIterator( seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > const & str, seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > it_tree )
			: m_string( str ),
			m_iter_underlying( seqan::begin( str.underlying() ) + cargo(*it_tree).physical_position ),
			m_iter_insertion( seqan::begin( str.insertion_string() ) + cargo(*it_tree).physical_position ),
			m_iter_tree( it_tree )
		{
			m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
		}

		inline void operator= ( Type other ){
			m_string = other.string();
			m_iter_underlying = other.iter_underlying();
			m_iter_insertion = other.iter_insertion();
			m_iter_tree = other.iter_tree();
			m_remaining_blocksize = other.remaining_blocksize();
		}

		inline bool operator== ( Type const & other ){
			return ( m_iter_underlying == other.iter_underlying() && m_iter_insertion == other.iter_insertion() && m_remaining_blocksize == other.remaining_blocksize() );
			//m_string == other.string() && && m_iter_tree == other.iter_tree() && m_remaining_blocksize == other.remaining_blocksize() );
		}

		inline bool operator!= ( Type const & other ){
			return !( operator==(other) );
		}

		inline TValue operator*(){
			if( !isInsertion( cargo(*m_iter_tree).op ) ){
				return *m_iter_underlying;
			}else{
				return *m_iter_insertion;
			}
		}

		inline Type & operator++(){
			if( m_remaining_blocksize > 0 ){
				++m_iter_insertion;
				++m_iter_underlying;
				--m_remaining_blocksize;
			}else{
				while( seqan::goNext( m_iter_tree ) && isDeletion( cargo(*m_iter_tree).op ) ) { }; //TODO: specify behavior for end of string reached -> is reset to tree root regardless of root-node-type atm.
				m_iter_insertion = seqan::begin( m_string.insertion_string() ) + cargo(*m_iter_tree).physical_position;
				m_iter_underlying = seqan::begin( m_string.underlying() ) + cargo(*m_iter_tree).physical_position;
				m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
			}
			return *this;
		}

		inline Type operator++( int ){
			Type temp_it = *this;
			++(*this);
			return temp_it;
		}

		inline Type & operator--(){
			minus_one();
			return *this;
		}

		inline Type operator--( int ){
			Type temp_it = *this;
			minus_one();
			return temp_it;
		}

		int remaining_blocksize() {
			return m_remaining_blocksize;
		}

		typename seqan::Iterator< String< TValue, TStringSpec > >::Type iter_underlying(){
			return m_iter_underlying;
		}

		typename seqan::Iterator< String< TValue > >::Type iter_insertion(){
			return m_iter_insertion;
		}

		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > iter_tree(){
			return m_iter_tree;
		}

		seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & string(){
			return m_string;
		}

	private:

		inline void minus_one(){
			if( m_remaining_blocksize < (int)cargo(*m_iter_tree).blocksize - 1 ){
				--m_iter_insertion;
				--m_iter_underlying;
				++m_remaining_blocksize;
			}else{
				while( seqan::goPrev( m_iter_tree ) && isDeletion( cargo(*m_iter_tree).op ) ) { }; //TODO: specify behavior for end of string reached -> is reset to tree root regardless of root-node-type atm.
				m_iter_insertion = seqan::begin( m_string.insertion_string() ) + cargo(*m_iter_tree).physical_position + cargo(*m_iter_tree).blocksize - 1;
				m_iter_underlying = seqan::begin( m_string.underlying() ) + cargo(*m_iter_tree).physical_position + cargo(*m_iter_tree).blocksize - 1;
				m_remaining_blocksize = 0;
			}
		}

		seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > const & m_string;
		typename seqan::Iterator< String< TValue, TStringSpec >  >::Type const m_iter_underlying;
		typename seqan::Iterator< String< TValue > >::Type const m_iter_insertion;
		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > const > m_iter_tree;
		int m_remaining_blocksize;
	};
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//meta functions

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
		struct Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > >, Standard >{
			typedef JournalIterator< TValue, TStringSpec, TJournalSpec > Type;
		};

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
		struct Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > >, Rooted >{
			typedef JournalIterator< TValue, TStringSpec, TJournalSpec > Type;
		};


#if 0
	template< typename TValue, typename TStringSpec, typename TJournalSpec >
		struct Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > const >{
			typedef ConstJournalIterator< TValue, TStringSpec, TJournalSpec > Type;
		};
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//begin functions
#if 1
	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TTagSpec>
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		begin(
				String< TValue, Journal< TStringSpec, TJournalSpec > > & str,
				Tag< TTagSpec > const )
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type begin_iter( str );
			return begin_iter;
		}

	template< typename TValue, typename TStringSpec, typename TJournalSpec>
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		begin(
				String< TValue, Journal< TStringSpec, TJournalSpec > > & str )
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type begin_iter( str );
			return begin_iter;
		}
#endif
	template <typename TValue, typename TStringSpec, typename TJournalSpec>
		inline typename Iterator< String< TValue, Journal < TStringSpec, TJournalSpec> > const >::Type
		begin (
				String<TValue, Journal<TStringSpec, TJournalSpec> > const & str	)
		{
			typename Iterator< String< TValue, Journal<TStringSpec, TJournalSpec> > const >::Type begin_iter( str );
			return begin_iter;
		}

	template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TTagSpec>
		inline typename Iterator< String< TValue, Journal < TStringSpec, TJournalSpec> > const, Tag<TTagSpec> const >::Type
		begin (
				String<TValue, Journal<TStringSpec, TJournalSpec> > const & str,
				Tag<TTagSpec> const )
		{
			typename Iterator< String< TValue, Journal<TStringSpec, TJournalSpec> > const, Tag<TTagSpec> const >::Type begin_iter( str );
			return begin_iter;
		}
#if 0
	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TSpec >
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		begin(
				String< TValue, Journal< TStringSpec, TJournalSpec > > const & str,
				Tag< TSpec > const)
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > const >::Type begin_iter( str );
			return begin_iter;
		}
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//end functions

	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TTagSpec >
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		end( String< TValue, Journal< TStringSpec, TJournalSpec > > & str,
				Tag<TTagSpec> const)
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type end_iter( str, seqan::end( journal_tree( str ) ) );
			return end_iter;
		}

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		end( String< TValue, Journal< TStringSpec, TJournalSpec > > & str )
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type end_iter( str, seqan::end( journal_tree( str ) ) );
			return end_iter;
		}

	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TTagSpec >
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		end( String< TValue, Journal< TStringSpec, TJournalSpec > > const & str,
				Tag<TTagSpec> const)
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > const >::Type end_iter( str, seqan::end( journal_tree( str ) ) );
			return end_iter;
		}

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		end( String< TValue, Journal< TStringSpec, TJournalSpec > > const & str )
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > const >::Type end_iter( str, seqan::end( journal_tree( str ) ) );
			return end_iter;
		}
#if 0
	template< typename TValue, typename TStringSpec, typename TJournalSpec, typename TTagSpec >
		inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type
		end( String< TValue, Journal< TStringSpec, TJournalSpec > > const & str,
				Tag<TTagSpec> const)
		{
			typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > const >::Type end_iter( str, seqan::end( journal_tree( str ) ) );
			return end_iter;
		}
#endif
};

#endif /* JOURNALITERATOR_H_ */
