/*
 * TreeIterator.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#ifndef TREEITERATOR_H_
#define TREEITERATOR_H_

namespace seqan{

		/*****************************************************************
		** Index and Offset Calculations for String-Stored Binary Tree  **
		*****************************************************************/

		inline unsigned int
		leftChildIndex( unsigned int parent_index )
		{
			return 2*parent_index + 1;
		}

		inline unsigned int
		rightChildIndex( unsigned int parent_index )
		{
			return 2*parent_index + 2;
		}

		inline unsigned int
		parentIndex( unsigned int child_index )
		{
			if( child_index % 2 == 0 ){
				return ( child_index - 2 ) / 2;
			}else{
				return ( child_index - 1 ) / 2;
			}
		}

		inline bool
		isLeftChild( unsigned int child_index )
		{
			return ( child_index % 2 == 1 );
		}

		inline bool
		isRightChild( unsigned int child_index )
		{
			return ( child_index != 0 ) && ( child_index % 2 == 0 );
		}

		inline unsigned int
		leftChildOffset( unsigned int parent_index ){
			return parent_index + 1;
		}

		inline unsigned int rightChildOffset( unsigned int parent_index ){
			return parent_index + 2;
		}

		//This Offset has to be subtracted rather than added to yield the correct parent_index
		inline unsigned int parentOffset( unsigned int child_index ){
			if( child_index % 2 == 0 ){
				return ( child_index + 2 ) / 2;
			}else{
				return ( child_index + 1 ) / 2;
			}
		}

		template< typename TNode, typename TSpec >
		class SeqTree;

		struct BinaryBasic;

		template< typename TTree >
		struct DefaultStringSpec{
			typedef seqan::Alloc<> Type;
		};

		template< typename TNode >
		struct DefaultStringSpec< SeqTree< TNode, BinaryBasic > >{
			typedef seqan::Alloc<> Type;
		};

		//template< typename TCargo, typename TNodeSpec = Basic, typename TTreeSpec = BinaryBasic>
		template< typename TTree >
		struct TreeIterator{
			typedef TreeIterator< TTree > Type;
			typedef typename Value< TTree >::Type TNode;
			typedef typename seqan::Iterator< seqan::String< TNode, typename DefaultStringSpec< TTree >::Type > >::Type TStorageIterator;

			TreeIterator( TTree & tree ) : m_iter( begin( storage( tree ) ) ), m_index( 0 ) { };
			TreeIterator( TStorageIterator iter, unsigned int index ) : m_iter( iter ), m_index( index ) { };
			TreeIterator( Type & other ) : m_iter( other.m_iter ), m_index( position( other ) ) { };
			TreeIterator( Type const & other ) : m_iter( other.m_iter ), m_index( position( other ) ) { };

			inline void operator= ( Type & other ){
				m_iter = other.m_iter;
				m_index = position( other );
			}

			inline void operator= ( Type const & other ){
				m_iter = other.m_iter;
				m_index = position( other );
			}

			inline Type operator+( unsigned int by ){
				return TreeIterator( m_iter + by, m_index + by );
			}

			inline Type operator+=( unsigned int by ){
				m_iter += by;
				m_index += by;
				return *this;
			}

			inline Type operator-=( unsigned int by ){
				m_iter -= by;
				m_index -= by;
				return *this;
			}

			inline Type operator++( int ){
				Type temp_it = *this;
				++m_iter;
				++m_index;
				return temp_it;
			}

			inline Type & operator++(){
				++m_iter;
				++m_index;
				return *this;
			}

			inline Type operator--( int ){
				Type temp_it = *this;
				--m_iter;
				--m_index;
				return temp_it;
			}

			inline Type & operator--(){
				--m_iter;
				--m_index;
				return *this;
			}

			inline TNode & operator*(){
				return *m_iter;
			}

			inline TNode operator*() const{
				return *m_iter;
			}

			TStorageIterator m_iter;
			unsigned int m_index;
		};

		template< typename TTree >
		inline size_t
		position( TreeIterator< TTree > & iter, TTree & tree ){
			return iter.m_index;
		}

		template< typename TTree >
		inline size_t
		position( TreeIterator< TTree > & iter ){
			return iter.m_index;
		}

		template< typename TTree >
		inline size_t
		position( TreeIterator< TTree > const & iter ){
			return iter.m_index;
		}

		template< typename TTree >
		inline bool goPrev( TreeIterator< TTree > & iter ){
			if( isLeaf( *iter ) ){
				return false;
			}else{
				if( goDownLeft( iter ) ){
					while( goDownRight( iter ) ) { };
					return true;
				}else if( isRightChild( position( iter ) ) ){
					return goUp( iter );
				}else{
					while( isLeftChild( position( iter ) ) ){
						goUp( iter );
					}
					return goUp( iter );
				}
			}
		}

		template< typename TTree >
		inline bool goNext( TreeIterator< TTree > & iter ){
			//std::cout << "At: " << info( *iter ) << std::endl;
			if( isLeaf( *iter ) ){
				return false;
			}else{
				if( goDownRight( iter ) ){
					while( goDownLeft( iter ) ) { };
					return true;
				}else if( isLeftChild( position( iter ) ) ){
					return goUp( iter );
				}else{
					while( isRightChild( position( iter ) ) ){
						goUp( iter );
					}
					return goUp( iter );
				}
			}
		}

		template< typename TTree >
		inline bool goDownLeft( TreeIterator< TTree > & iter ){
			if( isLeaf( *iter ) ){
				return false;
			}else if( !isLeaf( *( iter + leftChildOffset( position( iter ) ) ) ) ){
				iter += leftChildOffset( position( iter ) );
				return true;
			}else{
				return false;
			}
		}

		template< typename TTree >
		inline bool goDownRight( TreeIterator< TTree > & iter ){
			if( isLeaf( *iter ) ){
				return false;
			}else if( !isLeaf( *( iter + rightChildOffset( position( iter ) ) ) ) ){
				iter += rightChildOffset( position( iter ) );
				return true;
			}else{
				return false;
			}
		}

		template< typename TTree >
		inline bool goUp( TreeIterator< TTree > & iter ){
			if( position( iter ) != 0 ){ //Not in Root
				iter -= parentOffset( position( iter ) );
				return true;
			}else{
				return false;
			}
		}

		template< typename TTree >
		std::ostream & operator << ( std::ostream& os, TreeIterator< TTree > const & iter ) {
			os << info(*(iter.m_iter));
			os << "|";
			os << position( iter );
			os << "(" << &position( iter ) << ")";
			return os;
		}
};

#endif /* TREEITERATOR_H_ */
