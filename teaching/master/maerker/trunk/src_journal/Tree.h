/*
 * Tree.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#include <seqan.h>
#include <ostream>
#include "TreeIterator.h"

#ifndef TREE_H_
#define TREE_H_

namespace seqan{

		template< typename TNode, typename TSpec = BinaryBasic >
		class SeqTree{
		public:
			typedef seqan::String< TNode, typename DefaultStringSpec< SeqTree< TNode, TSpec > >::Type > TStorage;
			typedef typename seqan::Iterator< TStorage, Standard >::Type TStorageIterator;
			typedef TreeIterator< SeqTree< TNode, TSpec > > TTreeIterator;

			SeqTree() {
				seqan::appendValue( m_storage, TNode(), Generous() );
			};

			SeqTree( typename Value< TNode >::Type cargo ) {
			    seqan::appendValue( m_storage, TNode( cargo ), Generous() );
			    seqan::appendValue( m_storage, TNode(), Generous() );
    		    seqan::appendValue( m_storage, TNode(), Generous() );
			}

			SeqTree( TNode const & root ) {
				assert( isLeaf( root ) && " Single-Node-Tree must contain a Leaf!");
				seqan::appendValue( m_storage, root, Generous() );
			}

			TStorage m_storage;
		};

		template< typename TNode, typename TSpec >
		struct Value< SeqTree< TNode, TSpec > >{
			typedef TNode Type;
		};

		template< typename TNode, typename TSpec >
		struct Spec< SeqTree< TNode, TSpec > >{
			typedef TSpec Type;
		};

		template< typename TNode >
		inline seqan::String< TNode, typename DefaultStringSpec< SeqTree< TNode, BinaryBasic > >::Type > &
		storage( SeqTree< TNode, BinaryBasic > & tree )
		{
			return tree.m_storage;
		}

		template< typename TCargo, typename TNodeSpec, typename TComp >
		inline TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, BinaryBasic > >
		find( SeqTree< STNode< TCargo, TNodeSpec>, BinaryBasic > & tree, TComp const & pattern )
		{
			TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, BinaryBasic > > tree_it( tree );
			while( !isLeaf( *tree_it ) ){
				if( isLess( pattern, cargo( (*tree_it) ) ) ){
					tree_it += leftChildOffset( tree_it.m_index );
				}else if( isGreater( pattern, cargo( (*tree_it) ) ) ){
					tree_it += rightChildOffset( tree_it.m_index );
				}else{
					return tree_it;
				}
			}
			return tree_it;
		}

		template< typename TCargo, typename TNodeSpec, typename TTreeSpec >
		inline TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > >
		begin( SeqTree< STNode< TCargo, TNodeSpec>, TTreeSpec > & tree )
		{
			TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > > tree_it( tree );
			while( goDownLeft( tree_it ) ) { };
			return tree_it;
		}

		template< typename TCargo, typename TNodeSpec, typename TTreeSpec >
		inline TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > >
		end( SeqTree< STNode< TCargo, TNodeSpec>, TTreeSpec > & tree )
		{
			TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > > tree_it( tree );
			while( goDownRight( tree_it ) ) { };
			return tree_it;
		}

		template< typename TCargo, typename TNodeSpec, typename TTreeSpec>
		inline TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > >
		root( SeqTree< STNode< TCargo, TNodeSpec>, TTreeSpec > & tree )
		{
			return TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > >( tree );
		}

		template< typename TCargo, typename TNodeSpec, typename TTreeSpec, typename TInsert >
		inline bool
		insert( SeqTree< STNode< TCargo, TNodeSpec>, TTreeSpec > & tree, TInsert const & value )
		{
			TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > > tree_it = find( tree, value );
			if( isLeaf( *tree_it ) ){
				setLeaf( *tree_it, false );
				assign( *tree_it, value );
				unsigned int lci = leftChildIndex(tree_it.m_index);
				if( lci + 1 >= seqan::length( seqan::storage( tree ) ) ){
					resize( seqan::storage( tree ), lci + 2, Generous() ); //Resize will use the Dafault-Constructor on all new Elements which for a Node should always create a Leaf
				}
				return true;
			}else{
				return false;
			}
		}

		template< typename TCargo, typename TNodeSpec, typename TTreeSpec, typename TValue >
		inline bool
		erase( SeqTree< STNode< TCargo, TNodeSpec>, TTreeSpec > & tree, TValue const & value )
		{
			TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > > tree_it = find( tree, value );
			return erase( tree, tree_it );
		}

		template< typename TCargo, typename TNodeSpec, typename TTreeSpec >
		inline bool
		erase( SeqTree< STNode< TCargo, TNodeSpec>, TTreeSpec > & tree, TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > > tree_it )
		{
			if( !isLeaf( *tree_it ) ){ //we are at an internal node
				TreeIterator< SeqTree< STNode< TCargo, TNodeSpec >, TTreeSpec > > temp_it( tree_it );
				if( goDownRight( temp_it ) ){
					while( goDownLeft( temp_it ) ) { };
					cargo( (*tree_it) ) = cargo( (*temp_it) );
					erase( tree, temp_it );
					return true;
				}else if( goDownLeft( temp_it ) ){
					while( goDownRight( temp_it ) ) { };
					cargo( (*tree_it) ) = cargo( (*temp_it) );
					erase( tree, temp_it );
					return true;
				}else{ //no children, can be deleted
					makeLeaf( *tree_it ); //leafify
					return true;
				}
			}else{
				return false;
			}
		}

};

#endif /* TREE_H_ */
