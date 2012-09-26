/*
 * STNode.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#include <string> ///maybe?!
#include <sstream>
#include <seqan.h>


#ifndef STNODE_H_
#define STNODE_H_

namespace seqan{

//	namespace tree{

		struct Basic{
			Basic() : leaf( true ) { }; //Default-Constructor is a Leaf
			Basic( bool l ) : leaf( l ) { };
			bool leaf;
		};

		template< typename _TCargo, typename TSpec = Basic >
		class STNode{
		public:

			typedef _TCargo TCargo;
			typedef TCargo TValue;
			typedef TSpec TNodeSpec;

			STNode() : m_cargo(), m_spec( TSpec() ){	}; //Default-Constructor should always create a Leaf
			STNode( TCargo const & cargo ) : m_cargo( cargo ), m_spec( TSpec() ) { setLeaf( m_spec, false ); };
			STNode( TCargo & cargo ) : m_cargo( cargo ), m_spec( TSpec() ) { setLeaf( m_spec, false ); };
			STNode( TCargo const cargo, TSpec & spec ) : m_cargo( cargo ), m_spec( spec ) { };
			STNode( TCargo const cargo, TSpec const & spec ) : m_cargo( cargo ), m_spec( spec ) { };
			STNode( TCargo & cargo, TSpec & spec ) : m_cargo( cargo ), m_spec( spec ) { };

			TCargo m_cargo;
			TSpec  m_spec;
		};

		template< typename TCargo, typename TSpec >
		struct Value< STNode< TCargo, TSpec > >{
			typedef TCargo Type;
		};

		template< typename TCargo, typename TSpec >
		struct Spec< STNode< TCargo, TSpec > >{
			typedef TSpec Type;
		};

		template< typename TCargo, typename TSpec >
		inline typename Value< STNode< TCargo, TSpec > >::Type &
		cargo( STNode< TCargo, TSpec > & node ){
			return node.m_cargo;
		}

		template< typename TCargo, typename TSpec >
		inline typename Spec< STNode< TCargo, TSpec > >::Type &
		spec( STNode< TCargo, TSpec > & node ){
			return node.m_spec;
		}

		template< typename TCargo, typename TSpec >
		inline typename Value< STNode< TCargo, TSpec > >::Type const &
		cargo( STNode< TCargo, TSpec > const & node ){
			return node.m_cargo;
		}

		template< typename TCargo, typename TSpec >
		inline typename Spec< STNode< TCargo, TSpec > >::Type const &
		spec( STNode< TCargo, TSpec > const & node ){
			return node.m_spec;
		}

		template< typename TCargo, typename TSpec >
		inline void assignCargo( STNode< TCargo, TSpec > & node, TCargo const & value ){
			cargo( node ) = value;
		}

		template< typename TCargo >
		inline std::string stringify( TCargo const & cargo ){
			std::stringstream sstr;
			sstr << cargo;
			return sstr.str();
		}

		template< typename TValue, typename TSpec >
		inline std::string stringify( String< TValue, TSpec > const & cargo ){
			std::stringstream sstr;
			sstr << cargo;
			return sstr.str();
		}

		template<>
		inline std::string stringify( int const& cargo ){
			std::stringstream sstr;
			sstr << cargo;
			return sstr.str();
		}

		template< typename TCargo, typename TSpec >
		inline seqan::String< char > info( STNode< TCargo, TSpec > const & node ){
			String< char > out;
			if( isLeaf( node ) ){
				out += "Leaf<>";
			}else{
				out += "STNode< ";
				out += stringify( cargo( node ) );
				out += " />";
			}
			return out;
		}

		inline bool isLeaf( Basic & spec ){
			return spec.leaf;
		}

		inline bool isLeaf( Basic const & spec ){
			return spec.leaf;
		}

		inline bool setLeaf( Basic & spec, bool leaf ){
			spec.leaf = leaf;
			return true;
		}

		template< typename TCargo >
		inline bool isLeaf( STNode< TCargo, Basic > & node ){
			return isLeaf( spec( node ) );
		}

		template< typename TCargo >
		inline bool isLeaf( STNode< TCargo, Basic > const & node ){
			return isLeaf( spec( node ) );
		}

		template< typename TCargo >
		inline bool setLeaf( STNode< TCargo, Basic > & node, bool leaf ){
			return setLeaf( spec( node ), leaf );
		}

		template< typename TCargo >
		inline bool makeLeaf( STNode< TCargo, Basic > & node ){
			return setLeaf( spec( node ), true );
		}

//	};
};

#endif /* STNODE_H_ */
