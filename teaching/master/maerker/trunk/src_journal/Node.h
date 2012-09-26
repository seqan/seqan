/*
 * Node.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#include <string> ///maybe?!
#include <sstream>
#include <seqan.h>

#ifndef NODE
#define NODE

namespace seqan{

//	namespace tree{

		struct Basic{
			Basic() : leaf( true ) { }; //Default-Constructor is a Leaf
			Basic( bool l ) : leaf( l ) { };
			bool leaf;
		};

		template< typename _TCargo, typename TSpec = Basic >
		class Node{
		public:

			typedef _TCargo TCargo;
			typedef TCargo TValue;
			typedef TSpec TNodeSpec;

			Node() : m_cargo(), m_spec( TSpec() ){	}; //Default-Constructor should always create a Leaf
			Node( TCargo const & cargo ) : m_cargo( cargo ), m_spec( TSpec() ) { setLeaf( m_spec, false ); };
			Node( TCargo & cargo ) : m_cargo( cargo ), m_spec( TSpec() ) { setLeaf( m_spec, false ); };
			Node( TCargo const cargo, TSpec & spec ) : m_cargo( cargo ), m_spec( spec ) { };
			Node( TCargo const cargo, TSpec const & spec ) : m_cargo( cargo ), m_spec( spec ) { };
			Node( TCargo & cargo, TSpec & spec ) : m_cargo( cargo ), m_spec( spec ) { };

			TCargo const & cargo() const{
				return m_cargo;
			}

			TCargo & cargo(){
				return m_cargo;
			}

			TSpec const & spec() const{
				return m_spec;
			}

			TSpec & spec(){
				return m_spec;
			}
		private:
			TCargo m_cargo;
			TSpec  m_spec;
		};

		template< typename TCargo, typename TSpec >
		struct Value< Node< TCargo, TSpec > >{
			typedef TCargo Type;
		};

		template< typename TCargo, typename TSpec >
		struct Spec< Node< TCargo, TSpec > >{
			typedef TSpec Type;
		};

		template< typename TCargo, typename TSpec >
		inline void assign( Node< TCargo, TSpec > & node, TCargo const & value ){
			node.cargo() = value;
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
		inline seqan::String< char > info( Node< TCargo, TSpec > const & node ){
			String< char > out;
			if( isLeaf( node ) ){
				out += "Leaf<>";
			}else{
				out += "Node< ";
				out += stringify( node.cargo() );
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

		template< typename TCargo, typename TSpec >
		inline TSpec & spec( Node< TCargo, TSpec > & node ){
			return node.spec();
		}

		template< typename TCargo >
		inline bool isLeaf( Node< TCargo, Basic > & node ){
			return isLeaf( node.spec() );
		}

		template< typename TCargo >
		inline bool isLeaf( Node< TCargo, Basic > const & node ){
			return isLeaf( node.spec() );
		}

		template< typename TCargo >
		inline bool setLeaf( Node< TCargo, Basic > & node, bool leaf ){
			return setLeaf( spec( node ), leaf );
		}

		template< typename TCargo >
		inline bool makeLeaf( Node< TCargo, Basic > & node ){
			return setLeaf( spec( node ), true );
		}

//	};
};

#endif /* NODE_H_ */
