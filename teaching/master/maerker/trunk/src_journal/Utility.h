/*
 * Utility.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#ifndef UTILITY_H_
#define UTILITY_H_

namespace seqan{
	namespace util{
	    unsigned int randomNumber( int upper_bound )
	    {
	        int value = ( (float)rand() / RAND_MAX )*( upper_bound + 1 );

	        if( value > upper_bound ){
	            value = 0;
	        }

	        return value;
	    }

	    template< typename TSpec >
	    void generate_random_string( unsigned int length, seqan::String< char, TSpec > & the_string, seqan::String< char > const & alphabet_string ){

	        seqan::resize( the_string, length );
	        typename seqan::Iterator< seqan::String< char, TSpec > >::Type it = seqan::begin( the_string );
	        typename seqan::Iterator< seqan::String< char, TSpec > >::Type the_end = seqan::end( the_string );

	        while( it != the_end ){
	            *it = alphabet_string[ randomNumber( seqan::length( alphabet_string ) - 1 ) ];
	            ++it;
	        }

	    }
	};
};

#endif /* UTILITY_H_ */
