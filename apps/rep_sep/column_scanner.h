/*=========================================================================
  Copyright (C) 2009 by Stephan Aiche

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================
  $Id$
 ==========================================================================*/

#ifndef REPSEP_HEADER_COLUMN_SCANNER_H
#define REPSEP_HEADER_COLUMN_SCANNER_H

#include <seqan/store.h>
#include <map>

//#define REPSEP_DEBUG_COLUMN_SCANNER

// TODO: add DNP method

using namespace seqan;
using namespace std;

struct SimpleColumnConfig {

public:
    double error_probability;
    double test_level;

    SimpleColumnConfig(){
        error_probability = 0.1f;
        test_level = 0.05f;
    }
};

// simple column detection based on 1-column-statistics
struct SimpleColumn{
    SimpleColumnConfig parameters;
};

template <typename TSpec, typename TConfig, typename TId, typename TConsensusCharacter, typename TCandidateColumn>
bool isCandidate(FragmentStore<TSpec, TConfig> const& fragStore,
           TId const contigId,
           TConsensusCharacter const consensus,
           TCandidateColumn & candidate,
           SimpleColumn const algo_spec)
{
    // Get rid of unused variable warnings.
    (void)fragStore;
    (void)contigId;

	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	//typedef typename TFragmentStore::TReadPos TReadPos;
    
    // must be triplet
    typedef typename Value<TCandidateColumn>::Type TAssignedReadChar;
    typedef typename Value<TAssignedReadChar, 1>::Type TConsensusAlphabet;

	// All fragment store element types
	//typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	//typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
    

#ifdef REPSEP_DEBUG_COLUMN_SCANNER
    cout << "--- called SimpleScanner on candidate column with consensus [" << _sequenceCharacter(consensus) << "] --- " << endl;
    for(TSize i = 0 ; i < length(candidate) ; ++i) {
        cout << _sequenceCharacter(candidate[i]) << " ";
    }
    cout << endl;
#endif

    // cout deviations from consensus 
    std::map< TConsensusAlphabet, TSize > vote_map;
    TSize simple_deviations = 0;
    for(TSize i = 0 ; i < length(candidate) ; ++i) {
        // collect all deviations
        if( _sequenceCharacter(candidate[i]) != _sequenceCharacter(consensus) ) {
            ++simple_deviations;
        }
        // remember all characters and their cumulated occurrence
        ++vote_map[_sequenceCharacter(candidate[i])];
    }

    double lambda = length(candidate) * algo_spec.parameters.error_probability;
    double p = p_value(simple_deviations,lambda);
#ifdef REPSEP_DEBUG_COLUMN_SCANNER
    cout << "has p-value " << p << " with test level " << algo_spec.parameters.test_level << endl;
#endif
    return ( p <= algo_spec.parameters.test_level );
}

template <typename TSpec, typename TConfig, typename TId, typename TCandidateColumn>
void refineCandidates(FragmentStore<TSpec, TConfig> const& /*fragStore*/,
           TId const /*contigId*/,
           StringSet<TCandidateColumn> & /*candidate_set*/,
           SimpleColumn const)
{
#ifdef REPSEP_DEBUG_COLUMN_SCANNER
    cout << "refinement for SimpleScanner called --> no refinement will be made" << endl;
#endif
}

#endif
