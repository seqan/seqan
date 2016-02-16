 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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
 ==========================================================================*/

#ifndef SEQAN_HEADER_RAZERS_PARALLEL_H
#define SEQAN_HEADER_RAZERS_PARALLEL_H

#include "tbb/pipeline.h"
#include "tbb/spin_mutex.h"
#include "tbb/task_scheduler_init.h"

namespace seqan
{


	//////////////////////////////////////////////////////////////////////////////
	// Verification Token
	template <
		typename TGenome, 
		typename TReadSet, 
		typename TSpec >
	struct VerifierToken 
	{
		Segment<TGenome, InfixSegment> genomeInf;   // potential match genome region
		TReadSet *readSet;                          // q-gram index
		unsigned rseqNo;                            // read number
		unsigned gseqNo;							// genome number
		char orientation;							// genome strand F/R
	};


	//////////////////////////////////////////////////////////////////////////////
	// Filtration Pipe 
	//   Input:  <none>
	//   Output: verification tokens
	template <
		typename TGenome, 
		typename TSwiftFinder,
		typename TSwiftPattern, 
		typename TOptionSpec >
	class FiltrationPipe: public tbb::filter
	{
	public:
		typedef typename Host<TSwiftPattern>::Type				TReadIndex;
		typedef typename Host<TReadIndex>::Type					TReadSet;
		typedef VerifierToken<TGenome, TReadSet, TOptionSpec>	TToken;
		typedef RazerSOptions<TOptionSpec> const				TOptions;
	    
		static const size_t nTokens = 8;

		TGenome			&genome;
		TSwiftPattern	swiftPattern;
		TSwiftFinder	swiftFinder;
		TOptions const	&options;
	    
		TToken          tokens[nTokens];
		unsigned        nextToken;

		inline FiltrationPipe(TGenome &_genome, TSwiftPattern &_swiftPattern, unsigned gseqNo, char orientation, TOptions &_options) :
			tbb::filter(serial),
			genome(_genome),
			swiftPattern(_swiftPattern),
			swiftFinder(_genome, options.repeatLength, 1),
			options(_options),
			nextToken(0)
		{
			for (unsigned i = 0; i < nTokens; ++i)
			{
				tokens[i].readSet = &indexText(host(swiftPattern));
				tokens[i].orientation = orientation;
				tokens[i].gseqNo = gseqNo;
			}
		}

		void * operator () (void * /*_token*/)
		{
			if (find(swiftFinder, swiftPattern, options.errorRate, options._debugLevel))
			{
				TToken &token = tokens[nextToken];
				nextToken = (nextToken + 1) % nTokens;
				set(token.genomeInf, infix(swiftFinder));
				token.rseqNo = (*swiftFinder.curHit).ndlSeqNo;
				return &token;
			}
			else
				return NULL;
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// Verification Pipe 
	//   Input:  verification tokens
	//   Output: true matches
	template <
		typename TMatches, 
		typename TGenome, 
		typename TReadSet, 
		typename TVerificationPatterns,
		typename TOptionSpec,
		typename TSwiftSpec >
	class VerificationPipe: public tbb::filter
	{
	public:
		typedef typename Value<TMatches>::Type					TMatch;
		typedef typename Size<TGenome>::Type					TSize;
		typedef VerifierToken<TGenome, TReadSet, TOptionSpec>	TToken;
		typedef typename RazerSOptions<TOptionSpec>::TMutex		TMutex;
		typedef RazerSOptions<TOptionSpec>						TOptions;

		TMatches				&matches;
		TVerificationPatterns	&verificationPatterns;
		TOptions				&options;
		int64_t					FP;
		int64_t					TP;	 


		VerificationPipe(TMatches &_matches, TVerificationPatterns &_verificationPatterns, TOptions &_options):
			tbb::filter(parallel),
			matches(_matches),
			verificationPatterns(_verificationPatterns),
			options(_options),
			FP(0),
			TP(0) {}

		~VerificationPipe()
		{
			typename TMutex::scoped_lock lock(options.optionsMutex);
			options.FP += FP;
			options.TP += TP;
		}

		void * operator () (void *_token)
		{
			TToken &token = *reinterpret_cast<TToken*>(_token);
			typename TMutex::scoped_lock lock(options.patternMutex[token.rseqNo]);
/*
		::std::cout<<"Verify: "<<::std::endl;
		::std::cout<<"Genome: "<<token.genomeInf<<"\t" << beginPosition(token.genomeInf) << "," << endPosition(token.genomeInf) << ::std::endl;
		::std::cout<<"Read:   "<<(*token.readSet)[token.rseqNo]<<::std::endl;
*/	        
			TMatch m;
			if (matchVerify(m, token.genomeInf, token.rseqNo, *token.readSet, verificationPatterns, options, TSwiftSpec()))
			{
				// transform coordinates to the forward strand
				if (token.orientation == 'R') 
				{
					TSize gLength = length(host(token.genomeInf));
					TSize temp = m.gBegin;
					m.gBegin = gLength - m.gEnd;
					m.gEnd = gLength - temp;
				}
				m.rseqNo = token.rseqNo;
				m.gseqNo = token.gseqNo;
				m.orientation = token.orientation;

				if (!options.spec.DONT_DUMP_RESULTS)
				{
					typename TMutex::scoped_lock lock(options.matchMutex);
					appendValue(matches, m);
					if (length(matches) > options.compactThresh)
					{
						typename Size<TMatches>::Type oldSize = length(matches);
						maskDuplicates(matches);
						compactMatches(matches, options);
						options.compactThresh += (options.compactThresh >> 1);
						if (options._debugLevel >= 2)
							::std::cerr << '(' << oldSize - length(matches) << " matches removed)";
					}
				}

				++TP;
			} else
				++FP;

			return NULL;
		}
	};


//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TOptionSpec >
void findReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPattern,
	TVerifier &forwardPatterns,
	char orientation,				// q-gram index of reads
	RazerSOptions<TOptionSpec> &options)
{
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Value<TMatches>::Type					TMatch;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	tbb::pipeline pipeline;

	FiltrationPipe<
		TGenome,
        TSwiftFinder,
		TSwiftPattern,
		TOptionSpec > filtrationPipe(genome, swiftPattern, gseqNo, orientation, options);
    pipeline.add_filter(filtrationPipe);

	VerificationPipe<
		TMatches,
		TGenome,
		TReadSet,
		TVerifier,
		TOptionSpec,
		TSwiftSpec > verificationPipe(matches, forwardPatterns, options);
	pipeline.add_filter(verificationPipe);
    
	pipeline.run(filtrationPipe.nTokens);
	pipeline.clear();
}

}

#endif
