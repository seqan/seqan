///A tutorial about parsing a blast report.
#include <iostream>
#include <seqan/blast.h>

using namespace seqan;


template <typename TFile>
void read_blast_report(TFile & strm)
{

///First the type of Blast report needs to be specified.
///In this case, a BlastN report is parsed and for each alignment in the report (i.e. each HSP) we choose to parse 
///all the possible information delivered with the alignment (e.g. scores, percentage of gaps, orientation... see @Class.BlastHsp@).

	typedef BlastHsp<BlastN, FullInfo> TBlastHsp;

///The @Spec.StreamReport@ specialization determines that the alignments are parsed when iterating over them (i.e. only 
///one alignment is stored at a time, as opposed to the @Spec.StoreReport@ specialization).

	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;

	
///Our Blast report object
	TBlastReport blast;

///Counters
	unsigned hspcount = 0;
	unsigned hitcount = 0;
	unsigned highsignif = 0;

	while(!atEnd(strm,blast)) 
	{
///Get the current Blast report (there can be multiple reports in one file).
		read(strm,blast,Blast());
		std::cout << "Query: "<<getQueryName(blast) <<"\n";
		std::cout << "Database: "<<getDatabaseName(blast) <<"\n\n";

///Iterate over hits
		THitIterator hit_it(blast); 
		for(; !atEnd(strm,hit_it); goNext(strm,hit_it)) 
		{
			++hitcount;
			TBlastHit hit = getValue(strm,hit_it);
			std::cout << " Hit: " <<name(hit) <<"\n\n";

			/// iterate over alignments (HSPs)
			THspIterator hsp_it(hit);
			for(; !atEnd(strm,hsp_it); goNext(strm,hsp_it)) 
			{
				++hspcount;
 				TBlastHsp hsp = getValue(strm,hsp_it);
				
///Do something with the alignment, e.g.
///output score and length of alignment.
				std::cout << "  Score  = " << getBitScore(hsp) << "\n";
				std::cout << "  Length = " << length(hsp) << "\n\n";
///Count alignments with highly significant e-values.
				if(getEValue(hsp)<0.01)
					++highsignif;

			}
		}
	}
	std::cout <<"Total number of Hits: "<< hitcount<<std::endl;
	std::cout <<"Total number of HSPs: "<< hspcount<<std::endl;
	std::cout <<"Number of highly significant HSPs: "<< highsignif<<std::endl;

}


int main()
{
	std::fstream strm;
	strm.open("ecoln.out", std::ios_base::in | std::ios_base::binary);
	read_blast_report(strm);
	return 0;
}
