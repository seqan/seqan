// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_TEST_BLAST_PARSING_H
#define SEQAN_HEADER_TEST_BLAST_PARSING_H

using namespace std;
using namespace seqan;


/* IOREV
 *
 * uses some custom IO, test will have to be adapted anyways
 *
 */


			
//////////////////////////////////////////////////////////////////////////////
template<typename T>
void Test_BlastStoreReport() {
//IOREV


	typedef BlastHsp<BlastP,BasicInfo> TBlastHsp;
	typedef BlastReport<TBlastHsp,StoreReport<FullInfo> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;

	//typedef StringSet<String<AminoAcid>,Dependent<> > TStringSet;
 	//typedef Graph<Alignment<TStringSet> > TAliGraph;
	
	char buffer[1024];
	strcpy(buffer, SEQAN_PATH_TO_ROOT());
	strcat(buffer, "/extras/tests/blast/ecolip.out");

	std::fstream strm;
	strm.open(buffer, ios_base::in | ios_base::binary);

	TBlastReport blast;
	TBlastReport blast2;

	int alicount = 0;
	int hitcount = 0;
	int hspcount = 0;
	int repcount = 0;

	while(!atEnd(strm,blast)) 
	{
		getNext(strm,blast);   // complete report is now parsed and all the hits (and hsps) are stored
		SEQAN_ASSERT(queryName(blast)==getQueryName(blast));
		SEQAN_ASSERT(databaseName(blast)==getDatabaseName(blast));
		SEQAN_ASSERT(repcount!= 0 || (repcount==0 && databaseName(blast)=="ecoliKurz.aa "));
		SEQAN_ASSERT(repcount!= 0 || (repcount==0 && queryName(blast)=="gi|1786182|gb|AAC73112.1| (AE000111) thr operon leader peptide"));
		SEQAN_ASSERT(repcount!= 5 || (repcount==5 && queryName(blast)=="gi|1786187|gb|AAC73117.1| (AE000111) orf, hypothetical protein"));
		

		//std::cout << "Query  : " << queryName(blast)<<" = "<<	getQueryName(blast) <<"\n";
		//std::cout << "DB     : " << databaseName(blast)<<" = "<< getDatabaseName(blast)<<"\n";
		//std::cout << "NumHits: "<< numHits(blast)<<"  \nNumHsps: "<< numHsps(blast) <<"\n";
		hspcount+=numHsps(blast);
		hitcount+=numHits(blast);
		if(repcount==1) blast2 = blast;
		++repcount;
		SEQAN_ASSERT(repcount!=1 || (repcount==1 && length((blast.hits)[0])==21));
		SEQAN_ASSERT(repcount!=2 || (repcount==2 && getLength((blast.hits)[0])==820));
		SEQAN_ASSERT(repcount!=2 || (repcount==2 && length((blast.hits)[1])==342));

		THitIterator hit_it(blast); 
		SEQAN_ASSERT(atBegin(hit_it));
		++hit_it;
		SEQAN_ASSERT(!atBegin(hit_it));
		--hit_it;
		SEQAN_ASSERT(atBegin(hit_it));
		goEnd(hit_it);
		SEQAN_ASSERT(atEnd(hit_it));
		goBegin(hit_it);
		SEQAN_ASSERT(atBegin(hit_it));
		for(; !atEnd(hit_it); goNext(hit_it)) 
		{
			TBlastHit hit2 = getValue(hit_it);
			THitIterator hit_it2;
			hit_it2 = hit_it; 
			SEQAN_ASSERT(hit_it==hit_it2);
			TBlastHit hit = *hit_it2;
			SEQAN_ASSERT(getName(hit) == name(hit));
			SEQAN_ASSERT((repcount != 3 || !atBegin(hit_it)) || (repcount==3 && getName(hit) == "gb|AAC73114.1| (AE000111) homoserine kinase [Escherichia coli]"));
			
			THspIterator hsp_it2(hit);
			THspIterator hsp_it(hit);
			SEQAN_ASSERT(atBegin(hsp_it2));
			goNext(hsp_it2);
			SEQAN_ASSERT(!atBegin(hsp_it2) || atEnd(hsp_it2));
			goBegin(hsp_it2);
			SEQAN_ASSERT(atBegin(hsp_it2));
			SEQAN_ASSERT(hsp_it==hsp_it2);
			goEnd(hsp_it2);
			SEQAN_ASSERT(hsp_it!=hsp_it2);
			SEQAN_ASSERT(atEnd(hsp_it2));
			--hsp_it2;
			SEQAN_ASSERT(!atEnd(hsp_it2));
			hsp_it2 = hsp_it;
			for(; !atEnd(hsp_it) && !atEnd(hsp_it2); ++hsp_it, ++hsp_it2) 
			{

				TBlastHsp hsp = getValue(hsp_it);
				TBlastHsp hsp2 = *hsp_it2;
				SEQAN_ASSERT(getEValue(hsp) == eValue(hsp2));
				SEQAN_ASSERT(getQueryBegin(hsp) == queryBegin(hsp2));
				SEQAN_ASSERT(databaseBegin(hsp) == getDatabaseBegin(hsp2));
				SEQAN_ASSERT(getQueryEnd(hsp) == queryEnd(hsp2));
				SEQAN_ASSERT(databaseEnd(hsp) == getDatabaseEnd(hsp2));
				SEQAN_ASSERT(getQueryAlignmentString(hsp) == queryAlignmentString(hsp2));
				SEQAN_ASSERT(databaseAlignmentString(hsp) == getDatabaseAlignmentString(hsp2));
				if(eValue(hsp)<0.02) ++alicount;
				//
			}
		}
		if(atEnd(strm,blast))
		{
			SEQAN_ASSERT(eValueCutoff(blast)==getEValueCutoff(blast) && eValueCutoff(blast)==10);
			SEQAN_ASSERT(matrixName(blast)==getMatrixName(blast) && matrixName(blast)=="BLOSUM62");
	//		cout << "NumHsps insgesamt: "<<hspcount<<"\n";
	//		cout << "ECutoff: "<< eValueCutoff(blast)<< " = "<< getEValueCutoff(blast)<< "\n";
	//		cout << "Matrix : "<< matrixName(blast)<<" = "<< getMatrixName(blast)<< "\n";
			SEQAN_ASSERT(gapsAllowed(blast));
			SEQAN_ASSERT(gapOpen(blast)==getGapOpen(blast) && gapOpen(blast)==11);
			SEQAN_ASSERT(gapExtension(blast)==getGapExtension(blast) && gapExtension(blast)==1);
//			cout << "GapOpen: "<< gapOpen(blast) << " = " << getGapOpen(blast) <<"\n";
//			cout << "GapExt : "<< gapExtension(blast)<<" = "<<getGapExtension(blast)<<"\n";
			SEQAN_ASSERT(lambda(blast)==getLambda(blast) && lambda(blast)==0.311f);
			SEQAN_ASSERT(kappa(blast)==getKappa(blast) && kappa(blast)==0.134f);
			SEQAN_ASSERT(gappedLambda(blast)==getGappedLambda(blast) && gappedLambda(blast)==0.267f);
			SEQAN_ASSERT(gappedKappa(blast)==getGappedKappa(blast) && gappedKappa(blast)==0.0410f);
			SEQAN_ASSERT(entropy(blast)==getEntropy(blast) && entropy(blast)==0.378f);
			SEQAN_ASSERT(gappedEntropy(blast)==getGappedEntropy(blast) && gappedEntropy(blast)==0.140f);
			
			//cout << "Lambda : " <<lambda(blast)<<" = " <<getLambda(blast) << "\n";
			//cout << "K : "<< kappa(blast)<< " = "<<getKappa(blast)<<"  H : " << entropy(blast)<< " = "<<getEntropy(blast)<<"\n";
			//cout << "GLambda: " <<gappedLambda(blast)<<" = " <<getGappedLambda(blast) << "\n";
			//cout << "GK: "<< gappedKappa(blast)<< " = "<<getGappedKappa(blast)<<"  GH: " << gappedEntropy(blast)<< " = "<<getGappedEntropy(blast)<<"\n";
		}

	}

	SEQAN_ASSERT(hitcount==56);
	SEQAN_ASSERT(hspcount==58);
	SEQAN_ASSERT(alicount==15);
	SEQAN_ASSERT(repcount==16);

	TBlastReport blast1(blast);
	SEQAN_ASSERT(matrixName(blast)==getMatrixName(blast) && matrixName(blast)=="BLOSUM62");
	

}

//////////////////////////////////////////////////////////////////////////////
template<typename T>
void Test_BlastStoreReportBasic() {
//IOREV


	typedef BlastHsp<BlastP,BasicInfo> TBlastHsp;
	typedef BlastReport<TBlastHsp,StoreReport<BasicInfo> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;

	//typedef StringSet<String<AminoAcid>,Dependent<> > TStringSet;
 	//typedef Graph<Alignment<TStringSet> > TAliGraph;
	
	char buffer[1024];
	strcpy(buffer, SEQAN_PATH_TO_ROOT());
	strcat(buffer, "/extras/tests/blast/ecolip.out");
	
	std::fstream strm;
	strm.open(buffer, ios_base::in | ios_base::binary);
	
	TBlastReport blast;
	TBlastReport blast2;

	int alicount = 0;
	int hitcount = 0;
	int hspcount = 0;
	int repcount = 0;

	while(!atEnd(strm,blast)) 
	{
		read(strm,blast,Blast());   // complete report is now parsed and all the hits (and hsps) are stored
		SEQAN_ASSERT(queryName(blast)==getQueryName(blast));
		SEQAN_ASSERT(databaseName(blast)==getDatabaseName(blast));
		SEQAN_ASSERT(repcount!= 0 || (repcount==0 && databaseName(blast)=="ecoliKurz.aa "));
		SEQAN_ASSERT(repcount!= 0 || (repcount==0 && queryName(blast)=="gi|1786182|gb|AAC73112.1| (AE000111) thr operon leader peptide"));
		SEQAN_ASSERT(repcount!= 5 || (repcount==5 && queryName(blast)=="gi|1786187|gb|AAC73117.1| (AE000111) orf, hypothetical protein"));
		

		//std::cout << "Query  : " << queryName(blast)<<" = "<<	getQueryName(blast) <<"\n";
		//std::cout << "DB     : " << databaseName(blast)<<" = "<< getDatabaseName(blast)<<"\n";
		//std::cout << "NumHits: "<< numHits(blast)<<"  \nNumHsps: "<< numHsps(blast) <<"\n";
		hspcount+=numHsps(blast);
		hitcount+=numHits(blast);
		if(repcount==1) blast2 = blast;
		++repcount;

		THitIterator hit_it(blast); 
		THitIterator hit_it2(hit_it);
		++hit_it2;
		SEQAN_ASSERT(hit_it!=hit_it2);
		SEQAN_ASSERT(&hostReport(hit_it)==&blast && &hostReport(hit_it2)==&hostReport(hit_it) );
		for(; !atEnd(hit_it); goNext(hit_it)) 
		{
			TBlastHit hit = getValue(hit_it);
			TBlastHit hit2;
			hit2 = hit;
			THspIterator hsp_it(hit2);
			THspIterator hsp_it2(hsp_it);
			SEQAN_ASSERT(hsp_it==hsp_it2);
			for(; !atEnd(hsp_it); goNext(hsp_it)) 
			{
				TBlastHsp hsp = getValue(hsp_it);
				if(eValue(hsp)<0.02) ++alicount;
				//
			}
		}
	}

	SEQAN_ASSERT(hitcount==56);
	SEQAN_ASSERT(hspcount==58);
	SEQAN_ASSERT(alicount==15);
	SEQAN_ASSERT(repcount==16);

	SEQAN_ASSERT(numHits(blast2)==3);
	SEQAN_ASSERT(numHsps(blast2)==3);

	SEQAN_ASSERT(queryName(blast2)=="gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,");
	TBlastReport blast1(blast);
	SEQAN_ASSERT(queryName(blast1)==getQueryName(blast));
	

}


//////////////////////////////////////////////////////////////////////////////
template<typename T>
void Test_BlastParsing(BlastN) 
{
//IOREV

	typedef BlastHsp<BlastN, FullInfo> TBlastHsp;
	typedef ::std::fstream TFile;
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;
	//typedef typename Size<TBlastReport>::Type TSize;

	typedef StringSet<String<Dna>,Dependent<> > TStringSet;
 	typedef Graph<Alignment<TStringSet> > TAliGraph;

	
	char buffer[1024];
	strcpy(buffer, SEQAN_PATH_TO_ROOT());
	strcat(buffer, "/extras/tests/blast/ecoln.out");
	
	std::fstream strm2;
	strm2.open(buffer, ios_base::in | ios_base::binary);
	
		
	int alicount = 0;
	int hitcount = 0;
	int hspcount = 0;
	TBlastReport blast2;
	while(!atEnd(strm2,blast2)) 
	{
		read(strm2,blast2,Blast());
		SEQAN_ASSERT(databaseName(blast2) == "ecoli.nt ");
		SEQAN_ASSERT(queryName(blast2) == "Test");
		THitIterator hit_it(blast2); 
		for(; !atEnd(strm2,hit_it); goNext(strm2,hit_it)) 
		{
			++hitcount;
			TBlastHit hit = getValue(strm2,hit_it);
			SEQAN_ASSERT((hitcount!=8) || (hitcount==8 && name(hit) == "gb|AE000167.1|AE000167 Escherichia coli K-12 MG1655 section 57 of 400 of the complete genome"));
			THspIterator hsp_it(hit);
			for(; !atEnd(strm2,hsp_it); goNext(strm2,hsp_it)) 
			{
				++hspcount;
 				TBlastHsp hsp = getValue(strm2,hsp_it);
				if(eValue(hsp) < 0.02)
				{
					Align< String<Dna>, ArrayGaps> ali;
					getAlignment(hsp,ali,UnknownSource());
					TAliGraph ali_g;
					getAlignment(hsp,ali_g); //hit ID
					++alicount;
				}
				SEQAN_ASSERT((hspcount != 1) ||(hspcount == 1 && score(hsp)== 140));
				SEQAN_ASSERT((hspcount != 2) ||(hspcount == 2 && bitScore(hsp)== 182 && percentGaps(hsp)==0));
				SEQAN_ASSERT((hspcount != 3) ||(hspcount == 3 && percentIdentity(hsp)== 95 && getPercentGaps(hsp)==4));
				SEQAN_ASSERT((hspcount != 6) ||(hspcount == 6 && !databaseOrientationPlus(hsp)));
				SEQAN_ASSERT((hspcount != 7) ||(hspcount == 7 && eValue(hsp)== 0.52));
				SEQAN_ASSERT((hspcount != 12) ||(hspcount == 12 && getBlastMatchScore(hsp)== 15 && getQueryBegin(hsp)==472 && getQueryEnd(hsp)==486 && getQueryAlignmentString(hsp) == "cacctggtggcgatg" && getDatabaseAlignmentString(hsp) == "cacctggtggcgatg"));
				SEQAN_ASSERT((hspcount != 17) ||(hspcount == 17 && getBitScore(hsp)== 28.2f && queryBegin(hsp)==476 && queryEnd(hsp)==489));
				SEQAN_ASSERT((hspcount != 20) ||(hspcount == 20 && eValue(hsp)== 8.1 && getDatabaseBegin(hsp)==2092 && getDatabaseEnd(hsp)==2079));
				SEQAN_ASSERT((hspcount != 34) ||(hspcount == 34 && databaseBegin(hsp)==5787 && databaseEnd(hsp)==5770 && !databaseOrientationPlus(hsp)));

				if(hspcount == 8 || hspcount == 19)
				{
					TBlastHsp hsp2;
				    hsp2 = hsp;
					SEQAN_ASSERT(score(hsp) == getBlastMatchScore(hsp2));
					SEQAN_ASSERT(bitScore(hsp) == getBitScore(hsp2));
					SEQAN_ASSERT(getEValue(hsp) == eValue(hsp2));
					SEQAN_ASSERT(percentIdentity(hsp) == getPercentIdentity(hsp2));
					SEQAN_ASSERT(percentGaps(hsp) == getPercentGaps(hsp2));
					SEQAN_ASSERT(getNumGaps(hsp) == numGaps(hsp2));
					SEQAN_ASSERT(queryOrientationPlus(hsp) == queryOrientationPlus(hsp2));
					SEQAN_ASSERT(databaseOrientationPlus(hsp) == databaseOrientationPlus(hsp2));
					SEQAN_ASSERT(getQueryBegin(hsp) == queryBegin(hsp2));
					SEQAN_ASSERT(databaseBegin(hsp) == getDatabaseBegin(hsp2));
					SEQAN_ASSERT(getQueryEnd(hsp) == queryEnd(hsp2));
					SEQAN_ASSERT(databaseEnd(hsp) == getDatabaseEnd(hsp2));
					SEQAN_ASSERT(getQueryAlignmentString(hsp) == queryAlignmentString(hsp2));
					SEQAN_ASSERT(databaseAlignmentString(hsp) == getDatabaseAlignmentString(hsp2));
				}

			}
		}
	}



	SEQAN_ASSERT(hitcount==28);
	SEQAN_ASSERT(hspcount==34);
	SEQAN_ASSERT(alicount==5);

}


//////////////////////////////////////////////////////////////////////////////
template<typename T>
void Test_BlastParsing(BlastP) {
//IOREV


	typedef BlastHsp<BlastP, FullInfo> TBlastHsp;
	typedef ::std::fstream TFile;
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;
	//typedef typename Size<TBlastReport>::Type TSize;

	typedef StringSet<String<AminoAcid>,Dependent<> > TStringSet;
 	typedef Graph<Alignment<TStringSet> > TAliGraph;

	char buffer[1024];
	strcpy(buffer, SEQAN_PATH_TO_ROOT());
	strcat(buffer, "/extras/tests/blast/ecolip.out");
	
	std::fstream strm;
	strm.open(buffer, ios_base::in | ios_base::binary);
	
	int alicount = 0;
	int hitcount = 0;
	int hspcount = 0;
	int repcount = 0;

	TBlastReport blast;
	while(!atEnd(strm,blast)) 
	{
		read(strm,blast,Blast());
		SEQAN_ASSERT(getDatabaseName(blast) == "ecoliKurz.aa ");
		SEQAN_ASSERT(repcount!= 0 || (repcount==0 && getQueryName(blast) == "gi|1786182|gb|AAC73112.1| (AE000111) thr operon leader peptide"));
		THitIterator hit_it(blast); 
		for(; !atEnd(strm,hit_it); goNext(strm,hit_it)) 
		{
			++hitcount;
			TBlastHit hit = getValue(strm,hit_it);
			SEQAN_ASSERT(hitcount!=3 || (hitcount==3 && "gb|AAC76950.1| (AE000471) UDP-N-acetylenolpyruvoylglucosamine reductase"));
			THspIterator hsp_it(hit);
			for(; !atEnd(strm,hsp_it); goNext(strm,hsp_it)) 
			{
				++hspcount;
				TBlastHsp hsp = getValue(strm,hsp_it);
 				if(hspcount == 80 )
				{
					std::cout << "  HspQueryBegin: "<<getQueryBegin(hsp)<<"\n";
 				    std::cout << "  HspQueryEnd  : "<<getQueryEnd(hsp)<<"\n";
 					std::cout << "  HspDataBBegin: "<<getDatabaseBegin(hsp)<<"\n";
 					std::cout << "  HspDataBEnd  : "<<getDatabaseEnd(hsp)<<"\n";
 					std::cout << "         Score : "<<getBitScore(hsp)<<"\n";
 					std::cout << "        Expect : "<<getEValue(hsp)<<"\n";
 					std::cout << "      Identity : "<<getPercentIdentity(hsp)<<"\n";
 					std::cout << "     Positives : "<<getPercentPositives(hsp)<<"\n";
 					std::cout << "          Gaps : "<<getPercentGaps(hsp)<<"\n";
 					std::cout << "       AbsGaps : "<<getNumGaps(hsp)<<"\n";
					std::cout << "    QueryFrame : "<<queryFrame(hsp)<<"\n";
					std::cout << "QueryAliString : "<<queryAlignmentString(hsp) << "\n\n";
				}
				if(eValue(hsp) < 0.02)
				{
					Align< String<Dna>, ArrayGaps> ali;
					getAlignment(hsp,ali,UnknownSource());
					TAliGraph ali_g;
					getAlignment(hsp,ali_g); //hit ID
					++alicount;
				}
				SEQAN_ASSERT((hspcount != 1) ||(hspcount == 1 && getBitScore(hsp)== 20 && getDatabaseAlignmentString(hsp)=="MKRISTTITTTITITTGNGAG"));
				SEQAN_ASSERT((hspcount != 3) ||(hspcount == 3 && getPercentIdentity(hsp)== 34 && getPercentPositives(hsp)==46 && getQueryAlignmentString(hsp)=="LSYFGAKVLHPRTITPIAQFQIPCLIKNTGNP"));
				SEQAN_ASSERT((hspcount != 11) ||(hspcount == 11 && queryBegin(hsp)== 216 && databaseBegin(hsp)==151 && getNumGaps(hsp)==2));
				SEQAN_ASSERT((hspcount != 18) ||(hspcount == 18 && (databaseAlignmentString(hsp) == "PGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQD") && eValue(hsp) == 0.032));
				SEQAN_ASSERT((hspcount != 20) ||(hspcount == 20 && bitScore(hsp) == 19.6f && eValue(hsp) == 2.3));
				SEQAN_ASSERT((hspcount != 23) ||(hspcount == 23 && getBlastMatchScore(hsp)== 36 && percentPositives(hsp)==71 && getQueryBegin(hsp) == 82 && getDatabaseBegin(hsp) == 103 ));
				SEQAN_ASSERT((hspcount != 28) ||(hspcount == 28 && numGaps(hsp) == 4));
				SEQAN_ASSERT((hspcount != 36) ||(hspcount == 36 && (queryAlignmentString(hsp) == "MKQANQDRGTLLLALVAGLSINGTFAALFSSIVPFSVFPIISLVLTVYCLHQRYLNRTMPVGLPGLAAACFILGVLLYSTVVRAEYPDIGSNFFPAVLSVIMVFWIGAKMRNRKQEVAE") && databaseEnd(hsp) == 119 && queryEnd(hsp)== 119));
				SEQAN_ASSERT((hspcount != 41) ||(hspcount == 41 && percentIdentity(hsp)== 25 && score(hsp) == 37));
				SEQAN_ASSERT((hspcount != 50) ||(hspcount == 50 && queryFrame(hsp) == 1 && getEValue(hsp)== 0.15 && percentGaps(hsp)==14 && score(hsp) == 48));
				if(hspcount == 8 || hspcount == 19)
				{
					TBlastHsp hsp2;
				    hsp2 = hsp;
					SEQAN_ASSERT(percentPositives(hsp) == getPercentPositives(hsp2));
					SEQAN_ASSERT(databaseFrame(hsp) == getDatabaseFrame(hsp2));
					SEQAN_ASSERT(getQueryFrame(hsp) == queryFrame(hsp2));
					SEQAN_ASSERT(score(hsp) == getBlastMatchScore(hsp2));
					SEQAN_ASSERT(bitScore(hsp) == getBitScore(hsp2));
					SEQAN_ASSERT(getEValue(hsp) == eValue(hsp2));
					SEQAN_ASSERT(percentIdentity(hsp) == getPercentIdentity(hsp2));
					SEQAN_ASSERT(percentGaps(hsp) == getPercentGaps(hsp2));
					SEQAN_ASSERT(getNumGaps(hsp) == numGaps(hsp2));
					SEQAN_ASSERT(queryOrientationPlus(hsp) == queryOrientationPlus(hsp2));
					SEQAN_ASSERT(databaseOrientationPlus(hsp) == databaseOrientationPlus(hsp2));
					SEQAN_ASSERT(getQueryBegin(hsp) == queryBegin(hsp2));
					SEQAN_ASSERT(databaseBegin(hsp) == getDatabaseBegin(hsp2));
					SEQAN_ASSERT(getQueryEnd(hsp) == queryEnd(hsp2));
					SEQAN_ASSERT(databaseEnd(hsp) == getDatabaseEnd(hsp2));
					SEQAN_ASSERT(getQueryAlignmentString(hsp) == queryAlignmentString(hsp2));
					SEQAN_ASSERT(databaseAlignmentString(hsp) == getDatabaseAlignmentString(hsp2));
				}
				
			}
		}

		++repcount;
	}
	TBlastReport blast2(blast);
	SEQAN_ASSERT(getDatabaseName(blast) == getDatabaseName(blast2));


	SEQAN_ASSERT(hitcount==56);
	SEQAN_ASSERT(hspcount==58);
	SEQAN_ASSERT(alicount==15);
	SEQAN_ASSERT(repcount==16);


}



template<typename T>
void Test_BlastParsingBasic(BlastN) 
{
//IOREV

	typedef BlastHsp<BlastN, BasicInfo> TBlastHsp;
	typedef ::std::fstream TFile;
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;
	//typedef typename Size<TBlastReport>::Type TSize;

	typedef StringSet<String<Dna>,Dependent<> > TStringSet;
 	typedef Graph<Alignment<TStringSet> > TAliGraph;

	char buffer[1024];
	strcpy(buffer, SEQAN_PATH_TO_ROOT());
	strcat(buffer, "/extras/tests/blast/ecoln.out");
	
	std::fstream strm2;
	strm2.open(buffer, ios_base::in | ios_base::binary);
		
	int alicount = 0;
	int hitcount = 0;
	int hspcount = 0;
	TBlastReport blast2;
	while(!atEnd(strm2,blast2)) 
	{
		read(strm2,blast2,Blast());
		SEQAN_ASSERT(databaseName(blast2) == "ecoli.nt ");
		SEQAN_ASSERT(queryName(blast2) == "Test");
		THitIterator hit_it2(blast2); 
		for(; !atEnd(strm2,hit_it2); goNext(strm2,hit_it2)) 
		{
			THitIterator hit_it;
			hit_it = hit_it2;
			goNext(strm2,hit_it2);
			SEQAN_ASSERT(hit_it != hit_it2 || (atEnd(strm2,hit_it2)&&hit_it == hit_it2));
			hit_it2=hit_it;
			SEQAN_ASSERT(hit_it == hit_it2);
			++hitcount;
			TBlastHit hit = value(strm2,hit_it);
			SEQAN_ASSERT((hitcount!=8) || (hitcount==8 && name(hit) == "gb|AE000167.1|AE000167 Escherichia coli K-12 MG1655 section 57 of 400 of the complete genome"));
			THspIterator hsp_it2(hit);
			SEQAN_ASSERT(atBegin(strm2,hsp_it2));
			for(; !atEnd(strm2,hsp_it2); goNext(strm2,hsp_it2)) 
			{
				THspIterator hsp_it;
				hsp_it = hsp_it2;
				SEQAN_ASSERT(!(hsp_it!=hsp_it2));
				++hspcount;
 				TBlastHsp hsp = value(strm2,hsp_it);
				if(eValue(hsp) < 0.02)
				{
						Align< String<Dna>, ArrayGaps> ali;
						getAlignment(hsp,ali);
						TAliGraph ali_g;
						getAlignment(hsp,ali_g); //hit ID
						++alicount;
				}
				SEQAN_ASSERT((hspcount != 7) ||(hspcount == 7 && eValue(hsp)== 0.52));
				SEQAN_ASSERT((hspcount != 12) ||(hspcount == 12 && getQueryBegin(hsp)==472 && getQueryEnd(hsp)==486 && getQueryAlignmentString(hsp) == "cacctggtggcgatg" && getDatabaseAlignmentString(hsp) == "cacctggtggcgatg"));
				SEQAN_ASSERT((hspcount != 17) ||(hspcount == 17 && queryBegin(hsp)==476 && queryEnd(hsp)==489));
				SEQAN_ASSERT((hspcount != 20) ||(hspcount == 20 && eValue(hsp)== 8.1 && getDatabaseBegin(hsp)==2092 && getDatabaseEnd(hsp)==2079));
				SEQAN_ASSERT((hspcount != 34) ||(hspcount == 34 && databaseBegin(hsp)==5787 && databaseEnd(hsp)==5770 ));
				SEQAN_ASSERT(&hostReport(hit_it)==&blast2);
				SEQAN_ASSERT(&hostHit(hsp_it)==&hit);
				SEQAN_ASSERT(&hostHit(hsp_it2)==&hostHit(hsp_it));
//				SEQAN_ASSERT(&hostReport(hsp_it)==&blast2);
//				SEQAN_ASSERT(&hostReport(hsp_it2)==&blast2);
				if(hspcount == 8 || hspcount == 19)
				{
					TBlastHsp hsp2;
				    hsp2 = hsp;
					SEQAN_ASSERT(getEValue(hsp) == eValue(hsp2));
					SEQAN_ASSERT(getQueryBegin(hsp) == queryBegin(hsp2));
					SEQAN_ASSERT(databaseBegin(hsp) == getDatabaseBegin(hsp2));
					SEQAN_ASSERT(getQueryEnd(hsp) == queryEnd(hsp2));
					SEQAN_ASSERT(databaseEnd(hsp) == getDatabaseEnd(hsp2));
					SEQAN_ASSERT(getQueryAlignmentString(hsp) == queryAlignmentString(hsp2));
					SEQAN_ASSERT(databaseAlignmentString(hsp) == getDatabaseAlignmentString(hsp2));
				}

			}
		}
	}



	SEQAN_ASSERT(hitcount==28);
	SEQAN_ASSERT(hspcount==34);
	SEQAN_ASSERT(alicount==5);

}


//////////////////////////////////////////////////////////////////////////////
template<typename T>
void Test_BlastParsingBasic(BlastP) {
//IOREV


	typedef BlastHsp<BlastP, BasicInfo> TBlastHsp;
	typedef ::std::fstream TFile;
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;
	//typedef typename Size<TBlastReport>::Type TSize;

	typedef StringSet<String<AminoAcid>,Dependent<> > TStringSet;
 	typedef Graph<Alignment<TStringSet> > TAliGraph;

	char buffer[1024];
	strcpy(buffer, SEQAN_PATH_TO_ROOT());
	strcat(buffer, "/extras/tests/blast/ecolip.out");
	
	std::fstream strm;
	strm.open(buffer, ios_base::in | ios_base::binary);


	int alicount = 0;
	int hitcount = 0;
	int hspcount = 0;
	int repcount = 0;

	TBlastReport blast;
	while(!atEnd(strm,blast)) 
	{
		read(strm,blast,Blast());
		SEQAN_ASSERT(getDatabaseName(blast) == "ecoliKurz.aa ");
		SEQAN_ASSERT(repcount!= 0 || (repcount==0 && getQueryName(blast) == "gi|1786182|gb|AAC73112.1| (AE000111) thr operon leader peptide"));
		THitIterator hit_it2(blast);
		
		SEQAN_ASSERT(atBegin(strm,hit_it2));
		goNext(strm,hit_it2);
		SEQAN_ASSERT(!atBegin(strm,hit_it2) || atEnd(strm,hit_it2));
		goBegin(strm,hit_it2);
		SEQAN_ASSERT(atBegin(strm,hit_it2));

		for(; !atEnd(strm,hit_it2); goNext(strm,hit_it2)) 
		{		
			TBlastHit hit = getValue(strm,hit_it2);
			THitIterator hit_it;
			hit_it = hit_it2; 
			SEQAN_ASSERT(hit_it==hit_it2);
			++hitcount;
			SEQAN_ASSERT(hitcount!=3 || (hitcount==3 && "gb|AAC76950.1| (AE000471) UDP-N-acetylenolpyruvoylglucosamine reductase"));
			THspIterator hsp_it2(hit);
			SEQAN_ASSERT(atBegin(strm,hsp_it2));
			goNext(strm,hsp_it2);
			SEQAN_ASSERT(!atBegin(strm,hsp_it2) || atEnd(strm,hsp_it2));
			goBegin(strm,hsp_it2);
			SEQAN_ASSERT(atBegin(strm,hsp_it2));
			
			for(; !atEnd(strm,hsp_it2); goNext(strm,hsp_it2)) 
			{
				++hspcount;
				TBlastHsp hsp = getValue(strm,hsp_it2);

				THspIterator hsp_it;
				hsp_it = hsp_it2;
				SEQAN_ASSERT(hsp_it==hsp_it2);
 				if(hspcount == 80 )
				{
					std::cout << "  HspQueryBegin: "<<getQueryBegin(hsp)<<"\n";
 				    std::cout << "  HspQueryEnd  : "<<getQueryEnd(hsp)<<"\n";
 					std::cout << "  HspDataBBegin: "<<getDatabaseBegin(hsp)<<"\n";
 					std::cout << "  HspDataBEnd  : "<<getDatabaseEnd(hsp)<<"\n";
 					std::cout << "        Expect : "<<getEValue(hsp)<<"\n";
 					std::cout << "QueryAliString : "<<queryAlignmentString(hsp) << "\n\n";
					std::cout << "DataBAliString : "<<databaseAlignmentString(hsp) << "\n\n";
				}
				if(eValue(hsp) < 0.02)
				{
					if(hspcount==12)
					{
						String<AminoAcid> s1 = "MKKMQSIXXXXXXXXXXXXXXQAAEITLVPSVKLQIGDRDNRGYYXXXXXXXXXXXXKQHYEWRGNRW";
						String<AminoAcid> s2 = "MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDHGWWKQHYEWRGNRW";
						Align< String<AminoAcid>, ArrayGaps> ali;
						resize(rows(ali), 2);
						setSource(row(ali, 0), s1);
						setSource(row(ali, 1), s2);
						getAlignment(hsp,ali,KnownSource());
						//std::cout << ali;
						StringSet<String<AminoAcid>,Dependent<> > ali_str;
						assignValueById(ali_str,s1,0);
						assignValueById(ali_str,s2,1);
						TAliGraph ali_g(ali_str);
						getAlignment(hsp,ali_g,0,1); //hit ID
						//std::cout << ali_g;
						++alicount;
					}
					else
					{
					Align< String<Dna>, ArrayGaps> ali;
					getAlignment(hsp,ali,UnknownSource());
					TAliGraph ali_g;
					getAlignment(hsp,ali_g); //hit ID
					++alicount;
					}
				}
				SEQAN_ASSERT((hspcount != 1) ||(hspcount == 1 && getDatabaseAlignmentString(hsp)=="MKRISTTITTTITITTGNGAG"));
				SEQAN_ASSERT((hspcount != 3) ||(hspcount == 3 && getQueryAlignmentString(hsp)=="LSYFGAKVLHPRTITPIAQFQIPCLIKNTGNP"));
				SEQAN_ASSERT((hspcount != 11) ||(hspcount == 11 && queryBegin(hsp)== 216 && databaseBegin(hsp)==151 ));
				SEQAN_ASSERT((hspcount != 18) ||(hspcount == 18 && (databaseAlignmentString(hsp) == "PGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQD") && eValue(hsp) == 0.032));
				SEQAN_ASSERT((hspcount != 20) ||(hspcount == 20 && eValue(hsp) == 2.3));
				SEQAN_ASSERT((hspcount != 23) ||(hspcount == 23 && getQueryBegin(hsp) == 82 && getDatabaseBegin(hsp) == 103 ));
				SEQAN_ASSERT((hspcount != 36) ||(hspcount == 36 && (queryAlignmentString(hsp) == "MKQANQDRGTLLLALVAGLSINGTFAALFSSIVPFSVFPIISLVLTVYCLHQRYLNRTMPVGLPGLAAACFILGVLLYSTVVRAEYPDIGSNFFPAVLSVIMVFWIGAKMRNRKQEVAE") && databaseEnd(hsp) == 119 && queryEnd(hsp)== 119));
				SEQAN_ASSERT((hspcount != 50) ||(hspcount == 50 && getEValue(hsp)== 0.15));
				if(hspcount == 8 || hspcount == 19)
				{
					TBlastHsp hsp2;
				    hsp2 = hsp;
					SEQAN_ASSERT(getEValue(hsp) == eValue(hsp2));
					SEQAN_ASSERT(getQueryBegin(hsp) == queryBegin(hsp2));
					SEQAN_ASSERT(databaseBegin(hsp) == getDatabaseBegin(hsp2));
					SEQAN_ASSERT(getQueryEnd(hsp) == queryEnd(hsp2));
					SEQAN_ASSERT(databaseEnd(hsp) == getDatabaseEnd(hsp2));
					SEQAN_ASSERT(getQueryAlignmentString(hsp) == queryAlignmentString(hsp2));
					SEQAN_ASSERT(databaseAlignmentString(hsp) == getDatabaseAlignmentString(hsp2));
				}
				
			}
		}

		++repcount;
	}
	TBlastReport blast2(blast);
	SEQAN_ASSERT(getDatabaseName(blast) == getDatabaseName(blast2));


	SEQAN_ASSERT(hitcount==56);
	SEQAN_ASSERT(hspcount==58);
	SEQAN_ASSERT(alicount==15);
	SEQAN_ASSERT(repcount==16);

}


	
SEQAN_DEFINE_TEST(test_blast_store_report_int)
{
//IOREV
	Test_BlastStoreReport<int>();
}
	
SEQAN_DEFINE_TEST(test_blast_store_report_basic_int)
{
//IOREV
	Test_BlastStoreReportBasic<int>();
}

SEQAN_DEFINE_TEST(test_blast_parsing_int_blastn)
{
//IOREV
	Test_BlastParsing<int>(BlastN());
}

SEQAN_DEFINE_TEST(test_blast_parsing_int_blastp)
{
//IOREV
	Test_BlastParsing<int>(BlastP());
}


SEQAN_DEFINE_TEST(test_blast_parsing_basic_int_blastn)
{
//IOREV
	Test_BlastParsingBasic<int>(BlastN());
}

SEQAN_DEFINE_TEST(test_blast_parsing_basic_int_blastp)
{
//IOREV
	Test_BlastParsingBasic<int>(BlastP());
}





#endif

