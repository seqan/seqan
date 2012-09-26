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

#ifndef SEQAN_HEADER_BLAST_RUN_H
#define SEQAN_HEADER_BLAST_RUN_H

//SEQAN_NO_DDDOC: do not generate documentation for this file

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Different Blast Program Types
//////////////////////////////////////////////////////////////////////////////

// DNA alignments
// blastn compares a nucleotide query sequence against a nucleotide sequence database
struct TagRunBlastN_;
typedef Tag<TagRunBlastN_> const RunBlastN;
// megablast
struct TagRunMegaBlast_;
typedef Tag<TagRunMegaBlast_> const RunMegaBlast;
//blat
struct TagRunBlat_;
typedef Tag<TagRunBlat_> const RunBlat;

// Protein alignments
// blastp compares an amino acid query sequence against a protein sequence database
struct TagRunBlastP_;
typedef Tag<TagRunBlastP_> const RunBlastP;
// blastx compares a nucleotide query sequence translated in all reading frames against a protein sequence database
struct TagRunBlastX_;
typedef Tag<TagRunBlastX_> const RunBlastX;
// tblastn compares a protein query sequence against a nucleotide sequence database dynamically translated in all reading frames
struct TagRunTBlastN_;
typedef Tag<TagRunTBlastN_> const RunTBlastN;
// tblastx compares the six-frame translations of a nucleotide query sequence against the six-frame translations of a nucleotide sequence database. Please note that tblastx program cannot be used with the nr database on the BLAST Web page.
struct TagRunTBlastX_;
typedef Tag<TagRunTBlastX_> const RunTBlastX;

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
// Blast Program Calls
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//blastn
template<typename TString, typename TPath1, typename TPath, typename TParamString>
void
_runBlast(TPath1 blast_path,
	TPath db_path,
	TString db_name,
	TPath query_path,
	TString query_name,
	TPath out_path,
	TString outfile_name,
	Tag<TagRunBlastN_>,
	TParamString params)
{
SEQAN_CHECKPOINT


	std::stringstream blastcall;
	blastcall << blast_path << "blastall"; 
	blastcall << " -p blastn"; 
	blastcall << " -d " << db_path << db_name; 
	blastcall << " -i " << query_path << query_name; 
	blastcall << " -o " << out_path << outfile_name; 
	blastcall << " " << params; 
	//process params?

	std::cout <<"\n"<< blastcall.str()<<"\n";
	system((blastcall.str()).c_str());	


}
	
//run without saving the blast report file (saves temporarily)
template<typename TBlast, typename TBlastReport, typename TString, typename TParamString>
void
run(Tag<TBlast> tag,
	TString blast_path,
	TString db_path,
	TString db_name,
	TString query_path,
	TString query_name,
	TBlastReport & blastObj, //BlastReport<BlastHsp<BlastN,THspInfoSpec>,TStoreSpec,TInfoSpec >
	TParamString params)
{
SEQAN_CHECKPOINT

#ifdef PLATFORM_WINDOWS
	String<char> out_path(blast_path);
	String<char> out_path2("D:\\emde\\blast\\bin\\data\\");
#else
	String<char> out_path(blast_path);
#endif

	String<char> outfile_name = "tempseqanblast.out";

	_runBlast(blast_path,db_path,db_name,query_path,query_name,out_path,outfile_name,tag,params);

	std::fstream strm;
	std::stringstream s;
	s <<out_path<< outfile_name; 

	strm.open((s.str()).c_str(),::std::ios_base::in);
	read(strm,blastObj, Blast());

}
	
//run without saving the blast report file (saves temporarily)
template<typename TBlast, typename TBlastReport, typename TString, typename TParamString>
void
run(Tag<TBlast> tag,
	TString blast_path,
	TString db_path,
	TString db_name,
	TString query_path,
	TString query_name,
	TBlastReport & blastObj)
{
	SEQAN_CHECKPOINT;
	run(tag, blast_path, db_path, db_name, query_path, query_name, blastObj);
}



//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//megablast
template<typename TString, typename TPath1, typename TPath, typename TParamString>
void
_runBlast(TPath1 blast_path,
	TPath db_path,
	TString db_name,
	TPath query_path,
	TString query_name,
	TPath out_path,
	TString outfile_name,
	Tag<TagRunMegaBlast_>,
	TParamString params)
{
SEQAN_CHECKPOINT

	std::stringstream blastcall;
	blastcall << blast_path << "megablast"; 
	blastcall << " -d " << db_path << db_name; 
	blastcall << " -i " << query_path << query_name; 
	blastcall << " -o " << out_path << outfile_name; 
	blastcall << " " << params; 

	//process params

	system((blastcall.str()).c_str());	


}



//////////////////////////////////////////////////////////////////////////////
//blat
template<typename TString, typename TPath1, typename TPath, typename TParamString>
void
_runBlat(TPath1 blat_path,
	TPath db_path,
	TString db_name,
	TPath query_path,
	TString query_name,
	TPath out_path,
	TString outfile_name,
	Tag<TagRunBlat_>,
	TParamString params)
{
SEQAN_CHECKPOINT

	std::stringstream blatcall;
	blatcall << blat_path << "blat "; 
	blatcall << db_path << db_name<<" "; 
	blatcall << query_path << query_name<<" "; 
	blatcall << " " << params; 
	blatcall << " -out=blast " << out_path << outfile_name; 
	
	//process params

	std::cout << blatcall.str() << "\n";
	system((blatcall.str()).c_str());	


}



//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//blastp
template<typename TString, typename TPath1, typename TPath, typename TParamString>
void
_runBlast(TPath1 blast_path,
	TPath db_path,
	TString db_name,
	TPath query_path,
	TString query_name,
	TPath out_path,
	TString outfile_name,
	Tag<TagRunBlastP_>,
	TParamString params)
{
SEQAN_CHECKPOINT


	std::stringstream blastcall;
	blastcall << blast_path << "blastall"; 
	blastcall << " -p blastp"; 
	blastcall << " -d " << db_path << db_name; 
	blastcall << " -i " << query_path << query_name; 
	blastcall << " -o " << out_path << outfile_name; 
	blastcall << " " << params; 

	//process params

	system((blastcall.str()).c_str());	


}



//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//tblastn
template<typename TString, typename TPath1, typename TPath, typename TParamString>
void
_runBlast(TPath1 blast_path,
	TPath db_path,
	TString db_name,
	TPath query_path,
	TString query_name,
	TPath out_path,
	TString outfile_name,
	Tag<TagRunTBlastN_>,
	TParamString params)
{
SEQAN_CHECKPOINT


	std::stringstream blastcall;
	blastcall << blast_path << "blastall"; 
	blastcall << " -p tblastn"; 
	blastcall << " -d " << db_path << db_name; 
	blastcall << " -i " << query_path << query_name; 
	blastcall << " -o " << out_path << outfile_name; 
	blastcall << " " << params; 

	//process params

	system((blastcall.str()).c_str());	


}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//tblastx
template<typename TString, typename TPath1, typename TPath, typename TParamString>
void
_runBlast(TPath1 blast_path,
	TPath db_path,
	TString db_name,
	TPath query_path,
	TString query_name,
	TPath out_path,
	TString outfile_name,
	Tag<TagRunTBlastX_>,
	TParamString params)
{
SEQAN_CHECKPOINT


	std::stringstream blastcall;
	blastcall << blast_path << "blastall"; 
	blastcall << " -p tblastx"; 
	blastcall << " -d " << db_path << db_name; 
	blastcall << " -i " << query_path << query_name; 
	blastcall << " -o " << out_path << outfile_name; 
	blastcall << " " << params; 
	//process params

	system((blastcall.str()).c_str());	


}



//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
//blastx
template<typename TString, typename TPath1, typename TPath, typename TParamString>
void
_runBlast(TPath1 blast_path,
	TPath db_path,
	TString db_name,
	TPath query_path,
	TString query_name,
	TPath out_path,
	TString outfile_name,
	Tag<TagRunBlastX_>,
	TParamString params)
{
SEQAN_CHECKPOINT


	std::stringstream blastcall;
	blastcall << blast_path << "blastall"; 
	blastcall << " -p blastx"; 
	blastcall << " -d " << db_path << db_name; 
	blastcall << " -i " << query_path << query_name; 
	blastcall << " -o " << out_path << outfile_name; 
	blastcall << " " << params; 

	//process params

	system((blastcall.str()).c_str());	


}






//////////////////////////////////////////////////////////////////////////////
//formatdb
template<typename TStr,typename TString, typename TParamString>
void
_runFormatDB(TStr blast_path,
			 TStr db_path,
			 TString db_file,
			 TParamString params)
{
SEQAN_CHECKPOINT


	std::stringstream blastcall;
	blastcall << blast_path << "formatdb "; 
	blastcall << " -i " <<db_path << db_file; 
	blastcall << " " << params; 

	//process params

	system((blastcall.str()).c_str());	


}











}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
