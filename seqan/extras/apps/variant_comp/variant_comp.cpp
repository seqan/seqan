#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>


#include "seqan/platform.h"
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/refinement.h>

#include "variant_comp.h"

using namespace std;
using namespace seqan;




template<typename TChar>
inline bool
parseIsDigit(TChar const c)
{
	//return (((int ) c >=  48) && ((int) c <=  57));
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
		    (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}


template<typename TFile, typename TChar>
inline long double
parseReadDouble(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!parseIsDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

template<typename TFile, typename TChar>
inline void 
parseSkipWhitespace(TFile& file, TChar& c)
{
	if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) break;
	}
}

template<typename TFile, typename TChar, typename TString>
void
_parseReadWordUntilWhitespace(TFile& file, TString& str, TChar& c)
{
        append(str,c);
        if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
                c = _streamGet(file);
                return;
        }
        while (!_streamEOF(file)) {
                c = _streamGet(file);
                if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
                append(str, c);
        }
        return;
}   


template <typename TOptions>
bool loadGenomes(const char* fileName, 
		StringSet<Dna5String> &genomes,
		String<CharString> & genomeIDs, 
		::std::map<CharString,unsigned> &gIdStringToIdNumMap,
		TOptions &options)
{
	
	MultiFasta multiFasta;
	if (!open(multiFasta.concat,fileName,OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
	resize(genomeIDs,seqCount);
	if(options.sequenceContext || options.annotateRepeats > 0) 
		resize(genomes,seqCount);
	
	for(unsigned i = 0; i < seqCount; ++i)
	{
		CharString temp;
		if(options.sequenceContext || options.annotateRepeats > 0) 
			assignSeq(genomes[i], multiFasta[i], Fasta());
		assignSeqId(temp, multiFasta[i], Fasta());
		for (unsigned pos = 0; pos < length(temp); ++pos)
		{
			if(temp[pos]=='\t' || temp[pos]=='\b' || temp[pos]==' ')
			{
				resize(temp,pos);
				break;
			}
		}
		genomeIDs[i] = temp;
		gIdStringToIdNumMap.insert(::std::make_pair(temp,i)); 
	}
	return (seqCount > 0);
}

inline int
_parseReadNumber(CharString & file, unsigned & c)
{
	// Read number
	String<char> str(file[c]);
	while (c < length(file)-1) {
		++c;
		if (!_parseIsDigit(file[c])) break;
		append(str, file[c]);
	}
	return atoi(toCString(str));
}



inline void
_parseSkipWhitespace(CharString & file, unsigned & c)
{
	if ((unsigned) file[c] > 32) return;
	while (c < length(file)-1) {
		++c;
		if ((unsigned) file[c] > 32) break;
	}
}



inline void
_parseReadIdentifier(CharString & file, CharString & str, unsigned & c)
{
	// Read identifier
	append(str, file[c], Generous());
	while (c < length(file)-1) {
		++c;
		if (!_parseIsAlphanumericChar(file[c])) break;
		append(str, file[c], Generous());
	}
}

inline void
_parseReadWordUntilSemicolon(CharString & file, CharString& str, unsigned & c)
{
        clear(str);
        if (c == length(file)) {
                return;
        }
        while (c < length(file) && _parseIsAlphanumericChar(file[c]) && file[c] != ';'){//!(file[c]== ' ' || file[c] == '\t' || file[c] == ';'|| file[c] == '\r'|| file[c] == '\n')) {
                append(str, file[c]);
                ++c;
        }
        return;
}   

// parse last column of GFF file, containing id, tags,...
// add genotype information
template <typename TIndel, typename TMap, typename TOptions>
void
getInfoFromNinethCol(CharString &ninethCol, TIndel &indel, TMap & gIdStringToIdNumMap, TOptions &options )
{
	unsigned c = 0;
	CharString temp_str;
	_parseReadIdentifier(ninethCol,temp_str,c);
	if(options._debugLevel > 1)
		::std::cout << temp_str << "\t";
//      if(!(temp_str == "ID")) ::std::cout << "first feature field should be 'ID' but is " << temp_str<<::std::endl;

	// skip the "="
	++c;
	// read the ID
	clear(temp_str);
	clear(indel.idStr);
	CharString indelID;
	_parseReadIdentifier(ninethCol,indel.idStr,c);
	if(options._debugLevel > 1)
		::std::cout << "myID = "<< indel.idStr << "\n";

	// process tags in a loop
	CharString current_tag;
	clear(indel.insertionSeq);
	_parseSkipWhitespace(ninethCol,c);
	while(c < length(ninethCol))
	{
		// different tags are separated by ';'  
		while(ninethCol[c] != ';')
		{
			if(c == length(ninethCol)-1) // end of line
				break;
			++c;
		}
		if(c == length(ninethCol)-1) // end of line
			break;

		// get the current tag
		clear(current_tag);
		++c;
		_parseReadIdentifier(ninethCol,current_tag,c);
		if(options._debugLevel > 1)
			::std::cout << current_tag << " in features\n";
		if(current_tag=="size")
		{
			++c;
			indel.indelSize = _parseReadNumber(ninethCol,c);
			if(options._debugLevel > 1)
				::std::cout << indel.indelSize << " indel size\n";
            if(indel.type == INSERTION && indel.indelSize > 0)
                indel.indelSize = - indel.indelSize;
		}
		else
		{
			if(current_tag=="duplication") indel.duplication = true;
			else if(current_tag=="seq")
			{
				++c;
				while(c < length(ninethCol) && ninethCol[c] != ';')
				{
					appendValue(indel.insertionSeq,(Dna5)ninethCol[c]);
					++c;
				}
			}
		    if(current_tag=="geno") 
		    {
			    ++c;
      		    CharString temp_str;
                _parseReadWordUntilSemicolon(ninethCol,temp_str,c); 
                indel.genotype=temp_str;
       		    if(prefix(temp_str,3)=="hom")
	                indel.genotype=0; //hom
                else indel.genotype=1; //het
    	    }
			if(current_tag=="endChr") 
			{
				++c;
           		CharString temp_str;
		        _parseReadWordUntilSemicolon(ninethCol,temp_str,c); 
        		if(prefix(temp_str,3)=="chr")
		        	temp_str = suffix(temp_str,3);
		
        		//check if the genomeID is in our map of relevant genomeIDs, otherwise use fakeID
		        typename TMap::iterator it = gIdStringToIdNumMap.find(temp_str);
		        if(it != gIdStringToIdNumMap.end()) indel.secondGenomeId = it->second;
                else indel.secondGenomeId = 10000; //disable

			}
			if(current_tag=="endPos") 
			{
				++c;
		    	int temp = _parseReadNumber(ninethCol,c);
	    	    if(indel.type == TRANSLOCATION) 
                {
                    if(indel.genomeId == indel.secondGenomeId && abs(temp-(int)indel.originalPos) < MAX_DUPLICATION_SIZE)// duplication?
                    {
                        indel.type = INSERTION;
                        indel.duplication = true;
                        indel.indelSize = -abs(temp-(int)indel.originalPos);
                    }
                    else
                        indel.secondOriginalPos = temp;
                }
                if(indel.type == INVERSION)
                {
                    indel.indelSize = abs(temp-(int)indel.originalPos);
                    indel.originalPos = _min(temp,(int)indel.originalPos);
                }
			}
		}
	}
}


// parse last column of GFF file, containing id, tags,... for SNPs
template <typename TSnp, typename TMap, typename TOptions>
void
getSnpInfoFromNinethCol(CharString &ninethCol, TSnp &snp, TMap & , TOptions &options )
{
	unsigned c = 0;
	CharString temp_str;
	_parseReadIdentifier(ninethCol,temp_str,c);
	if(options._debugLevel > 1)
		::std::cout << temp_str << "\t";
//      if(!(temp_str == "ID")) ::std::cout << "first feature field should be 'ID' but is " << temp_str<<::std::endl;

	// skip the "="
	++c;
	// read the ID
	clear(temp_str);
	clear(snp.idStr);
	CharString snpID;
	_parseReadIdentifier(ninethCol,snp.idStr,c);
	if(options._debugLevel > 1)
		::std::cout << "myID = "<< snp.idStr << "\n";

	// process tags in a loop
	CharString current_tag;
	_parseSkipWhitespace(ninethCol,c);
	while(c < length(ninethCol))
	{
		// different tags are separated by ';'  
		while(ninethCol[c] != ';')
		{
			if(c == length(ninethCol)-1) // end of line
				break;
			++c;
		}
		if(c == length(ninethCol)-1) // end of line
			break;

		// get the current tag
		clear(current_tag);
		++c;
		_parseReadIdentifier(ninethCol,current_tag,c);
		if(options._debugLevel > 1)
		    ::std::cout << current_tag << " in features\n";
		if(current_tag=="geno") 
		{
			++c;
      		CharString temp_str;
            _parseReadWordUntilSemicolon(ninethCol,temp_str,c); 
            snp.genotype=temp_str;
       		if(prefix(temp_str,3)=="hom")
	            snp.genotype=0; //hom
            else snp.genotype=1; //het
    	}
    }
	return;
}

template<typename TSnp, typename TFile, typename TChar, typename TMap, typename TOptions>
void
_parseSnp(TSnp & snp, TFile & file, TChar & c, TMap & gIdStringToIdNumMap, TOptions & options)
{
    typedef int	TContigPos;
	
    // skip whitespaces and read entry in column 4  --> genomic begin position
	_parseSkipWhitespace(file, c);
	snp.originalPos = (TContigPos) _parseReadNumber(file,c) - 1;
	if(options._debugLevel > 1) 
		::std::cout << snp.originalPos << "\t";
		
	// skip whitespaces and read entry in column 5  --> genomic end position // not needed here
	_parseSkipWhitespace(file, c);
    snp.simPos = _parseReadNumber(file,c) - 1; // this is not the sim pos, we are just using this field to avoid unused variable warnings
	SEQAN_ASSERT_EQ(snp.simPos,snp.originalPos); // make sure begin and end pos are the same, this is a SNP!

    // skip whitespaces and read entry in column 6  --> score (percent identity or mapping quality) or a '.'
	int readSupport = 1000; //  --> no information about read support (reference indel)
	_parseSkipWhitespace(file, c);
	if(c=='.')
		c = _streamGet(file);               // 
	else 
		readSupport = (TContigPos) _parseReadDouble(file,c); // number of supporting reads
			
	if(options._debugLevel > 1) 
		::std::cout << readSupport << "\t";
	
    snp.quality = readSupport;	

	// skip whitespaces and read entry in column 7  --> strand information: '+' or '-' // not needed here
	_parseSkipWhitespace(file, c);
	c = _streamGet(file);
		
	// skip whitespaces and read entry in column 8  --> always '.' here
	_parseSkipWhitespace(file, c);
	c = _streamGet(file);
		
	// skip whitespaces and read entry in column 9  --> tags, extra information. first tag is always "ID"
	_parseSkipWhitespace(file, c);
		
	//remember this position and store the whole 9th column
	//typename std::fstream::pos_type idtagstart = file.tellg();
	clear(snp.ninethCol);
	while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n'))) // while in same line
	{
		append(snp.ninethCol,c);
		c = _streamGet(file);
	}
	if(options._debugLevel > 1) std::cout << "9thcol=" <<  snp.ninethCol << "\n\n";
		
	getSnpInfoFromNinethCol(snp.ninethCol,snp,gIdStringToIdNumMap,options);
  
	_parseSkipWhitespace(file, c);

}

/////////////////////////////////////////////////////////////
// read Gff input file containing indels and SNPs
template <
	typename TIndelSet,
	typename TSnpSet,
	typename TGenomeMap,
	typename TOptions
>
int readGFF(
	const char*				&filename,
	TIndelSet 				&indelSet,
	TSnpSet 				&snpSet,
	TGenomeMap				&gIdStringToIdNumMap,
	TOptions				&options)
{
	typedef typename Value<TIndelSet>::Type	TIndel;
	typedef typename Value<TSnpSet>::Type	TSnp;
	typedef int				TId;
	typedef int				TContigPos;
	
	
	::std::ifstream file;
	file.open(filename, ::std::ios_base::in | ::std::ios_base::binary);
	if (!file.is_open()) return 1;
		
	
	clear(indelSet);
	char c = _streamGet(file);
	while (!_streamEOF(file))
	{
		TIndel indel = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		TSnp snp = {0,0,0,0,0,0,0,0,0};
		
		if(c == '#')
			_parseSkipLine(file,c);	
				
		// skip whitespaces just in case (actually there shouldnt be a whitespace at the beginning of a line)
		_parseSkipWhitespace(file, c);
	
		if(c == '#')
			_parseSkipLine(file,c);	
		// and read entry in column 1  --> genomeID
		CharString temp_str;
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(prefix(temp_str,3)=="chr")
			temp_str = suffix(temp_str,3);
		
		TId contigId;
		//check if the genomeID is in our map of relevant genomeIDs, otherwise skip match
		typename TGenomeMap::iterator it = gIdStringToIdNumMap.find(temp_str);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if(it != gIdStringToIdNumMap.end()) contigId = it->second;
		else
		{
		//	std::cout  << "No1\n";
			_parseSkipLine(file,c);
			continue;
		}
		indel.genomeId = contigId;
		snp.genomeId = contigId;

		// skip whitespaces and read entry in column 2
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		indel.field2 = temp_str;
		snp.field2 = temp_str;

		// skip whitespaces and read entry in column 3
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		indel.field3 = temp_str;
        bool skip = true;
        // skip everything that is not an insertion or deletion
        if(temp_str == "snp" || snp.field2 == "snp")
        {
            snp.field3 = temp_str;
            _parseSnp(snp, file, c, gIdStringToIdNumMap, options);
            appendValue(snpSet,snp);
            continue;
        }
        else if(indel.field3 == "insertion" || indel.field2 == "insertion")
        {
            indel.type = INSERTION;
            skip = false;
        }
        else if(indel.field3 == "deletion" || indel.field2 == "deletion")
        {
            indel.type = DELETION;
            skip = false;
        }
        else if(indel.field3 == "intron" || indel.field2 == "intron")
        {
            indel.type = DELETION;
            skip = false;
        }
        else if(indel.field3 == "inversion" || indel.field2 == "inversion")
        {
            indel.type = INVERSION;
            skip = false;
        }
        else if(indel.field3 == "translocation" || indel.field2 == "translocation")
        {
            indel.type = TRANSLOCATION;
            skip = false;
        }
        else if(indel.field3 == "duplication" || indel.field2 == "duplication")
        {
            indel.type = INSERTION;
            indel.duplication = true;
            skip = false;
        }

        if(skip)
        {
            _parseSkipLine(file,c);
            continue;
        }
		// skip whitespaces and read entry in column 4  --> genomic begin position
		_parseSkipWhitespace(file, c);
		indel.originalPos = (TContigPos) _parseReadNumber(file,c) - 1;
		if(options._debugLevel > 1) 
			::std::cout << indel.originalPos << "\t";
		
		// skip whitespaces and read entry in column 5  --> genomic end position // not needed here
		_parseSkipWhitespace(file, c);
		indel.indelSize = _parseReadNumber(file,c) - indel.originalPos;
		if(indel.indelSize > 0 && indel.type == INSERTION)
        {
            indel.indelSize = -indel.indelSize;
        }
		// skip whitespaces and read entry in column 6  --> score (percent identity or mapping quality) or a '.'
		double readSupport = 1000.0; //  --> no information about read support (reference indel)
		_parseSkipWhitespace(file, c);
		if(c=='.')
			c = _streamGet(file);               // 
		else 
			readSupport = _parseReadDouble(file,c); // number of supporting reads
			
        indel.quality = readSupport;
        snp.quality = (int)readSupport;

		if(options._debugLevel > 1) 
			::std::cout << readSupport << "\t";
		
		// skip whitespaces and read entry in column 7  --> strand information: '+' or '-' // not needed here
		_parseSkipWhitespace(file, c);
		c = _streamGet(file);
		
		// skip whitespaces and read entry in column 8  --> always '.' here
		_parseSkipWhitespace(file, c);
		c = _streamGet(file);
		
		// skip whitespaces and read entry in column 9  --> tags, extra information. first tag is always "ID"
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		
		//remember this position and store the whole 9th column
		//typename std::fstream::pos_type idtagstart = file.tellg();
		clear(indel.ninethCol);
		while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n'))) // while in same line
		{
			append(indel.ninethCol,c);
			c = _streamGet(file);
		}
		if(options._debugLevel > 1) std::cout << "9thcol=" <<  indel.ninethCol << "\n\n";
		
		getInfoFromNinethCol(indel.ninethCol,indel,gIdStringToIdNumMap,options);
    
        if(indel.type != TRANSLOCATION) 
            appendValue(indelSet,indel);

		_parseSkipWhitespace(file, c);

/*		_parseReadIdentifier(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\n";
		if(temp_str!="ID") ::std::cout << "first feature field should be 'ID'"<<::std::endl;
		
		// skip the "="
		c = _streamGet(file);
		
		// read the ID
		clear(temp_str);
		clear(indel.idStr);
		CharString indelID;
		_parseReadIdentifier(file,indel.idStr,c);
		if(options._debugLevel > 1) 
			::std::cout << "myID = "<< indel.idStr << "\n";
		
		// process tags in a loop
		CharString current_tag;
		_parseSkipWhitespace(file,c); 
		while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n'))) // while in same line
		{
			// different tags are separated by ';'  
			while(c != ';')
			{
				if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
					break;
				c = _streamGet(file);
			}
			if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
				break;
			
			// get the current tag
			clear(current_tag);
			c = _streamGet(file);
			_parseReadIdentifier(file,current_tag,c);
			if(options._debugLevel > 1) 
				::std::cout << current_tag << " in features\n";
			if(current_tag=="size")
			{
		//		if(c == '-') 
					c = _streamGet(file);
				indel.indelSize = _parseReadNumber(file,c); 
				if(options._debugLevel > 1) 
					::std::cout << indel.indelSize << " indel size\n";
			}
			else
			{
				if(current_tag=="duplication") indel.duplication = true;
				if(current_tag=="seq")
				{
					c = _streamGet(file);
					while(!_streamEOF(file) && c != ';' && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')))
					{
						 appendValue(indel.insertionSeq,(Dna5)c);
						c = _streamGet(file);
					}
				}
		//		else _parseSkipLine(file,c);
			}
		}
		appendValue(indelSet,indel);
		_parseSkipWhitespace(file, c);*/
		
	}

	file.close();
	
	if(options._debugLevel > 0) ::std::cout << "Parsed "<<length(indelSet)<<" indels." << ::std::endl;
	
	return 0;
}



//////////////////////////////////////////////////////////////////////////////
// Print usage
template<typename TOptions>
void printHelp(int, const char *[],TOptions &, bool longHelp = false) 
{
	cerr << "***********************" << endl;
	cerr << "*** compareVariants ***" << endl;
	cerr << "***********************" << endl << endl;
	cerr << "Usage: compareVariants [OPTIONS]... <SOURCE SEQUENCE FILE>" << endl;
	cerr << "\n";
	if (longHelp) {
		cerr << "  -ip,  --input-predicted FILE     \t" << "input gff file containing predicted indels" << endl;
		cerr << "  -ir,  --input-reference FILE     \t" << "input gff file containing reference indels" << endl;
		cerr << "  -o,   --output FILE              \t" << "indel output filename" << endl;
		cerr << "  -on,  --outputFN FILE            \t" << "indel output filename for false negatives (unmatched reference indels)" << endl;
		cerr << "  -os,  --outputSnp FILE           \t" << "snp output filename" << endl;
		cerr << "  -osn, --outputSnpFN FILE         \t" << "snp output filename for false negatives (unmatched reference snps) " << endl;
		cerr << "  -gta, --genotype-aware           \t" << "genotype aware (het/hom) matching of snps" << endl;
		cerr << "  -pt,  --position-tolerance NUM   \t" << "position tolerance in bp" << endl;
		cerr << "  -st,  --size-tolerance NUM       \t" << "size tolerance in percent" << endl;
		cerr << "  -sc,  --sequence-context         \t" << "switch on sequence-context mode" << endl;
		cerr << "  -at,  --attach-tag STR           \t" << "string to attach as tag to overlapped indels" << endl;
		cerr << "  -v,   --verbose                  \t" << "verbose mode" << endl;
//		cerr << "  -vv,  --very-verbose             \t" << "very verbose mode" << endl;
		cerr << "  -h,   --help                     \t" << "print this help" << endl;
	} else {
		cerr << "Try 'compareVariants --help' for more information." << endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// Main part
int main(int argc, const char *argv[])
{
	srand(time(NULL));

	unsigned fnameCount = 0;
	const char *fname[1] = {""};
	IndelCompareOptions<> options;

	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option

			if (strcmp(argv[arg], "-ip") == 0 || strcmp(argv[arg], "--input-predicted") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.inputPredicted = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-ir") == 0 || strcmp(argv[arg], "--input-reference") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.inputReference = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-pt") == 0 || strcmp(argv[arg], "--position-tolerance") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.positionTolerance;
					if (!istr.fail())
					{
						if (options.positionTolerance < 0)
							cerr << "PositionTolerance must be a positive integer value" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-st") == 0 || strcmp(argv[arg], "--size-tolerance") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.sizeTolerance;
					if (!istr.fail())
					{
						if (options.sizeTolerance < 0 || options.sizeTolerance > 100)
							cerr << "SizeTolerance must be a value between 0 and 100" << endl << endl;
						else{
							options.sizeTolerance /= 100;
							continue;
						}
					}
				}
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-r") == 0 || strcmp(argv[arg], "--ranges") == 0) {
				if (arg + 1 < argc) {
					++arg;
					fstream file;
					clear(options.ranges);
					file.open(argv[arg],ios_base::in | ios_base::binary);
					char c = _streamGet(file);
					while (!_streamEOF(file))
					{
						parseSkipWhitespace(file,c);
						int rangeBegin = static_cast<int>(parseReadDouble(file,c));
						parseSkipWhitespace(file,c);
						int rangeEnd = static_cast<int>(parseReadDouble(file,c));
						appendValue(options.ranges,Pair<int,int>(rangeBegin,rangeEnd));
					}
					file.close();
					continue;
				}
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-at") == 0 || strcmp(argv[arg], "--attach-tag") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.attachTag = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-on") == 0 || strcmp(argv[arg], "--outputFN") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.outputFN = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.output = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-osn") == 0 || strcmp(argv[arg], "--outputSnpFN") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.outputSnpFN = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-os") == 0 || strcmp(argv[arg], "--outputSnp") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.outputSnp = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-gta") == 0 || strcmp(argv[arg], "--genotype-aware") == 0) {
				options.genotypeAware = true;
				continue;
			}
			if (strcmp(argv[arg], "-sc") == 0 || strcmp(argv[arg], "--sequence-context") == 0) {
				options.sequenceContext = true;
				continue;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--very-verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 2);
				continue;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, true);
				return 0;
			}
		}
		else {
			// parse file name
			if (fnameCount == 1) {
				printHelp(argc, argv, options);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 1) {
		printHelp(argc, argv, options);
		return 0;
	}
	
	StringSet<IndelInfo>		refIndels, predictedIndels;
	StringSet<SnpInfo>		    refSnps, predictedSnps;
	
	::std::map<CharString,unsigned> gIdStringToIdNumMap;
	String<CharString> genomeIDs;
	StringSet<Dna5String> genomes;
	
	loadGenomes(fname[0],genomes,genomeIDs,gIdStringToIdNumMap,options);
	
	if (readGFF(options.inputReference, refIndels, refSnps, gIdStringToIdNumMap, options) > 0) 
	{
		cerr << "Reference variants " << options.inputReference << " can't be loaded." << endl;
		return 0;
	}
	if (readGFF(options.inputPredicted, predictedIndels, predictedSnps, gIdStringToIdNumMap, options) > 0) 
	{
		cerr << "Predicted variants " << options.inputPredicted << " can't be loaded." << endl;
		return 0;
	}
	
	if(options._debugLevel > 0 )
	{
		::std::cout << endl << "Number of reference indels: " << length(refIndels) << endl;
		::std::cout << "Number of predicted indels: " << length(predictedIndels) << endl;
		::std::cout << endl << "Number of reference SNPs: " << length(refSnps) << endl;
		::std::cout << "Number of predicted SNPs: " << length(predictedSnps) << endl << endl;
	}
	
	int result = compareIndels(refIndels,predictedIndels,genomes,genomeIDs,options);
	if(result > 0)
	{
		cerr << "Something went wrong.. Exiting..\n";
		return 1;
	}
    if(!empty(predictedSnps) || !empty(refSnps))
        result = compareSnps(refSnps,predictedSnps,genomes,genomeIDs,options);
	if(result > 0)
	{
		cerr << "Something went wrong.. Exiting..\n";
		return 1;
	}

	return 0;
}
