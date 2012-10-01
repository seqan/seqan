/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2007-2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: David Weese <david.weese@fu-berlin.de>
  ===========================================================================
  Swiss Army Knife tool... It slices and dices and makes the laundry!

  This tool allows to cut sequences and parts of sequences out of sequence
  files.  It supports all formats supported by the AutoSeqFormat class from
  SeqAn, including FASTA, FASTQ and QSeq (Illumina format).


  Usage: sak [OPTION]... <SOURCE SEQUENCE FILE>

  Main Options:
    -o,  --output FILE              set output filename (default: use stdout)
    -q,  --qual                     enable Fastq output (default: Fasta)
    -h,  --help                     print this help

  Extract Options:
    -s,  --sequence NUM             select a single sequence by index
    -sn, --sequence-name NAME       select a single sequence by name
    -ss, --sequences START END      select sequences (default: select all)
    -i,  --infix START END          extract infix
    -rc, --revcomp                  reverse complement
    -l,  --max-length               maximal number of sequence characters to
                                    write out
		-ll, --line-length              maximal characters per output line
  ===========================================================================*/

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
//#include "../library/apps/razers/mmap_fasta.h"

using namespace std;
using namespace seqan;

//____________________________________________________________________________
// Global Parameters

	int			optionSeqStart = 0;
	int			optionSeqEnd = MaxValue<int>::VALUE;
	int			optionInfStart = 0;
	int			optionInfEnd = -1;
	bool		optionRevComp = false;
	const char	*optionOutput = NULL;
	bool		optionFastQ = false;
    bool        optionSeqNameSet = false;
    CharString  optionSeqName = "";
    int         optionMaxLength = -1;
    unsigned    optionLineLength = 200;

	typedef Dna5 TAlphabet;
	typedef String<TAlphabet> TSeqString;

//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TSeqSet, typename TQuals, typename TIDs>
bool loadSeqs(TSeqSet &seqs, TQuals &quals, TIDs &ids, const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);

	int seqCount = length(multiFasta);

	if (optionSeqEnd > seqCount) optionSeqEnd = seqCount;
	if (optionSeqStart < 0) optionSeqStart = 0;

    if (optionSeqNameSet) {
        clear(seqs);
        clear(quals);
        clear(ids);
        for (unsigned i = 0; i < length(multiFasta); ++i) {
            CharString thisName;
            assignSeqId(thisName, multiFasta[i], format);
            if (optionSeqName == infix(thisName, 0, length(optionSeqName))) {
                resize(seqs, 1, Exact());
                assignSeq(seqs[0], multiFasta[i], format);
                resize(ids, 1, Exact());
                assignSeqId(ids[0], multiFasta[i], format);
                resize(quals, 1, Exact());
                resize(quals[0], length(seqs[0]), 33 + 80);
                assignQual(quals[0], multiFasta[i], format);
                return true;
            }
        }
    } else {
        if (optionSeqStart < optionSeqEnd)
        {
            resize(seqs, optionSeqEnd - optionSeqStart, Exact());
            resize(ids, optionSeqEnd - optionSeqStart, Exact());
            resize(quals, optionSeqEnd - optionSeqStart, Exact());
            for(int i = 0, j = optionSeqStart; j < optionSeqEnd; ++i, ++j)
            {
                assignSeq(seqs[i], multiFasta[j], format);		// read Genome sequence
                assignSeqId(ids[i], multiFasta[j], format);		// read Genome ids
                resize(quals[i], length(seqs[i]), 33 + 80);
                assignQual(quals[i], multiFasta[j], format);	// read qualities
            }
        }
    }

	return (seqCount > 0);
}

template < 
	typename TStream, 
	typename TId >
void dumpFastaId(
	TStream &out,
	TId &id)
{
  /*
	unsigned size = length(id);
	for (unsigned i = 0; i < size; ++i)
		if (id[i] == ' ')
		{
			size = i;
			break;
		}
	out << infix(id, 0, size) << endl;
	*/
	out << id << endl;
}

template < 
	typename TStream, 
	typename TSeq >
void dumpFastaSeq(
	TStream &out,
    int n,
	TSeq const &seq)
{
    if (n == -1) {
        out << seq << endl;
    } else {
        SEQAN_ASSERT_GEQ(n, 0);
        unsigned size = _min(length(seq), static_cast<unsigned>(n));
        //n -= size;
        out << infix(seq, 0, size) << endl;
    }
}


template < 
	typename TReadSet, 
	typename TReadIDs >
void saveFasta(
	TReadSet const &readSet,		// generated read sequences
	TReadIDs const &readIDs)		// corresponding Fasta ids
{
	ostream *out = &cout;
	ofstream file;
	
	if (optionOutput != NULL)
	{
		file.open(optionOutput, ios_base::out | ios_base::trunc);
		if (!file.is_open()) {
			cerr << "Failed to open output file" << endl;
			return;
		}
		else
			cout << "Writing sequences to " << optionOutput << "\n";
		out = &file;
	}

	unsigned reads = length(readSet);
	for(unsigned i = 0; i < reads && optionMaxLength != 0; ++i)
	{
		(*out) << '>';
		dumpFastaId(*out, readIDs[i]);
		unsigned len = length(readSet[i]);
		if (optionMaxLength != -1 && len > (unsigned)optionMaxLength)
			len = optionMaxLength;
		for(unsigned j = 0; j < len; j += optionLineLength)
			dumpFastaSeq(*out, min(optionLineLength, len - j), suffix(readSet[i],j));
	}
	
	file.close();
}

template < 
	typename TReadSet,
	typename TQualSet,
	typename TReadIDs >
void saveFastq(
	TReadSet const &readSet,		// generated read sequences
	TQualSet const &qualSet,		// qualities
	TReadIDs const &readIDs)		// corresponding Fasta ids
{
	ostream *out = &cout;
	ofstream file;
	
	if (optionOutput != NULL)
	{
		file.open(optionOutput, ios_base::out | ios_base::trunc);
		if (!file.is_open()) {
			cerr << "Failed to open output file" << endl;
			return;
		}
		else
			cout << "Writing reads to " << optionOutput << "\n";
		out = &file;
	}

	unsigned reads = length(readSet);
	for(unsigned i = 0; i < reads && optionMaxLength != 0; ++i)
	{
        int tmp = optionMaxLength;
		(*out) << '@';
		dumpFastaId(*out, readIDs[i]);
		dumpFastaSeq(*out, optionMaxLength, readSet[i]);
		(*out) << '+';
		dumpFastaId(*out, readIDs[i]);
		dumpFastaSeq(*out, tmp, qualSet[i]);
	}
	
	file.close();
}

//////////////////////////////////////////////////////////////////////////////
// Print usage
void printHelp(int, const char *[], bool longHelp = false) 
{
	cerr << "************************" << endl;
	cerr << "*** Swiss Army Knife ***" << endl;
	cerr << "************************" << endl << endl;
	cerr << "Usage: sak [OPTION]... <SOURCE SEQUENCE FILE>" << endl;
	cerr << "\n";
	if (longHelp) {
		cerr << endl << "Main Options:" << endl;
		cerr << "  -o,  --output FILE            \t" << "set output filename (default: use stdout)" << endl;
		cerr << "  -q,  --qual                   \t" << "enable Fastq output (default: Fasta)" << endl;
		cerr << "  -h,  --help                   \t" << "print this help" << endl;
		cerr << endl << "Extract Options:" << endl;
		cerr << "  -s,  --sequence NUM           \t" << "select a single sequence by index" << endl;
        cerr << "  -sn, --sequence-name NAME     \t" << "select a single sequence by name" << endl;
		cerr << "  -ss, --sequences START END    \t" << "select sequences (default: select all)" << endl;
		cerr << "  -i,  --infix START END        \t" << "extract infix" << endl;
		cerr << "  -rc, --revcomp                \t" << "reverse complement" << endl;
		cerr << "  -l,  --max-length             \t" << "maximal number of sequence characters" << endl;
		cerr << "  -ll, --line-length            \t" << "maximal characters per output line" << endl;
	} else {
		cerr << "Try 'sak --help' for more information." << endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// Main part
int main(int argc, const char *argv[])
{
	unsigned fnameCount = 0;
	const char *fname[2] = { "" , "" };
	
	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option

			if (strcmp(argv[arg], "-s") == 0 || strcmp(argv[arg], "--sequence") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionSeqStart;
					if (!istr.fail())
					{
						if (optionSeqStart < 0)
							cerr << "sequence number must be a value >=0" << endl << endl;
						else
						{
							optionSeqEnd = optionSeqStart + 1;
							continue;
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

            if (strcmp(argv[arg], "-sn") == 0 || strcmp(argv[arg], "--sequence-name") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
                    optionSeqNameSet = true;
                    optionSeqName = argv[arg];
                    continue;
                }
            }

			if (strcmp(argv[arg], "-ss") == 0 || strcmp(argv[arg], "--sequences") == 0) {
				if (arg + 2 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionSeqStart;
					++arg;
					if (!istr.fail())
					{
						if (optionSeqStart < 0)
							cerr << "first sequence number must be a value >=0" << endl << endl;
						else
						{
							istringstream istr(argv[arg]);
							if (!istr.fail())
							{
								istr >> optionSeqEnd;
								if (optionSeqEnd < 0)
									cerr << "last sequence number must be a value >=0" << endl << endl;
								else
								{
									++optionSeqEnd;
									continue;
								}
							}
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-i") == 0 || strcmp(argv[arg], "--infix") == 0) {
				if (arg + 2 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionInfStart;
					++arg;
					if (!istr.fail())
					{
						if (optionInfStart < 0)
							cerr << "infix start a value >=0" << endl << endl;
						else
						{
							istringstream istr(argv[arg]);
							if (!istr.fail())
							{
								istr >> optionInfEnd;
								if (optionInfEnd < 0)
									cerr << "infix end a value >=0" << endl << endl;
								else
									continue;
							}
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-rc") == 0 || strcmp(argv[arg], "--revcomp") == 0) {
				optionRevComp = true;
				continue;
			}

			if (strcmp(argv[arg], "-q") == 0 || strcmp(argv[arg], "--qual") == 0) {
				optionFastQ = true;
				continue;
			}

			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv);
					return 0;
				}
				++arg;
				optionOutput = argv[arg];
				continue;
			}

			if (strcmp(argv[arg], "-ll") == 0 || strcmp(argv[arg], "--line-length") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
					istringstream istr(argv[arg]);
					istr >> optionLineLength;
                    continue;
                }
			}
 
            if (strcmp(argv[arg], "-l") == 0 || strcmp(argv[arg], "--max-length") == 0) {
                if (arg + 1 < argc) {
                    ++arg;
					istringstream istr(argv[arg]);
					istr >> optionMaxLength;
                    continue;
                }
			}
            
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, true);
				return 0;
			}
		}
		else {
			// parse file name
			if (fnameCount == 1) {
				printHelp(argc, argv);
				return 1;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 1) {
		printHelp(argc, argv);
		return 0;
	}
	
//____________________________________________________________________________
// input

	StringSet<TSeqString>	seqsIn, seqsOut;
	StringSet<CharString>	qualsIn, qualsOut;
	StringSet<CharString>	seqNamesIn, seqNamesOut;	// genome names, taken from the Fasta file

	if (!loadSeqs(seqsIn, qualsIn, seqNamesIn, fname[0]))
	{
		cerr << "Failed to open file" << fname[0] << endl;
		return 0;
	}
//	cout << lengthSum(seqsIn) << " bps of " << length(seqsIn) << " source sequence loaded." << endl;

//____________________________________________________________________________
// data processing

	resize(seqsOut, length(seqsIn));
	resize(qualsOut, length(qualsIn));
	for (unsigned i = 0; i < length(seqsIn); ++i)
	{
		unsigned end = optionInfEnd;
		if (end > length(seqsIn[i]))
			end = length(seqsIn[i]);

		if (optionInfStart < (int)length(seqsIn[i]))
		{
			seqsOut[i] = infix(seqsIn[i], optionInfStart, end);
			if (optionRevComp)
				reverseComplement(seqsOut[i]);
            qualsOut[i] = infix(qualsIn[i], optionInfStart, end);
			if (optionRevComp)
				reverse(qualsOut[i]);
		}
	}
	seqNamesOut = seqNamesIn;

//____________________________________________________________________________
// output

	if (optionFastQ)
		saveFastq(seqsOut, qualsOut, seqNamesOut);
	else
		saveFasta(seqsOut, seqNamesOut);
	
//____________________________________________________________________________

	return 0;
}
