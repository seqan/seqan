#include <iostream>
#include <fstream>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace std;
using namespace seqan;


template <typename TIndex>
void findMUMs(TIndex &esa, unsigned minLen)
{
	typename Iterator<TIndex, Mums>::Type it(esa, minLen);  // set minimum MUM length
	String< typename SAValue<TIndex>::Type > occs;          // temp. string storing the hit positions

	cout << resetiosflags(ios::left);
	while (!atEnd(it)) 
	{
		occs = getOccurrences(it);                          // gives hit positions (seqNo,seqOfs)
		orderOccurrences(occs);                             // order them by seqNo
		
		for(unsigned i = 0; i < length(occs); ++i)
			cout << setw(8)
					<< getValueI2(occs[i]) + 1              // output them in MUMmer's output format
					<< "  ";

		cout << setw(8) 
				<< repLength(it)
				<< endl;

		++it;
	}
	cout << setiosflags(ios::left) << endl;
}


template <typename TSpec>
int runMummy(int argc, const char *argv[], unsigned seqCount, unsigned minLen)
{
	typedef String<Dna5, TSpec> TString;
	typedef StringSet<TString>  TStringSet;
	typedef Index<TStringSet>   TIndex;

	TIndex	index;

	// import sequences
    StringSet<CharString> meta;
	for(int arg = 1, seq = 0; arg < argc; ++arg)
	{
		// skip two argument options
		if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--profile") == 0 ||
			strcmp(argv[arg], "-l") == 0 || strcmp(argv[arg], "--minlen") == 0) 
		{
			++arg;
			continue;
		}

		if (argv[arg][0] != '-')
		{
		    try
            {
                SeqFileIn file(argv[arg]);
                readRecords(meta, indexText(index), file);
            }
            catch (FileOpenError const & err)
            {
                std::cerr << "Problem opening file " << argv[arg] << ": "
                          << err.what() << "\n";
                return 1;
            }
            catch (ParseError const & err)
            {
                std::cerr << "Problem reading file " << argv[arg] << ": "
                          << err.what() << "\n";
                return 1;
            }
            clear(meta);
			++seq;
		}
	}
	cout << lengthSum(indexText(index)) << " bps sequence imported." << endl;

	findMUMs(index, minLen);

	return 0;
}


void printHelp(int, const char *[], bool longHelp = false) 
{
	cout << "***************************************" << endl;
	cout << "***        Simple MUM finder        ***" << endl;
	cout << "*** written by David Weese (c) 2007 ***" << endl;
	cout << "***************************************" << endl << endl;
	cout << "Usage: mummy [OPTION]... <SEQUENCE FILE> ... <SEQUENCE FILE>" << endl;
	if (longHelp) {
		cout << endl << "Options:" << endl;
		cout << "  -e, --extern          \t" << "use external memory (for large datasets)" << endl;
		cout << "  -l, --minlen          \t" << "set minimum MUM length" << endl;
		cout << "                        \t" << "if not set, default value is 20" << endl;
		cout << "  -h, --help            \t" << "print this help" << endl;
	} else {
		cout << "Try 'mummy --help' for more information." << endl;
	}
}


int main(int argc, const char *argv[])
{
	if (argc < 2) {
		printHelp(argc, argv);
		return 0;
	}

	unsigned optMinLen = 20;
	bool	 optExtern = false;

	unsigned seqCount = 0;
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option

			if (strcmp(argv[arg], "-e") == 0 || strcmp(argv[arg], "--extern") == 0) {
				// use external memory algorithms
				optExtern = true;
				continue;
			}
			if (strcmp(argv[arg], "-l") == 0 || strcmp(argv[arg], "--minlen") == 0) {
				// minimum match length
				if (arg + 1 == argc) {
					printHelp(argc, argv);
					return 0;
				}
				++arg;
				optMinLen = atoi(argv[arg]);
				continue;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, true);
				return 0;
			}
		} else {
			// parse sequence file
			++seqCount;
		}
	}

	if (optExtern)
		return runMummy<External<> >(argc, argv, seqCount, optMinLen);
	else
		return runMummy<Alloc<> >(argc, argv, seqCount, optMinLen);
}
