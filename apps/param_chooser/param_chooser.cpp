 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by Anne-Katrin Emde

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

#define USE_LOGVALUES		// this is recommended when using probability values
#define RUN_RAZERS
#define RAZERS_CONCATREADS		// use <ConcatDirect> StringSet to store reads
#define RAZERS_MASK_READS		// remove matches with max-hits optimal hits on-the-fly
#define RAZERS_MEMOPT			// optimize memory consumption
#define RAZERS_PRUNE_QGRAM_INDEX
#define SEQAN_PROFILE
#define NON_REDUNDANT
//#define LOSSRATE_VALIDATION	//generates output for loss rate validation (empirical vs computed)

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "../splazers/paramChooser.h"



using namespace seqan;
using namespace std;

typedef float TFloat;





//////////////////////////////////////////////////////////////////////////////
// Print usage
void printHelp(int, const char *[], ParamChooserOptions & pm_options, bool longHelp = false) 
{
	cerr << "*******************************************************************" << endl;
	cerr << "********** Compute filter parameters for given loss rate **********" << endl;
	cerr << "*******************************************************************" << endl << endl;
	cerr << "Usage: paramChooser [OPTIONS]... " << endl;
	if (longHelp) {
		cerr << endl << "Options:" << endl;
		cerr << std::endl << " Parameter choosing options: " << std::endl;
		cerr << "  -n,   --length NUM            \t" << "sequence length ("<<pm_options.totalN<<")" << endl;
		cerr << "  -i,   --percent-identity NUM  \t" << "set the percent identity threshold (95)" << endl;
		cerr << "  -rr,  --recognition-rate NUM  \t" << "set the percent recognition rate (99.0)" << endl;
		cerr << "  -ug,  --ungapped              \t" << "only consider ungapped shapes (off)" << endl;
		cerr << "  -og,  --one-gapped            \t" << "only consider shapes containing at most one gap (of arbitrary length)" << endl;
		cerr << std::endl << " Parameter computing options: " << std::endl;
		cerr << "  -pf,  --param-folder STR      \t" << "output directory where filter parameter files will be stored" << endl;
		cerr << "  -sf,  --shape-file STR        \t" << "file containing additional shapes (optional)" << endl;
		cerr << "  -ap,  --append                \t" << "append parameters to previously computed files (default: off)" << endl;
		cerr << "  -d,   --error-distribution    \t" << "file containing mismatch probabilities (must contain at least n values, one value per line)" << endl;
		cerr << "  -pi,  --prob-insert           \t" << "probability of an insertion (" << pm_options.optionProbINSERT << ")" << endl;
		cerr << "  -pd,  --prob-delete           \t" << "probability of a deletion (" << pm_options.optionProbDELETE << ")" << endl;
		cerr << "                                \t" << "(for hamming-only filters use -pi 0 -pd 0)" << endl;
		cerr << "  -prbf,--prb-folder STR        \t" << "directory of [_prb.txt|.fastq|fastqint] files containing qualitiy values (optional)" << endl;
		cerr << "  -pq,  --phred-qualities       \t" << "fastq files contain Phred qualities (default: Solexa qualities)" << endl;
		cerr << "  -h,   --help                  \t" << "print this help" << endl << endl;
	} else {
		cerr << "Try 'chooseFilteringParameters --help' for more information." << endl;
	}
}

int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line
//	static const TFloat epsilon = 0.0000001;	
	static const TFloat epsilon = (TFloat)0.0000000000001;	

        RazerSOptions<> r_options;
        ParamChooserOptions pm_options;

	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option
			if (strcmp(argv[arg], "-n") == 0 || strcmp(argv[arg], "--length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.totalN;
					if (!istr.fail())
					{
						if (pm_options.totalN < 1 || pm_options.totalN > 100)
							cerr << "sequence length must be a value between 1 and 100" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-i") == 0 || strcmp(argv[arg], "--percent-identity") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionErrorRate;
					pm_options.optionErrorRate = (100.0 - pm_options.optionErrorRate) / 100.0;
					if (!istr.fail())
					{
						if (pm_options.optionErrorRate < 0 || pm_options.optionErrorRate > 0.1)
							cerr << "Percent identity threshold must be a value between 90 and 100" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-rr") == 0 || strcmp(argv[arg], "--recognition-rate") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionLossRate;
					if (!istr.fail())
					{
						if (pm_options.optionLossRate < 80.0 || pm_options.optionLossRate > 100.0)
							cerr << "Loss rate must be a value between 0 and 100" << endl << endl;
						else
						{
							pm_options.optionLossRate = 100.0 - pm_options.optionLossRate;
							pm_options.optionLossRate /= 100.0;
							continue;
						}
					}
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-pi") == 0 || strcmp(argv[arg], "--prob-insert") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionProbINSERT;
					if (!istr.fail())
					{
						if (pm_options.optionProbINSERT < 0 || pm_options.optionProbINSERT > 1)
							cerr << "Insert probability must be a value between 0 and 1" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}

			if (strcmp(argv[arg], "-pd") == 0 || strcmp(argv[arg], "--prob-delete") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionProbDELETE;
					if (!istr.fail())
					{
						if (pm_options.optionProbDELETE < 0 || pm_options.optionProbDELETE > 1)
							cerr << "Delete probability must be a value between 0 and 1" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-prbf") == 0 || strcmp(argv[arg], "--prb-folder") == 0) {
				if (arg + 1 < argc) {
					++arg;
					pm_options.fnameCount0 = true;
					fstream file;
					pm_options.fname[0] = argv[arg];
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-pf") == 0 || strcmp(argv[arg], "--param-folder") == 0) { //prefix for previously computed param files
				if (arg + 1 < argc) {
					++arg;
					pm_options.paramFolder = argv[arg];
//					cout << "Session id prefix specified\n";
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-sf") == 0 || strcmp(argv[arg], "--shape-file") == 0) { //prefix for previously computed param files
				if (arg + 1 < argc) {
					++arg;
					pm_options.shapeFile = argv[arg];
					pm_options.useDefaultShapes = false;
//					cout << "Session id prefix specified\n";
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--prefix") == 0) { //prefix for previously computed param files
				if (arg + 1 < argc) {
					++arg;
					pm_options.prefixCount = true;
					fstream file;
					pm_options.fprefix[0] = argv[arg];
//					cout << "Session id prefix specified\n";
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-d") == 0 || strcmp(argv[arg], "--distribution-file") == 0) { //should also support fastq files
				if (arg + 1 < argc) {
					++arg;
					pm_options.fnameCount1 = true;
					fstream file;
					pm_options.fname[1] = argv[arg];
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-ha") == 0 || strcmp(argv[arg], "--hamming") == 0) {
				pm_options.optionHammingOnly = true;
				continue;
			}
			if (strcmp(argv[arg], "-ap") == 0 || strcmp(argv[arg], "--append") == 0) {
				pm_options.appendToPrevious = true;
				continue;
			}

			if (strcmp(argv[arg], "-pq") == 0 || strcmp(argv[arg], "--phred-qualities") == 0) {
				pm_options.solexaQual = false;
				continue;
			}

			if (strcmp(argv[arg], "-og") == 0 || strcmp(argv[arg], "--one-gapped") == 0) {
				pm_options.chooseOneGappedOnly = true;       //optionChooseOneGappedOnly chooses shape with at most one gap
                                pm_options.chooseUngappedOnly = false;
				continue;
			}

			if (strcmp(argv[arg], "-ug") == 0 || strcmp(argv[arg], "--ungapped") == 0) {
                                if(pm_options.chooseOneGappedOnly) continue;     //if both ungapped and onegapped specified --> optionChooseOneGappedOnly 
                                pm_options.chooseUngappedOnly = true;
				continue;
			}
			if (strcmp(argv[arg], "-mt") == 0 || strcmp(argv[arg], "--min-threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.minThreshold;
					if (!istr.fail())
					{
						if (pm_options.minThreshold < 1 || pm_options.minThreshold > 20)
							cerr << "minimum threshold should be a value between 1 and 20" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-mq") == 0 || strcmp(argv[arg], "--max-weight") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.maxWeight;
					if (!istr.fail())
					{
						if (pm_options.maxWeight < 6 || pm_options.maxWeight > 14)
							cerr << "maximum weight should be a value between 6 and 14" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}

			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, pm_options, true);
				return 0;
			}
		}
	}
	pm_options.optionErrorRate += epsilon;
	pm_options.optionLossRate += epsilon;
	
//	pm_options.verbose = true;
	r_options._debugLevel = 1;
	r_options.errorRate = pm_options.optionErrorRate;

	if(length(pm_options.paramFolder) == 0) 
	{
		string appFolder = argv[0];
		size_t lastPos = appFolder.find_last_of('/') + 1;
		if (lastPos == appFolder.npos + 1) lastPos = appFolder.find_last_of('\\') + 1;
		if (lastPos == appFolder.npos + 1) lastPos = 0;
		appFolder.erase(lastPos);
		pm_options.paramFolderPath = appFolder;
	}

	pm_options.verbose = true;

	
	if(pm_options.optionProbINSERT <= epsilon && pm_options.optionProbDELETE <= epsilon)
		pm_options.optionHammingOnly=true;

	r_options.hammingOnly = pm_options.optionHammingOnly;


        chooseParams(r_options, pm_options);

	return 0;
}
