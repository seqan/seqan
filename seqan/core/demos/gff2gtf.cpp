///A simple annotation converter. Convert a Gff to Gtf or vice versa.
#include <fstream>
#include <iostream>
#include <string>

#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

int main(int argc, const char *argv[])
{
	typedef FragmentStore<> TFragStore;
	
	//////////////////////////////////////////////////////////////////////////////
	// Define options
	CommandLineParser parser;
	addUsageLine(parser, "[OPTION]... <infile> <outfile>");
	
	addOption(parser, CommandLineOption("gff",  "",    "write annotation in Gff format", OptionType::Bool));
	addOption(parser, CommandLineOption("gtf",  "",    "write annotation in Gtf format (default)", OptionType::Bool));
	requiredArguments(parser, 2);

	bool stop = !parse(parser, argc, argv, std::cerr);
	if (stop) return 0;

	//////////////////////////////////////////////////////////////////////////////
	// Extract and check options	
	TFragStore store;
	std::ifstream inFile(toCString(getArgumentValue(parser, 0)), std::ios_base::in);
	std::ofstream outFile(toCString(getArgumentValue(parser, 1)), std::ios_base::out);

	if (!inFile.is_open() && (stop = true))
		std::cerr << "Failed to open annotation infile for reading." << std::endl;
	else
		read(inFile, store, Gff());

	if (!stop)
	{
		if (!outFile.is_open() && (stop = true))
			std::cerr << "Failed to open annotation outfile for writing." << std::endl ;
		else
		{
			if (isSetShort(parser, "gff"))
				write(outFile, store, Gff());
			else
				write(outFile, store, Gtf());
		}
	}
	
	if (stop)
	{
		std::cerr << "Exiting ..." << std::endl;
		return 1;
	}

	return 0;
}
