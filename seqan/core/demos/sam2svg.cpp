///A simple Sam viewer. Convert a Sam alignment file into a SVG vector graphic.
#include <fstream>
#include <iostream>
#include <string>

#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/misc/misc_svg.h>

using namespace seqan;

int main(int argc, const char *argv[])
{
	typedef FragmentStore<> TFragStore;

		typedef  TFragStore::TContigStore					TContigStore;
        typedef  Value<TContigStore>::Type						TContig;
		typedef TFragStore::TContigPos TContigPos;

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	CommandLineParser parser;
	addUsageLine(parser, "[OPTION]... <Sam file> <SVG file>");
	addUsageLine(parser, "[OPTION]... <Sam file> <GENOME file> <SVG file>");

	addOption(parser, CommandLineOption("c",  "contig",    "display only contig #NUM (default: show all contigs)", OptionType::Int | OptionType::Label | OptionType::List));
	addOption(parser, CommandLineOption("p",  "pos",    2, "set begin and end position (default: show whole strands)", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("l",  "lines",  2, "set first and last line of the alignment (default: show whole alignment)", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("a",  "ascii",     "output alignment in ASCII format instead of SVG", OptionType::Bool | OptionType::Label));
	addOption(parser, CommandLineOption("gs", "gap-space", "begin and end position (-p) are given in gap space", OptionType::Bool | OptionType::Label));
	requiredArguments(parser, 2);

	bool stop = !parse(parser, argc, argv, std::cerr);
	if (stop) return 0;

	//////////////////////////////////////////////////////////////////////////////
	// Extract and check options
	String<unsigned> contigs;
	TContigPos left = 0;
	TContigPos right = MaxValue<TContigPos>::VALUE;
	unsigned firstLine = 0;
	unsigned lastLine = MaxValue<unsigned>::VALUE;
	bool inASCII = false;

	if (isSetLong(parser, "pos"))
	{
		__int64 l = 0, r = 0;
		getOptionValueLong(parser, "pos", 0, l);
		getOptionValueLong(parser, "pos", 1, r);
		if ((l >= r) && (stop = true))
			std::cerr << "Begin position must be less than end position." << std::endl;
		left = l;
		right = r;
	}

	if (isSetLong(parser, "lines"))
	{
		getOptionValueLong(parser, "lines", 0, firstLine);
		getOptionValueLong(parser, "lines", 1, lastLine);
		if ((firstLine >= lastLine) && (stop = true))
			std::cerr << "First line must be less or equal than last line." << std::endl;
	}

	TFragStore store;
	std::ifstream samFile(toCString(getArgumentValue(parser, 0)), std::ios_base::in | std::ios_base::binary);
	std::ofstream ascii;
	SVGFile svg;

	//////////////////////////////////////////////////////////////////////////////
	// Optionally load genome file
	unsigned outArgNo = 1;
	if (argumentCount(parser) > 2)
	{
		if (!loadContigs(store, getArgumentValue(parser, 1)) && (stop = true))
			std::cerr << "Failed to load genome." << std::endl;
		++outArgNo;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Load Sam file
	if (!stop) read(samFile, store, Sam());

	//////////////////////////////////////////////////////////////////////////////
	// Choose contigs
	if (isSetLong(parser, "contig"))
	{
		resize(contigs, length(getOptionValuesLong(parser, "contig")));
		for (unsigned i = 0; i < length(contigs); ++i)
			getOptionValueLong(parser, "contig", i, contigs[i]);
	} else {
		resize(contigs, length(store.contigStore));
		for (unsigned i = 0; i < length(contigs); ++i)
			contigs[i] = i;
	}

	if (isSetLong(parser, "ascii"))
		inASCII = true;

	//////////////////////////////////////////////////////////////////////////////
	// Optionally load genome and open SVG file for writing
	if (!stop)
	{
		if (inASCII)
		{
			ascii.open(toCString(getArgumentValue(parser, outArgNo)), std::ios_base::out | std::ios_base::trunc);
			if (!ascii.is_open()) stop = true;
		} else
			if (!open(svg, toCString(getArgumentValue(parser, outArgNo)))) stop = true;

		if (stop) std::cerr << "Failed to open output file for writing." << std::endl;
	}

	// something went wrong
	if (stop)
	{
		std::cerr << "Exiting ..." << std::endl;
		return 1;
	}

	for(unsigned o=0;o<length(store.contigStore);++o)
	std::cerr<<store.contigNameStore[o]<<std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Output alignment
	AlignedReadLayout layout;
	std::cout << "Layouting reads ... " << std::flush;
	layoutAlignment(layout, store);
	std::cout << "done" << std::endl;

	for (unsigned i = 0; i < length(contigs); ++i)
		if (i < length(store.contigStore))
		{
			std::cout << "Writing contig " << contigs[i] << " ... " << std::flush;

			__int64 l = left;
			__int64 r = right;

			if (!isSetLong(parser, "gap-space"))
			{
				typedef Gaps<Nothing, AnchorGaps< TContig::TGapAnchors> >	TContigGaps;
				TContigGaps	contigGaps(store.contigStore[i].gaps);
				l = positionSeqToGap(contigGaps, l);
				if (r != MaxValue<TContigPos>::VALUE)
					r = positionSeqToGap(contigGaps, r);
			}

			if (r == MaxValue<TContigPos>::VALUE)
			{
				r = 0;
				for (unsigned j = 0; j < length(layout.contigRows[i]); ++j)
				{
					unsigned id = back(layout.contigRows[i][j]);
					if (r < store.alignedReadStore[id].beginPos)
						r = store.alignedReadStore[id].beginPos;
					if (r < store.alignedReadStore[id].endPos)
						r = store.alignedReadStore[id].endPos;
				}
			}

			std::cout <<l<<'\t'<<r<<'\t'<<firstLine<<'\t'<<lastLine<<std::endl;
			if (inASCII)
				printAlignment(ascii, Raw(), layout, store, contigs[i], l, r, firstLine, lastLine);
			else
				printAlignment(svg, Raw(), layout, store, contigs[i], l, r, firstLine, lastLine);

			std::cout << "done" << std::endl;
		}

	return 0;
}
