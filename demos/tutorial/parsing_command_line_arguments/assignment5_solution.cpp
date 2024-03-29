    #include <iostream>

    #include <seqan/arg_parse.h>

    struct ModifyStringOptions
    {
	unsigned period;
	unsigned rangeBegin, rangeEnd;
	bool toUppercase;
	bool toLowercase;
	seqan2::CharString inputFileName;

	ModifyStringOptions() :
	    period(1), rangeBegin(0), rangeEnd(0),toUppercase(false),
	    toLowercase(false)
	{}
    };

    seqan2::ArgumentParser::ParseResult
    parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
    {
	// Setup ArgumentParser.
	seqan2::ArgumentParser parser("modify_string");

	// Define Options
	addOption(parser, seqan2::ArgParseOption(
	    "I", "input-file",
	    "A text file that will printed with the modifications applied.",
	    seqan2::ArgParseArgument::INPUT_FILE));
	setValidValues(parser, "input-file", "txt");
	setRequired(parser, "input-file");

	addOption(parser, seqan2::ArgParseOption(
	    "i", "period", "Period to use for the index.",
	    seqan2::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "period", "1");
	setDefaultValue(parser, "period", "1");
	addOption(parser, seqan2::ArgParseOption(
	    "U", "uppercase", "Select to-uppercase as operation."));
	addOption(parser, seqan2::ArgParseOption(
	    "L", "lowercase", "Select to-lowercase as operation."));

	// Parse command line.
	seqan2::ArgumentParser::ParseResult res = seqan2::parse(parser, argc, argv);

	// Only extract  options if the program will continue after parseCommandLine()
	if (res != seqan2::ArgumentParser::PARSE_OK)
	    return res;

	// Extract option values.
	getOptionValue(options.period, parser, "period");
	options.toUppercase = isSet(parser, "uppercase");
	options.toLowercase = isSet(parser, "lowercase");
	getOptionValue(options.inputFileName, parser, "input-file");

	// If both to-uppercase and to-lowercase were selected then this is an error.
	if (options.toUppercase && options.toLowercase)
	{
	    std::cerr << "ERROR: You cannot specify both to-uppercase and to-lowercase!\n";
	    return seqan2::ArgumentParser::PARSE_ERROR;
	}

	return seqan2::ArgumentParser::PARSE_OK;
    }

    seqan2::CharString modifyString(seqan2::CharString const & text,
				   ModifyStringOptions const & options)
    {
	seqan2::CharString result;

	if (options.toLowercase)
	{
	    for (unsigned i = 0; i < length(text); ++i)
	    {
		if (i % options.period == 0u)
		    appendValue(result, tolower(text[i]));
		else
		    appendValue(result, text[i]);
	    }
	}
	else
	{
	    for (unsigned i = 0; i < length(text); ++i)
	    {
		if (i % options.period == 0u)
		    appendValue(result, toupper(text[i]));
		else
		    appendValue(result, text[i]);
	    }
	}

	return result;
    }

    int main(int argc, char const ** argv)
    {
	// Parse the command line.
	ModifyStringOptions options;
	seqan2::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan2::ArgumentParser::PARSE_OK)
	    return res == seqan2::ArgumentParser::PARSE_ERROR;

	std::fstream inFile(toCString(options.inputFileName), std::ios::binary | std::ios::in);
	if (inFile.good())
	{
	    std::cerr << "ERROR: Could not open input file " << options.inputFileName << '\n';
	    return 1;
	}
	seqan2::CharString text;
	while (inFile.good())
	{
	    char c = inFile.get();
	    if (inFile.good())
		appendValue(text, c);
	}
	std::cout << modifyString(text, options);

	return 0;
    }


