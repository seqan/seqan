#include <fstream>
#include <iostream>
#include <string>

#include <seqan/arg_parse.h>
#include <seqan/store.h>
#include "razers.h"

using namespace seqan;


int main(int argc, const char * argv[])
{

    //////////////////////////////////////////////////////////////////////////////
    // Define options
    ArgumentParser parser;
    RazerSOptions<> options;
    double delta = 0;

    addUsageLine(parser, "[OPTION]... <reads.fasta>");
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS FILE"));
    addOption(parser, ArgParseOption("mr", "mutation-rate", "set the percent mutation rate", ArgParseOption::DOUBLE));
    addOption(parser, ArgParseOption("qd", "quality-delta", "add a delta value to qualities", ArgParseOption::DOUBLE));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    FragmentStore<> store;
    getOptionValue(options.mutationRate, parser, "mutation-rate");
    getOptionValue(delta, parser, "quality-delta");
    options.mutationRate /= 100;

    //////////////////////////////////////////////////////////////////////////////
    // Load reads
    CharString readsFilename;
    getArgumentValue(readsFilename, parser, 0);

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(readsFilename)) || !loadReads(store, seqFileIn, options))
    {
        std::cerr << "Failed to load reads." << std::endl;
        return 1;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Output error distribution

    for (unsigned i = 0; i < length(options.avrgQuality); ++i)
    {
        double sequencingError = exp((options.avrgQuality[i] + delta) * log(10.0) / -10.0);
        double errorProb = 1.0 - (1.0 - sequencingError) * (1.0 - options.mutationRate);

//        std::cout << options.avrgQuality[i] << std::endl;
        std::cout << i << '\t' << errorProb << std::endl;
    }
    return 0;
}
