#include <fstream>
#include <iostream>
#include <string>

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/store.h>
#include "razers.h"

using namespace seqan;


int main(int argc, const char * argv[])
{

    //////////////////////////////////////////////////////////////////////////////
    // Define options
    CommandLineParser parser;
    RazerSOptions<> options;
    double delta = 0;

    addUsageLine(parser, "[OPTION]... <reads.fasta>");
    addOption(parser, CommandLineOption("mr", "mutation-rate", "set the percent mutation rate", OptionType::Double | OptionType::Label, 100.0 * options.mutationRate));
    addOption(parser, CommandLineOption("qd", "quality-delta", "add a delta value to qualities", OptionType::Double | OptionType::Label, delta));

    requiredArguments(parser, 1);

    bool stop = !parse(parser, argc, argv, std::cerr);
    if (stop)
        return 0;

    FragmentStore<> store;
    getOptionValueLong(parser, "mutation-rate", options.mutationRate);
    getOptionValueLong(parser, "quality-delta", delta);
    options.mutationRate /= 100;

    //////////////////////////////////////////////////////////////////////////////
    // Load reads

    if (!loadReads(store, toCString(getArgumentValue(parser, 0)), options) && (stop = true))
        std::cerr << "Failed to load reads." << std::endl;

    // something went wrong
    if (stop)
    {
        std::cerr << "Exiting ..." << std::endl;
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
