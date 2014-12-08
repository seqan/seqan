// NOTE(esiragusa): this tool is commented out in the cmake file.

#include <fstream>
#include <iostream>
#include <string>

// NOTE(esiragusa): this tool uses the old cmdparser.
#include <seqan/misc/misc_cmdparser.h>
#include "razers.h"

using namespace seqan;

int main(int argc, const char * argv[])
{

    //////////////////////////////////////////////////////////////////////////////
    // Define options
    typedef RazerSOptions<> TOptions;
    CommandLineParser parser;
    TOptions options;

    addUsageLine(parser, "[OPTION]... <error_dist file>");

    unsigned maxOverlap = 10;
    unsigned maxErrors = 10;

    addOption(parser, CommandLineOption("d", "delta", 2, "delta range", OptionType::Int | OptionType::Label));
    addOption(parser, CommandLineOption("mo", "max-overlap", "estimate for overlaps 0,1,...,max-overlap", OptionType::Int | OptionType::Label, maxOverlap));
    addOption(parser, CommandLineOption("me", "max-errors", "estimate for errors 0,1,...,max-errors", OptionType::Int | OptionType::Label, maxErrors));
    requiredArguments(parser, 1);

    bool stop = !parse(parser, argc, argv, std::cerr);
    if (stop)
        return 0;

    //////////////////////////////////////////////////////////////////////////////
    // Extract and check options
    getOptionValueLong(parser, "max-overlap", maxOverlap);
    getOptionValueLong(parser, "max-errors", maxErrors);

    int deltaMin = 11;
    int deltaMax = 31;
    if (isSetLong(parser, "delta"))
    {
        getOptionValueLong(parser, "delta", 0, deltaMin);
        getOptionValueLong(parser, "delta", 1, deltaMax);
    }

    std::ifstream errorDistFile(toCString(getArgumentValue(parser, 0)), std::ios_base::in | std::ios_base::binary);
    double prob;
//    char buf[256];

    while (true)
    {
//        int pos, dummy;
//		errorDistFile >> pos;
//		errorDistFile >> dummy;
//		errorDistFile >> dummy;
//		errorDistFile >> dummy;
        errorDistFile >> prob;
        if (!errorDistFile.good())
            break;
//        errorDistFile.getline(buf, 256);
//        std::cout<< errorProb<<std::endl;
        appendValue(options.errorProb, prob);
    }

    if (empty(options.errorProb))
    {
        std::cerr << "Couldn't read error distribution file." << std::endl;
        return 1;
    }

//    clear(options.errorProb);
//    resize(options.errorProb, 22,0.4);

    resize(options.readLengths, length(options.errorProb), 0);
    appendValue(options.readLengths, 1);
    options.lossRate = 1.0;

//	for(unsigned k=0;k<length(options.errorProb);++k)
//		std::cout<<"  "<<options.errorProb[k];

    //////////////////////////////////////////////////////////////////////////////
    // Output losses to file


    std::cout << "# ==========================================================" << std::endl;
    std::cout << "# Estimated losses for " << getArgumentValue(parser, 0) << std::endl;
    std::cout << "# Format Description" << std::endl;
    std::cout << "# <q>\t<overlap>\t<real time>";
    for (unsigned i = 0; i <= maxErrors; ++i)
        std::cout << "\t<" << i << "-error loss>";
    std::cout << std::endl;
    std::cout << "# ==========================================================" << std::endl;

    //////////////////////////////////////////////////////////////////////////////
    // Compute losses for different errors, overlaps

    Shape<Dna, GenericShape> dummyShape;
    String<unsigned> deltas;

    for (int delta = deltaMax; delta >= deltaMin; --delta)
    {
        clear(deltas);
        resize(deltas, maxOverlap + 1, delta);
        unsigned maxLength = length(options.readLengths) - 1;
        options.errorRate = (double)maxErrors / (double)maxLength + 0.00001;

        typedef TOptions::TProb TFloat;
        String<TFloat> estLosses;
        estimatePigeonholeLosses(estLosses, deltas, options);

//            for (unsigned e = 0; e <= maxErrors; ++e)
//                std::cerr << "prob to see " <<e<< " errors: " << 1000000*estLosses[2 + e] << std::endl;


        for (int idx = length(estLosses) - (3 + maxErrors); idx > 0; idx -= 3 + maxErrors)
        {
            unsigned ol = estLosses[idx];
            unsigned q  = estLosses[idx + 1];

            std::cout << q - ol << '\t' << ol << "\t0:00.00";

            for (unsigned e = 0; e <= maxErrors; ++e)
                std::cout << '\t' << estLosses[idx + 2 + e] / estLosses[2 + e];
            std::cout << std::endl;

//            for (unsigned e = 0; e <= maxErrors; ++e)
//                std::cerr << '\t' << estLosses[idx + 2 + e];
//            std::cerr << std::endl;
        }
    }

    return 0;
}
