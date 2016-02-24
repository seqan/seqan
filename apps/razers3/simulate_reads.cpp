#include <fstream>
#include <iostream>
#include <string>
#include <random>

#include <seqan/arg_parse.h>
#include <seqan/store.h>
#include <seqan/parallel.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace seqan;

template <typename TReads, typename TErrorDist, typename TStore>
void simulateReads(TReads & reads, TErrorDist const & errorDist, TStore const & store)
{
    String<std::mt19937 > rng;
#ifdef _OPENMP
    resize(rng, omp_get_max_threads());
#else
    resize(rng, 1);
#endif

    std::uniform_int_distribution<unsigned> distContig(0, length(store.contigStore) - 1);
    std::uniform_int_distribution<unsigned> distOrientation(0, 1);
    std::uniform_real_distribution<double> distFrac(0.0, 1.0);
    std::uniform_int_distribution<unsigned> distSubstitute(1, 3);
    unsigned readLen = length(errorDist);

    SEQAN_OMP_PRAGMA(parallel for)
    for (int i = 0; i < (int)length(reads); ++i)
    {
        unsigned id = 0;
#ifdef _OPENMP
        id = omp_get_thread_num();
#endif
        unsigned contigId = distContig(rng[id]);
        unsigned pos;
        bool ok;

        do
        {
            std::uniform_int_distribution<unsigned> distPos(0, length(store.contigStore[contigId].seq) - readLen);
            pos = distPos(rng[id]);
            reads[i] = infix(store.contigStore[contigId].seq, pos, pos + readLen);

            ok = true;
            for (unsigned j = 0; j < readLen; ++j)
                if (ordValue(reads[i][j]) == 4)
                {
                    ok = false;
                    break;
                }
        }
        while (!ok);

        if (distOrientation(rng[id]) == 1)
            reverseComplement(reads[i]);

        for (unsigned j = 0; j < readLen; ++j)
        {
            if (distFrac(rng[id]) < errorDist[j])
            {
                reads[i][j] = (Dna5)((ordValue(reads[i][j]) + distSubstitute(rng[id])) & 3);
                assignQualityValue(reads[i][j], 0);
            }
            else
                assignQualityValue(reads[i][j], 40);
        }
    }
}

int main(int argc, const char * argv[])
{

    //////////////////////////////////////////////////////////////////////////////
    // Define options
    ArgumentParser parser;

    addUsageLine(parser, "[OPTION]... <reference.fasta> <error_dist file>");

    CharString outputFilename = "reads.fa";
    unsigned numReads = 1000000;
    String<double> errorProb;

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "REFERENCE FILE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "ERROR DIST FILE"));
    addOption(parser, ArgParseOption("n", "num-reads", "number of reads", ArgParseOption::INTEGER));
    addOption(parser, ArgParseOption("o", "output", "output read filename", ArgParseOption::STRING));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    //////////////////////////////////////////////////////////////////////////////
    // Extract and check options
    getOptionValue(numReads, parser, "num-reads");
    getOptionValue(outputFilename, parser, "output");

    CharString distFile;
    getArgumentValue(distFile, parser, 1);
    std::ifstream errorDistFile(toCString(distFile), std::ios_base::in | std::ios_base::binary);

    while (true)
    {
        double error;
        errorDistFile >> error;
        if (!errorDistFile.good())
            break;
        appendValue(errorProb, error);
    }

    FragmentStore<> store;

    CharString referenceFile;
    getArgumentValue(distFile, parser, 0);
    if (!loadContigs(store, toCString(referenceFile)))
    {
        std::cerr << "Failed to load genome." << std::endl;
        return 1;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Simulate reads

    String<String<Dna5Q> > reads;
    resize(reads, numReads);
    simulateReads(reads, errorProb, store);

    //////////////////////////////////////////////////////////////////////////////
    // Output losses to file

    std::ofstream readFile(toCString(outputFilename));
    for (unsigned i = 0; i < length(reads); ++i)
    {
        readFile << "@read" << i << std::endl;
        readFile << reads[i] << std::endl;
        readFile << "+" << std::endl;
        for (unsigned j = 0; j < length(reads[i]); ++j)
            readFile << (char)('!' + getQualityValue(reads[i][j]));
        readFile << std::endl;
    }

    return 0;
}
