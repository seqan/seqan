#include <fstream>
#include <iostream>
#include <string>

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/store.h>
#include <seqan/random.h>
#include <seqan/parallel.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace seqan;

template <typename TReads, typename TErrorDist, typename TStore>
void simulateReads(TReads & reads, TErrorDist const & errorDist, TStore const & store)
{
    String<Rng<MersenneTwister> > rng;
#ifdef _OPENMP
    resize(rng, omp_get_max_threads());
#else
    resize(rng, 1);
#endif

    Pdf<Uniform<unsigned> > pdfContig(0, length(store.contigStore) - 1);
    Pdf<Uniform<unsigned> > pdfOrientation(0, 1);
    Pdf<Uniform<double> >   pdfFrac(0.0, 1.0);
    Pdf<Uniform<unsigned> > pdfSubstitute(1, 3);
    unsigned readLen = length(errorDist);

    SEQAN_OMP_PRAGMA(parallel for)
    for (int i = 0; i < (int)length(reads); ++i)
    {
        unsigned id = 0;
#ifdef _OPENMP
        id = omp_get_thread_num();
#endif
        unsigned contigId = pickRandomNumber(rng[id], pdfContig);
        unsigned pos;
        bool ok;

        do
        {
            Pdf<Uniform<unsigned> > pdfPos(0, length(store.contigStore[contigId].seq) - readLen);
            pos = pickRandomNumber(rng[id], pdfPos);
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

        if (pickRandomNumber(rng[id], pdfOrientation) == 1)
            reverseComplement(reads[i]);

        for (unsigned j = 0; j < readLen; ++j)
        {
            if (pickRandomNumber(rng[id], pdfFrac) < errorDist[j])
            {
                reads[i][j] = (Dna5)((ordValue(reads[i][j]) + pickRandomNumber(rng[id], pdfSubstitute)) & 3);
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
    CommandLineParser parser;

    addUsageLine(parser, "[OPTION]... <reference.fasta> <error_dist file>");

    CharString outputFilename = "reads.fa";
    unsigned numReads = 1000000;
    String<double> errorProb;

    addOption(parser, CommandLineOption("n", "num-reads", "number of reads", OptionType::Int | OptionType::Label, numReads));
    addOption(parser, CommandLineOption("o", "output", "output read filename", OptionType::String | OptionType::Label, outputFilename));
    requiredArguments(parser, 2);

    bool stop = !parse(parser, argc, argv, std::cerr);
    if (stop)
        return 0;

    //////////////////////////////////////////////////////////////////////////////
    // Extract and check options
    getOptionValueLong(parser, "num-reads", numReads);
    getOptionValueLong(parser, "output", outputFilename);

    std::ifstream errorDistFile(toCString(getArgumentValue(parser, 1)), std::ios_base::in | std::ios_base::binary);

    while (true)
    {
        double error;
        errorDistFile >> error;
        if (!errorDistFile.good())
            break;
        appendValue(errorProb, error);
    }

    FragmentStore<> store;

    if (!loadContigs(store, getArgumentValue(parser, 0)) && (stop = true))
        std::cerr << "Failed to load genome." << std::endl;

    // something went wrong
    if (stop)
    {
        std::cerr << "Exiting ..." << std::endl;
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
