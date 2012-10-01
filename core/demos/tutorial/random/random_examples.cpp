// FRAGMENT(header)
#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/random.h>

using namespace seqan;

const int SEED = 42;

int main()
{
// FRAGMENT(random-number-generation-raw)
    Rng<MersenneTwister> rng(SEED);
    std::cout << "pickRandomNumber(rng) == " << pickRandomNumber(rng) << std::endl;

// FRAGMENT(random-number-generation-metafunction-value)
    typedef Value<Rng<MersenneTwister> >::Type TMTValue;
    TMTValue value = pickRandomNumber(rng);

// FRAGMENT(ignore-variable-value)
    // This fragment is not included from the tutorial.
    (void) value;  // Suppress unused variable warning.

// FRAGMENT(random-number-generation-pdf)
    Pdf<Uniform<double> > uniformDouble(0, 1);
    std::cout << "pickRandomNumber(rng, uniformDouble) == " << pickRandomNumber(rng, uniformDouble) << std::endl;
    Pdf<Uniform<int> > uniformInt(0, 42);
    std::cout << "pickRandomNumber(rng, uniformInt) == " << pickRandomNumber(rng, uniformInt) << std::endl;
    Pdf<Normal> normal(0, 1);
    std::cout << "pickRandomNumber(rng, normal) == " << pickRandomNumber(rng, normal) << std::endl;
    Pdf<LogNormal> logNormal(0, 1);
    std::cout << "pickRandomNumber(rng, logNormal) == " << pickRandomNumber(rng, logNormal) << std::endl;

// FRAGMENT(random-number-generation-log-normal)
    Pdf<LogNormal> logNormal2(0, 1, MuSigma());
    std::cout << "pickRandomNumber(rng, logNormal2) == " << pickRandomNumber(rng, logNormal2) << std::endl;
    Pdf<LogNormal> logNormal3(0.1, 1, MeanStdDev());
    std::cout << "pickRandomNumber(rng, logNormal3) == " << pickRandomNumber(rng, logNormal3) << std::endl;

// FRAGMENT(shuffling)
    CharString container = "Hello World!";
    shuffle(container, rng);
    std::cout << "shuffle(\"Hello World!\") == " << container << std::endl;

    return 0;
}
