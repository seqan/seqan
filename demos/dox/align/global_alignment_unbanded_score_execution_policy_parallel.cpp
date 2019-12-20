#include <iostream>

#include <seqan/align_parallel.h>
#include <seqan/stream.h>  // for printint strings

int main()
{
    using TSequence = seqan::String<seqan::Dna>;
    using TThreadModel = seqan::Parallel;
    using TVectorSpec = seqan::Vectorial;
    using TExecPolicy = seqan::ExecutionPolicy<TThreadModel, TVectorSpec>;

    // dummy sequences
    TSequence seqH = "CGATT";
    TSequence seqV = "CGAAATT";

    seqan::StringSet<TSequence> seqs1;
    seqan::StringSet<TSequence> seqs2;

    for (size_t i = 0; i < 100; ++i)
    {
        appendValue(seqs1, seqH);
        appendValue(seqs2, seqV);
    }

    TExecPolicy execPolicy;
    setNumThreads(execPolicy, 4);

    seqan::Score<int16_t, seqan::Simple> scoreAffine(2, -2, -1, -4);

    seqan::String<int16_t> scores = seqan::globalAlignmentScore(execPolicy, seqs1, seqs2, scoreAffine);

    for (int16_t score : scores)
        std::cout << "Score: " << score << "\n";

    return EXIT_SUCCESS;
}
