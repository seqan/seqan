#include <iostream>

#include <seqan/align_parallel.h>
#include <seqan/stream.h>  // for printint strings

int main()
{
    using TSequence = seqan2::String<seqan2::Dna>;
    using TThreadModel = seqan2::Parallel;
    using TVectorSpec = seqan2::Vectorial;
    using TExecPolicy = seqan2::ExecutionPolicy<TThreadModel, TVectorSpec>;

    // dummy sequences
    TSequence seqH = "CGATT";
    TSequence seqV = "CGAAATT";

    seqan2::StringSet<TSequence> seqs1;
    seqan2::StringSet<TSequence> seqs2;

    for (size_t i = 0; i < 100; ++i)
    {
        appendValue(seqs1, seqH);
        appendValue(seqs2, seqV);
    }

    TExecPolicy execPolicy;
    setNumThreads(execPolicy, 4);

    seqan2::Score<int16_t, seqan2::Simple> scoreAffine(2, -2, -1, -4);

    seqan2::String<int16_t> scores = seqan2::globalAlignmentScore(execPolicy, seqs1, seqs2, scoreAffine);

    for (int16_t score : scores)
        std::cout << "Score: " << score << "\n";

    return EXIT_SUCCESS;
}
