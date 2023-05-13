#include <iostream>

#include <seqan/align_parallel.h>
#include <seqan/stream.h>  // for printint strings

int main()
{
    using TSequence = seqan::String<seqan::Dna>;
    using TThreadModel = seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>;
    using TVectorSpec = seqan::Vectorial;
    using TExecPolicy = seqan::ExecutionPolicy<TThreadModel, TVectorSpec>;

    // dummy sequences
    TSequence seqH;
    TSequence seqV;

    for (size_t i = 0; i < 10000; ++i)
    {
        seqan::appendValue(seqH, 'A');
        seqan::appendValue(seqV, 'A');
    }

    seqan::StringSet<TSequence> seqs1;
    seqan::StringSet<TSequence> seqs2;

    for (size_t i = 0; i < 100; ++i)

    {
        seqan::appendValue(seqs1, seqH);
        seqan::appendValue(seqs2, seqV);
    }

    TExecPolicy execPolicy;
    seqan::setBlockSize(execPolicy, 500); // Sets the size of blocks used internally to partition the alignment matrix.
    seqan::setParallelAlignments(execPolicy, 10); // Compute ten alignments at the same time.
    seqan::setNumThreads(execPolicy, 4); // Use four threads to compute the actual alignment.

    seqan::Score<int32_t, seqan::Simple> scoreAffine(2, -2, -1, -4);

    seqan::String<int32_t> scores = seqan::globalAlignmentScore(execPolicy, seqs1, seqs2, scoreAffine);

    for (int32_t score : scores)
        std::cout << "Score: " << score << "\n";

    return EXIT_SUCCESS;
}
