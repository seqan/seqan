#include <iostream>

#include <seqan/align_parallel.h>
#include <seqan/stream.h>  // for printint strings

int main()
{
    using TSequence = seqan2::String<seqan2::Dna>;
    using TThreadModel = seqan2::WavefrontAlignment<seqan2::BlockOffsetOptimization>;
    using TVectorSpec = seqan2::Vectorial;
    using TExecPolicy = seqan2::ExecutionPolicy<TThreadModel, TVectorSpec>;

    // dummy sequences
    TSequence seqH;
    TSequence seqV;

    for (size_t i = 0; i < 10000; ++i)
    {
        seqan2::appendValue(seqH, 'A');
        seqan2::appendValue(seqV, 'A');
    }

    seqan2::StringSet<TSequence> seqs1;
    seqan2::StringSet<TSequence> seqs2;

    for (size_t i = 0; i < 100; ++i)

    {
        seqan2::appendValue(seqs1, seqH);
        seqan2::appendValue(seqs2, seqV);
    }

    TExecPolicy execPolicy;
    seqan2::setBlockSize(execPolicy, 500); // Sets the size of blocks used internally to partition the alignment matrix.
    seqan2::setParallelAlignments(execPolicy, 10); // Compute ten alignments at the same time.
    seqan2::setNumThreads(execPolicy, 4); // Use four threads to compute the actual alignment.

    seqan2::Score<int32_t, seqan2::Simple> scoreAffine(2, -2, -1, -4);

    seqan2::String<int32_t> scores = seqan2::globalAlignmentScore(execPolicy, seqs1, seqs2, scoreAffine);

    for (int32_t score : scores)
        std::cout << "Score: " << score << "\n";

    return EXIT_SUCCESS;
}
