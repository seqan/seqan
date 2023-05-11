#include <iostream>

#include <seqan/align_parallel.h>
#include <seqan/stream.h>  // for printing strings

int main()
{
    using TSequence = seqan2::String<seqan2::Dna>;
    using TAlignedSequence = seqan2::Gaps<TSequence>;
    using TThreadModel = seqan2::Parallel;
    using TVectorSpec = seqan2::Vectorial;
    using TExecPolicy = seqan2::ExecutionPolicy<TThreadModel, TVectorSpec>;

    // dummy sequences
    TSequence seqH = "CGATT";
    TSequence seqV = "CGAAATT";

    seqan2::StringSet<TAlignedSequence> seqs1;
    seqan2::StringSet<TAlignedSequence> seqs2;

    for (size_t i = 0; i < 100; ++i) // Create a data set of 100 dummy sequences
    {
        appendValue(seqs1, TAlignedSequence(seqH));
        appendValue(seqs2, TAlignedSequence(seqV));
    }

    TExecPolicy execPolicy;
    setNumThreads(execPolicy, 4);

    seqan2::Score<int16_t, seqan2::Simple> scoreAffine(2, -2, -1, -4);

    seqan2::String<int16_t> scores = seqan2::globalAlignment(execPolicy, seqs1, seqs2, scoreAffine);

    for (int16_t score : scores)
        std::cout << "Score: " << score << "\n";

    for (size_t pos = 0; pos < seqan2::length(seqs1); ++pos) // print out alignments
    {
        std::cout << seqs1[pos] << "\n";
        std::cout << seqs2[pos] << "\n\n";
    }

    return EXIT_SUCCESS;
}
