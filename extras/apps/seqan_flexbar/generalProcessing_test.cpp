#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "generalProcessing.h"

using namespace seqan;

SEQAN_DEFINE_TEST(findN_test)
{
    StringSet<String<Dna5> >seqs;
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(seqs, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(seqs, "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN");
    appendValue(seqs, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
    appendValue(seqs, "ACTGTACGTGATCG.AATGCTGACTGACTGAC");
    appendValue(seqs, ".ACTGTACGTGATCG.AATGCTGACTGACTGAC");
    appendValue(seqs, "ACT.GTACGTGATCG.AATGCTGACTGA.CTGAC");
    appendValue(seqs, "ACT.GTACGTGATCG.AATGC.TGACTGACTG.AC");
    appendValue(seqs, ".......");
    appendValue(seqs, "");
    unsigned allowed = 3;
    StringSet<int> exspectedNoSub;
    appendValue(exspectedNoSub, 0);
    appendValue(exspectedNoSub, 1);
    appendValue(exspectedNoSub, 2);
    appendValue(exspectedNoSub, 3);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, 1);
    appendValue(exspectedNoSub, 2);
    appendValue(exspectedNoSub, 3);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, 0);
    StringSet<int> res;
    resize(res, length(seqs));
    for (unsigned i = 0; i < length(seqs); ++i)
    {
        res[i] = findN(seqs[i], allowed);           //Test WITHOUT substitutions
    }
    for (unsigned i = 0; i < length(exspectedNoSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], res[i]);
    }

    StringSet<String<Dna5> > seqs2 = seqs;
    seqs2[1][4] = 'A';
    seqs2[2][0] = 'A';
    seqs2[2][31] = 'A';
    seqs2[3][0] = 'A';
    seqs2[3][12] = 'A';
    seqs2[3][31] = 'A';
    seqs2[4][0] = 'A';
    seqs2[4][12] = 'A';
    seqs2[4][20] = 'A';
    seqs2[5][0] = 'A';
    seqs2[5][1] = 'A';
    seqs2[5][2] = 'A';

     for (unsigned i = 0; i < length(seqs); ++i)
    {
        res[i] = findN(seqs[i], allowed, 'A');      //Test WITH substitutions
    }
    for (unsigned i = 0; i < 5; ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], res[i]);
    }
}

SEQAN_DEFINE_TEST(processN_test)
{
    GeneralStats stats;

    StringSet<String<Dna5> >seqs;
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(seqs, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(seqs, "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN");
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    StringSet<String<Dna5> >seqs2 = seqs;
    StringSet<String<char> > ids;
    appendValue(ids, "Null");
    appendValue(ids, "Eins");
    appendValue(ids, "Zwei");
    appendValue(ids, "Drei");
    appendValue(ids, "Vier");
    appendValue(ids, "Null2");
    StringSet<String<char> > ids2 = ids;
    unsigned allowed = 3;
    StringSet<String<Dna5> > exspectedNoSub;
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");

    processN(seqs2, ids2, allowed, stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedSeqs, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(exspectedNoSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], seqs2[i]);
    }

    stats.removedSeqs = 0;
    stats.uncalledBases = 0;
    seqs2 = seqs;
    ids2 = ids;
    Dna substitute = 'A';
    StringSet<String<Dna5> > exspectedSub = exspectedNoSub;
    exspectedSub[1][4] = substitute;
    exspectedSub[2][0] = substitute;
    exspectedSub[2][31] = substitute;
    exspectedSub[3][0] = substitute;
    exspectedSub[3][12] = substitute;
    exspectedSub[3][31] = substitute;

    processN(seqs2, ids2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedSeqs, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(exspectedSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedSub[i], seqs2[i]);
    }
}

SEQAN_DEFINE_TEST(processN_paired_test)
{
    GeneralStats stats;
    
    StringSet<String<Dna5> >seqs;
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(seqs, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(seqs, "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN");
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    StringSet<String<Dna5> >seqsRev;
    appendValue(seqsRev, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqsRev, "ATGNGNGGGTANCACGTGATCGTNACGTAGCANGC");
    appendValue(seqsRev, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(seqsRev, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(seqsRev, "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN");
    appendValue(seqsRev, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    StringSet<String<Dna5> >seqs2 = seqs;
    StringSet<String<Dna5> >seqsRev2 = seqsRev;
    StringSet<String<char> > ids;
    appendValue(ids, "Null");
    appendValue(ids, "Eins");
    appendValue(ids, "Zwei");
    appendValue(ids, "Drei");
    appendValue(ids, "Vier");
    appendValue(ids, "Null2");
    StringSet<String<char> > exspectedIds;
    appendValue(exspectedIds, "Null");
    appendValue(exspectedIds, "Null2"); //Due to swap action!
    appendValue(exspectedIds, "Zwei");
    appendValue(exspectedIds, "Drei");
   
    StringSet<String<char> > ids2 = ids;
    StringSet<String<char> > idsRev2 = ids;
    unsigned allowed = 3;
    StringSet<String<Dna5> > exspectedNoSub;
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC"); //Due to swap action!
    appendValue(exspectedNoSub, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
   
    processN(seqs2, ids2, seqsRev2, idsRev2, allowed, stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedSeqs, 4u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(exspectedNoSub); ++i)
    {
        SEQAN_ASSERT_EQ(ids2[i], exspectedIds[i]);
        SEQAN_ASSERT_EQ(exspectedNoSub[i], seqs2[i]);
        SEQAN_ASSERT_EQ(idsRev2[i], exspectedIds[i]);
    }

    stats.removedSeqs = 0;
    stats.uncalledBases = 0;
    seqs2 = seqs;
    seqsRev2 = seqsRev;
    ids2 = ids;
    idsRev2 = ids;
    Dna substitute = 'A';
    StringSet<String<Dna5> > exspectedSub = exspectedNoSub;
    exspectedSub[2][0] = substitute;
    exspectedSub[2][31] = substitute;
    exspectedSub[3][0] = substitute;
    exspectedSub[3][12] = substitute;
    exspectedSub[3][31] = substitute;

    processN(seqs2, ids2, seqsRev2, idsRev2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedSeqs, 4u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(exspectedSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedSub[i], seqs2[i]);
    }
}

SEQAN_DEFINE_TEST(processN_multiplex_test)
{
    GeneralStats stats;

    StringSet<String<Dna5> >multiplex;
    appendValue(multiplex, "ACTGTA");
    appendValue(multiplex, "TGACGT");
    appendValue(multiplex, "GTACGA");
    appendValue(multiplex, "GTACTG");
    appendValue(multiplex, "AAAAAA");
    appendValue(multiplex, "GGGTAC");
    StringSet<String<Dna5> >multiplex2 = multiplex;
    
    StringSet<String<Dna5> > exspectedMultiplex;
    appendValue(exspectedMultiplex, "ACTGTA");
    appendValue(exspectedMultiplex, "TGACGT");
    appendValue(exspectedMultiplex, "GTACGA");
    appendValue(exspectedMultiplex, "GTACTG");
    appendValue(exspectedMultiplex, "GGGTAC");

    StringSet<String<Dna5> >seqs;
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(seqs, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(seqs, "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN");
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    StringSet<String<Dna5> >seqs2 = seqs;
    StringSet<String<char> > ids;
    appendValue(ids, "Null");
    appendValue(ids, "Eins");
    appendValue(ids, "Zwei");
    appendValue(ids, "Drei");
    appendValue(ids, "Vier");
    appendValue(ids, "Null2");
    StringSet<String<char> > ids2 = ids;
    unsigned allowed = 3;
    StringSet<String<Dna5> > exspectedNoSub;
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");

    processN(seqs2, ids2, multiplex2, allowed, stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedSeqs, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(exspectedNoSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], seqs2[i]);
        SEQAN_ASSERT_EQ(exspectedMultiplex[i], multiplex2[i]);
    }

    stats.removedSeqs = 0;
    stats.uncalledBases = 0;
    seqs2 = seqs;
    ids2 = ids;
    multiplex2 = multiplex;
    Dna substitute = 'A';
    StringSet<String<Dna5> > exspectedSub = exspectedNoSub;
    exspectedSub[1][4] = substitute;
    exspectedSub[2][0] = substitute;
    exspectedSub[2][31] = substitute;
    exspectedSub[3][0] = substitute;
    exspectedSub[3][12] = substitute;
    exspectedSub[3][31] = substitute;

    processN(seqs2, ids2, multiplex2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedSeqs, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(exspectedSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedSub[i], seqs2[i]);
        SEQAN_ASSERT_EQ(exspectedMultiplex[i], multiplex2[i]);
    }
}

SEQAN_DEFINE_TEST(processN_paired_multiplex_test)
{
    GeneralStats stats;

    StringSet<String<Dna5Q> >multiplex;
    appendValue(multiplex, "ACTGTA");
    appendValue(multiplex, "TGACGT");
    appendValue(multiplex, "GTACGA");
    appendValue(multiplex, "GTACTG");
    appendValue(multiplex, "AAAAAA");
    appendValue(multiplex, "GGGTAC");
    StringSet<String<Dna5Q> >multiplex2 = multiplex;

    StringSet<String<Dna5Q> > exspectedMultiplex;
    appendValue(exspectedMultiplex, "ACTGTA");
    appendValue(exspectedMultiplex, "GGGTAC");
    appendValue(exspectedMultiplex, "GTACGA");
    appendValue(exspectedMultiplex, "GTACTG");

    StringSet<String<Dna5> >seqs;
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqs, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(seqs, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(seqs, "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN");
    appendValue(seqs, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    StringSet<String<Dna5> >seqsRev;
    appendValue(seqsRev, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(seqsRev, "ATGNGNGGGTANCACGTGATCGTNACGTAGCANGC");
    appendValue(seqsRev, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(seqsRev, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(seqsRev, "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN");
    appendValue(seqsRev, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    StringSet<String<Dna5> >seqs2 = seqs;
    StringSet<String<Dna5> >seqsRev2 = seqsRev;
    StringSet<String<char> > ids;
    appendValue(ids, "Null");
    appendValue(ids, "Eins");
    appendValue(ids, "Zwei");
    appendValue(ids, "Drei");
    appendValue(ids, "Vier");
    appendValue(ids, "Null2");
    StringSet<String<char> > exspectedIds;
    appendValue(exspectedIds, "Null");
    appendValue(exspectedIds, "Null2"); //Due to swap action!
    appendValue(exspectedIds, "Zwei");
    appendValue(exspectedIds, "Drei");
    StringSet<String<char> > ids2 = ids;
    StringSet<String<char> > idsRev2 = ids;
    unsigned allowed = 3;
    StringSet<String<Dna5> > exspectedNoSub;
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC"); //Due to swap action!
    appendValue(exspectedNoSub, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");

    processN(seqs2, ids2, seqsRev2, idsRev2, multiplex2, allowed, stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedSeqs, 4u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(exspectedNoSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], seqs2[i]);
        SEQAN_ASSERT_EQ(ids2[i], exspectedIds[i]);
        SEQAN_ASSERT_EQ(idsRev2[i], exspectedIds[i]);
        SEQAN_ASSERT_EQ(multiplex2[i], exspectedMultiplex[i]);
    }

    stats.removedSeqs = 0;
    stats.uncalledBases = 0;
    seqs2 = seqs;
    seqsRev2 = seqsRev;
    ids2 = ids;
    idsRev2 = ids;
    multiplex2 = multiplex;
    Dna substitute = 'A';
    StringSet<String<Dna5> > exspectedSub = exspectedNoSub;
    exspectedSub[2][0] = substitute;
    exspectedSub[2][31] = substitute;
    exspectedSub[3][0] = substitute;
    exspectedSub[3][12] = substitute;
    exspectedSub[3][31] = substitute;

    processN(seqs2, ids2, seqsRev2, idsRev2, multiplex2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedSeqs, 4u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(exspectedSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedSub[i], seqs2[i]);
        SEQAN_ASSERT_EQ(multiplex2[i], exspectedMultiplex[i]);
    }
}

SEQAN_DEFINE_TEST(preTrim_test)
{
    GeneralStats stats;

    StringSet<String<Dna> > seqs;
    appendValue(seqs, "ACGTAACTGA");
    appendValue(seqs, "AAAAAACTTTTT");
    appendValue(seqs, "AAAAAAG");
    appendValue(seqs, "TACGG");
    appendValue(seqs, "TAAAAAA");
    appendValue(seqs, "");
    appendValue(seqs, "GATTACAGATTACA");
    StringSet<String<Dna> > seqs2 = seqs;

    StringSet<String<char> > ids;
    appendValue(ids, "loeschenTrim");
    appendValue(ids, "Head/Tail");
    appendValue(ids, "Head");
    appendValue(ids, "loeschenNone");
    appendValue(ids, "Tail");
    appendValue(ids, "loeschenLeer");
    appendValue(ids, "None");
    StringSet<String<char> > ids2 = ids;

    preTrim(seqs2, ids2, 3, 3, 4, stats);
    SEQAN_ASSERT_EQ(seqs2[0], "TAAC");
    SEQAN_ASSERT_EQ(ids2[0], "loeschenTrim");
    SEQAN_ASSERT_EQ(seqs2[1], "AAACTT");
    SEQAN_ASSERT_EQ(ids2[1], "Head/Tail");
    SEQAN_ASSERT_EQ(seqs2[2], "TACAGATT");
    SEQAN_ASSERT_EQ(ids2[2], "None");
    SEQAN_ASSERT_EQ(length(ids2), 3u);
    SEQAN_ASSERT_EQ(length(seqs2), 3u);

    seqs2 = seqs;
    ids2 = ids;
    preTrim(seqs2, ids2, 6, 5, 1, stats);
    SEQAN_ASSERT_EQ(seqs2[1], "C");
    SEQAN_ASSERT_EQ(ids2[1], "Head/Tail");
    SEQAN_ASSERT_EQ(seqs2[0], "AGA");
    SEQAN_ASSERT_EQ(ids2[0], "None");
    SEQAN_ASSERT_EQ(length(ids2), 2u);
    SEQAN_ASSERT_EQ(length(seqs2), 2u);

    seqs2 = seqs;
    ids2 = ids;
    preTrim(seqs2, ids2, 6, 0, 1, stats);
    SEQAN_ASSERT_EQ(seqs2[2], "G");
    SEQAN_ASSERT_EQ(ids2[2], "Head");

    seqs2 = seqs;
    ids2 = ids;
    preTrim(seqs2, ids2, 0, 6, 1, stats);
    SEQAN_ASSERT_EQ(seqs2[4], "T");
    SEQAN_ASSERT_EQ(ids2[4], "Tail");

    seqs2 = seqs;
    ids2 = ids;
    preTrim(seqs2, ids2, 0, 0, 6, stats);
    SEQAN_ASSERT_EQ(seqs2[4], "TAAAAAA");
    SEQAN_ASSERT_EQ(ids2[4], "Tail");
    SEQAN_ASSERT_EQ(seqs2[3], "GATTACAGATTACA");
    SEQAN_ASSERT_EQ(ids2[3], "None");
    SEQAN_ASSERT_EQ(length(seqs2), 5u);
    SEQAN_ASSERT_EQ(length(ids2), 5u);
}

SEQAN_DEFINE_TEST(preTrim_paired_test)
{
    GeneralStats stats;

    StringSet<String<Dna> > seqs;
    appendValue(seqs, "ACGTAACTGA");
    appendValue(seqs, "AAAAAACTTTTT");
    appendValue(seqs, "AAAAAAG");
    appendValue(seqs, "TACGG");
    appendValue(seqs, "TAAAAAA");
    appendValue(seqs, "");
    appendValue(seqs, "GATTACAGATTACA");
    StringSet<String<Dna> > seqs2 = seqs;
    StringSet<String<Dna> > seqsRev = seqs;
    seqsRev[6] = "GTCA";
    StringSet<String<Dna> > seqsRev2 = seqsRev;

    StringSet<String<char> > ids;
    appendValue(ids, "loeschenTrim");
    appendValue(ids, "Head/Tail");
    appendValue(ids, "Head");
    appendValue(ids, "loeschenNone");
    appendValue(ids, "Tail");
    appendValue(ids, "loeschenLeer");
    appendValue(ids, "loeschenRev");
    StringSet<String<char> > ids2 = ids;
    StringSet<String<char> > idsRev = ids;
    StringSet<String<char> > idsRev2 = idsRev;

    preTrim(seqs2, ids2, seqsRev2, idsRev2, 3, 3, 4, stats);
    SEQAN_ASSERT_EQ(seqs2[0], "TAAC");
    SEQAN_ASSERT_EQ(ids2[0], "loeschenTrim");
    SEQAN_ASSERT_EQ(seqs2[1], "AAACTT");
    SEQAN_ASSERT_EQ(ids2[1], "Head/Tail");
    SEQAN_ASSERT_EQ(length(ids2), 2u);
    SEQAN_ASSERT_EQ(length(seqs2), 2u);

    seqs2 = seqs;
    seqsRev2 = seqsRev;
    ids2 = ids;
    idsRev2 = idsRev;
    preTrim(seqs2, ids2, seqsRev2, idsRev2, 6, 5, 1, stats);
    SEQAN_ASSERT_EQ(seqs2[0], "C");
    SEQAN_ASSERT_EQ(ids2[0], "Head/Tail");
    SEQAN_ASSERT_EQ(length(ids2), 1u);
    SEQAN_ASSERT_EQ(length(seqs2), 1u);

    seqs2 = seqs;
    seqsRev2 = seqsRev;
    ids2 = ids;
    idsRev2 = idsRev;
    preTrim(seqs2, ids2, seqsRev2, idsRev2, 6, 0, 1, stats);
    SEQAN_ASSERT_EQ(seqs2[2], "G");
    SEQAN_ASSERT_EQ(ids2[2], "Head");

    seqs2 = seqs;
    seqsRev2 = seqsRev;
    ids2 = ids;
    idsRev2 = idsRev;
    preTrim(seqs2, ids2, seqsRev2, idsRev2, 0, 6, 1, stats);
    SEQAN_ASSERT_EQ(seqs2[3], "T");
    SEQAN_ASSERT_EQ(ids2[3], "Tail");

    seqs2 = seqs;
    seqsRev2 = seqsRev;
    ids2 = ids;
    idsRev2 = idsRev;
    preTrim(seqs2, ids2, seqsRev2, idsRev2, 0, 0, 6, stats);
    SEQAN_ASSERT_EQ(seqs2[3], "TAAAAAA");
    SEQAN_ASSERT_EQ(ids2[3], "Tail");
    SEQAN_ASSERT_EQ(length(seqs2), 4u);
    SEQAN_ASSERT_EQ(length(ids2), 4u);
}

SEQAN_DEFINE_TEST(trimTo_test)
{
    GeneralStats stats;
    
    StringSet<String<char> > seqs;
    appendValue(seqs, "123456789");
    appendValue(seqs, "123456");
    appendValue(seqs, "1234567");
    appendValue(seqs, "");

    StringSet<String<char> > ids;
    appendValue(ids, "neun");
    appendValue(ids, "sechs");
    appendValue(ids, "sieben");
    appendValue(ids, "null");

    trimTo(seqs, ids, 7, stats);

    SEQAN_ASSERT_EQ(seqs[0], "1234567");
    SEQAN_ASSERT_EQ(seqs[1], "1234567");
    SEQAN_ASSERT_EQ(ids[0], "neun");
    SEQAN_ASSERT_EQ(ids[1], "sieben");
    SEQAN_ASSERT_EQ(length(seqs), 2u);
    SEQAN_ASSERT_EQ(length(ids), 2u);
}

SEQAN_DEFINE_TEST(trimTo_paired_test)
{
    GeneralStats stats;

    StringSet<String<char> > seqs;
    appendValue(seqs, "123456789");
    appendValue(seqs, "123456");
    appendValue(seqs, "1234567");
    appendValue(seqs, "");

    StringSet<String<char> > seqsRev = seqs;
    seqsRev[0] = "123456";

    StringSet<String<char> > ids;
    appendValue(ids, "neun");
    appendValue(ids, "fuenf");
    appendValue(ids, "sieben");
    appendValue(ids, "null");
    StringSet<String<char> > idsRev = ids;

    trimTo(seqs, ids, seqsRev, idsRev, 7, stats);

    SEQAN_ASSERT_EQ(seqs[0], "1234567");
    SEQAN_ASSERT_EQ(ids[0], "sieben");
    SEQAN_ASSERT_EQ(length(seqs), 1u);
    SEQAN_ASSERT_EQ(length(ids), 1u);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    int tnum = 1;
#ifdef _OPENMP
	omp_set_num_threads(8);
	tnum = omp_get_max_threads();
#endif
	std::cout<<"\nRunning Tests using " << tnum << " thread(s).\n\n";
    SEQAN_CALL_TEST(findN_test);
    SEQAN_CALL_TEST(processN_test);
    SEQAN_CALL_TEST(processN_paired_test);
    SEQAN_CALL_TEST(processN_multiplex_test);
    SEQAN_CALL_TEST(processN_paired_multiplex_test);
    SEQAN_CALL_TEST(preTrim_test);
    SEQAN_CALL_TEST(preTrim_paired_test);
    SEQAN_CALL_TEST(trimTo_test);
    SEQAN_CALL_TEST(trimTo_paired_test);
}
SEQAN_END_TESTSUITE
