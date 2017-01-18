#include <iostream>
#include <seqan/basic.h>
#ifndef STDLIB_VS
#include <seqan/blast.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
      std::cerr << "USAGE: FILE_OUT\n";
      return 0;
    }

    typedef String<AminoAcid>               TSequence;
    typedef std::string                     TId;
    typedef Gaps<TSequence, ArrayGaps>      TGaps;
    typedef BlastMatch<TGaps, TGaps>        TBlastMatch;
    typedef BlastRecord<TBlastMatch>        TBlastRecord;
    typedef BlastIOContext<Blosum62>        TContext;

    StringSet<TSequence>    queries;
    StringSet<TSequence>    subjects;
    StringSet<TId>          qIds;
    StringSet<TId>          sIds;

    appendValue(queries, "VAYAQPRKLCYP");
    appendValue(queries, "NAYPRUTEIN");
    appendValue(queries, "AVITSFTQ");

    appendValue(subjects,
        "SSITEEKHIPHKEQDKDAEFLSKEALKTHMTENVLQMDRRAVQDPSTSFLQLLKAKGLLG"
        "LPDYEVNLADVNSPGFRKVAYAQTKPRRLCFPNGGTRRGSFIMDTAVVVMVSLRYVNIGK"
        "VIFPGATDVSEGEDEFWAGLPQAYGCLATEFLCIHIAIYSWIHVQSSRYDDMNASVIRAK"
        "LNLAVITSWTQLIQAEKETI");

    appendValue(subjects,
        "GATRDSKGNAVITSFTQARLRVYADLLGPYWIILHVIELTGVGNTGQKCTLNHMGTYAVF"
        "DLKQPPATNDLGLPKPCFIGFDIQNELAIGTVGHSEAVIAAFTQRDRLEERAESKQSLAR"
        "PVISPKLIAEVSTVLESALNQMYSSLGFYRVERAEDYAQPRKLCVVKKKSFNCLNADIWL"
        "EYRMEDQKSVPKVFKIMMDD");

    appendValue(qIds, "Query_Numero_Uno with args");
    appendValue(qIds, "Query_Numero_Dos with args");
    appendValue(qIds, "Query_Numero_Tres with args");

    appendValue(sIds, "Subject_Numero_Uno");
    appendValue(sIds, "Subject_Numero_Dos");

    BlastTabularFileOut<TContext> outfile(argv[1]);
    String<TBlastRecord> records;

    // protein vs protein search is BLASTP
    context(outfile).blastProgram = BlastProgram::BLASTP;

    // set gap parameters in blast notation
    setScoreGapOpenBlast(context(outfile).scoringScheme, -11);
    setScoreGapExtend(context(outfile).scoringScheme, -1);
    SEQAN_ASSERT(isValid(context(outfile).scoringScheme));

    // set the database properties in the context
    context(outfile).dbName = "The Foo Database";
    context(outfile).dbTotalLength = length(concat(subjects));
    context(outfile).dbNumberOfSeqs = length(subjects);

    writeHeader(outfile); // write file header

    for (unsigned q = 0; q < length(queries); ++q)
    {
        appendValue(records, TBlastRecord(qIds[q]));
        TBlastRecord & r = back(records);

        r.qLength = length(queries[q]);

        for (unsigned s = 0; s < length(subjects); ++s)
        {
            appendValue(r.matches, TBlastMatch(sIds[s]));
            TBlastMatch & m = back(records[q].matches);

            assignSource(m.alignRow0, queries[q]);
            assignSource(m.alignRow1, subjects[s]);

            localAlignment(m.alignRow0, m.alignRow1, seqanScheme(context(outfile).scoringScheme));

            m.qStart = beginPosition(m.alignRow0);
            m.qEnd   = endPosition(m.alignRow0);
            m.sStart = beginPosition(m.alignRow1);
            m.sEnd   = endPosition(m.alignRow1);

            m.sLength = length(subjects[s]);

            computeAlignmentStats(m, context(outfile));
            computeBitScore(m, context(outfile));
            computeEValue(m, r.qLength, context(outfile));

            if (m.eValue > 1)
                eraseBack(records[q].matches);
        }

        r.matches.sort(); // sort by bitscore

        writeRecord(outfile, r);
    }

    writeFooter(outfile);

    return 0;
}
#else
int main()
{
    std::cerr << "USAGE: FILE_OUT\n";
    return 0;
}
#endif
