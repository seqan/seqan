#include <iostream>
#include <seqan/basic.h>
#ifdef SEQAN_CXX11_COMPLETE
#include <seqan/blast.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    if (argc != 2)
        return 1;

    typedef String<AminoAcid>               TSequence;
    typedef std::string                     TId;
    typedef Align<TSequence, ArrayGaps>     TAlign;
    typedef BlastMatch<TAlign>              TBlastMatch;
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
            appendValue(r.matches, TBlastMatch(qIds[q], sIds[s]));
            TBlastMatch & m = back(records[q].matches);

            resize(rows(m.align), 2);
            assignSource(row(m.align, 0), queries[q]);
            assignSource(row(m.align, 1), subjects[s]);

            localAlignment(m.align, seqanScheme(context(outfile).scoringScheme));

            m.qStart = beginPosition(row(m.align, 0));
            m.qEnd   = endPosition(row(m.align, 0));
            m.sStart = beginPosition(row(m.align, 1));
            m.sEnd   = endPosition(row(m.align, 1));

            m.qLength = length(queries[q]);
            m.sLength = length(subjects[s]);

            computeAlignmentStats(m, context(outfile));
            computeBitScore(m, context(outfile));
            computeEValue(m, context(outfile));

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
    std::cerr << "Demo not run, because you don't have full C++11 support.\n";
    return 0;
}
#endif
