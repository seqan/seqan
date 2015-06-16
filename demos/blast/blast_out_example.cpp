#include <iostream>
#include <seqan/basic.h>
#ifdef SEQAN_CXX11_COMPLETE
#include <seqan/blast.h>

using namespace seqan;

int main()
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<TAlign> TBlastMatch;
    typedef BlastRecord<TBlastMatch> TBlastRecord;
    typedef BlastIOContext<Blosum62> TContext;

    StringSet<String<AminoAcid>, Owner<ConcatDirect<>>> queries;
    StringSet<CharString, Owner<ConcatDirect<>>> qIds;
    StringSet<String<AminoAcid>> subjects;
    StringSet<CharString> sIds;

    resize(subjects, 2);
    resize(sIds, 2);

    subjects[0] =
    "SSITEEKHIPHKEQDKDAEFLSKEALKTHMTENVLQMDRRAVQDPSTSFLQLLKAKGLLG"
    "LPDYEVNLADVNSPGFRKVAYAQTKPRRLCFPNGGTRRGSFIMDTAVVVMVSLRYVNIGK"
    "VIFPGATDVSEGEDEFWAGLPQAYGCLATEFLCIHIAIYSWIHVQSSRYDDMNASVIRAK"
    "LNLAVITSWTQLIQAEKETI";

    subjects[1] =
    "GATRDSKGNAVITSFTQARLRVYADLLGPYWIILHVIELTGVGNTGQKCTLNHMGTYAVF"
    "DLKQPPATNDLGLPKPCFIGFDIQNELAIGTVGHSEAVIAAFTQRDRLEERAESKQSLAR"
    "PVISPKLIAEVSTVLESALNQMYSSLGFYRVERAEDYAQPRKLCVVKKKSFNCLNADIWL"
    "EYRMEDQKSVPKVFKIMMDD";

    sIds[0] = "Subject_Numero_Uno";
    sIds[1] = "Subject_Numero_Dos";

    appendValue(queries, "VAYAQPRKLCYP");
    appendValue(queries, "NAYPRUTEIN");
    appendValue(queries, "AVITSFTQ");

    appendValue(qIds, "Query_Numero_Uno with args");
    appendValue(qIds, "Query_Numero_Dos with args");
    appendValue(qIds, "Query_Numero_Tres with args");

    String<TBlastRecord> records;
    resize(records, length(queries)); // always one record for every query!

//     BlastTabularFileOut<TContext> outfile("/tmp/output.m9");
    BlastReportFileOut<TContext> outfile("/tmp/output.m0");

    // set gap parameters in blast notation
    setScoreGapOpenBlast(context(outfile).scoringScheme, -11);
    setScoreGapExtend(context(outfile).scoringScheme, -1);

    // protein vs protein search is BLASTP
    context(outfile).blastProgram = BlastProgram::BLASTP;

    // set the database properties in the context
    context(outfile).dbName = "The Foo Database";
    context(outfile).dbTotalLength = length(concat(subjects));
    context(outfile).dbNumberOfSeqs = length(subjects);

    writeHeader(outfile); // write file header

    for (unsigned q = 0; q < length(queries); ++q)
    {
        records[q].qId = qIds[q];
        records[q].qLength = length(queries[q]);

        for (unsigned s = 0; s < length(subjects); ++s)
        {
            records[q].matches.emplace_back(qIds[q], sIds[s]);
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

            if (m.eValue > 1) // or other cutoff
            {
                records[q].matches.pop_back();
                continue;
            }
        }

        records[q].matches.sort(); // sort by bitscore

        writeRecord(outfile, records[q]);
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
