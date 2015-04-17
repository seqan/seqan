#include <iostream>
#include <seqan/blast.h>

using namespace seqan;

int main()
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;
    typedef BlastRecord<CharString, CharString, uint32_t, TAlign> TBlastRecord;
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

//     BlastTabularOut<TContext> outfile("/tmp/output.blast");
    BlastReportOut<TContext> outfile("/tmp/output.blast");

    Blosum62 scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);

    // upon assigning, this is also converted to SeqAn's scoring behaviour
    setBlastScoringScheme(context(outfile), scheme);

    // protein vs protein search is BLASTP
    setBlastProgram(context(outfile), BlastProgram::BLASTP);

    // set the database properties in the context
    context(outfile).dbName = "The Foo Database";
    context(outfile).dbTotalLength = length(concat(subjects));
    context(outfile).dbNumberOfSeqs = length(subjects);

    writeHeader(outfile); // write file header

    for (int q = 0; q < length(queries); ++q)
    {
        records[q].qId = qIds[q];
        records[q].qLength = length(queries[q]);

        for (int s = 0; s < length(subjects); ++s)
        {
            records[q].matches.emplace_back(qIds[q], sIds[s]);
            TBlastMatch & m = records[q].matches.back();

            resize(rows(m.align), 2);
            assignSource(row(m.align, 0), queries[q]);
            assignSource(row(m.align, 1), subjects[s]);

            localAlignment(m.align, context(outfile).scoringAdapter.scheme);

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