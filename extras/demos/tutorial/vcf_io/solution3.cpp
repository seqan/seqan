#include <seqan/basic.h>
#include <seqan/vcf_io.h>

#include <sstream>

int main()
{
    seqan::VcfStream out("-", seqan::VcfStream::WRITE);

    // Fill sequence names.
    appendValue(out.header.sequenceNames, "20");

    // Fill sample names.
    appendValue(out.header.sampleNames, "NA00001");
    appendValue(out.header.sampleNames, "NA00002");
    appendValue(out.header.sampleNames, "NA00002");

    // Write out headers.
    //
    // This is somewhat tedious.
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("fileDate", "20090805"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("source", "myImputationProgramV3.1"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("reference", "file:///seq/references/1000GenomesPilot-NCBI36.fasta"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("contig", "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("phasing", "partial"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("FILTER", "<ID=q10,Description=\"Quality below 10\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("FILTER", "<ID=s50,Description=\"Less than 50% of samples have data\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("ID", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("ID", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("ID", "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"));
    appendValue(out.header.headerRecords, seqan::VcfHeaderRecord("ID", "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"));

    // Write out the records.
    seqan::VcfRecord record;

    record.rID = 0;
    record.beginPos = 14369;
    record.id = "rs6054257";
    record.ref = "G";
    record.alt = "A";
    record.qual = 29;
    record.filter = "PASS";
    record.info = "NS=3;DP=14;AF=0.5;DB;H2";
    record.format = "GT:GQ:DP:HQ";
    appendValue(record.genotypeInfos, "0|0:48:1:51,51");
    appendValue(record.genotypeInfos, "1|0:48:8:51,51");
    appendValue(record.genotypeInfos, "1/1:43:5:.,.");
    if (writeRecord(out, record) != 0)
        std::cerr << "ERROR: Problem writing output file.";
    
    return 0;
}
