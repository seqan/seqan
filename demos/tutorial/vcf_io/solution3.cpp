#include <seqan/vcf_io.h>

using namespace seqan;

int main()
{
    // Open output file.
    VcfFileOut out(std::cout, Vcf());

    // Fill sequence names.
    appendValue(contigNames(context(out)), "20");

    // Fill sample names.
    appendValue(sampleNames(context(out)), "NA00001");
    appendValue(sampleNames(context(out)), "NA00002");
    appendValue(sampleNames(context(out)), "NA00002");

    // Fill and write out headers - This is somewhat tedious.
    VcfHeader header;
    appendValue(header, VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(header, VcfHeaderRecord("fileDate", "20090805"));
    appendValue(header, VcfHeaderRecord("source", "myImputationProgramV3.1"));
    appendValue(header, VcfHeaderRecord("reference", "file:///seq/references/1000GenomesPilot-NCBI36.fasta"));
    appendValue(header, VcfHeaderRecord("contig", "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"));
    appendValue(header, VcfHeaderRecord("phasing", "partial"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"));
    appendValue(header, VcfHeaderRecord("FILTER", "<ID=q10,Description=\"Quality below 10\">"));
    appendValue(header, VcfHeaderRecord("FILTER", "<ID=s50,Description=\"Less than 50% of samples have data\">"));
    appendValue(header, VcfHeaderRecord("ID", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(header, VcfHeaderRecord("ID", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
    appendValue(header, VcfHeaderRecord("ID", "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"));
    appendValue(header, VcfHeaderRecord("ID", "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"));
    writeHeader(out, header);

    // Fill and write out the record.
    VcfRecord record;
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
    writeRecord(out, record);

    return 0;
}
