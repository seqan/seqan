// ==========================================================================
//                         Mason - A Read Simulator
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "mason_options.h"

// ----------------------------------------------------------------------------
// Function getYesNoStr()
// ----------------------------------------------------------------------------

char const * getYesNoStr(bool b)
{
    return b ? "YES" : "NO";
}

// ----------------------------------------------------------------------------
// Function getVerbosityStr()
// ----------------------------------------------------------------------------

char const * getVerbosityStr(int verbosity)
{
    switch (verbosity)
    {
        case 0:
            return "QUIET";
        case 1:
            return "NORMAL";
        case 2:
            return "VERBOSE";
        case 3:
            return "VERY VERBOSE";
        default:
            return "<invalid>";
    }
}

// ----------------------------------------------------------------------------
// Function getFragmentSizeModelStr()
// ----------------------------------------------------------------------------

char const * getFragmentSizeModelStr(FragmentSamplerOptions::FragmentSizeModel model)
{
    switch (model)
    {
        case FragmentSamplerOptions::UNIFORM:
            return "UNIFORM";
        case FragmentSamplerOptions::NORMAL:
            return "NORMAL";
        default:
            return "<invalid>";
    }
}

// ----------------------------------------------------------------------------
// Function getMateOrientationStr()
// ----------------------------------------------------------------------------

char const * getMateOrientationStr(SequencingOptions::MateOrientation orientation)
{
    switch (orientation)
    {
        case SequencingOptions::FORWARD_REVERSE:
            return "FORWARD-REVERSE (R1 --> <-- R2)";
        case SequencingOptions::REVERSE_FORWARD:
            return "REVERSE-FORWARD (R1 <-- --> R2)";
        case SequencingOptions::FORWARD_FORWARD:
            return "FORWARD-FORWARD (R1 --> --> R2)";
        case SequencingOptions::FORWARD_FORWARD2:
            return "FORWARD-FORWARD2 (R2 --> --> R1)";
        default:
            return "<invalid>";
    }
}

// ----------------------------------------------------------------------------
// Function getSourceStrandsStr()
// ----------------------------------------------------------------------------

char const * getSourceStrandsStr(SequencingOptions::SourceStrands strands)
{
    switch (strands)
    {
        case SequencingOptions::BOTH:
            return "BOTH";
        case SequencingOptions::FORWARD:
            return "FORWARD";
        case SequencingOptions::REVERSE:
            return "REVERSE";
        default:
            return "<invalid>";
    }
}

// ----------------------------------------------------------------------------
// Function getSequencingTechnologyStr()
// ----------------------------------------------------------------------------

char const * getSequencingTechnologyStr(SequencingOptions::SequencingTechnology technology)
{
    switch (technology)
    {
        case SequencingOptions::ILLUMINA:
            return "ILLUMINA";
        case SequencingOptions::ROCHE_454:
            return "ROCHE 454";
        case SequencingOptions::SANGER:
            return "SANGER";
        default:
            return "<invalid>";
    }
}

// ----------------------------------------------------------------------------
// Function getFragmentSizeModelStr()
// ----------------------------------------------------------------------------

char const * getFragmentSizeModelStr(Roche454SequencingOptions::ReadLengthModel model)
{
    switch (model)
    {
        case Roche454SequencingOptions::UNIFORM:
            return "UNIFORM";
        case Roche454SequencingOptions::NORMAL:
            return "NORMAL";
        default:
            return "<invalid>";
    }
}

// ----------------------------------------------------------------------------
// Function getBSSeqProtocolStr()
// ----------------------------------------------------------------------------

char const * getBSSeqProtocolStr(BSSeqOptions::BSProtocol protocol)
{
    switch (protocol)
    {
        case BSSeqOptions::DIRECTIONAL:
            return "DIRECTIONAL";
        case BSSeqOptions::UNDIRECTIONAL:
            return "UNDIRECTIONAL";
        default:
            return "<invalid>";
    }
}

// ----------------------------------------------------------------------------
// Function MethylationLevelSimulatorOptions::addOptions()
// ----------------------------------------------------------------------------

void MethylationLevelSimulatorOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Methylation Level Simulation");

    addOption(parser, seqan::ArgParseOption("", "methylation-levels", "Enable methylation level simulation."));

    addOption(parser, seqan::ArgParseOption("", "meth-cg-mu", "Median of beta distribution for methylation "
                                            "level of CpG loci.", seqan::ArgParseOption::DOUBLE, "MU"));
    setMinValue(parser, "meth-cg-mu", "0");
    setMaxValue(parser, "meth-cg-mu", "1");
    setDefaultValue(parser, "meth-cg-mu", "0.6");

    addOption(parser, seqan::ArgParseOption("", "meth-cg-sigma", "Standard deviation of beta distribution for "
                                            "methylation level of CpG loci.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "meth-cg-sigma", "0");
    setMaxValue(parser, "meth-cg-sigma", "1");
    setDefaultValue(parser, "meth-cg-sigma", "0.03");

    addOption(parser, seqan::ArgParseOption("", "meth-chg-mu", "Median of beta distribution for methylation "
                                            "level of CHG loci.", seqan::ArgParseOption::DOUBLE, "MU"));
    setMinValue(parser, "meth-chg-mu", "0");
    setMaxValue(parser, "meth-chg-mu", "1");
    setDefaultValue(parser, "meth-chg-mu", "0.08");

    addOption(parser, seqan::ArgParseOption("", "meth-chg-sigma", "Standard deviation of beta distribution for "
                                            "methylation level of CHG loci.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "meth-chg-sigma", "0");
    setMaxValue(parser, "meth-chg-sigma", "1");
    setDefaultValue(parser, "meth-chg-sigma", "0.008");

    addOption(parser, seqan::ArgParseOption("", "meth-chh-mu", "Median of beta distribution for methylation "
                                            "level of CHH loci.", seqan::ArgParseOption::DOUBLE, "MU"));
    setMinValue(parser, "meth-chh-mu", "0");
    setMaxValue(parser, "meth-chh-mu", "1");
    setDefaultValue(parser, "meth-chh-mu", "0.05");

    addOption(parser, seqan::ArgParseOption("", "meth-chh-sigma", "Standard deviation of beta distribution for "
                                            "methylation level of CHH loci.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "meth-chh-sigma", "0");
    setMaxValue(parser, "meth-chh-sigma", "1");
    setDefaultValue(parser, "meth-chh-sigma", "0.005");
}

// ----------------------------------------------------------------------------
// Function MethylationLevelSimulatorOptions::addTextSections()
// ----------------------------------------------------------------------------

void MethylationLevelSimulatorOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    (void)parser;
}

// ----------------------------------------------------------------------------
// Function MethylationLevelSimulatorOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MethylationLevelSimulatorOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    getOptionValue(simulateMethylationLevels, parser, "methylation-levels");
    getOptionValue(methMuCG, parser, "meth-cg-mu");
    getOptionValue(methSigmaCG, parser, "meth-cg-sigma");
    getOptionValue(methMuCHG, parser, "meth-chg-mu");
    getOptionValue(methSigmaCHG, parser, "meth-chg-sigma");
    getOptionValue(methMuCHH, parser, "meth-chh-mu");
    getOptionValue(methSigmaCHH, parser, "meth-chh-sigma");
}

// ----------------------------------------------------------------------------
// Function MethylationLevelSimulatorOptions::print()
// ----------------------------------------------------------------------------

void MethylationLevelSimulatorOptions::print(std::ostream & out) const
{
    out << "METHYLATION LEVELS OPTIONS\n"
        << "  VERBOSITY      \t" << getVerbosityStr(verbosity) << "\n"
        << "\n"
        << "  ENABLED\t" << getYesNoStr(simulateMethylationLevels) << "\n"
        << "\n"
        << "  MEDIAN CG\t" << methMuCG << "\n"
        << "  STDDEV CG\t" << methSigmaCG << "\n"
        << "  MEDIAN CHG\t" << methMuCHG << "\n"
        << "  STDDEV CHG\t" << methSigmaCHG << "\n"
        << "  MEDIAN CHH\t" << methMuCHH << "\n"
        << "  STDDEV CHH\t" << methSigmaCHH << "\n";
}

// ----------------------------------------------------------------------------
// Function BSSeqOptions::addOptions()
// ----------------------------------------------------------------------------

void BSSeqOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "BS-Seq Options");

    addOption(parser, seqan::ArgParseOption("", "enable-bs-seq", "Enable BS-seq simulation."));

    addOption(parser, seqan::ArgParseOption("", "bs-seq-protocol", "Protocol to use for BS-Seq simulation.",
                                            seqan::ArgParseOption::STRING, "PROTOCOL"));
    setDefaultValue(parser, "bs-seq-protocol", "directional");
    setValidValues(parser, "bs-seq-protocol", "directional undirectional");

    addOption(parser, seqan::ArgParseOption("", "bs-seq-conversion-rate", "Conversion rate for unmethylated Cs to "
                                            "become Ts.", seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "bs-seq-conversion-rate", "0");
    setMaxValue(parser, "bs-seq-conversion-rate", "1");
    setDefaultValue(parser, "bs-seq-conversion-rate", "0.99");
}

// ----------------------------------------------------------------------------
// Function BSSeqOptions::addTextSections()
// ----------------------------------------------------------------------------

void BSSeqOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    (void)parser;
}

// ----------------------------------------------------------------------------
// Function BSSeqOptions::getOptionValues()
// ----------------------------------------------------------------------------

void BSSeqOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    getOptionValue(bsSimEnabled, parser, "enable-bs-seq");
    getOptionValue(bsConversionRate, parser, "bs-seq-conversion-rate");
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "bs-seq-protocol");
    bsProtocol = (tmp == "undirectional") ? UNDIRECTIONAL : DIRECTIONAL;
}

// ----------------------------------------------------------------------------
// Function BSSeqOptions::print()
// ----------------------------------------------------------------------------

void BSSeqOptions::print(std::ostream & out) const
{
    out << "BS-SEQ OPTIONS\n"
        << "  VERBOSITY      \t" << getVerbosityStr(verbosity) << "\n"
        << "\n"
        << "  ENABLED\t" << getYesNoStr(bsSimEnabled) << "\n"
        << "  PROTOCOL\t" << getBSSeqProtocolStr(bsProtocol) << "\n"
        << "  CONVERSION RATE\t" << bsConversionRate << "\n";
}

// ----------------------------------------------------------------------------
// Function MaterializerOptions::addOptions()
// ----------------------------------------------------------------------------

void MaterializerOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Apply VCF Variants to Reference");

    addOption(parser, seqan::ArgParseOption("ir", "input-reference", "Path to FASTA file to read the reference from.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN.fa"));
    setValidValues(parser, "input-reference", seqan::SeqFileIn::getFileExtensions());
    setRequired(parser, "input-reference");

    addOption(parser, seqan::ArgParseOption("iv", "input-vcf", "Path to the VCF file with variants to apply.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN.vcf"));
    setValidValues(parser, "input-vcf", seqan::VcfFileIn::getFileExtensions());
}

// ----------------------------------------------------------------------------
// Function MaterializerOptions::addTextSections()
// ----------------------------------------------------------------------------

void MaterializerOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    addTextSection(parser, "VCF Variant Notes");

    addText(parser,
            "If the option \\fB--input-vcf\\fP/\\fB-iv\\fP is given then the given VCF file is read and the variants "
            "are applied to the input reference file.  If it is not given then the input reference file is taken "
            "verbatimly for simulating reads.");
    addText(parser,
            "There are some restrictions on the VCF file and the application of the variants to the reference will "
            "fail if the VCF file is non-conforming.  VCF files from the \\fImason_variator\\fP program are "
            "guaranteed to be read.");
     addText(parser,
             "Only the haplotypes of the first individual will be generated.");
}

// ----------------------------------------------------------------------------
// Function MaterializerOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MaterializerOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    getOptionValue(fastaFileName, parser, "input-reference");
    getOptionValue(vcfFileName, parser, "input-vcf");
}

// ----------------------------------------------------------------------------
// Function MaterializerOptions::print()
// ----------------------------------------------------------------------------

void MaterializerOptions::print(std::ostream & out) const
{
    out << "MATERIALIZER OPTIONS\n"
        << "  VERBOSITY         \t" << getVerbosityStr(verbosity) << "\n"
        << "  REFERENCE FASTA   \t" << fastaFileName << "\n"
        << "  VARIANTS VCF      \t" << vcfFileName << "\n";
}

// ----------------------------------------------------------------------------
// Function FragmentSamplerOptions::addOptions()
// ----------------------------------------------------------------------------

void FragmentSamplerOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Fragment Size (Insert Size) Options");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-model", "The model to use for the fragment size simulation.",
                                            seqan::ArgParseOption::STRING, "MODEL"));
    setDefaultValue(parser, "fragment-size-model", "normal");
    setValidValues(parser, "fragment-size-model", "normal uniform");

    addOption(parser, seqan::ArgParseOption("", "fragment-min-size", "Smallest fragment size to use when using "
                                            "uniform fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-min-size", "100");
    setMinValue(parser, "fragment-min-size", "1");

    addOption(parser, seqan::ArgParseOption("", "fragment-max-size", "Largest fragment size to use when using "
                                            "uniform fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-max-size", "400");
    setMinValue(parser, "fragment-max-size", "1");

    addOption(parser, seqan::ArgParseOption("", "fragment-mean-size", "Mean fragment size for normally distributed "
                                            "fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-mean-size", "300");
    setMinValue(parser, "fragment-mean-size", "1");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-std-dev", "Fragment size standard deviation when using "
                                            "normally distributed fragment size simulation.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setDefaultValue(parser, "fragment-size-std-dev", "30");
    setMinValue(parser, "fragment-size-std-dev", "1");
}

// ----------------------------------------------------------------------------
// Function FragmentSamplerOptions::addTextSections()
// ----------------------------------------------------------------------------

void FragmentSamplerOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    addTextSection(parser, "Fragment Size (Insert Size) Simulation");

    addText(parser,
            "You can choose between a normal and a uniform distribution of fragment lengths.  When sequencing "
            "these fragments from both sides in a paired protocol, the fragment size will become the insert "
            "size.");
}

// ----------------------------------------------------------------------------
// Function FragmentSamplerOptions::getOptionValues()
// ----------------------------------------------------------------------------

void FragmentSamplerOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "fragment-size-model");
    model = (tmp == "normal") ? FragmentSamplerOptions::NORMAL : FragmentSamplerOptions::UNIFORM;

    getOptionValue(minFragmentSize, parser, "fragment-min-size");
    getOptionValue(maxFragmentSize, parser, "fragment-max-size");
    getOptionValue(meanFragmentSize, parser, "fragment-mean-size");
    getOptionValue(stdDevFragmentSize, parser, "fragment-size-std-dev");
}

// ----------------------------------------------------------------------------
// Function FragmentSamplerOptions::print()
// ----------------------------------------------------------------------------

void FragmentSamplerOptions::print(std::ostream & out) const
{
    out << "FRAGMENT SAMPLING\n"
        << "  SIZE MODEL     \t" << getFragmentSizeModelStr(model) << "\n"
        << "  MIN SIZE       \t" << minFragmentSize << "\n"
        << "  MAX SIZE       \t" << maxFragmentSize << "\n"
        << "  MEAN SIZE      \t" << meanFragmentSize << "\n"
        << "  STD DEV SIZE   \t" << stdDevFragmentSize << "\n";
}

// ----------------------------------------------------------------------------
// Function SequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void SequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Global Read Simulation Options");

    addOption(parser, seqan::ArgParseOption("", "seq-technology", "Set sequencing technology to simulate.",
                                            seqan::ArgParseOption::STRING, "TECH"));
    setValidValues(parser, "seq-technology", "illumina 454 sanger");
    setDefaultValue(parser, "seq-technology", "illumina");

    addOption(parser, seqan::ArgParseOption("", "seq-mate-orientation", "Orientation for paired reads.  See section Read "
                                            "Orientation below.", seqan::ArgParseOption::STRING, "ORIENTATION"));
    setValidValues(parser, "seq-mate-orientation", "FR RF FF FF2");
    setDefaultValue(parser, "seq-mate-orientation", "FR");

    addOption(parser, seqan::ArgParseOption("", "seq-strands", "Strands to simulate from, only applicable to paired "
                                            "sequencing simulation.", seqan::ArgParseOption::STRING, "STRAND"));
    setValidValues(parser, "seq-strands", "forward reverse both");
    setDefaultValue(parser, "seq-strands", "both");

    addOption(parser, seqan::ArgParseOption("", "embed-read-info", "Whether or not to embed read information."));

    addOption(parser, seqan::ArgParseOption("", "read-name-prefix", "Read names will have this prefix.",
                                            seqan::ArgParseOption::STRING, "STR"));
    setDefaultValue(parser, "read-name-prefix", "simulated.");

    // Add options for nested options structs.
    bsSeqOptions.addOptions(parser);
}

// ----------------------------------------------------------------------------
// Function SequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void SequencingOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    addTextSection(parser, "Sequencing Simulation");

    addText(parser,
            "Simulation of base qualities is disabled when writing out FASTA files.  Simulation of paired-end "
            "sequencing is enabled when specifying two output files.");

    addTextSection(parser, "Read Orientation");
    addText(parser,
            "You can use the \\fB--mate-orientation\\fP to set the relative orientation when doing paired-end "
            "sequencing.  The valid values are given in the following.");
    addListItem(parser, "FR",
                "Reads are inward-facing, the same as Illumina paired-end reads: R1 --> <-- R2.");
    addListItem(parser, "RF",
                "Reads are outward-facing, the same as Illumina mate-pair reads: R1 <-- --> R2.");
    addListItem(parser, "FF",
                "Reads are on the same strand: R1 --> --> R2.");
    addListItem(parser, "FF2",
                "Reads are on the same strand but the \"right\" reads are sequenced to the left of the \"left\" reads, "
                "same as 454 paired: R2 --> --> R1.");

    // Add text sections for nested options structs.
    bsSeqOptions.addTextSections(parser);
}

// ----------------------------------------------------------------------------
// Function SequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void SequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;

    getOptionValue(tmp, parser, "seq-mate-orientation");
    if (tmp == "FR")
        mateOrientation = FORWARD_REVERSE;
    else if (tmp == "RF")
        mateOrientation = REVERSE_FORWARD;
    else if (tmp == "FF")
        mateOrientation = FORWARD_FORWARD;
    else if (tmp == "FF2")
        mateOrientation = FORWARD_FORWARD2;

    getOptionValue(tmp, parser, "seq-strands");
    if (tmp == "forward")
        strands = FORWARD;
    else if (tmp == "reverse")
        strands = REVERSE;
    else
        strands = BOTH;

    getOptionValue(tmp, parser, "seq-technology");
    if (tmp == "illumina")
        sequencingTechnology = SequencingOptions::ILLUMINA;
    else if (tmp == "454")
        sequencingTechnology = SequencingOptions::ROCHE_454;
    else
        sequencingTechnology = SequencingOptions::SANGER;

    getOptionValue(embedReadInfo, parser, "embed-read-info");
    getOptionValue(readNamePrefix, parser, "read-name-prefix");

    // Get option values for nested options.
    bsSeqOptions.getOptionValues(parser);
}

// ----------------------------------------------------------------------------
// Function SequencingOptions::print()
// ----------------------------------------------------------------------------

void SequencingOptions::print(std::ostream & out) const
{
    out << "SEQUENCING\n"
        << "  SIMULATE QUALITIES \t" << getYesNoStr(simulateQualities) << "\n"
        << "  SIMULATE MATE PAIRS\t" << getYesNoStr(simulateMatePairs) << "\n"
        << "  MATE ORIENTATION   \t" << getMateOrientationStr(mateOrientation) << "\n"
        << "  SOURCE STRANDS     \t" << getSourceStrandsStr(strands) << "\n"
        << "  SEQUENCING TECH    \t" << getSequencingTechnologyStr(sequencingTechnology) << "\n"
        << "\n";

    // Print options for nested options.
    bsSeqOptions.print(out);
}

// ----------------------------------------------------------------------------
// Function IlluminaSequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void IlluminaSequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Illumina Options");

    addOption(parser, seqan::ArgParseOption("", "illumina-read-length", "Read length for Illumina simulation.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "illumina-read-length", "1");
    setDefaultValue(parser, "illumina-read-length", "100");

    addOption(parser, seqan::ArgParseOption("", "illumina-error-profile-file",
                                            "Path to file with Illumina error profile.  The file must be a text file "
                                            "with floating point numbers separated by space, each giving a positional "
                                            "error rate.",
                                            seqan::ArgParseOption::INPUT_FILE, "FILE"));
    setValidValues(parser, "illumina-error-profile-file", "txt");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-insert",
                                            "Insert per-base probability for insertion in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-insert", "0");
    setMaxValue(parser, "illumina-prob-insert", "1");
    setDefaultValue(parser, "illumina-prob-insert", "0.00005");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-deletion",
                                            "Insert per-base probability for deletion in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-deletion", "0");
    setMaxValue(parser, "illumina-prob-deletion", "1");
    setDefaultValue(parser, "illumina-prob-deletion", "0.00005");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch-scale",
                                            "Scaling factor for Illumina mismatch probability.",
                                            seqan::ArgParseOption::DOUBLE, "FACTOR"));
    setMinValue(parser, "illumina-prob-mismatch-scale", "0");
    setDefaultValue(parser, "illumina-prob-mismatch-scale", "1.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch",
                                            "Average per-base mismatch probability in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-mismatch", "0.0");
    setMaxValue(parser, "illumina-prob-mismatch", "1.0");
    setDefaultValue(parser, "illumina-prob-mismatch", "0.004");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch-begin",
                                            "Per-base mismatch probability of first base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-mismatch-begin", "0.0");
    setMaxValue(parser, "illumina-prob-mismatch-begin", "1.0");
    setDefaultValue(parser, "illumina-prob-mismatch-begin", "0.002");

    addOption(parser, seqan::ArgParseOption("", "illumina-prob-mismatch-end",
                                            "Per-base mismatch probability of last base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "illumina-prob-mismatch-end", "0.0");
    setMaxValue(parser, "illumina-prob-mismatch-end", "1.0");
    setDefaultValue(parser, "illumina-prob-mismatch-end", "0.012");

    addOption(parser, seqan::ArgParseOption("", "illumina-position-raise",
                                            "Point where the error curve raises in relation to read length.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "illumina-position-raise", "0.0");
    setMaxValue(parser, "illumina-position-raise", "1.0");
    setDefaultValue(parser, "illumina-position-raise", "0.66");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-mean-begin",
                                            "Mean PHRED quality for non-mismatch bases of first base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-mean-begin", "40.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-mean-end",
                                            "Mean PHRED quality for non-mismatch bases of last base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-mean-end", "39.5");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-stddev-begin",
                                            "Standard deviation of PHRED quality for non-mismatch bases of first base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-stddev-begin", "0.05");

    addOption(parser, seqan::ArgParseOption("", "illumina-quality-stddev-end",
                                            "Standard deviation of PHRED quality for non-mismatch bases of last base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-quality-stddev-end", "10.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-mean-begin",
                                            "Mean PHRED quality for mismatch bases of first base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-mean-begin", "40.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-mean-end",
                                            "Mean PHRED quality for mismatch bases of last base in Illumina sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-mean-end", "30.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-stddev-begin",
                                            "Standard deviation of PHRED quality for mismatch bases of first base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-stddev-begin", "3.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-mismatch-quality-stddev-end",
                                            "Standard deviation of PHRED quality for mismatch bases of last base "
                                            "in Illumina sequencing.", seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "illumina-mismatch-quality-stddev-end", "15.0");

    addOption(parser, seqan::ArgParseOption("", "illumina-left-template-fastq",
                                            "FASTQ file to use for a template for left-end reads.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN.fq"));
    setValidValues(parser, "illumina-left-template-fastq", seqan::SeqFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("", "illumina-right-template-fastq",
                                            "FASTQ file to use for a template for right-end reads.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN.fq"));
    setValidValues(parser, "illumina-right-template-fastq", seqan::SeqFileIn::getFileExtensions());
}

// ----------------------------------------------------------------------------
// Function IlluminaSequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void IlluminaSequencingOptions::addTextSections(seqan::ArgumentParser & /*parser*/) const
{}

// ----------------------------------------------------------------------------
// Function IlluminaSequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void IlluminaSequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    getOptionValue(readLength, parser, "illumina-read-length");

    getOptionValue(probabilityMismatchFile, parser, "illumina-error-profile-file");

    getOptionValue(probabilityInsert, parser, "illumina-prob-insert");
    getOptionValue(probabilityDelete, parser, "illumina-prob-deletion");

    getOptionValue(probabilityMismatchScale, parser, "illumina-prob-mismatch-scale");

    getOptionValue(probabilityMismatch, parser, "illumina-prob-mismatch");
    getOptionValue(probabilityMismatchBegin, parser, "illumina-prob-mismatch-begin");
    getOptionValue(probabilityMismatchEnd, parser, "illumina-prob-mismatch-end");
    getOptionValue(positionRaise, parser, "illumina-position-raise");

    getOptionValue(meanQualityBegin, parser, "illumina-quality-mean-begin");
    getOptionValue(meanQualityEnd, parser, "illumina-quality-mean-end");
    getOptionValue(stdDevQualityBegin, parser, "illumina-quality-stddev-begin");
    getOptionValue(stdDevQualityEnd, parser, "illumina-quality-stddev-end");

    getOptionValue(meanMismatchQualityBegin, parser, "illumina-mismatch-quality-mean-begin");
    getOptionValue(meanMismatchQualityEnd, parser, "illumina-mismatch-quality-mean-end");
    getOptionValue(stdDevMismatchQualityBegin, parser, "illumina-mismatch-quality-stddev-begin");
    getOptionValue(stdDevMismatchQualityEnd, parser, "illumina-mismatch-quality-stddev-end");

    getOptionValue(leftTemplateFastq, parser, "illumina-left-template-fastq");
    getOptionValue(rightTemplateFastq, parser, "illumina-right-template-fastq");
}
// ----------------------------------------------------------------------------
// Function IlluminaSequencingOptions::print()
// ----------------------------------------------------------------------------

void IlluminaSequencingOptions::print(std::ostream & out) const
{
    out << "ILLUMINA SEQUENCING\n"
        << "  READ LENGTH                  \t" << readLength << "\n"
        << "  ERROR PROFILE FILE           \t" << probabilityMismatchFile << "\n"
        << "\n"
        << "  PROBABILITY INSERT           \t" << probabilityInsert << "\n"
        << "  PROBABILITY DELETE           \t" << probabilityDelete << "\n"
        << "\n"
        << "  PROBABILITY MISMATCH SCALE   \t" << probabilityMismatchScale << "\n"
        << "  PROBABILITY MISMATCH AVG     \t" << probabilityMismatch << "\n"
        << "  PROBABILITY MISMATCH BEGIN   \t" << probabilityMismatchBegin << "\n"
        << "  PROBABILITY MISMATCH END     \t" << probabilityMismatchEnd << "\n"
        << "  MISMATCH RAISE POINT         \t" << positionRaise << "\n"
        << "\n"
        << "  MEAN QUALITY BEGIN           \t" << meanQualityBegin << "\n"
        << "  MEAN QUALITY END             \t" << meanQualityEnd << "\n"
        << "  STD DEV QUALITY BEGIN        \t" << stdDevQualityBegin << "\n"
        << "  STD DEV QUALITY END          \t" << stdDevQualityEnd << "\n"
        << "\n"
        << "  MEAN MISMATCH QUALITY BEGIN  \t" << meanMismatchQualityBegin << "\n"
        << "  MEAN MISMATCH QUALITY END    \t" << meanMismatchQualityEnd << "\n"
        << "  STDDEV MISMATCH QUALITY BEGIN\t" << stdDevMismatchQualityBegin << "\n"
        << "  STDDEV MISMATCH QUALITY END  \t" << stdDevMismatchQualityEnd << "\n"
        << "\n"
        << "  LEFT TEMPLATE FASTQ          \t" << leftTemplateFastq << "\n"
        << "  RIGHT TEMPLATE FASTQ         \t" << rightTemplateFastq << "\n";
}

// ----------------------------------------------------------------------------
// Function SangerSequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void SangerSequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "Sanger Sequencing Options");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-model", "The model to use for sampling the Sanger "
                                            "read length.", seqan::ArgParseOption::STRING, "MODEL"));
    setValidValues(parser, "sanger-read-length-model", "normal uniform");
    setDefaultValue(parser, "sanger-read-length-model", "normal");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-min", "The minimal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "sanger-read-length-min", "0");
    setDefaultValue(parser, "sanger-read-length-min", "400");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-max", "The maximal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "sanger-read-length-max", "0");
    setDefaultValue(parser, "sanger-read-length-max", "600");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-mean", "The mean read length when the read length is "
                                            "sampled with normal distribution.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "sanger-read-length-mean", "0");
    setDefaultValue(parser, "sanger-read-length-mean", "400");

    addOption(parser, seqan::ArgParseOption("", "sanger-read-length-error", "The read length standard deviation when the "
                                            "read length is sampled uniformly.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "sanger-read-length-error", "0");
    setDefaultValue(parser, "sanger-read-length-error", "40");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-mismatch-scale",
                                            "Scaling factor for Sanger mismatch probability.",
                                            seqan::ArgParseOption::DOUBLE, "FACTOR"));
    setMinValue(parser, "sanger-prob-mismatch-scale", "0");
    setDefaultValue(parser, "sanger-prob-mismatch-scale", "1.0");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-mismatch-begin",
                                            "Per-base mismatch probability of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-mismatch-begin", "0.0");
    setMaxValue(parser, "sanger-prob-mismatch-begin", "1.0");
    setDefaultValue(parser, "sanger-prob-mismatch-begin", "0.005");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-mismatch-end",
                                            "Per-base mismatch probability of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-mismatch-end", "0.0");
    setMaxValue(parser, "sanger-prob-mismatch-end", "1.0");
    setDefaultValue(parser, "sanger-prob-mismatch-end", "0.001");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-insertion-begin",
                                            "Per-base insertion probability of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-insertion-begin", "0.0");
    setMaxValue(parser, "sanger-prob-insertion-begin", "1.0");
    setDefaultValue(parser, "sanger-prob-insertion-begin", "0.0025");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-insertion-end",
                                            "Per-base insertion probability of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-insertion-end", "0.0");
    setMaxValue(parser, "sanger-prob-insertion-end", "1.0");
    setDefaultValue(parser, "sanger-prob-insertion-end", "0.005");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-deletion-begin",
                                            "Per-base deletion probability of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-deletion-begin", "0.0");
    setMaxValue(parser, "sanger-prob-deletion-begin", "1.0");
    setDefaultValue(parser, "sanger-prob-deletion-begin", "0.0025");

    addOption(parser, seqan::ArgParseOption("", "sanger-prob-deletion-end",
                                            "Per-base deletion probability of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "PROB"));
    setMinValue(parser, "sanger-prob-deletion-end", "0.0");
    setMaxValue(parser, "sanger-prob-deletion-end", "1.0");
    setDefaultValue(parser, "sanger-prob-deletion-end", "0.005");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-start-mean",
                                            "Mean PHRED quality for non-mismatch bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-start-mean", "40.0");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-end-mean",
                                            "Mean PHRED quality for non-mismatch bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-end-mean", "39.5");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-start-stddev",
                                            "Mean PHRED quality for non-mismatch bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-start-stddev", "0.1");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-match-end-stddev",
                                            "Mean PHRED quality for non-mismatch bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-match-end-stddev", "2");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-start-mean",
                                            "Mean PHRED quality for erroneous bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-start-mean", "30");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-end-mean",
                                            "Mean PHRED quality for erroneous bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-end-mean", "20");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-start-stddev",
                                            "Mean PHRED quality for erroneous bases of first base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-start-stddev", "2");

    addOption(parser, seqan::ArgParseOption("", "sanger-quality-error-end-stddev",
                                            "Mean PHRED quality for erroneous bases of last base in Sanger sequencing.",
                                            seqan::ArgParseOption::DOUBLE, "QUAL"));
    setDefaultValue(parser, "sanger-quality-error-end-stddev", "5");
}

// ----------------------------------------------------------------------------
// Function SangerSequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void SangerSequencingOptions::addTextSections(seqan::ArgumentParser & /*parser*/) const
{}

// ----------------------------------------------------------------------------
// Function SangerSequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void SangerSequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "sanger-read-length-model");
    readLengthIsUniform = (tmp == "uniform");

    getOptionValue(readLengthMean, parser, "sanger-read-length-mean");
    getOptionValue(readLengthError, parser, "sanger-read-length-error");
    getOptionValue(readLengthMin, parser, "sanger-read-length-min");
    getOptionValue(readLengthMax, parser, "sanger-read-length-max");

    getOptionValue(probabilityMismatchBegin, parser, "sanger-prob-mismatch-begin");
    getOptionValue(probabilityMismatchEnd, parser, "sanger-prob-mismatch-end");
    getOptionValue(probabilityInsertBegin, parser, "sanger-prob-insertion-begin");
    getOptionValue(probabilityInsertEnd, parser, "sanger-prob-insertion-end");
    getOptionValue(probabilityDeleteBegin, parser, "sanger-prob-deletion-begin");
    getOptionValue(probabilityDeleteEnd, parser, "sanger-prob-deletion-end");

    getOptionValue(qualityMatchStartMean, parser, "sanger-quality-match-start-mean");
    getOptionValue(qualityMatchEndMean, parser, "sanger-quality-match-end-mean");
    getOptionValue(qualityMatchStartStdDev, parser, "sanger-quality-match-start-stddev");
    getOptionValue(qualityMatchEndStdDev, parser, "sanger-quality-match-end-stddev");
    getOptionValue(qualityErrorStartMean, parser, "sanger-quality-error-start-mean");
    getOptionValue(qualityErrorEndMean, parser, "sanger-quality-error-end-mean");
    getOptionValue(qualityErrorStartStdDev, parser, "sanger-quality-error-start-stddev");
    getOptionValue(qualityErrorEndStdDev, parser, "sanger-quality-error-end-stddev");
}

// ----------------------------------------------------------------------------
// Function SangerSequencingOptions::print()
// ----------------------------------------------------------------------------

void SangerSequencingOptions::print(std::ostream & out) const
{
    out << "SANGER SEQUENCING\n"
        << "  UNIFORM READ LENGTH       \t" << getYesNoStr(readLengthIsUniform) << "\n"
        << "  READ LENGTH MEAN          \t" << readLengthMean << "\n"
        << "  READ LENGTH ERROR         \t" << readLengthError << "\n"
        << "  READ LENGTH MIN           \t" << readLengthMin << "\n"
        << "  READ LENGTH MAX           \t" << readLengthMax << "\n"
        << "\n"
        << "  PROBABILITY MISMATCH BEGIN\t" << probabilityMismatchBegin << "\n"
        << "  PROBABILITY MISMATCH END  \t" << probabilityMismatchEnd << "\n"
        << "  PROBABILITY INSERT BEGIN  \t" << probabilityInsertBegin << "\n"
        << "  PROBABILITY INSERT END    \t" << probabilityInsertEnd << "\n"
        << "  PROBABILITY DELETE BEGIN  \t" << probabilityDeleteBegin << "\n"
        << "  PROBABILITY DELETE END    \t" << probabilityDeleteEnd << "\n"
        << "\n"
        << "  QUALITY MATCH START MEAN  \t" << qualityMatchStartMean << "\n"
        << "  QUALITY MATCH END MEAN    \t" << qualityMatchEndMean << "\n"
        << "  QUALITY MATCH START STDDEV\t" << qualityMatchStartStdDev << "\n"
        << "  QUALITY MATCH END STDDEV  \t" << qualityMatchEndStdDev << "\n"
        << "\n"
        << "  QUALITY ERROR START MEAN  \t" << qualityErrorStartMean << "\n"
        << "  QUALITY ERROR END MEAN    \t" << qualityErrorEndMean << "\n"
        << "  QUALITY ERROR START STDDEV\t" << qualityErrorStartStdDev << "\n"
        << "  QUALITY ERROR END STDDEV  \t" << qualityErrorEndStdDev << "\n";
}

// ----------------------------------------------------------------------------
// Function Roche454SequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void Roche454SequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
    addSection(parser, "454 Sequencing Options");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-model", "The model to use for sampling the 454 read length.",
                                            seqan::ArgParseOption::STRING, "MODEL"));
    setValidValues(parser, "454-read-length-model", "normal uniform");
    setDefaultValue(parser, "454-read-length-model", "normal");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-min", "The minimal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "454-read-length-min", "0");
    setDefaultValue(parser, "454-read-length-min", "10");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-max", "The maximal read length when the read length is "
                                            "sampled uniformly.", seqan::ArgParseOption::INTEGER, "LEN"));
    setMinValue(parser, "454-read-length-max", "0");
    setDefaultValue(parser, "454-read-length-max", "600");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-mean", "The mean read length when the read length is "
                                            "sampled with normal distribution.", seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "454-read-length-mean", "0");
    setDefaultValue(parser, "454-read-length-mean", "400");

    addOption(parser, seqan::ArgParseOption("", "454-read-length-stddev", "The read length standard deviation when the "
                                            "read length is sampled with normal distribution.",
                                            seqan::ArgParseOption::DOUBLE, "LEN"));
    setMinValue(parser, "454-read-length-stddev", "0");
    setDefaultValue(parser, "454-read-length-stddev", "40");

    addOption(parser, seqan::ArgParseOption("", "454-no-sqrt-in-std-dev", "For error model, if set then "
                                            "(sigma = k * r)) is used, otherwise (sigma = k * sqrt(r))."));

    addOption(parser, seqan::ArgParseOption("", "454-proportionality-factor", "Proportionality factor for calculating the standard "
                                            "deviation proportional to the read length.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "454-proportionality-factor", "0");
    setDefaultValue(parser, "454-proportionality-factor", "0.15");

    addOption(parser, seqan::ArgParseOption("", "454-background-noise-mean", "Mean of lognormal distribution to use for "
                                            "the noise.", seqan::ArgParseOption::DOUBLE, "MEAN"));
    setMinValue(parser, "454-background-noise-mean", "0");
    setDefaultValue(parser, "454-background-noise-mean", "0.23");

    addOption(parser, seqan::ArgParseOption("", "454-background-noise-stddev", "Standard deviation of lognormal "
                                            "distribution to use for the noise.", seqan::ArgParseOption::DOUBLE,
                                            "SIGMA"));
    setMinValue(parser, "454-background-noise-stddev", "0");
    setDefaultValue(parser, "454-background-noise-stddev", "0.15");
}

// ----------------------------------------------------------------------------
// Function Roche454SequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void Roche454SequencingOptions::addTextSections(seqan::ArgumentParser & /*parser*/) const
{}

// ----------------------------------------------------------------------------
// Function Roche454SequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void Roche454SequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    seqan::CharString tmp;
    getOptionValue(tmp, parser, "454-read-length-model");
    lengthModel = (tmp == "uniform") ? Roche454SequencingOptions::UNIFORM : Roche454SequencingOptions::NORMAL;

    getOptionValue(minReadLength, parser, "454-read-length-min");
    getOptionValue(maxReadLength, parser, "454-read-length-max");
    getOptionValue(meanReadLength, parser, "454-read-length-mean");
    getOptionValue(stdDevReadLength, parser, "454-read-length-stddev");
    sqrtInStdDev = !isSet(parser, "454-no-sqrt-in-std-dev");
    getOptionValue(k, parser, "454-proportionality-factor");
    getOptionValue(backgroundNoiseMean, parser, "454-background-noise-mean");
    getOptionValue(backgroundNoiseStdDev, parser, "454-background-noise-stddev");
}

// ----------------------------------------------------------------------------
// Function Roche454SequencingOptions::print()
// ----------------------------------------------------------------------------

void Roche454SequencingOptions::print(std::ostream & out) const
{
    out << "ROCHE 454 SEQUENCING\n"
        << "  READ LENGTH MODEL \t" << getFragmentSizeModelStr(lengthModel) << "\n"
        << "  MIN READ LENGTH   \t" << minReadLength << "\n"
        << "  MAX READ LENGTH   \t" << maxReadLength << "\n"
        << "  MEAN READ LENGTH  \t" << meanReadLength << "\n"
        << "  STDDEV READ LENGTH\t" << stdDevReadLength << "\n";
}

// ----------------------------------------------------------------------------
// Function MasonSimulatorOptions::addOptions()
// ----------------------------------------------------------------------------

void MasonSimulatorOptions::addOptions(seqan::ArgumentParser & parser) const
{
    // Add top-level options.

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Low verbosity."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    addOption(parser, seqan::ArgParseOption("", "seed", "Seed to use for random number generator.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("", "meth-seed", "Seed to use for methylation level random number "
                                            "generator.", seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "meth-seed", "0");

    addOption(parser, seqan::ArgParseOption("", "seed-spacing", "Offset for seeds to use when multi-threading.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "seed-spacing", "2048");

    addOption(parser, seqan::ArgParseOption("", "num-threads", "Number of threads to use.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "num-threads", "1");
    setDefaultValue(parser, "num-threads", "1");

    addOption(parser, seqan::ArgParseOption("", "force-single-end", "Force single-end simulation although --out-right "
                                            "file is given."));

    addOption(parser, seqan::ArgParseOption("", "chunk-size", "Number of fragments to simulate in one batch.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "chunk-size", "65536");
    setDefaultValue(parser, "chunk-size", "65536");

    addOption(parser, seqan::ArgParseOption("n", "num-fragments", "Number of reads/pairs to simulate.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setRequired(parser, "num-fragments");
    setMinValue(parser, "num-fragments", "1");

    addOption(parser, seqan::ArgParseOption("", "meth-fasta-in", "FASTA file with methylation levels of the input file.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN"));
    setValidValues(parser, "meth-fasta-in", seqan::SeqFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("o", "out", "Output of single-end/left end reads.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
    setRequired(parser, "out");
    setValidValues(parser, "out", seqan::SeqFileOut::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("or", "out-right", "Output of right reads.  Giving this options enables "
                                            "paired-end simulation.", seqan::ArgParseOption::OUTPUT_FILE, "OUT2"));
    setValidValues(parser, "out-right", seqan::SeqFileOut::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("oa", "out-alignment", "SAM/BAM file with alignments.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
    setValidValues(parser, "out-alignment", seqan::BamFileOut::getFileExtensions());

    // Add options of the component options.
    matOptions.addOptions(parser);
    methOptions.addOptions(parser);
    fragSamplerOptions.addOptions(parser);
    seqOptions.addOptions(parser);
    illuminaOptions.addOptions(parser);
    sangerOptions.addOptions(parser);
    rocheOptions.addOptions(parser);
}

// ----------------------------------------------------------------------------
// Function MasonSimulatorOptions::addTextSections()
// ----------------------------------------------------------------------------

void MasonSimulatorOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    // Add top-level text sections.
    addTextSection(parser, "Simulation Overview");

    addText(parser,
            "The first step is the application of VCF variants to the input reference file.");

    addText(parser,
            "After the generation of the haplotypes, fragments are sampled from the sequence.  These fragments "
            "correspond to the fragments in the real preparation step for the sequencing.  They are later sequenced "
            "from one or both sides depending on whether a single-end or a paired protocol is used.");

    addTextSection(parser, "Important Parameters");

    addText(parser, "For most users, the following options are most important.");

    addListItem(parser, "Paired-End Simulation",
                "Use --fragment-length-model to switch between normally and uniformly distributed insert sizes. "
                "Use the --fragment-* options for configuring the insert size simulation.");

    addTextSection(parser, "Multi-Threading");

    addText(parser,
            "When using multi-threading, each thread gets its own random number generator (RNG).  The RNG of thread "
            "i is initialized with the value of \\fB--seed\\fP plus i.");

    addTextSection(parser, "BAM/SAM Tags");
    addText(parser,
            "Mason can write out a BAM or SAM file with alignments of the reads against the reference.  The records "
            "have tags that give information about the simulated reads.  Below is a list of the tags and their meaning.");

    addListItem(parser, "NM", "Edit distance when aligned to the reference (i).");
    addListItem(parser, "MD", "String for mismatching positions (Z).");

    addListItem(parser, "oR", "Name of \\fBo\\fPriginal \\fBr\\fPeference, (Z).");
    addListItem(parser, "oH", "Number of the \\fBo\\fPriginal \\fBh\\fPhaplotype (1-based), (i).");
    addListItem(parser, "oP", "\\fBo\\fPriginal \\fBp\\fPosition on the original reference (i).");
    addListItem(parser, "oS", "\\fBo\\fPriginal \\fBs\\fPtrand, \\fIF/R\\fP for forward and reverse strand (A).");
    addListItem(parser, "uR",
                "Reason for being unaligned, \\fII/B\\fP for being in insertion or spanning over breakpoint.");

    addListItem(parser, "XE", "Number of sequencing \\fIe\\fPrrors in the read (i).");
    addListItem(parser, "XS", "Number of \\fIS\\fPNPs in the read alignment (i).");
    addListItem(parser, "XI", "Number of small \\fIi\\fPndels in the read alignment (i).");

    // Add text sections of the component options.
    matOptions.addTextSections(parser);
    methOptions.addTextSections(parser);
    fragSamplerOptions.addTextSections(parser);
    seqOptions.addTextSections(parser);
    illuminaOptions.addTextSections(parser);
    sangerOptions.addTextSections(parser);
    rocheOptions.addTextSections(parser);
}

// ----------------------------------------------------------------------------
// Function MasonSimulatorOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MasonSimulatorOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    // Get top-level options.
    if (isSet(parser, "quiet"))
        verbosity = 0;
    if (isSet(parser, "verbose"))
        verbosity = 2;
    if (isSet(parser, "very-verbose"))
        verbosity = 3;
    getOptionValue(seed, parser, "seed");
    getOptionValue(methSeed, parser, "meth-seed");
    getOptionValue(seedSpacing, parser, "seed-spacing");
#if SEQAN_HAS_OPENMP
    getOptionValue(numThreads, parser, "num-threads");
#else  // #if SEQAN_HAS_OPENMP
    numThreads = 1;
#endif  // #if SEQAN_HAS_OPENMP
    getOptionValue(chunkSize, parser, "chunk-size");
    getOptionValue(numFragments, parser, "num-fragments");
    getOptionValue(forceSingleEnd, parser, "force-single-end");
    getOptionValue(methFastaInFile, parser, "meth-fasta-in");
    getOptionValue(outFileNameLeft, parser, "out");
    getOptionValue(outFileNameRight, parser, "out-right");
    getOptionValue(outFileNameSam, parser, "out-alignment");

    // Get options for the other components that we use.
    methOptions.getOptionValues(parser);
    matOptions.getOptionValues(parser);
    fragSamplerOptions.getOptionValues(parser);
    seqOptions.getOptionValues(parser);
    illuminaOptions.getOptionValues(parser);
    sangerOptions.getOptionValues(parser);
    rocheOptions.getOptionValues(parser);

    // Copy in the verbosity flag into the component options.
    matOptions.verbosity = verbosity;
    methOptions.verbosity = verbosity;
    fragSamplerOptions.verbosity = verbosity;
    seqOptions.verbosity = verbosity;
    illuminaOptions.verbosity = verbosity;
    sangerOptions.verbosity = verbosity;
    rocheOptions.verbosity = verbosity;

    // Configure simulation of pairs and mates depending on output files.
    seqOptions.simulateQualities = (endsWith(outFileNameLeft, ".fastq") || endsWith(outFileNameLeft, ".fq"));
    seqOptions.simulateMatePairs = !forceSingleEnd && !empty(outFileNameRight);
    methOptions.simulateMethylationLevels = !empty(methFastaInFile);
}

// ----------------------------------------------------------------------------
// Function MasonSimulatorOptions::print()
// ----------------------------------------------------------------------------

void MasonSimulatorOptions::print(std::ostream & out) const
{
    out << "MASON OPTIONS\n"
        << "-------------\n"
        << "\n"
        << "VERBOSITY\t" << getVerbosityStr(verbosity) << "\n"
        << "\n"
        << "SEED\t" << seed << "\n"
        << "METHYLATION SEED\t" << methSeed << "\n"
        << "SEED SPACING\t" << seedSpacing << "\n"
        << "\n"
        << "FORCE SINGLE END\t" << getYesNoStr(forceSingleEnd) << "\n"
        << "NUM FRAGMENTS\t" << numFragments << "\n"
        << "\n"
        << "NUM THREADS\t" << numThreads << "\n"
        << "CHUNK SIZE\t" << chunkSize << "\n"
        << "\n"
        << "METHYLATION FASTA IN\t" << methFastaInFile << "\n"
        << "OUTPUT FILE LEFT\t" << outFileNameLeft << "\n"
        << "OUTPUT FILE RIGHT\t" << outFileNameRight << "\n"
        << "PAIRED END SIMULATION\t" << getYesNoStr(!forceSingleEnd && !empty(outFileNameRight)) << "\n"
        << "\n";
    matOptions.print(out);
    out << "\n";
    methOptions.print(out);
    out << "\n";
    fragSamplerOptions.print(out);
    out << "\n";
    seqOptions.print(out);
    out << "\n";
    illuminaOptions.print(out);
    out << "\n";
    sangerOptions.print(out);
    out << "\n";
    rocheOptions.print(out);
}

// ----------------------------------------------------------------------------
// Function MasonMaterializerOptions::addOptions()
// ----------------------------------------------------------------------------

void MasonMaterializerOptions::addOptions(seqan::ArgumentParser & parser) const
{
    // Add top-level options.

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Low verbosity."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    addOption(parser, seqan::ArgParseOption("", "seed", "Seed for random number generation.",
                                            seqan::ArgParseOption::INTEGER, "Int"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("", "meth-seed", "Seed for methylation simulation random number generation.",
                                            seqan::ArgParseOption::INTEGER, "Int"));
    setDefaultValue(parser, "meth-seed", "0");

    addOption(parser, seqan::ArgParseOption("o", "out", "Output of materialized contigs.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
    setRequired(parser, "out");
    setValidValues(parser, "out", seqan::SeqFileOut::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("", "out-breakpoints", "TSV file to write breakpoints in variants to.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "TSV"));
    setValidValues(parser, "out-breakpoints", "tsv txt");

    addOption(parser, seqan::ArgParseOption("", "haplotype-name-sep",
                                            "String separating contig name from haplotype number.",
                                            seqan::ArgParseOption::STRING, "SEP"));
    setDefaultValue(parser, "haplotype-name-sep", "/");

    addOption(parser, seqan::ArgParseOption("", "meth-fasta-in", "FASTA file with methylation levels of the input file.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN"));
    setValidValues(parser, "meth-fasta-in", seqan::SeqFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("", "meth-fasta-out", "FASTA file with methylation levels of the output file.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
    setValidValues(parser, "meth-fasta-out", seqan::SeqFileOut::getFileExtensions());

    // Add options of the component options.
    matOptions.addOptions(parser);
    methOptions.addOptions(parser);
}

// ----------------------------------------------------------------------------
// Function MasonMaterializerOptions::addTextSections()
// ----------------------------------------------------------------------------

void MasonMaterializerOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    // Add text sections of the component options.
    matOptions.addTextSections(parser);
    methOptions.addTextSections(parser);
}

// ----------------------------------------------------------------------------
// Function MasonMaterializerOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MasonMaterializerOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    // Get top-level options.
    if (isSet(parser, "quiet"))
        verbosity = 0;
    if (isSet(parser, "verbose"))
        verbosity = 2;
    if (isSet(parser, "very-verbose"))
        verbosity = 3;
    getOptionValue(seed, parser, "seed");
    getOptionValue(methSeed, parser, "meth-seed");
    getOptionValue(outputFileName, parser, "out");
    getOptionValue(outputBreakpointFile, parser, "out-breakpoints");
    getOptionValue(haplotypeNameSep, parser, "haplotype-name-sep");
    getOptionValue(methFastaInFile, parser, "meth-fasta-in");
    getOptionValue(methFastaOutFile, parser, "meth-fasta-out");

    // Get options for the other components that we use.
    matOptions.getOptionValues(parser);
    methOptions.getOptionValues(parser);

    // Copy in the verbosity flag into the component options.
    matOptions.verbosity = verbosity;
    methOptions.verbosity = verbosity;

    // Set methylation simulation enabled flag.
    methOptions.simulateMethylationLevels = !empty(methFastaInFile) && !empty(methFastaOutFile);
}

// ----------------------------------------------------------------------------
// Function MasonMaterializerOptions::print()
// ----------------------------------------------------------------------------

void MasonMaterializerOptions::print(std::ostream & out) const
{
    out << "MASON MATERIALIZER OPTIONS\n"
        << "--------------------------\n"
        << "\n"
        << "VERBOSITY               \t" << getVerbosityStr(verbosity) << "\n"
        << "\n"
        << "SEED                    \t" << seed << "\n"
        << "METHYLATION SEED        \t" << methSeed << "\n"
        << "\n"
        << "OUTPUT FILE             \t" << outputFileName << "\n"
        << "BREAKPOINT TSV OUT      \t" << outputBreakpointFile << "\n"
        << "METHYLATION LEVEL INPUT \t" << methFastaInFile << "\n"
        << "METHYLATION LEVEL OUTPUT\t" << methFastaOutFile << "\n"
        << "\n"
        << "HAPLOTYPE NAME SEP      \t" << haplotypeNameSep << "\n"
        << "\n";
    matOptions.print(out);
    out << "\n";
    methOptions.print(out);
    out << "\n";
}

// ----------------------------------------------------------------------------
// Function MasonSplicingOptions::addOptions()
// ----------------------------------------------------------------------------

void MasonSplicingOptions::addOptions(seqan::ArgumentParser & parser) const
{
    // Add top-level options.

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Low verbosity."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    addOption(parser, seqan::ArgParseOption("", "seed", "Seed for random number generation.",
                                            seqan::ArgParseOption::INTEGER, "Int"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("o", "out", "Output of materialized contigs.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
    setRequired(parser, "out");
    setValidValues(parser, "out", seqan::SeqFileOut::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("", "haplotype-name-sep",
                                            "String separating contig name from haplotype number.",
                                            seqan::ArgParseOption::STRING, "SEP"));
    setDefaultValue(parser, "haplotype-name-sep", "/");

    addOption(parser, seqan::ArgParseOption("ig", "in-gff", "Path to input GFF or GTF file, must be "
                                            "sorted by reference name.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN.gff"));
    setValidValues(parser, "in-gff", seqan::GffFileIn::getFileExtensions());
    setRequired(parser, "in-gff");

    addOption(parser, seqan::ArgParseOption("", "gff-type", "Splicing will filter to the records that have this type.",
                                            seqan::ArgParseOption::INPUT_FILE, "TYPE"));
    setDefaultValue(parser, "gff-type", "exon");

    addOption(parser, seqan::ArgParseOption("", "gff-group-by", "Assign features to their parent using the tag "
                                            "with this name.", seqan::ArgParseOption::INPUT_FILE, "KEY"));
    setDefaultValue(parser, "gff-group-by", "Parent");

    // Add options of the component options.
    matOptions.addOptions(parser);
}

// ----------------------------------------------------------------------------
// Function MasonSplicingOptions::addTextSections()
// ----------------------------------------------------------------------------

void MasonSplicingOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    // Add text sections of the component options.
    matOptions.addTextSections(parser);
}

// ----------------------------------------------------------------------------
// Function MasonSplicingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MasonSplicingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    // Get top-level options.
    if (isSet(parser, "quiet"))
        verbosity = 0;
    if (isSet(parser, "verbose"))
        verbosity = 2;
    if (isSet(parser, "very-verbose"))
        verbosity = 3;
    getOptionValue(seed, parser, "seed");
    getOptionValue(outputFileName, parser, "out");
    getOptionValue(haplotypeNameSep, parser, "haplotype-name-sep");
    getOptionValue(inputGffFile, parser, "in-gff");
    getOptionValue(gffType, parser, "gff-type");
    getOptionValue(gffGroupBy, parser, "gff-group-by");

    // Get options for the other components that we use.
    matOptions.getOptionValues(parser);

    // Copy in the verbosity flag into the component options.
    matOptions.verbosity = verbosity;
}

// ----------------------------------------------------------------------------
// Function MasonSplicingOptions::print()
// ----------------------------------------------------------------------------

void MasonSplicingOptions::print(std::ostream & out) const
{
    out << "MASON SPLICING OPTIONS\n"
        << "-----------------------\n"
        << "\n"
        << "VERBOSITY               \t" << getVerbosityStr(verbosity) << "\n"
        << "\n"
        << "SEED                    \t" << seed << "\n"
        << "\n"
        << "OUTPUT FILE             \t" << outputFileName << "\n"
        << "\n"
        << "HAPLOTYPE NAME SEP      \t" << haplotypeNameSep << "\n"
        << "\n"
        << "INPUT GFF FILE          \t" << inputGffFile << "\n"
        << "GFF TYPE                \t" << gffType << "\n"
        << "GFF GROUP BY            \t" << gffGroupBy << "\n"
        << "\n";
    matOptions.print(out);
    out << "\n";
}

// ----------------------------------------------------------------------------
// Function MasonFragmentSequencingOptions::addOptions()
// ----------------------------------------------------------------------------

void MasonFragmentSequencingOptions::addOptions(seqan::ArgumentParser & parser) const
{
    // Add top-level options.

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Low verbosity."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    addOption(parser, seqan::ArgParseOption("", "seed", "Seed to use for random number generator.",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("i", "in", "Path to input file.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN"));
    setRequired(parser, "in");
    setValidValues(parser, "in", seqan::SeqFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("o", "out", "Output of single-end/left end reads.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
    setRequired(parser, "out");
    setValidValues(parser, "out", seqan::SeqFileOut::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("or", "out-right", "Output of right reads.  Giving this options enables "
                                            "paired-end simulation.", seqan::ArgParseOption::OUTPUT_FILE, "OUT2"));
    setValidValues(parser, "out-right", seqan::SeqFileOut::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("", "force-single-end", "Force single-end simulation although --out-right "
                                            "is given."));

    // Add options of the component options.
    seqOptions.addOptions(parser);
    illuminaOptions.addOptions(parser);
    sangerOptions.addOptions(parser);
    rocheOptions.addOptions(parser);
}

// ----------------------------------------------------------------------------
// Function MasonFragmentSequencingOptions::addTextSections()
// ----------------------------------------------------------------------------

void MasonFragmentSequencingOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    // Add text sections of the component options.
    seqOptions.addTextSections(parser);
    illuminaOptions.addTextSections(parser);
    sangerOptions.addTextSections(parser);
    rocheOptions.addTextSections(parser);
}

// ----------------------------------------------------------------------------
// Function MasonFragmentSequencingOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MasonFragmentSequencingOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    // Get top-level options.
    if (isSet(parser, "quiet"))
        verbosity = 0;
    if (isSet(parser, "verbose"))
        verbosity = 2;
    if (isSet(parser, "very-verbose"))
        verbosity = 3;
    getOptionValue(seed, parser, "seed");
    getOptionValue(inputFileName, parser, "in");
    getOptionValue(outFileNameLeft, parser, "out");
    getOptionValue(outFileNameRight, parser, "out-right");
    getOptionValue(forceSingleEnd, parser, "force-single-end");

    // Get options for the other components that we use.
    seqOptions.getOptionValues(parser);
    illuminaOptions.getOptionValues(parser);
    sangerOptions.getOptionValues(parser);
    rocheOptions.getOptionValues(parser);

    // Copy in the verbosity flag into the component options.
    seqOptions.verbosity = verbosity;
    illuminaOptions.verbosity = verbosity;
    sangerOptions.verbosity = verbosity;
    rocheOptions.verbosity = verbosity;

    // Configure simulation of pairs and mates depending on output files.
    seqOptions.simulateQualities = (endsWith(outFileNameLeft, ".fastq") || endsWith(outFileNameLeft, ".fq"));
    seqOptions.simulateMatePairs = !forceSingleEnd && !empty(outFileNameRight);
}

// ----------------------------------------------------------------------------
// Function MasonFragmentSequencingOptions::print()
// ----------------------------------------------------------------------------

void MasonFragmentSequencingOptions::print(std::ostream & out) const
{
    out << "MASON FRAGMENT SEQUENCING OPTIONS\n"
        << "---------------------------------\n"
        << "\n"
        << "VERBOSITY        \t" << getVerbosityStr(verbosity) << "\n"
        << "\n"
        << "SEED             \t" << seed << "\n"
        << "\n"
        << "INPUT FILE NAME  \t" << inputFileName << "\n"
        << "\n"
        << "OUTPUT FILE LEFT \t" << outFileNameLeft << "\n"
        << "OUTPUT FILE RIGHT\t" << outFileNameRight << "\n"
        << "\n";
    seqOptions.print(out);
    out << "\n";
    illuminaOptions.print(out);
    out << "\n";
    sangerOptions.print(out);
    out << "\n";
    rocheOptions.print(out);
}

// ----------------------------------------------------------------------------
// Function MasonMethylationOptions::addOptions()
// ----------------------------------------------------------------------------

void MasonMethylationOptions::addOptions(seqan::ArgumentParser & parser) const
{
    // Add top-level options.

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Low verbosity."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Higher verbosity."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Highest verbosity."));

    addOption(parser, seqan::ArgParseOption("", "seed", "Seed for RNG.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "seed", "0");

    addOption(parser, seqan::ArgParseOption("i", "in", "Input FASTA file with genome.",
                                            seqan::ArgParseOption::INPUT_FILE, "IN.fa"));
    setRequired(parser, "in");
    setValidValues(parser, "in", seqan::SeqFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("o", "out", "Input FASTA file with genome.",
                                            seqan::ArgParseOption::INPUT_FILE, "OUT.fa"));
    setRequired(parser, "out");
    setValidValues(parser, "out", seqan::SeqFileOut::getFileExtensions());

    // Add options of the component options.
    methOptions.addOptions(parser);
}

// ----------------------------------------------------------------------------
// Function MasonMethylationOptions::addTextSections()
// ----------------------------------------------------------------------------

void MasonMethylationOptions::addTextSections(seqan::ArgumentParser & parser) const
{
    // Add text sections of the component options.
    methOptions.addTextSections(parser);
}

// ----------------------------------------------------------------------------
// Function MasonMethylationOptions::getOptionValues()
// ----------------------------------------------------------------------------

void MasonMethylationOptions::getOptionValues(seqan::ArgumentParser const & parser)
{
    // Get top-level options.
    if (isSet(parser, "quiet"))
        verbosity = 0;
    if (isSet(parser, "verbose"))
        verbosity = 2;
    if (isSet(parser, "very-verbose"))
        verbosity = 3;

    getOptionValue(seed, parser, "seed");
    getOptionValue(fastaInFile, parser, "in");
    getOptionValue(methFastaOutFile, parser, "out");

    // Get options for the other components that we use.
    methOptions.getOptionValues(parser);

    // Copy in the verbosity flag into the component options.
    methOptions.verbosity = verbosity;
}

// ----------------------------------------------------------------------------
// Function MasonMethylationOptions::print()
// ----------------------------------------------------------------------------

void MasonMethylationOptions::print(std::ostream & out) const
{
    out << "MASON METHYLATION OPTIONS\n"
        << "-------------------------\n"
        << "\n"
        << "VERBOSITY          \t" << getVerbosityStr(verbosity) << "\n"
        << "\n"
        << "SEED               \t" << seed << "\n"
        << "\n"
        << "FASTA IN FILE      \t" << fastaInFile << "\n"
        << "METH FASTA OUT FILE\t" << methFastaOutFile << "\n"
        << "\n";
    methOptions.print(out);
    out << "\n";
}
