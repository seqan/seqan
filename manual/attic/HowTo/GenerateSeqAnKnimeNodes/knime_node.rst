::

    #cpp
    // ==========================================================================
    //                                 knime_node
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
    // Author: Your Name <your.email@example.net>
    // ==========================================================================

    #include <seqan/basic.h>
    #include <seqan/sequence.h>
    #include <seqan/seq_io.h>

    #include <seqan/arg_parse.h>

    // ==========================================================================
    // Classes
    // ==========================================================================

    // --------------------------------------------------------------------------
    // Class KnimeNodeOptions
    // --------------------------------------------------------------------------

    // This struct stores the options from the command line.
    //
    // You might want to rename this to reflect the name of your app.

    struct KnimeNodeOptions
    {
        // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
        int verbosity;

        // The arguments of the program are stored here.
        seqan::CharString inputFile;
        seqan::CharString outputFile;

        KnimeNodeOptions() :
            verbosity(1)
        {}
    };

    // ==========================================================================
    // Functions
    // ==========================================================================

    // --------------------------------------------------------------------------
    // Function parseCommandLine()
    // --------------------------------------------------------------------------

    seqan::ArgumentParser::ParseResult
    parseCommandLine(KnimeNodeOptions & options, int argc, char const ** argv)
    {
        // Setup ArgumentParser.
        seqan::ArgumentParser parser("knime_node");
        // Set short description, version, and date.
        setShortDescription(parser, "This is a very simple KNIME node providing an input and output port.");
        setVersion(parser, "0.1");
        setDate(parser, "Sep 2013");

        // Define usage line and long description.
        addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
        addDescription(parser, "This is a very simple KNIME node providing an input and output port. The code should be modified such that the node does something useful");

        // We require one argument.
        addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "IN"));
        setValidValues(parser, 0, "fastq fq");

        addOption(parser, seqan::ArgParseOption("o", "outputFile", "Name of the multi-FASTA output.", seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
        setValidValues(parser, "outputFile", "fastq fq");
        setDefaultValue(parser, "outputFile", "result.fastq");
        setRequired(parser, "o");

        // The verbosity option should be used to help debugging
        addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
        addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
        addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

        // Add Examples Section.
        addTextSection(parser, "Examples");
        addListItem(parser, "\\fBknime_node\\fP \\fB-v\\fP \\fItext\\fP",
                    "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

        // Parse command line.
        seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

        std::cout << <u>LINE</u> << std::endl;

        // Only extract  options if the program will continue after parseCommandLine()
        if (res != seqan::ArgumentParser::PARSE_OK)
            return res;

        std::cout << <u>LINE</u> << std::endl;
        // Extract option values.
        if (isSet(parser, "quiet"))
            options.verbosity = 0;
        if (isSet(parser, "verbose"))
            options.verbosity = 2;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 3;

        std::cout << <u>LINE</u> << std::endl;
        seqan::getArgumentValue(options.inputFile, parser, 0);

        // Get output file name from command line if set.  Otherwise, autogenerate from input file name.
        if (isSet(parser, "outputFile"))
        {
            std::cout << <u>LINE</u> << std::endl;
            seqan::getOptionValue(options.outputFile, parser, "outputFile");
            std::cout << options.inputFile << std::endl;
            std::cout << options.outputFile << std::endl;
        }
        else
        {
            std::cout << <u>LINE</u> << std::endl;
            options.outputFile = options.inputFile;
            std::cout << options.outputFile << std::endl;
            seqan::append(options.outputFile, ".fastq");
            std::cout << options.outputFile << std::endl;
        }

        std::cout << <u>LINE</u> << std::endl;
        return seqan::ArgumentParser::PARSE_OK;
    }

    // --------------------------------------------------------------------------
    // Function main()
    // --------------------------------------------------------------------------

    // Program entry point.

    int main(int argc, char const ** argv)
    {
        // Parse the command line.
        seqan::ArgumentParser parser;
        KnimeNodeOptions options;
        seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

        // If there was an error parsing or built-in argument parser functionality
        // was triggered then we exit the program.  The return code is 1 if there
        // were errors and 0 if there were none.
        if (res != seqan::ArgumentParser::PARSE_OK)
            return res == seqan::ArgumentParser::PARSE_ERROR;

        std::cout << "EXAMPLE PROGRAM\n"
                  << "===============\n\n";

        // Print the command line arguments back to the user.
        if (options.verbosity > 0)
        {
            std::cout << "<u>OPTIONS</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u><u>_</u>_\n"
                      << '\n'
                      << "VERBOSITY\t" << options.verbosity << '\n'
                      << "INPUT_FILE\t" << options.inputFile << "\n\n"
                      << "OUTPUT_FILE\t" << options.outputFile << "\n\n";
        }

        // Reading the input
        seqan::SequenceStream seqIn(seqan::toCString(options.inputFile));
        seqan::StringSet<seqan::CharString> ids;
        seqan::StringSet<seqan::Dna5String> seqs;
        seqan::StringSet<seqan::CharString> quals;

        std::cout << <u>LINE</u> << std::endl;

        if (atEnd(seqIn))
        {
            std::cout << "ERROR: File does not contain any sequences!\n";
            return 1;
        }
        if(seqan::readAll(ids, seqs, quals, seqIn) != 0)
        {
            std::cout << "ERROR: Could not read first record!\n";
            return 1;
        }

        std::cout << <u>LINE</u> << std::endl;
        // DO SOMETHING HERE
        // *
        // *
        // *

        seqan::SequenceStream seqOut(seqan::toCString(options.outputFile), seqan::SequenceStream::WRITE);
        std::cout << options.outputFile << std::endl;
        std::cout << "id: " << ids[0] << " " << length(ids) << std::endl;
        std::cout << seqs[0] << " " << length(seqs) << std::endl;
        std::cout << quals[0] << " " << length(quals) << std::endl;
        if (seqan::writeAll(seqOut, ids, seqs, quals) != 0)
        {
            std::cerr << "ERROR: Could not write records!\n";
            return 1;
        }
        std::cout << <u>LINE</u> << std::endl;

        return 0;
    }

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
