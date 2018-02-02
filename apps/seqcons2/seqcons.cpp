// ==========================================================================
//                                 SeqCons
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

#include <stdexcept>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "seqcons_app.h"
#include "seqcons_options.h"

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    SeqConsOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Print the command line arguments back to the user.
    if (options.verbosity >= 1)
    {
        std::cerr << "SeqCons\n"
                  << "=======\n"
                  << "\n"
                  << "__OPTIONS____________________________________________________________________\n"
                  << '\n';
        options.print(std::cerr);
    }

    // Perform additional consistency checking of options.
    try
    {
        options.checkConsistency();
    }
    catch (std::runtime_error & e)
    {
        std::cerr << "\nERROR: Inconsistent command line options:\n"
                  << "  " << e.what() << "\n";
        return 1;
    }

    // Run aplication.
    double startTime = seqan::sysTime();
    try
    {
        SeqConsApp app(options);
        app.run();
    }
    catch (std::runtime_error & e)
    {
        std::cerr << "\nERROR: An error occurred during the program's execution:\n"
                  << "  " << e.what() << "\n";
        return 1;
    }

    if (options.verbosity >= 1)
    {
        std::cerr << "\nOverall time: " << seqan::sysTime() - startTime << " s\n"
                  << "\n"
                  << "Done.  Have a nice day.\n";
    }
    return 0;
}
