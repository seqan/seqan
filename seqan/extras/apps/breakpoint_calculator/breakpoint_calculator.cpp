// ==========================================================================
//                           Breakpoint Calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#include <seqan/arg_parse.h>

#include "breakpoint_calculator.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{
    // Setup argument parser.
    ArgumentParser parser("breakpoint_calculator");
    Options options;
    setupCommandLineParser(parser);
    
    // Then, parse the command line.
    ArgumentParser::ParseResult res = parseArgumentsAndCheck(options, parser, argc, argv);
    
    // Finally, launch the program.
    if (res == ArgumentParser::PARSE_OK)
        return mainWithOptions(options);
    else
        return res == ArgumentParser::PARSE_ERROR; // 1 - error, 0 - otherwise
}
