// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Globally shared code for the rabema tool.
// ==========================================================================

#ifndef SEQAN_APPS_RABEMA_RABEMA_H_
#define SEQAN_APPS_RABEMA_RABEMA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

// Tag for global options.
struct Global_ {};
typedef Tag<Global_> Global;

// Tag for the build golden standard subprogram, also used for
// specializing the Options class for the subprogram.
struct BuildGoldStandard_ {};
typedef Tag<BuildGoldStandard_> BuildGoldStandard;

// Tag for the evaluation subprogram, also used for specializing the
// Options class for the subprogram.
struct EvaluateResults_ {};
typedef Tag<EvaluateResults_> EvaluateResults;

// Enum for selecting a subprogram.
enum SelectedCommand
{
    COMMAND_NONE,
    COMMAND_BUILD_STANDARD,
    COMMAND_EVALUATE
};

// Class for storing options.
template <typename TSpec>
class Options;

// Global options used in rabema.cpp.
template <>
class Options<Global>
{
public:
    // True iff verbose mode is enabled.
    bool verbose;
    // True iff very verbose mode is enabled.
    bool veryVerbose;

    // The selected command to execute.
    SelectedCommand selectedCommand;

    // Path to the file to write the output to.  "-" for stdout.
    CharString outputFile;

    Options() :
        verbose(false), veryVerbose(false), selectedCommand(COMMAND_NONE), outputFile("-")
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Ceil away from Zero.
//
// ceilAwayFromZero(-1.5) == -2
// ceilAwayFromZero(1.5) == 2
template <typename T>
inline T ceilAwayFromZero(const T & x)
{
    if (x < 0)
        return floor(x);

    return ceil(x);
}

void setUpCommandLineParser(CommandLineParser & parser)
{
    addVersionLine(parser, "0.1");

    addTitleLine(parser, "*************************************");
    addTitleLine(parser, "* RABEMA - Read Alignment Benchmark *");
    addTitleLine(parser, "*************************************");
    addTitleLine(parser, "");
    addTitleLine(parser, "(c) 2010 by Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>");
}

#endif  // #ifndef SEQAN_APPS_RABEMA_RABEMA_H_
