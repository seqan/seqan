// ==========================================================================
//                   NGS: Regions of Interest Analysis
// ==========================================================================
// Copyright (c) 2012-2018, Bernd Jagla, Institut Pasteur
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
// Create ROI coverage plot thumbnail grid.
// ==========================================================================

#include <fstream>
#include <iostream>
#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/roi_io.h>

#include "png_canvas.h"

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

// Options for the bam2roi Application.

struct Options
{
    // Verbosity of output: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // -----------------------------------------------------------------------
    // Input / Output Options
    // -----------------------------------------------------------------------

    // Paths to input ROI file.
    seqan::CharString inputFileName;

    // Output prefix, will append "%d.png".
    seqan::CharString outputFileName;

    // -----------------------------------------------------------------------
    // Plot Options
    // -----------------------------------------------------------------------

    // Number of columns and rows in grid.
    int numCols;
    int numRows;

    // Width and height of one plot in px.
    int plotWidth;
    int plotHeight;

    // Border and space between plots.
    int borderWidth;
    int space;

    // Maximal number of plots to create, 0 for infinity.
    int maxRois;

    // Largest value to print in plots.  In case of going beyond, we'll make a little red dot at the top.  A value of 0
    // indicates no limit.
    int maxValue;

    Options() :
            verbosity(1), numCols(10), numRows(10), plotWidth(10), plotHeight(10), borderWidth(0), space(0),
            maxRois(0), maxValue(0)
    {}
};

// --------------------------------------------------------------------------
// Class Plotter
// --------------------------------------------------------------------------

class Plotter
{
public:
    // The options to use.
    Options options;

    // The index of the current grid and plot on grid.
    int gridNo;
    int plotNo;

    // The current PNG canvas.
    PngCanvas canvas;

    // Colors.
    PngColor color;
    PngColor bgColor;

    Plotter(Options const & options) :
            options(options), gridNo(0), plotNo(0),
            canvas(options.numCols * options.plotWidth + 2 * options.borderWidth + (options.numCols - 1) * options.space,
                   options.numRows * options.plotHeight + 2 * options.borderWidth + (options.numRows - 1) * options.space),
            color(0, 0, 0, 0xff), bgColor(0, 0, 0, 10)
    {}

    ~Plotter()
    {
        // Write PngCanvas if there is any plot on it.
        if (plotNo)
            writeGrid();
    }

    // Plot grid to canvas, advance gridNo, plotNo, and write out grid if necessary.
    void plotGrid(seqan::String<unsigned> const & points)
    {
        _doPlot(points);

        // Advance plotNo and gridNo if necessary, write out plot to file.
        plotNo += 1;
        if (plotNo == options.numCols * options.numRows)
        {
            writeGrid();
            clearGrid();
            plotNo = 0;
            gridNo += 1;
        }
    }

    std::pair<int, int> plotStart(int idx)
    {
        int row = idx / options.numCols;
        int col = idx % options.numCols;
        int x = options.borderWidth + col * (options.plotWidth + options.space);
        int y = options.borderWidth + row * (options.plotHeight + options.space);
        return std::make_pair(x, y);
    }

    void _doPlot(seqan::String<unsigned> const & counts)
    {
        std::pair<int, int> startPos = plotStart(plotNo);
        // Compute polyline coordinates.
        int numPoints = length(counts);
        double h = 0.0;
        double l = 0.0;
        for (unsigned i = 0; i < length(counts); ++i)
            h = std::max(h, (double)counts[i]);
        if (h == 0.0)
            h = 1.0;
        if (h == l)  // Center line if highest == lowest.
            h *= 2;

        if (options.maxValue != 0)
            h = options.maxValue;

        canvas.color = bgColor;
        canvas.filledRectangle(startPos.first, startPos.second,
                                    startPos.first + options.plotWidth,
                                    startPos.second + options.plotHeight);

        canvas.color = color;
        std::vector<std::pair<int, int> > points;
        for (unsigned i = 0; i < length(counts); ++i)
        {
            int x = static_cast<int>(startPos.first + (double)i / numPoints * options.plotWidth);
            unsigned val = counts[i];
            if (options.maxValue != 0 && (int)val >= options.maxValue)
                val = options.maxValue;
            int y = static_cast<int>(startPos.second + options.plotHeight - val / h * options.plotHeight);
            points.push_back(std::make_pair(x, y));
        }
        canvas.polyline(points);

        // Print polyline instructions if --very-verbose.
        if (options.verbosity >= 3)
        {
            std::cerr << "polyline\n";
            for (unsigned i = 0; i < points.size(); ++i)
                std::cerr << "  " << points[i].first << ", " << points[i].second << "\n";
        }

        PngColor red(0xff, 0, 0, 0xff);
        canvas.color = red;
        for (unsigned i = 0; i < points.size(); ++i)
            if (options.maxValue != 0 && (int)counts[i] >= options.maxValue)
                canvas.point(points[i].first, points[i].second);
    }

    // Write grid on PngCanvas to output.
    void writeGrid()
    {
        std::stringstream ss;
        ss << options.outputFileName << gridNo << ".png";
        canvas.write(ss.str().c_str());
    }

    // Clear grid.
    void clearGrid()
    {
        canvas.color = PngColor::WHITE();
        canvas.filledRectangle(0, 0, canvas.width - 1, canvas.height - 1);
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function yesNo()                                                 [Options]
// --------------------------------------------------------------------------

char const * yesNo(bool b)
{
    return b ? "YES" : "NO";
}

// --------------------------------------------------------------------------
// Function print()                                                 [Options]
// --------------------------------------------------------------------------

void print(std::ostream & out, Options const & options)
{
    out << "__OPTIONS_____________________________________________________________________\n"
        << "\n"
        << "INPUT FILE       \t" << options.inputFileName << "\n"
        << "OUTPUT FILE      \t" << options.outputFileName << "\n"
        << "\n"
        << "NUM COLS         \t" << options.numCols << "\n"
        << "NUM ROWS         \t" << options.numRows << "\n"
        << "\n"
        << "PLOT WIDTH       \t" << options.plotWidth << "\n"
        << "PLOT HEIGHT      \t" << options.plotHeight << "\n"
        << "\n"
        << "BORDER WIDTH     \t" << options.borderWidth << "\n"
        << "SPACE            \t" << options.space << "\n"
        << "\n"
        << "MAX ROIS         \t" << options.maxRois << "\n"
        << "MAX VALUE        \t" << options.maxValue << "\n"
        << "\n";
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

// Parse command line with ArgumentParser class.

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("roi_plot_thumbnails");
    setCategory(parser, "NGS ROI Analysis");

    // Set short description, version, and date.
    setShortDescription(parser, "Create plot grids for ROI file.");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "\\fB-if\\fP \\fIIN.roi\\fP \\fB-o\\fP \\fIOUT\\fP");

    addDescription(parser,
                   "Create PNG images with plot grids to \\fIOUT${i}.png\\fP from ROI records "
                   "in \\fIIN.roi\\fP.");

    // -----------------------------------------------------------------------
    // General Options
    // -----------------------------------------------------------------------

    addSection(parser, "General Options");

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose mode."));
    addOption(parser, seqan::ArgParseOption("vv", "vverbose", "Very verbose mode."));

    // -----------------------------------------------------------------------
    // Input / Output Options
    // -----------------------------------------------------------------------

    addSection(parser, "Input / Output Parameters");

    addOption(parser, seqan::ArgParseOption("if", "input-file", "ROI to plot.",
                                            seqan::ArgParseOption::INPUT_FILE));
    setValidValues(parser, "input-file", seqan::RoiFileIn::getFileExtensions());
    setRequired(parser, "input-file");

    addOption(parser, seqan::ArgParseOption("o", "output-prefix", "Prefix of output file.",
                                            seqan::ArgParseOption::OUTPUT_FILE));
    setRequired(parser, "output-prefix");

    // -----------------------------------------------------------------------
    // Plot Configuration
    // -----------------------------------------------------------------------

    addSection(parser, "PlotConfiguration");

	addOption(parser, seqan::ArgParseOption("", "max-rois", "Maximal number of plots to create (0 for all).",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "max-rois", "0");
    setDefaultValue(parser, "max-rois", "20000");

	addOption(parser, seqan::ArgParseOption("", "max-value", "Fix maximal y value.  0 for no limit..",
                                            seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "max-value", "0");
    setDefaultValue(parser, "max-value", "0");

	addOption(parser, seqan::ArgParseOption("", "num-cols", "Number of plots in one row.",
                                            seqan::ArgParseOption::INTEGER, "COLS"));
    setMinValue(parser, "num-cols", "1");
    setDefaultValue(parser, "num-cols", "40");

	addOption(parser, seqan::ArgParseOption("", "num-rows", "Number of plots in one column.",
                                            seqan::ArgParseOption::INTEGER, "ROWS"));
    setMinValue(parser, "num-rows", "1");
    setDefaultValue(parser, "num-rows", "50");

	addOption(parser, seqan::ArgParseOption("", "plot-height", "Height of one plot in px.",
                                            seqan::ArgParseOption::INTEGER, "HEIGHT"));
    setMinValue(parser, "plot-height", "1");
    setDefaultValue(parser, "plot-height", "30");

	addOption(parser, seqan::ArgParseOption("", "plot-width", "Width of one plot in px.",
                                            seqan::ArgParseOption::INTEGER, "WIDTH"));
    setMinValue(parser, "plot-width", "1");
    setDefaultValue(parser, "plot-width", "30");

	addOption(parser, seqan::ArgParseOption("", "border-width", "Border width.",
                                            seqan::ArgParseOption::INTEGER, "WIDTH"));
    setMinValue(parser, "border-width", "0");
    setDefaultValue(parser, "border-width", "0");

	addOption(parser, seqan::ArgParseOption("", "spacing", "Space between plots.",
                                            seqan::ArgParseOption::INTEGER, "WIDTH"));
    setDefaultValue(parser, "spacing", "2");

    // -----------------------------------------------------------------------
    // Parsing and Value Extraction
    // -----------------------------------------------------------------------

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.

    options.verbosity = 1;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "vverbose"))
        options.verbosity = 3;

    getOptionValue(options.inputFileName, parser, "input-file");
    getOptionValue(options.outputFileName, parser, "output-prefix");

    getOptionValue(options.numCols, parser, "num-cols");
    getOptionValue(options.numRows, parser, "num-rows");
    getOptionValue(options.plotWidth, parser, "plot-width");
    getOptionValue(options.plotHeight, parser, "plot-height");
    getOptionValue(options.borderWidth, parser, "border-width");
    getOptionValue(options.space, parser, "spacing");
    getOptionValue(options.maxRois, parser, "max-rois");
    getOptionValue(options.maxValue, parser, "max-value");

	return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
	// Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0 (e.g. help
    // was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Print program name and options.

    if (options.verbosity >= 1)
    {
        std::cerr << "ROI PLOT GRIDS\n"
                  << "==============\n"
                  << "\n";
        print(std::cerr, options);
    }

    // -----------------------------------------------------------------------
    // Open input and output file
    // -----------------------------------------------------------------------

    if (options.verbosity >= 1)
        std::cerr << "__OPENING FILES_______________________________________________________________\n"
                  << "\n";

    if (options.verbosity >= 1)
         std::cerr << "Opening " << options.inputFileName << " ...";
    seqan::RoiFileIn roiFileIn;
    if (!open(roiFileIn, toCString(options.inputFileName)))
	{
		std::cerr << "\nERROR: Could not open " << options.inputFileName << "\n";
        return 1;
	}
    if (options.verbosity >= 1)
        std::cerr << " OK\n";

    // Read header (actually, skip it).
    seqan::RoiHeader roiHeader;
    readHeader(roiHeader, roiFileIn);

    // -----------------------------------------------------------------------
    // Create plot grids.
    // -----------------------------------------------------------------------

    if (options.verbosity >= 1)
        std::cerr << "\n__PLOTTING____________________________________________________________________\n"
                  << "\n"
                  << "WORKING...";

    Plotter plotter(options);

    seqan::RoiRecord record;
    for (int i = 0; (options.maxRois == 0 || i < options.maxRois) && !atEnd(roiFileIn); ++i)
    {
        if (options.verbosity >= 2)
            std::cerr << " " << i << "(" << length(record.count) << ", " << plotter.gridNo << ")";
        readRecord(record, roiFileIn);

        plotter.plotGrid(record.count);
    }

    if (options.verbosity >= 1)
        std::cerr << " OK\n\nDone.\n";

    return 0;
}
