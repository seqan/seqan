#include <fstream>
#include <iostream>
#include <string>

#include <seqan/store.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/svg.h>

using namespace seqan;

// Struct for storing command line options.
struct Options
{
    int contigID;
    int beginPos;
    int endPos;
    int beginLine;
    int endLine;
    bool writeAscii;
    bool gapSpace;

    std::string fileAliIn;
    std::string fileOut;
    std::string fileRefIn;

    Options() :
        contigID(-1), beginPos(-1), endPos(-1), beginLine(-1), endLine(-1), writeAscii(false), gapSpace(false)
    {}
};

// Parse command line, write results to options, return PARSE_OK if everything went well.
ArgumentParser::ParseResult parseCommandLine(Options & options, int argc, char const ** argv)
{
    ArgumentParser parser;
    setShortDescription(parser, "draw SAM/BAM file as SVG vector graphics");
    addUsageLine(parser, "[OPTION] -a <SAM/BAM file> [-r <GENOME file>] -o <SVG/ASCII file>");

    addOption(parser, ArgParseOption("a", "alignment", "SAM/BAM file to load.", ArgParseOption::INPUT_FILE));
    setRequired(parser, "alignment", true);
    setValidValues(parser, "alignment", BamFileIn::getFileExtensions());

    addOption(parser, ArgParseOption("r", "reference", "FASTA file to use as the reference.", ArgParseOption::INPUT_FILE));
    setValidValues(parser, "reference", SeqFileIn::getFileExtensions());

    addOption(parser, ArgParseOption("o", "out", "SVG/ASCII file to write.", ArgParseOption::OUTPUT_FILE));
    setRequired(parser, "out", true);
    setValidValues(parser, "out", ".txt .svg");

    addOption(parser, ArgParseOption("", "contig", "Display contig with this numeric ID only (zero-based), default is to show all contigs.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "contig", -1);

    addOption(parser, ArgParseOption("", "begin-pos", "Begin position of the region to show, default is to show all.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "begin-pos", -1);

    addOption(parser, ArgParseOption("", "end-pos", "End position of the region to show, default is to show all.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "end-pos", -1);

    addOption(parser, ArgParseOption("", "begin-line", "First line to show, zero-based", ArgParseOption::INTEGER));
    setDefaultValue(parser, "begin-line", -1);

    addOption(parser, ArgParseOption("", "end-line", "Last line to show, zero-based, default is to show all", ArgParseOption::INTEGER));
    setDefaultValue(parser, "end-line", -1);

    addOption(parser, ArgParseOption("", "gap-space", "begin and end position are given in gap space instead of in sequence space"));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.contigID, parser, "contig");
    getOptionValue(options.beginPos, parser, "begin-pos");
    getOptionValue(options.endPos, parser, "end-pos");
    getOptionValue(options.beginLine, parser, "begin-line");
    getOptionValue(options.endLine, parser, "end-line");

    getOptionValue(options.fileAliIn, parser, "alignment");
    getOptionValue(options.fileRefIn, parser, "reference");
    getOptionValue(options.fileOut, parser, "out");

    options.writeAscii = endsWith(options.fileOut, ".txt");
    options.gapSpace = isSet(parser, "gap-space");

    // Begin position cannot be greater than end position.
    if (options.beginPos != -1 && options.endPos != -1 && options.beginPos > options.endPos)
    {
        std::cerr << "ERROR: begin position cannot be greater than end position!\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // First line number cannot be greater than last line number.
    if (options.beginLine != -1 && options.endLine != -1 && options.beginLine < options.endLine)
    {
        std::cerr << "ERROR: first line cannot be greater than end position!\n";
        return ArgumentParser::PARSE_ERROR;
    }

    return ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
    // -----------------------------------------------------------------------
    // Parse Options
    // -----------------------------------------------------------------------
    Options options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // -----------------------------------------------------------------------
    // Load Files
    // -----------------------------------------------------------------------

    typedef FragmentStore<> TFragStore;

    typedef TFragStore::TContigStore TContigStore;
    typedef Value<TContigStore>::Type TContig;
    typedef TFragStore::TContigPos TContigPos;

    TFragStore store;
    BamFileIn samFile(options.fileAliIn.c_str());
    std::ofstream ascii;
    SVGFile svg;

    // Optionally load genome file
    if (!options.fileRefIn.empty())
    {
        if (!loadContigs(store, options.fileRefIn.c_str()))
        {
            std::cerr << "Failed to load genome.\n";
            return 1;
        }
    }

    // Load Sam file
    readRecords(store, samFile);

    // Choose contigs
    std::vector<unsigned> contigs;
    if (options.contigID == -1)
        for (unsigned i = 0; i < length(store.contigStore); ++i)
            contigs.push_back(i);
    else
        contigs.push_back(options.contigID);

    // Optionally load genome and open SVG file for writing
    if (options.writeAscii)
    {
        ascii.open(options.fileOut.c_str(), std::ios_base::out | std::ios_base::trunc);
        if (!ascii.is_open())
        {
            std::cerr << "ERROR: could not open output file for writing.\n";
            return 1;
        }
    }
    else if (!open(svg, options.fileOut.c_str()))
    {
        std::cerr << "ERROR: could not open output file for writing.\n";
        return 1;
    }

    // -----------------------------------------------------------------------
    // Write Alignment
    // -----------------------------------------------------------------------

    AlignedReadLayout layout;
    std::cerr << "Layouting reads ... ";
    layoutAlignment(layout, store);
    std::cerr << "done\n";

    std::cerr << "Writing " << contigs.size() << " contigs...\n";

    int beginLine = (options.beginLine == -1) ? 0 : options.beginLine;
    int endLine = (options.endLine == -1) ? std::numeric_limits<int>::max() : options.endLine;

    for (unsigned i = 0; i < contigs.size(); ++i)
        if (contigs[i] < length(store.contigStore))
        {
            std::cerr << "Writing contig " << store.contigNameStore[contigs[i]] << " ... ";

            int64_t l = (options.beginPos == -1) ? 0 : options.beginPos;
            int64_t r = (options.endPos == -1) ? std::numeric_limits<TContigPos>::max() : options.endPos;

            if (!options.gapSpace)
            {
                typedef Gaps<Nothing, AnchorGaps<TContig::TGapAnchors> >   TContigGaps;
                TContigGaps contigGaps(store.contigStore[i].gaps);
                l = positionSeqToGap(contigGaps, l);
                if (r != std::numeric_limits<TContigPos>::max())
                    r = positionSeqToGap(contigGaps, r);
            }

            if (r == std::numeric_limits<TContigPos>::max())
            {
                r = 0;
                for (unsigned j = 0; j < length(layout.contigRows[i]); ++j)
                {
                    unsigned id = back(layout.contigRows[i][j]);
                    if (r < store.alignedReadStore[id].beginPos)
                        r = store.alignedReadStore[id].beginPos;
                    if (r < store.alignedReadStore[id].endPos)
                        r = store.alignedReadStore[id].endPos;
                }
            }

            std::cerr << l << '\t' << r << '\t' << beginLine << '\t' << endLine << "\n";
            if (options.writeAscii)
                printAlignment(ascii, layout, store, contigs[i], l, r, beginLine, endLine);
            else
                printAlignment(svg, layout, store, contigs[i], l, r, beginLine, endLine);

            std::cout << "done\n";
        }

    return 0;
}
