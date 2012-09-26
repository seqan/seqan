// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Bjoern Kahlert <Bjoern.Kahlert@fu-berlin.de>
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_CTD_SUPPORT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_CTD_SUPPORT_H_

#include <seqan/sequence.h>

#include <seqan/arg_parse/xml_support.h>
#include <seqan/arg_parse/argument_parser.h>
#include <seqan/arg_parse/arg_parse_doc.h>

#include <fstream>

namespace seqan {

// ----------------------------------------------------------------------------
// Function _join()
// ----------------------------------------------------------------------------
template <typename TSequence>
inline TSequence
_toHTML(TSequence const & sequence)
{
    HtmlToolDocPrinter_ docPrinter;
    return docPrinter._toHtml(sequence);
}

// ----------------------------------------------------------------------------
// Function _join()
// ----------------------------------------------------------------------------

/**
 * joins all elements of the the passed StringSet into a single CharString
 * the provided delimiter is used to separate the single entries in the
 * resulting CharString
 */
template <typename TValue>
inline std::string
_join(std::vector<TValue> const & v, std::string const & delimiter)
{
    typedef typename std::vector<TValue>::const_iterator TStringSetIterator;

    std::stringstream joined;
    for (TStringSetIterator it = v.begin(); it != v.end(); ++it)
    {
        if (it != v.begin())
            joined << delimiter;
        joined << *it;
    }
    return joined.str();
}

// ----------------------------------------------------------------------------
// Function _getPrefixedOptionName()
// ----------------------------------------------------------------------------

inline std::string
_getPrefixedOptionName(ArgParseOption const & opt)
{
    std::string optName = "";
    if (!empty(opt.longName))
        optName = "--" + opt.longName;
    else
        optName = "-" + opt.shortName;

    return optName;
}

// ----------------------------------------------------------------------------
// Function _getOptionName()
// ----------------------------------------------------------------------------

inline std::string
_getOptionName(ArgParseOption const & opt)
{
    if (!empty(opt.longName))
        return opt.longName;
    else
        return opt.shortName;
}

// ----------------------------------------------------------------------------
// Function _addMinMaxRestrictions()
// ----------------------------------------------------------------------------

inline void
_addMinMaxRestrictions(std::vector<std::string> & restrictions, ArgParseArgument const & opt)
{

    std::string minMaxRestriction = "";
    if (opt.minValue != "")
    {
        append(minMaxRestriction, opt.minValue);
        append(minMaxRestriction, ":");
    }
    if (opt.maxValue != "")
    {
        if (minMaxRestriction == "")
            append(minMaxRestriction, ":");
        append(minMaxRestriction, opt.maxValue);
    }

    if (minMaxRestriction != "")
        appendValue(restrictions, minMaxRestriction);

}

// ----------------------------------------------------------------------------
// Function _addValidValuesRestrictions()
// ----------------------------------------------------------------------------

inline void
_addValidValuesRestrictions(std::vector<std::string> & restrictions, ArgParseArgument const & opt)
{
    if (length(opt.validValues) != 0)
    {
        for (std::vector<std::string>::const_iterator valid = opt.validValues.begin();
             valid != opt.validValues.end();
             ++valid)
        {
            // for files we set *.(Name of the format)
            if (isOutputFileArgument(opt) || isInputFileArgument(opt))
            {
                std::string filetype = "*.";
                append(filetype, *valid);
                appendValue(restrictions, filetype);
            }
            else
            {
                appendValue(restrictions, *valid);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function _includeInCTD()
// ----------------------------------------------------------------------------

/*
 * returns true if this option should be included in the ctd
 */
inline bool
_includeInCTD(ArgParseOption const & opt)
{
    return !(opt.longName == "help" || opt.longName == "version" || opt.longName == "write-ctd" || opt.longName == "export-help" || (opt.shortName == "" && opt.longName == ""));
}

// ----------------------------------------------------------------------------
// Function _indent()
// ----------------------------------------------------------------------------

std::string _indent(const int currentIndent)
{
    std::string indent = "";
    for (int i = 0; i < currentIndent; ++i)
        indent += "\t";
    return indent;
}

void _writeCLIElement(std::ofstream & ctdfile, int currentIndent, std::string const & optionIdentifier, std::string const & ref_name, bool isList)
{
    ctdfile << _indent(currentIndent)
            << "<clielement optionIdentifier=\"" << optionIdentifier
            << "\" isList=\"" << (isList ? "true" : "false") << "\">\n";

    ctdfile << _indent(currentIndent + 1) << "<mapping ref_name=\"" << ref_name << "\" />\n";

    ctdfile << _indent(currentIndent) << "</clielement>\n";
}

// ----------------------------------------------------------------------------
// Function writeCTD()
// ----------------------------------------------------------------------------

/**
.Function.writeCTD
..summary:Exports the app's interface description to a .ctd file.
..cat:Miscellaneous
..signature:writeCTD(parser)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..returns:$true$ if the ctd file could be created correctly, $false$ otherwise.
..include:seqan/arg_parse.h
*/

inline bool
writeCTD(ArgumentParser const & me)
{
    typedef ArgumentParser::TOptionMap::const_iterator   TOptionMapIterator;
    typedef ArgumentParser::TArgumentMap::const_iterator TArgumentMapIterator;
    typedef ArgumentParser::TArgumentMapSize TArgumentMapSize;

    // create file [appname].ctd in working directory
    std::string ctdfilename;
    getOptionValue(ctdfilename, me, "write-ctd");

    std::ofstream ctdfile;
    ctdfile.open(toCString(ctdfilename));

    if (!ctdfile.is_open())
    {
        std::cerr << getAppName(me) << ": Unable to create ctd file: " << ctdfilename << std::endl;
        return false;
    }

    ctdfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ctdfile << "<tool status=\"external\">\n";

    int currentIndent = 1;

    std::string toolname(toCString(xmlEscape(getAppName(me))));

    ctdfile << _indent(currentIndent) << "<name>" << toolname << "</name>\n";
    ctdfile << _indent(currentIndent) << "<version>" << xmlEscape(getVersion(me)) << "</version>\n";
    ctdfile << _indent(currentIndent) << "<description><![CDATA[" << xmlEscape(getShortDescription(me)) << ".]]></description>\n";
    ctdfile << _indent(currentIndent) << "<manual><![CDATA[" << xmlEscape(getAppName(me)) << ".]]></manual>\n"; // TODO(aiche): as soon as we have a more sophisticated documentation embedded into the CmdParser, we should at this here
    ctdfile << _indent(currentIndent) << "<docurl>Direct links in docs</docurl>\n";
    ctdfile << _indent(currentIndent) << "<category>" << xmlEscape(getCategory(me)) << "</category>\n";
    ctdfile << _indent(currentIndent++) << "<cli>\n";

    // the unix way 1st the options
    for (TOptionMapIterator optionMapIterator = me.optionMap.begin();
         optionMapIterator != me.optionMap.end();
         ++optionMapIterator)
    {
        ArgParseOption const & opt = *optionMapIterator;
        std::string optionIdentifier = _getPrefixedOptionName(opt);
        std::string refName = toolname + "." + _getOptionName(opt);

        if (_includeInCTD(opt))
        {
            _writeCLIElement(ctdfile, currentIndent, optionIdentifier, refName, isListArgument(opt));
        }
    }

    // add a warning to the CTD that arguments are hard to interpret by the users
    if (me.argumentList.size() > 0)
    {
        ctdfile << _indent(currentIndent)
                << "<!-- Following clielements are arguments."
                << " You should consider providing a help text to ease understanding. -->\n";
    }
    // then the arguments
    for (TArgumentMapSize argIdx = 0; argIdx != me.argumentList.size(); ++argIdx)
    {
        // arguments do not have an option identifier
        std::string optionIdentifier = "";
        std::stringstream refName;
        refName << toolname << "." << "argument-" << argIdx;
        _writeCLIElement(ctdfile, currentIndent, optionIdentifier, refName.str(), isListArgument(me.argumentList[argIdx]));
    }

    ctdfile << _indent(--currentIndent) << "</cli>\n";
    ctdfile << _indent(currentIndent++) << "<PARAMETERS version=\"1.3\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/Param_1_3.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    ctdfile << _indent(currentIndent++) << "<NODE name=\"" << toolname << "\" description=\"???\">\n";

    for (TOptionMapIterator optionMapIterator = me.optionMap.begin();
         optionMapIterator != me.optionMap.end();
         ++optionMapIterator)
    {
        ArgParseOption const & opt = *optionMapIterator;

        // exclude help, version, etc.
        if (!_includeInCTD(opt))
            continue;

        // prefer short name for options
        std::string optionName = _getOptionName(opt);

        std::string type;

        if (isStringArgument(opt) || isBooleanOption(opt))
            type = "string";
        else if (isIntegerArgument(opt))
            type = "int";
        else if (isDoubleArgument(opt))
            type = "double";

        // set up tags
        std::vector<std::string> tags;
        if (isInputFileArgument(opt))
        {
            appendValue(tags, "input file");
        }
        if (isOutputFileArgument(opt))
        {
            appendValue(tags, "output file");
        }
        if (isRequired(opt))
        {
            appendValue(tags, "required");
        }

        // set up restrictions
        std::vector<std::string> restrictions;
        _addValidValuesRestrictions(restrictions, opt);
        _addMinMaxRestrictions(restrictions, opt);


        if (isListArgument(opt))
        {
            ctdfile << _indent(currentIndent)
                    << "<ITEMLIST " << "name=\"" << xmlEscape(optionName) << "\" "
                    << "type=\"" << type << "\" "
                    << "description=\"" << xmlEscape(_toHTML(opt._helpText)) << "\" "
                    << "tags=\"" << xmlEscape(_join(tags, ",")) << "\" "
                    << "restrictions=\"" << xmlEscape(_join(restrictions, ",")) << "\""
                    << ">\n";

            for (size_t i = 0; i < opt.defaultValue.size(); ++i)
            {
                ctdfile << _indent(currentIndent + 1) << "<LISTITEM value=\"" << xmlEscape(opt.defaultValue[i]) << "\"/>\n";
            }
            ctdfile << _indent(currentIndent) << "</ITEMLIST>\n";
        }
        else
        {
            ctdfile << _indent(currentIndent)
                    << "<ITEM " << "name=\"" << xmlEscape(optionName) << "\" "
                    << "value=\"" << xmlEscape(_join(opt.defaultValue, ",")) << "\" "
                    << "type=\"" << type << "\" "
                    << "description=\"" << xmlEscape(_toHTML(opt._helpText)) << "\" "
                    << "tags=\"" << xmlEscape(_join(tags, ",")) << "\" "
                    << "restrictions=\"" << xmlEscape(_join(restrictions, ",")) << "\""
                    << "/>\n";
        }
    }

    for (TArgumentMapSize argIdx = 0; argIdx != me.argumentList.size(); ++argIdx)
    {
        ArgParseArgument arg = me.argumentList[argIdx];

        // prefer short name for options
        std::stringstream argumentNameStream;
        argumentNameStream << "argument-" << argIdx;
        std::string optionName = argumentNameStream.str();

        std::string type;

        if (isStringArgument(arg))
            type = "string";
        else if (isIntegerArgument(arg))
            type = "int";
        else if (isDoubleArgument(arg))
            type = "double";

        // set up tags
        std::vector<std::string> tags;
        appendValue(tags, "required");
        if (isInputFileArgument(arg))
        {
            appendValue(tags, "input file");
        }
        if (isOutputFileArgument(arg))
        {
            appendValue(tags, "output file");
        }

        // set up restrictions
        std::vector<std::string> restrictions;
        _addValidValuesRestrictions(restrictions, arg);
        _addMinMaxRestrictions(restrictions, arg);


        ctdfile << _indent(currentIndent)
                << "<ITEM" << (isListArgument(arg) ? "LIST" : "") << " name=\"" << xmlEscape(optionName) << "\" "
                << (isListArgument(arg) ? " " : "value=\"\" ")
                << "type=\"" << type << "\" "
                << "description=\"" << xmlEscape(_toHTML(arg._helpText)) << "\" " // it will be "" in most cases but we try
                << "tags=\"" << xmlEscape(_join(tags, ",")) << "\" "
                << "restrictions=\"" << xmlEscape(_join(restrictions, ",")) << "\""
                << "/>\n";
    }

    ctdfile << _indent(--currentIndent) << "</NODE>\n";
    ctdfile << _indent(--currentIndent) << "</PARAMETERS>\n";
    ctdfile << "</tool>" << std::endl;

    ctdfile.close();
    return true;
}

} // namespace seqan

#endif // SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_CTD_SUPPORT_H_
