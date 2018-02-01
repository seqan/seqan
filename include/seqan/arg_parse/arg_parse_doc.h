// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_
#define SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_

#include <seqan/arg_parse/tool_doc.h>
#include <seqan/arg_parse/argument_parser.h>

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function getAppName()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getAppName
 * @brief Return program name of ArgumentParser.
 *
 * @signature TCharStringRef getAppName(parser);
 *
 * @param[in] parser The ArgumentParser to get the app name for.
 *
 * @return TCharStringRef The app name, const-ref to @link CharString @endlink.
 */

inline CharString const & getAppName(ArgumentParser const & parser)
{
    return getName(parser._toolDoc);
}

// ----------------------------------------------------------------------------
// Helper Function _parseAppName()
// ----------------------------------------------------------------------------

inline void _parseAppName(ArgumentParser & parser, std::string const & candidate)
{
    //IOREV _notio_ irrelevant for io-revision
    int i = length(candidate) - 1;

    for (; i >= 0; --i)
        if (candidate[i] == '\\' || candidate[i] == '/')
            break;

    setName(parser._toolDoc, candidate.substr(i + 1));
}

// ----------------------------------------------------------------------------
// Helper Function _addLine()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addLine
 * @brief Adds a line of text to the help output of the ArgumentParser.
 *
 * The line of text will be added to the block of the options.
 *
 * @signature void addLine(parser, line);
 *
 * @param[in,out] parser The ArgumentParser to add the line to.
 * @param[in]     line   The line of text to add, @link StringConcept @endlink of <tt>char</tt>.
 */

template <typename TString>
inline void addLine(ArgumentParser & me, TString const & line)
{
    addOption(me, ArgParseOption("", "", line));
}

// ----------------------------------------------------------------------------
// Function addSection()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addSection
 * @brief Begins a new section of the option block of the ArgumentParser help output.
 *
 * @signature void addSection(parser, title);
 *
 * @param[in,out] parser The ArgumentParser to add the line to.
 * @param[in]     title  The title to add, @link StringConcept @endlink of <tt>char</tt>.
 *
 * @code{.cpp}
 * ArgumentParser parser;
 *
 * [...] // init parser
 *
 * addSection(parser, "In-/Output-Options");
 * addOption("i", ... );
 * addOption("o", ... );
 *
 * addSection(parser, "Other Options");
 * addOption("x", ... );
 * @endcode
 */

template <typename TString>
inline void addSection(ArgumentParser & me, TString const & line)
{
    addLine(me, "");
    addLine(me, line);
}

// ----------------------------------------------------------------------------
// Function addUsageLine()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addUseLine
 * @brief Adds a line of text to the usage output of the ArgumentParser.
 *
 * @signature void addUsageLine(parser, line);
 *
 * @param[in,out] parser The ArgumentParser to add the line to.
 * @param[in]     line   The line to add, a <tt>std::string</tt>.
 */

inline void addUsageLine(ArgumentParser & me, std::string const & line)
{
    me._usageText.push_back(line);
}

// ----------------------------------------------------------------------------
// Helper Function _addUsage()
// ----------------------------------------------------------------------------

inline void _addUsage(ToolDoc & toolDoc, ArgumentParser const & me)
{
    for (unsigned i = 0; i < length(me._usageText); ++i)
    {
        std::string text = "\\fB";
        append(text, getAppName(me));
        append(text, "\\fP ");
        append(text, me._usageText[i]);
        addText(toolDoc, text, false);
    }
}

// ----------------------------------------------------------------------------
// Function addDescription()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addDescription
 * @brief Appends a description paragraph to the ArgumentParser documentation.
 *
 * @signature void addDescription(parser, description);
 *
 * @param[in,out] parser      The ArgumentParser to add the line to.
 * @param[in]     description The description text, a <tt>std::string</tt>.
 */

inline void addDescription(ArgumentParser & me, std::string const & description)
{
    me._description.push_back(description);
}

// ----------------------------------------------------------------------------
// Function setAppName()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setAppName
 * @brief Sets application name of ArgumentParser.
 *
 * @signature void setAppName(parser, name);
 *
 * @param[in,out] parser The ArgumentParser to set the name of.
 * @param[in]     name   The application name, <tt>std::string</tt>.
 */

inline void setAppName(ArgumentParser & me, std::string const & name)
{
    setName(me._toolDoc, name);
}

// ----------------------------------------------------------------------------
// Function setShortDescription()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setShortDescription
 * @brief Sets shortDescription of ArgumentParser.
 *
 * @signature void setShortDescription(parser, desc);
 *
 * @param[in,out] parser The ArgumentParser to set the short description of.
 * @param[in]     desc   The short description, <tt>std::string</tt>.
 */

inline void setShortDescription(ArgumentParser & me, std::string const & description)
{
    setShortDescription(me._toolDoc, description);
}

// ----------------------------------------------------------------------------
// Function getShortDescription()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getShortDescription
 * @brief Returns the short description.
 *
 * @signature CharString getShortDescription(parser);
 *
 * @param[in,out] parser The ArgumentParser to get short description for.
 *
 * @return CharString A @link CharString @endlink with the short description.
 */

inline CharString getShortDescription(ArgumentParser const & me)
{
    return getShortDescription(me._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setUrl()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setUrl
 * @brief Sets url of ArgumentParser.
 *
 * @signature void setUrl(parser, url);
 *
 * @param[in,out] parser  The ArgumentParser to set the url of.
 * @param[in]     url     The url string to set, @link CharString @endlink.
 */

inline void setUrl(ArgumentParser & me, CharString const & urlString)
{
    setUrl(me._toolDoc, urlString);
}

// --------------------------------------------------------------------------
// Function getUrl()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getUrl
 * @brief Returns the url string.
 *
 * @signature TCharStringRef getUrl(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the url string from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the url string.
 */

inline CharString const & getUrl(ArgumentParser const & me)
{
    return getUrl(me._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setVersion()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setVersion
 * @brief Sets version of ArgumentParser.
 *
 * @signature void setVersion(parser, version);
 *
 * @param[in,out] parser  The ArgumentParser to set the version of.
 * @param[in]     version The version string to set, <tt>std::string</tt>.
 */

inline void setVersion(ArgumentParser & me, std::string const & versionString)
{
    setVersion(me._toolDoc, versionString);
    if (!hasOption(me, "version"))
        addOption(me, ArgParseOption("", "version", "Display version information."));
}

// --------------------------------------------------------------------------
// Function getVersion()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getVersion
 * @brief Returns the version string.
 *
 * @signature TCharStringRef getVersion(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the version string from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the version string.
 */

inline CharString const & getVersion(ArgumentParser const & me)
{
    return getVersion(me._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setShortCopyright()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setShortCopyright
 * @brief Sets short copyright of ArgumentParser.
 *
 * @signature void setShortCopyright(parser, short copyright);
 *
 * @param[in,out] parser  The ArgumentParser to set the short copyright of.
 * @param[in]     short copyright The short copyright string to set, <tt>std::string</tt>.
 */

inline void setShortCopyright(ArgumentParser & me, CharString const & shortCopyrightString)
{
    setShortCopyright(me._toolDoc, shortCopyrightString);
}

// --------------------------------------------------------------------------
// Function getShortCopyright()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getShortCopyright
 * @brief Returns the short copyright string.
 *
 * @signature TCharStringRef getShortCopyright(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the short copyright string from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the short copyright string.
 */

inline CharString const & getShortCopyright(ArgumentParser const & me)
{
    return getShortCopyright(me._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setLongCopyright()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setLongCopyright
 * @brief Sets long copyright of ArgumentParser.
 *
 * @signature void setLongCopyright(parser, long copyright);
 *
 * @param[in,out] parser  The ArgumentParser to set the long copyright of.
 * @param[in]     long copyright The long copyright string to set, <tt>std::string</tt>.
 */

inline void setLongCopyright(ArgumentParser & me, CharString const & longCopyrightString)
{
    setLongCopyright(me._toolDoc, longCopyrightString);
    if (!hasOption(me, "copyright"))
        addOption(me, ArgParseOption("", "copyright", "Display long copyright information."));
}

// --------------------------------------------------------------------------
// Function getLongCopyright()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getLongCopyright
 * @brief Returns the long copyright string.
 *
 * @signature TCharStringRef getLongCopyright(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the long copyright string from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the long copyright string.
 */

inline CharString const & getLongCopyright(ArgumentParser const & me)
{
    return getLongCopyright(me._toolDoc);
}


// ----------------------------------------------------------------------------
// Function setCitation()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setCitation
 * @brief Sets citation of ArgumentParser.
 *
 * @signature void setCitation(parser, citation);
 *
 * @param[in,out] parser  The ArgumentParser to set the citation of.
 * @param[in]     citation The citation string to set, <tt>std::string</tt>.
 */

inline void setCitation(ArgumentParser & me, CharString const & citationString)
{
    setCitation(me._toolDoc, citationString);
}

// --------------------------------------------------------------------------
// Function getCitation()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getCitation
 * @brief Returns the citation string.
 *
 * @signature TCharStringRef getCitation(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the citation string from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the citation string.
 */

inline CharString const & getCitation(ArgumentParser const & me)
{
    return getCitation(me._toolDoc);
}

// --------------------------------------------------------------------------
// Function setCategory()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setCategory
 * @brief Sets category of ArgumentParser.
 *
 * @signature void setCategory(parser, category);
 *
 * @param[in,out] parser  The ArgumentParser to set the category of.
 * @param[in]     category The category to set, <tt>std::string</tt>.
 */

inline void setCategory(ArgumentParser & parser, CharString const & category)
{
    setCategory(parser._toolDoc, category);
}

// --------------------------------------------------------------------------
// Function getCategory()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getCategory
 * @brief Returns the category.
 *
 * @signature TCharStringRef getCategory(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the category from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the category.
 */

inline CharString const & getCategory(ArgumentParser const & parser)
{
    return getCategory(parser._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setDate()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setDate
 * @brief Sets date string of ArgumentParser.
 *
 * @signature void setDate(parser, date);
 *
 * @param[in,out] parser The ArgumentParser to set the date string of.
 * @param[in]     date   The date string to set, <tt>std::string</tt>.
 */

inline void setDate(ArgumentParser & me, std::string const & date)
{
    setDate(me._toolDoc, date);
}

// ----------------------------------------------------------------------------
// Function addTextSection()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addTextSection
 * @brief Add a text section to the ArgumentParser.
 *
 * @signature void addTextSection(parser, title);
 *
 * @param[in,out] parser The ArgumentParser to add the text section title to.
 * @param[in]     title  The section title to add, <tt>std::string</tt>.
 */

inline void addTextSection(ArgumentParser & me, std::string const & title)
{
    addSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addTextSubSection()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addTextSubSection
 * @brief Add a text sub section to the ArgumentParser.
 *
 * @signature void addTextSubSection(parser, title);
 *
 * @param[in,out] parser The ArgumentParser add the subsection title to of.
 * @param[in]     title  The sub section title to add, <tt>std::string</tt>.
 */

inline void addTextSubSection(ArgumentParser & me, std::string const & title)
{
    addSubSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addText()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addText
 * @brief Add text to an ArgumentParser.
 *
 * @signature void addText(parser, text);
 *
 * @param[in,out] parser ArgumentParser to add text to.
 * @param[in]     text   The <tt>std::string</tt> to add to the parser.
 */

inline void addText(ArgumentParser & me, std::string const & text)
{
    addText(me._toolDoc, text);
}

// ----------------------------------------------------------------------------
// Function addListItem()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addListItem
 * @brief Appends a list item to the ArgumentParser
 *
 * @signature void addListItem(parser, item, description);
 *
 * @param[in,out] parser      The ArgumentParser to add the list item to.
 * @param[in]     item        The item to add, <tt>std::string</tt>.
 * @param[in]     description The item to add, <tt>std::string</tt>.
 */

inline void addListItem(ArgumentParser & me, std::string const & item, std::string const & description)
{
    addListItem(me._toolDoc, item, description);
}

// ----------------------------------------------------------------------------
// Function printShortHelp()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#printShortHelp
 * @brief Prints a short help message for the parser to a stream.
 *
 * @signature void printShortHelp(parser, out);
 *
 * @param[in,out] parser The ArgumentParser to print help for.
 * @param[in,out] out    The <tt>std::ostream</tt> to print help to.
 */

inline void printShortHelp(ArgumentParser const & me, std::ostream & stream)
{
    // TODO: maybe we can get this a bit prettier
    ToolDoc shortDoc(me._toolDoc);
    clearEntries(shortDoc);

    _addUsage(shortDoc, me);

    std::stringstream shortHelp;
    shortHelp << "Try '" << getAppName(me) << " --help' for more information.\n";
    addText(shortDoc, shortHelp.str());

    print(stream, shortDoc, "txt");
}

inline void printShortHelp(ArgumentParser const & me)
{
    printShortHelp(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function printVersion()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#printVersion
 * @brief Prints the version information of the parser to a stream.
 *
 * @signature void printVersion(parser, stream);
 *
 * @param[in,out] parser The ArgumenParser to print for.
 * @param[in,out] stream The <tt>std::ostream</tt> to print to.
 */

inline void printVersion(ArgumentParser const & me, std::ostream & stream)
{
    stream << getAppName(me) << " version: " << getVersion(me) << std::endl;
    stream << "SeqAn version: " << SEQAN_VERSION_MAJOR << '.' <<  SEQAN_VERSION_MINOR << '.'
           << SEQAN_VERSION_PATCH;
    if (SEQAN_VERSION_PRE_RELEASE != 0)
        stream << "-pre" << SEQAN_VERSION_PRE_RELEASE;
    stream << "\n";
}

inline void printVersion(ArgumentParser const & me)
{
    printVersion(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function printLongCopyright()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#printLongCopyright
 * @brief Prints the long copyright information of the parser to a stream.
 *
 * @signature void printLongCopyright(parser, stream);
 *
 * @param[in,out] parser The ArgumenParser to print for.
 * @param[in,out] stream The <tt>std::ostream</tt> to print to.
 */

inline void printLongCopyright(ArgumentParser const & me, std::ostream & stream)
{
    stream << "=============================================================================" << std::endl
           << "Copyright information for " << getAppName(me) << ":" << std::endl
           << "-----------------------------------------------------------------------------" << std::endl
           << me._toolDoc._longCopyright << std::endl << std::endl
           << "=============================================================================" << std::endl
           << "This program contains SeqAn code licensed under the following terms:" << std::endl
           << "-----------------------------------------------------------------------------" << std::endl
           << " Copyright (c) 2006-2018, Knut Reinert, FU Berlin" << std::endl
           << " All rights reserved." << std::endl
           << "" << std::endl
           << " Redistribution and use in source and binary forms, with or without" << std::endl
           << " modification, are permitted provided that the following conditions are met:" << std::endl
           << "" << std::endl
           << "     * Redistributions of source code must retain the above copyright" << std::endl
           << "       notice, this list of conditions and the following disclaimer." << std::endl
           << "     * Redistributions in binary form must reproduce the above copyright" << std::endl
           << "       notice, this list of conditions and the following disclaimer in the" << std::endl
           << "       documentation and/or other materials provided with the distribution." << std::endl
           << "     * Neither the name of Knut Reinert or the FU Berlin nor the names of" << std::endl
           << "       its contributors may be used to endorse or promote products derived" << std::endl
           << "       from this software without specific prior written permission." << std::endl
           << "" << std::endl
           << " THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"" << std::endl
           << " AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE" << std::endl
           << " IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE" << std::endl
           << " ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE" << std::endl
           << " FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL" << std::endl
           << " DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR" << std::endl
           << " SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER" << std::endl
           << " CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT" << std::endl
           << " LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY" << std::endl
           << " OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH" << std::endl
           << " DAMAGE." << std::endl;
}

inline void printLongCopyright(ArgumentParser const & me)
{
    printLongCopyright(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function _addNumericalRestriction()
// ----------------------------------------------------------------------------

inline void _addNumericalRestriction(std::string & text, ArgParseArgument const & arg)
{
    // expand min/max restrictions
    if (!empty(arg.minValue) || !empty(arg.maxValue))
    {
        append(text, " In range [");

        if (!empty(arg.minValue))
            append(text, arg.minValue);
        else
            append(text, "-inf");

        append(text, "..");

        if (!empty(arg.maxValue))
            append(text, arg.maxValue);
        else
            append(text, "inf");

        append(text, "].");
    }
}

// ----------------------------------------------------------------------------
// Function _expandList()
// ----------------------------------------------------------------------------

// expands the given vector as text in the form v1, v2, and v3, while respecting
// the size with respect to the used commas and "and"s
inline void _expandList(std::string & text, std::vector<std::string> const & list)
{
    for (std::vector<std::string>::size_type i = 0; i < list.size(); ++i)
    {
        if (i + 1 == list.size() && list.size() == 2u)
            append(text, " and ");
        else if (i + 1 == list.size()  && list.size() > 2u)
            append(text, ", and ");
        else if (i != 0)
            append(text, ", ");

        append(text, "\\fI");
        append(text, list[i]);
        append(text, "\\fP");

    }
}

// ----------------------------------------------------------------------------
// Function _addDefaultValues()
// ----------------------------------------------------------------------------

inline void _addDefaultValues(std::string & text, ArgParseArgument const & arg)
{
    if (!empty(arg.defaultValue))
    {
        append(text, " Default: ");
        _expandList(text, arg.defaultValue);
        append(text, ".");
    }
}

inline void _addDefaultValues(std::string & text, ArgParseOption const & arg)
{
    if (!isFlagOption(arg))
        _addDefaultValues(text, static_cast<ArgParseArgument>(arg));
}

// ----------------------------------------------------------------------------
// Function _seperateExtensionsForPrettyPrinting()
// ----------------------------------------------------------------------------

inline void _seperateExtensionsForPrettyPrinting(std::vector<std::string> & file_ext,
                                                 std::vector<std::string> & comp_ext,
                                                 std::vector<std::string> const & validValues)
{
    // seperate file extensions and compression extensions
    for (std::vector<std::string>::size_type i = 0; i < validValues.size(); ++i)
    {
        std::regex rgx("^(\\.)?([A-z0-9]+)(\\.)?([A-z0-9]+)?");
        std::smatch result;

        std::regex_search(validValues[i], result, rgx);

        if (!result[4].str().empty())
        {
            comp_ext.push_back(result[4].str());
            file_ext.push_back("." + result[2].str() + "[.*]");
        }
        else
        {
            file_ext.push_back("." + result[2].str());
        }
    }

    std::sort(file_ext.rbegin(), file_ext.rend()); // sort extensions in reverse order such that '.fa[.x]'
    std::sort(comp_ext.rbegin(), comp_ext.rend()); // comes before '.fa' and will be chosen by std::unique()

    comp_ext.erase(std::unique(comp_ext.begin(), comp_ext.end()), comp_ext.end()); // remove duplicates
    file_ext.erase(std::unique(file_ext.begin(), file_ext.end(),
                    [&](auto& lhs, auto& rhs)
                    {
                        return lhs.substr(0, lhs.find('[')) == rhs.substr(0, rhs.find('['));
                    }), file_ext.end());
}

// ----------------------------------------------------------------------------
// Function _addValidValuesRestrictions()
// ----------------------------------------------------------------------------

inline void _addValidValuesRestrictions(std::string & text, ArgParseArgument const & arg)
{
    if (!empty(arg.validValues))
    {
        if (isInputFileArgument(arg) || isOutputFileArgument(arg))
        {
            std::vector<std::string> file_extensions;
            std::vector<std::string> compresssion_extensions;

            _seperateExtensionsForPrettyPrinting(file_extensions, compresssion_extensions, arg.validValues);

            append(text, " Valid filetype");

            if (file_extensions.size() > 1)
                append(text, "s are: ");
            else
                append(text, " is: ");

            _expandList(text, file_extensions);

            if (compresssion_extensions.size() != 0)
            {
                append(text, ", where * is any of the following extensions: ");
                _expandList(text, compresssion_extensions);
                append(text, " for transparent (de)compression");
            }
        }
        else
        {
            append(text, " One of ");
            _expandList(text, arg.validValues);
        }

        append(text, ".");
    }
}

inline void _addValidValuesRestrictions(std::string & text, ArgParseOption const & opt)
{
    if (!isFlagOption(opt))
        _addValidValuesRestrictions(text, static_cast<ArgParseArgument>(opt));
}

// ----------------------------------------------------------------------------
// Function _addTypeAndListInfo()
// ----------------------------------------------------------------------------

inline void _addTypeAndListInfo(std::string & text, ArgParseArgument const & arg)
{
    std::string type = getArgumentTypeAsString(arg);
    for (auto & c: type)
         c = toupper(c);

    // Write arguments to term line -> only exception, boolean flags
    if (!empty(type))
    {
        append(text, " ");

        if (isListArgument(arg))
            append(text, "List of ");

        if (arg._numberOfValues != 1)
            append(text,  std::to_string(arg._numberOfValues) + " ");

        append(text, "\\fI");
        append(text, type);
        append(text, "\\fP");

        if (isListArgument(arg) || arg._numberOfValues != 1)
            append(text, "'s");
    }
}

// ----------------------------------------------------------------------------
// Function printHelp()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Parameter order.

/*!
 * @fn ArgumentParser#printHelp
 * @brief Prints the help message for the parser.
 *
 * @signature void printHelp(parser, out, format, showAdvancedOptions);
 *
 * @param[in,out] parser                The ArgumentParser print the help for.
 * @param[out]    out                   The output stream to print to (<tt>std::ostream</tt>).
 * @param[in]     format                The format to print, one of "html", "man", and "txt".
 * @param[in]     showAdvancedOptions   Also show advanced options to user (default = false).
 */

inline void printHelp(ArgumentParser const & me,
                      std::ostream & stream,
                      CharString const & format,
                      bool const showAdvancedOptions)
{
    ToolDoc toolDoc(me._toolDoc);
    clearEntries(toolDoc);  // We will append me._toolDoc later.

    // Build synopsis section.
    addSection(toolDoc, "Synopsis");
    _addUsage(toolDoc, me);

    // Add description to tool documentation.
    addSection(toolDoc, "Description");
    for (unsigned i = 0; i < me._description.size(); ++i)
        addText(toolDoc, me._description[i]);

    // Add arguments to arguments section
    if (length(me.argumentList) != 0)
        addSection(toolDoc, "Required Arguments");

    for (unsigned i = 0; i < length(me.argumentList); ++i)
    {
        ArgParseArgument const & arg = me.argumentList[i];

        // Build list item term.
        std::string term;
        if (!empty(arg._argumentLabel))
        {
            std::regex space(" ");
            term = "\\fB";
            append(term, std::regex_replace(arg._argumentLabel, space,"_"));
            append(term, "\\fP");
        }
        else
        {
            term = "\\fBARGUMENT ";
            append(term, std::to_string(i));
            append(term, "\\fP");
        }

        // expand type, list and numValues information
        _addTypeAndListInfo(term, arg);

        std::string helpText = arg._helpText;

        // expand min/max restrictions
        _addNumericalRestriction(helpText, arg);

        // expand validValues restrictions
        _addValidValuesRestrictions(helpText, arg);

        // expand defaultValue
        _addDefaultValues(helpText, arg);

        // Add list item.
        addListItem(toolDoc, term, helpText);
    }

    // Add options to options section.
    if (length(me.optionMap) != 0)
        addSection(toolDoc, "Options");

    for (unsigned i = 0; i < length(me.optionMap); ++i)
    {
        ArgParseOption const & opt = me.optionMap[i];
        if (empty(opt.shortName) && empty(opt.longName))  // this is not an option but a text line
        {
            if (empty(opt._helpText))  // TODO(holtgrew): Should go away in future.
                continue;  // Skip empty lines.

            // Is command line parser section, maps to ToolDoc subsection.
            for (unsigned j = i + 1; j < length(me.optionMap); ++j)
            {
                ArgParseOption const & nextopt = me.optionMap[j];
                if (empty(nextopt.shortName) && empty(nextopt.longName))
                    break;
                // has visible children
                if (!isHidden(nextopt) && (!isAdvanced(nextopt) || showAdvancedOptions))
                {
                    std::string title = opt._helpText;
                    append(title, ":");
                    addSubSection(toolDoc, title);
                    break;
                }
            }
        }
        else if (!isHidden(opt) && (!isAdvanced(opt) || showAdvancedOptions))
        {
            // Build list item term.
            std::string term;
            if (!empty(opt.shortName))
            {
                term = "\\fB-";
                append(term, opt.shortName);
                append(term, "\\fP");
            }
            if (!empty(opt.shortName) && !empty(opt.longName))
                append(term, ", ");
            if (!empty(opt.longName))
            {
                append(term, "\\fB--");
                append(term, opt.longName);
                append(term, "\\fP");
            }

            // expand type, list and numValues information
            if (!opt._isFlag)
                _addTypeAndListInfo(term, opt);

            std::string helpText = opt._helpText;

            // expand min/max restrictions
            _addNumericalRestriction(helpText, opt);

            // expand validValues restrictions
            _addValidValuesRestrictions(helpText, opt);

            // expand defaultValue
            _addDefaultValues(helpText, opt);

            // Add list item.
            addListItem(toolDoc, term, helpText);
        }
    }

    append(toolDoc, me._toolDoc);
    print(stream, toolDoc, format);
}

inline void printHelp(ArgumentParser const & me, std::ostream & stream, CharString const & format)
{
    printHelp(me, stream, format, false);
}

inline void printHelp(ArgumentParser const & me, std::ostream & stream)
{
    printHelp(me, stream, "txt", false);
}

inline void printHelp(ArgumentParser const & me)
{
    printHelp(me, std::cerr, "txt", false);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_
